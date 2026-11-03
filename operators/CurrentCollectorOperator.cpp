#include "operators/CurrentCollectorOperator.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>

// BATREE_DIR is set by the makefile to the project folder (-DBATREE_DIR=...).
// When built without the makefile, fall back to paths relative to where batree is run.
#ifndef BATREE_DIR
#define BATREE_DIR "."
#endif

namespace
{
/// Domain (tab-free) attribute of both collector meshes
const int DOMAIN_ATTRIBUTE = 1;

/// Load a collector mesh from <project folder>/mesh and refine it uniformly. Refinement splits each
/// element into children numbered consecutively, so the domain elements stay first and in the
/// same order in both collectors.
mfem::Mesh
LoadCollectorMesh(const std::string & file, int refinements)
{
  mfem::Mesh mesh(std::string(BATREE_DIR) + "/mesh/" + file, 1, 0);
  for (int r = 0; r < refinements; r++)
    mesh.UniformRefinement();
  return mesh;
}

/// Partition the positive collector with METIS, as ParMesh does by default.
mfem::Array<int>
PartitionPositive(mfem::Mesh & pos)
{
  int * p = pos.GeneratePartitioning(mfem::Mpi::WorldSize());
  mfem::Array<int> partitioning(pos.GetNE());
  for (int e = 0; e < pos.GetNE(); e++)
    partitioning[e] = p[e];
  delete[] p;
  return partitioning;
}

/// Partition the negative collector to match the positive one: each domain element goes to the
/// rank of the positive domain element at the same place, each tab element to the rank of the
/// nearest positive domain element.
mfem::Array<int>
PartitionNegative(mfem::Mesh & pos, const mfem::Array<int> & pos_partitioning, mfem::Mesh & neg)
{
  mfem::Array<int> pos_domain;
  for (int e = 0; e < pos.GetNE(); e++)
    if (pos.GetAttribute(e) == DOMAIN_ATTRIBUTE)
      pos_domain.Append(e);

  mfem::Array<int> partitioning(neg.GetNE());
  mfem::Vector c(neg.Dimension()), cp(pos.Dimension());
  for (int e = 0; e < neg.GetNE(); e++)
  {
    neg.GetElementCenter(e, c);
    int nearest = -1;
    mfem::real_t dist = mfem::infinity();
    for (int ep : pos_domain)
    {
      pos.GetElementCenter(ep, cp);
      if (c.DistanceTo(cp) < dist)
      {
        dist = c.DistanceTo(cp);
        nearest = ep;
      }
    }
    if (neg.GetAttribute(e) == DOMAIN_ATTRIBUTE && dist > 1e-8 * pos.GetElementSize(nearest))
      mfem::mfem_error("CurrentCollectorOperator: the domains of the positive and negative "
                       "collector meshes do not match.");
    partitioning[e] = pos_partitioning[nearest];
  }
  return partitioning;
}

/// The coupling blocks C (ParMixedBilinearForm) pair element e of the current space with element e
/// of an H1 space. Check that they are the same element on this rank.
void
CheckSameElements(mfem::ParMesh & currents, mfem::ParMesh & collector)
{
  mfem::Vector c(currents.Dimension()), cc(collector.Dimension());
  for (int e = 0; e < currents.GetNE(); e++)
  {
    bool same = e < collector.GetNE();
    if (same)
    {
      currents.GetElementCenter(e, c);
      collector.GetElementCenter(e, cc);
      same = c.DistanceTo(cc) <= 1e-8 * currents.GetElementSize(e);
    }
    if (!same)
      mfem::mfem_error("CurrentCollectorOperator: the current space and a collector space do not "
                       "share their element order on this rank.");
  }
}
}

CurrentCollectorOperator::CurrentCollectorOperator(int order,
                                                   int refinements,
                                                   mfem::BlockVector & x)
  : _pos_smesh(LoadCollectorMesh("rectangle_current_collector_positive_scaled.msh", refinements)),
    _neg_smesh(LoadCollectorMesh("rectangle_current_collector_negative_scaled.msh", refinements)),
    _pos_partitioning(PartitionPositive(_pos_smesh)),
    _neg_partitioning(PartitionNegative(_pos_smesh, _pos_partitioning, _neg_smesh)),
    _pos_pmesh(MPI_COMM_WORLD, _pos_smesh, _pos_partitioning.GetData()),
    _neg_pmesh(MPI_COMM_WORLD, _neg_smesh, _neg_partitioning.GetData()),
    _dom_pmesh(
        mfem::ParSubMesh::CreateFromDomain(_pos_pmesh, mfem::Array<int>({DOMAIN_ATTRIBUTE}))),
    _h1_coll(order, _pos_smesh.Dimension()),
    _l2_coll(0, _pos_smesh.Dimension()),
    h1_fespace_current_collector_positive(&_pos_pmesh, &_h1_coll),
    h1_fespace_current_collector_negative(&_neg_pmesh, &_h1_coll),
    l2_fespace_currents(&_dom_pmesh, &_l2_coll),
    _Solver(MPI_COMM_WORLD),
    _x(x)
{
  CheckSameElements(_dom_pmesh, _pos_pmesh);
  CheckSameElements(_dom_pmesh, _neg_pmesh);

  // GMRES must solve tightly: Newton on the element currents stops at 1e-10, and at 1e-8 the
  // total current is only right to about 2e-6.
  const mfem::real_t rel_tol = 1e-12;

  _Solver.iterative_mode = false;
  _Solver.SetRelTol(rel_tol);
  _Solver.SetAbsTol(0.0);
  _Solver.SetMaxIter(500);
  _Solver.SetPrintLevel(0); // only warn if GMRES does not converge
  _Solver.SetKDim(50);
  _Solver.SetPreconditioner(_I_prec);

  _block_trueOffsets.SetSize(settings::NMSMDEQS + 1);
  _block_trueOffsets[0] = 0;
  _block_trueOffsets[1] = h1_fespace_current_collector_positive.GetTrueVSize();
  _block_trueOffsets[2] = h1_fespace_current_collector_negative.GetTrueVSize();
  _block_trueOffsets[3] = l2_fespace_currents.GetTrueVSize();
  // _block_trueOffsets[4] = 1;

  // std::cout << "size of _block_trueOffsets = " << _block_trueOffsets[0] << ", " <<
  // _block_trueOffsets[1] << ", " << _block_trueOffsets[2] << ", "
  //           << _block_trueOffsets[3] << ", " << _block_trueOffsets[4] << std::endl;
  std::cout << "size of _block_trueOffsets = " << _block_trueOffsets[0] << ", "
            << _block_trueOffsets[1] << ", " << _block_trueOffsets[2] << ", "
            << _block_trueOffsets[3] << std::endl;

  _block_trueOffsets.PartialSum();

  _x.Update(_block_trueOffsets);
  _f.Update(_block_trueOffsets);
  // How to distinguish between the positive and negative current collector?
  _cbp = new ChargeBalance(h1_fespace_current_collector_positive, l2_fespace_currents, PE);
  _cbn = new ChargeBalance(h1_fespace_current_collector_negative, l2_fespace_currents, NE);
  _cc = new CurrentConstraint(l2_fespace_currents);

  _Bcc = nullptr;

  std::cout << "CurrentCollectorOperator constructed" << std::endl;
  _cbp->Update();
  _cbn->Update();
  mfem::ConstantCoefficient voltage(4.2 / phi_scale);
  // Update assembles the L2 mass matrix M used by SolveCoupled. Its RHS (M · 4.2/phi_scale)
  // is no longer used: block-row 2 of the RHS now comes from the local model.
  _cc->Update(voltage);
  // Solve();
  // std::cout << "Solve completed" << std::endl;
  // The solve is now called from main through SolveCoupled(), which sets the (2,2) block
  // and block-row 2 of the RHS before each linear solve.
}
// replace Mult with Solve, does not take any parameters!
bool
CurrentCollectorOperator::Solve()
{
  _Bcc(0, 0) = &_cbp->GetK();
  delete _Bcc(0, 2);
  _Bcc(0, 2) = _cbp->GetC().Transpose();

  _f.GetBlock(0) = _cbp->GetZ();

  _Bcc(1, 1) = &_cbn->GetK();
  delete _Bcc(1, 2);
  _Bcc(1, 2) = _cbn->GetC().Transpose();

  _f.GetBlock(1) = _cbn->GetZ();

  _Bcc(2, 0) = &_cbp->GetC();
  _Bcc(2, 1) = &_cbn->GetC();
  // 3x3 system: the reference potential s (block 3) and its row Pᵀ i = I are removed.
  // The total current is now set by the tab flux alone: Pᵀ i = ΣZp.
  // _Bcc(2, 3) = &_cc->GetProduct();

  // _f.GetBlock(2) = _cc->GetZ();
  // Block-row 2 of the RHS, M (V̄ + R̄ ī), and _Bcc(2, 2) = D are set by SolveCoupled().

  //  delete _Bcc(3, 2);
  // _Bcc(3, 2) = _cc->GetProduct().Transpose();

  // _f.GetBlock(3) = I;

  // ============================================================
  // Construct the global block matrix
  // ============================================================

  delete _Acc;
  _Acc = mfem::HypreParMatrixFromBlocks(_Bcc);

#ifdef MFEM_USE_SUPERLU
  // Exact solve with SuperLU_DIST when MFEM has it: works at any mesh refinement.
  mfem::SuperLURowLocMatrix A_superlu(*_Acc);
  mfem::SuperLUSolver direct_solver(_Acc->GetComm());
  direct_solver.SetPrintStatistics(false);
  direct_solver.SetOperator(A_superlu);
  direct_solver.Mult(_f, _x);
  return true;
#else
  // Otherwise GMRES. With this preconditioner it converges on the default collector mesh and one
  // refinement of it, but not on finer meshes.
  _Solver.SetOperator(*_Acc);
  _Solver.Mult(_f, _x);
  return _Solver.GetConverged();
#endif
}

void
CurrentCollectorOperator::SolveCoupled(const LocalModel & model,
                                       const std::string & slope,
                                       mfem::real_t omega,
                                       mfem::real_t tol,
                                       int max_it)
{
  if (slope != "exact" && slope != "frozen" && slope != "none")
    mfem::mfem_error("SolveCoupled: slope must be 'exact', 'frozen' or 'none'.");

  MPI_Comm comm = l2_fespace_currents.GetComm();
  const int n = l2_fespace_currents.GetTrueVSize();

  // L2 (P0) mass matrix: diagonal, with M(e,e) = area of element e
  const mfem::HypreParMatrix & M = _cc->GetM();

  auto global = [comm](mfem::real_t x, MPI_Op op)
  {
    MPI_Allreduce(MPI_IN_PLACE, &x, 1, MFEM_MPI_REAL_T, op, comm);
    return x;
  };

  // Element areas: area = M 1
  mfem::Vector ones(n), area(n);
  ones = 1.0;
  M.Mult(ones, area);

  // Initial guess: the tab current spread evenly, ī = ΣZp / A
  if (_ibar.Size() != n)
  {
    const mfem::real_t I_tab = global(_cbp->GetZ().Sum(), MPI_SUM);
    const mfem::real_t A = global(area.Sum(), MPI_SUM);
    _ibar.SetSize(n);
    _ibar = I_tab / A;
  }

  mfem::Vector Vbar, Rbar, Rused(n), w(n), ihat(n);
  Rused = 0.0;

  if (mfem::Mpi::Root())
  {
    std::printf(
        "\n===== Collector-cell coupling: slope = %s, omega = %g =====\n", slope.c_str(), omega);
    std::printf(" pass  max|i_hat - i_bar|   integral i_hat     min i_hat     max i_hat\n");
  }

  bool converged = false;
  bool linear_solve_failed = false;
  for (int it = 0; it < max_it; it++)
  {
    // 1. Local cell models at the current guess: V̄ₑ(īₑ) and R̄ₑ = −dV̄ₑ/dī
    model.Eval(_ibar, Vbar, Rbar);

    // 2. Slope used in the matrix. "none" leaves it at 0 (the report's Option 1).
    if (slope == "exact" || (slope == "frozen" && it == 0))
      Rused = Rbar;

    // 3. (2,2) block D = M diag(R_used). M is diagonal, so scaling its rows gives exactly
    //    M diag(R_used), and D (and the whole 3x3 matrix) stays symmetric.
    delete _D;
    _D = new mfem::HypreParMatrix(M);
    _D->ScaleRows(Rused);
    _Bcc(2, 2) = _D;

    // 4. Block-row 2 of the RHS: M (V̄ + R_used ī). The R_used terms cancel when î = ī,
    //    so the converged answer has φ⁺ − φ⁻ = V̄(i) on every element, whatever slope is used.
    for (int e = 0; e < n; e++)
      w(e) = Vbar(e) + Rused(e) * _ibar(e);
    M.Mult(w, _f.GetBlock(2));

    // 5. Linear solve of the 3x3 system (diagnostics printed on the first pass only)
    _verbose = (it == 0);
    const bool solved = Solve();

    // 6. Convergence check on the element currents
    ihat = _x.GetBlock(2);
    mfem::real_t err = 0.0;
    for (int e = 0; e < n; e++)
      err = std::max(err, std::abs(ihat(e) - _ibar(e)));
    err = global(err, MPI_MAX);
    const mfem::real_t total = mfem::InnerProduct(comm, area, ihat);
    const mfem::real_t imin = global(n ? ihat.Min() : mfem::infinity(), MPI_MIN);
    const mfem::real_t imax = global(n ? ihat.Max() : -mfem::infinity(), MPI_MAX);

    if (mfem::Mpi::Root())
      std::printf("%5d  %18.3e  %15.8f  %12.6f  %12.6f\n", it + 1, err, total, imin, imax);

    // An unconverged linear solve gives a wrong î, so Newton cannot go on from it
    if (!solved)
    {
      linear_solve_failed = true;
      break;
    }

    // 7. Update: ī ← (1 − ω) ī + ω î
    for (int e = 0; e < n; e++)
      _ibar(e) = (1.0 - omega) * _ibar(e) + omega * ihat(e);

    if (err < tol)
    {
      converged = true;
      break;
    }
    if (!std::isfinite(err) || err > 1e8) // diverging (expected for slope = "none")
      break;
  }

  // Residual of the coupled equations, block-row 2: Cp φ⁺ + Cn φ⁻ − M V̄(î).
  // This is the true coupled problem, with no R̄ terms, so it checks the answer itself.
  mfem::Vector r(n);
  model.Eval(ihat, Vbar, Rbar);
  M.Mult(Vbar, r);
  const mfem::real_t scale = global(r.Normlinf(), MPI_MAX);
  _cbp->GetC().Mult(1.0, _x.GetBlock(0), -1.0, r);
  _cbn->GetC().Mult(1.0, _x.GetBlock(1), 1.0, r);
  const mfem::real_t res = global(r.Normlinf(), MPI_MAX) / scale;
  const mfem::real_t rmin = global(n ? Rbar.Min() : mfem::infinity(), MPI_MIN);
  const mfem::real_t rmax = global(n ? Rbar.Max() : -mfem::infinity(), MPI_MAX);

  if (mfem::Mpi::Root())
  {
    std::printf("Coupling %s; relative coupled residual = %.3e\n",
                converged             ? "converged"
                : linear_solve_failed ? "did NOT converge (the linear solver failed)"
                                      : "did NOT converge",
                res);
    std::printf("Local slope R at the answer: %g .. %g\n", rmin, rmax);
  }

  std::cout << "\n===== Converged solution =====\n";
  for (int i = 0; i < 3; i++)
  {
    std::cout << "\n===== Converged solution block " << i << " =====\n";
    _x.GetBlock(i).Print(std::cout);
  }
  std::cout << "Solve completed" << std::endl;
}

CurrentCollectorOperator::~CurrentCollectorOperator()
{
  // for (int i = 0; i < _Bcc.NumRows(); i++)
  //   for (int j = 0; j < _Bcc.NumCols(); j++)
  //     delete _Bcc(i, j);
  // Only some blocks are owned here. (0,0), (1,1), (2,0) and (2,1) point at Kmat/Cmat inside
  // _cbp/_cbn, which free them; deleting them here freed memory that was never allocated
  // on its own. (2,2) is _D, deleted below.
  delete _Bcc(0, 2);
  delete _Bcc(1, 2);
  delete _D;

  delete _Acc;

  delete _cbp;
  delete _cbn;
  delete _cc;
}
