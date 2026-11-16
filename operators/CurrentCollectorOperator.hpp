#pragma once

#include "equations/ChargeBalance.hpp"
#include "equations/CurrentConstraint.hpp"
#include "operators/LocalModel.hpp"

#include <memory>

/// Summary of one SolveCoupled call
struct CouplingInfo
{
  int passes = 0;
  bool converged = false;
  mfem::real_t residual = 0.0;         // relative residual of the coupled equations at the answer
  mfem::real_t imin = 0.0, imax = 0.0; // range of the element currents
  mfem::real_t rmin = 0.0, rmax = 0.0; // range of the local slopes R̄ the model last returned
};

// remove inheritence from Operator class
//  possibly rename bcc by acc
class CurrentCollectorOperator
{
protected:
  /// 2d serial meshes of the positive and negative current collectors (domain plus tab)
  mfem::Mesh _pos_smesh;
  mfem::Mesh _neg_smesh;

  /// Element-to-rank maps of the two serial meshes. Each domain element goes to the same rank in
  /// both collectors, so that on every rank element e of the current space is element e of both
  /// H1 spaces (the coupling blocks C pair elements by their local number).
  mfem::Array<int> _pos_partitioning;
  mfem::Array<int> _neg_partitioning;

  /// Parallel meshes defined by a partitioning of the respective serial mesh
  mfem::ParMesh _pos_pmesh;
  mfem::ParMesh _neg_pmesh;

  /// Collector domain without the tab (attribute 1), where the element currents live
  mfem::ParSubMesh _dom_pmesh;

  /// FE collections and spaces: H1 potentials in each collector, L2 element currents
  mfem::H1_FECollection _h1_coll;
  mfem::L2_FECollection _l2_coll;
  mfem::ParFiniteElementSpace h1_fespace_current_collector_positive;
  mfem::ParFiniteElementSpace h1_fespace_current_collector_negative;
  mfem::ParFiniteElementSpace l2_fespace_currents;

  ChargeBalance * _cbp = nullptr; // charge balance for positive current collector
  ChargeBalance * _cbn = nullptr; // charge balance for negative current collector
  CurrentConstraint * _cc = nullptr;

  mfem::real_t current_dt;

  mfem::CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
  mfem::HypreSmoother M_prec; // Preconditioner for the mass matrix M

  // CGSolver _Solver;      // Implicit solver for T = M + dt K
  mfem::GMRESSolver _Solver;   // GMRES solver for the block system
  mfem::HypreSmoother _I_prec; // Preconditioner for the implicit solver

  mutable mfem::Vector b; // auxiliary vector

  /// For solution true vector (4 MSMD equations)
  mfem::Array<int> _block_trueOffsets;

  /// System matrices for concentration and potential eqs
  mfem::HypreParMatrix * _Acc = nullptr;

  /// Auxiliary rhs vectors for MSMD equations
  mutable mfem::BlockVector _f;

  // mutable Array2D<const HypreParMatrix *> _Bcc{4, 4};
  mutable mfem::Array2D<const mfem::HypreParMatrix *> _Bcc{3, 3}; // reference potential s removed

  /// (2,2) block D = M diag(R̄) of the coupled system. Owned here; _Bcc(2,2) points to it.
  mfem::HypreParMatrix * _D = nullptr;

  /// Element currents ī fed to the local cell models (true dofs of the L2 space).
  /// Kept between calls, so a later time step can start from the previous answer.
  mfem::Vector _ibar;

  /// Print the block, RHS and solution diagnostics in Solve()
  bool _verbose = true;

  /// Reference to solution true dof vector
  mfem::BlockVector & _x;

  /// ParaView output, created on the first SaveParaView call: the grid functions for φ⁺, φ⁻
  /// and i, and one data collection per mesh
  std::unique_ptr<mfem::ParGridFunction> _phi_pos_gf, _phi_neg_gf, _i_gf;
  std::unique_ptr<mfem::ParaViewDataCollection> _pv_pos, _pv_neg, _pv_dom;

public:
  /// order: FE order of the collector potentials. refinements: uniform refinements of both
  /// collector meshes (0 uses the meshes as read).
  CurrentCollectorOperator(int order, int refinements, mfem::BlockVector & x);

  /// One linear solve of the 3x3 system. _Bcc(2,2) and block 2 of the RHS must be set
  /// first (SolveCoupled does this). Returns false if the iterative solver (used when MFEM has
  /// no SuperLU) did not converge.
  virtual bool Solve();

  /// Couple the collectors to a local cell model V̄(ī), at one instant in time.
  /// Each pass solves the 3x3 system with D = M diag(R̄) in the (2,2) block and
  /// M (V̄ + R̄ ī) in block-row 2 of the RHS, then sets ī ← (1 − ω) ī + ω î.
  ///   slope = "exact"  : R̄ recomputed every pass (Newton)
  ///           "frozen" : R̄ from the first pass only (chord)
  ///           "none"   : R̄ = 0, the report's Option 1 (lagged Picard)
  /// The model is asked for R̄ only on the passes that use it. report = false prints nothing.
  CouplingInfo SolveCoupled(CellModel & model,
                            const std::string & slope,
                            mfem::real_t omega = 1.0,
                            mfem::real_t tol = 1e-10,
                            int max_it = 40,
                            bool report = true);

  /// Terminal voltage (scaled) of the last solution: the mean of φ⁺ over the positive tab, where
  /// the current enters uniformly. φ⁻ is 0 on the negative tab.
  mfem::real_t TerminalVoltage() const;

  /// Number of local collector elements, and their areas (scaled)
  int NumElements() const { return l2_fespace_currents.GetTrueVSize(); }
  void ElementAreas(mfem::Vector & area) const;

  /// Write φ⁺ and φ⁻ (H1, on the positive and negative collector meshes) and i (L2, on the
  /// tab-free domain) to ParaView/, one collection per mesh. cycle and time label this output
  /// in the time series.
  void SaveParaView(int cycle, mfem::real_t time);

  virtual ~CurrentCollectorOperator();
};
