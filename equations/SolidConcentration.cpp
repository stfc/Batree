#include "equations/SolidConcentration.hpp"

void
SolidConcentration::Update(const mfem::Coefficient & j)
{
  MFEM_ASSERT(particle_region == NE || particle_region == PE, "Particle not in electrode!");

  const mfem::real_t Rpar = particle_region == NE ? RN : RP;
  const mfem::real_t D = particle_region == NE ? DN : DP;
  const mfem::real_t t_scale = particle_region == NE ? tn_scale : tp_scale;

  mfem::FunctionCoefficient r2([](const mfem::Vector & r) { return r(0) * r(0); });
  mfem::ProductCoefficient dr2(D / Rpar / Rpar, r2);
  mfem::ProductCoefficient jjr2(const_cast<mfem::Coefficient &>(j), r2);
  mfem::ProductCoefficient jr2(-1. / Rpar / t_scale, jjr2);

  if (!M)
  {
    M = new mfem::ParBilinearForm(&fespace);
    M->AddDomainIntegrator(new mfem::MassIntegrator(r2));
    M->Assemble(0); // keep sparsity pattern of M and K the same
    M->FormSystemMatrix(ess_tdof_list, Mmat);
  }

  if (!K) // revisit if nonlinear
  {
    K = new mfem::ParBilinearForm(&fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(dr2));
    K->Assemble(0); // keep sparsity pattern of M and K the same
    K->FormSystemMatrix(ess_tdof_list, Kmat);
  }

  if (!Q)
  {
    Q = new mfem::ParLinearForm(&fespace);
    Q->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(jr2),
                             const_cast<mfem::Array<int> &>(surface_bdr));
  }
  Q->Assemble();
  Q->ParallelAssemble(b);
}

mfem::real_t
SolidConcentration::SurfaceConcentration(const mfem::BlockVector & x)
{
  mfem::real_t csurf = IsSurfaceOwned() ? x.GetBlock(SC + particle_id)[surface_tdof] : 0;

  if (GetParticleRank() == GetSurfaceRank())
    return csurf;

  if (IsSurfaceOwned())
    MPI_Send(&csurf, 1, MFEM_MPI_REAL_T, particle_rank, particle_id, MPI_COMM_WORLD);
  if (IsParticleOwned())
    MPI_Recv(
        &csurf, 1, MFEM_MPI_REAL_T, surface_rank, particle_id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  return csurf;
}

mfem::real_t
SolidConcentration::AverageConcentration(const mfem::BlockVector & x)
{
  mfem::ParGridFunction sc_gf(&fespace);
  sc_gf.SetFromTrueDofs(x.GetBlock(SC + particle_id));
  mfem::GridFunctionCoefficient sc_gfc(&sc_gf);
  mfem::FunctionCoefficient r2([=](const mfem::Vector & x) { return x(0) * x(0); });
  mfem::ProductCoefficient scr2(sc_gfc, r2);

  mfem::QuadratureSpace x_qspace(fespace.GetParMesh(), fespace.FEColl()->GetOrder() + 2);
  return 3 * x_qspace.Integrate(scr2);
}

int
SolidConcentration::FindSurfaceTrueDof()
{
  mfem::Array<int> nat_dofs;
  fespace.GetEssentialTrueDofs(surface_bdr, nat_dofs);
  return nat_dofs.IsEmpty() ? -1 : nat_dofs[0];
}

int
SolidConcentration::FindSurfaceRank()
{
  mfem::Array<bool> is_surface_rank(mfem::Mpi::WorldSize());
  MPI_Allgather(
      &surface_owned, 1, MPI_CXX_BOOL, is_surface_rank.GetData(), 1, MPI_CXX_BOOL, MPI_COMM_WORLD);
  return std::distance(is_surface_rank.begin(),
                       std::find(is_surface_rank.begin(), is_surface_rank.end(), true));
}
