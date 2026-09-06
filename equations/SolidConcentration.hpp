#pragma once

#include "equations/Equation.hpp"

class SolidConcentration : public Equation
{
private:
  const unsigned particle_id;
  const int particle_rank;
  const bool particle_owned;
  const int particle_dof;
  const Region particle_region;
  const mfem::Array<int> surface_bdr;
  const int surface_tdof;
  const bool surface_owned;
  const int surface_rank;

public:
  SolidConcentration(mfem::ParFiniteElementSpace & f,
                     const unsigned & id,
                     const int & rank,
                     const int & dof,
                     const Region & region)
    : Equation(f),
      particle_id(id),
      particle_rank(rank),
      particle_owned(particle_rank == mfem::Mpi::WorldRank()),
      particle_dof(dof),
      particle_region(region),
      surface_bdr({0, 1}),
      surface_tdof(FindSurfaceTrueDof()),
      surface_owned(surface_tdof != -1),
      surface_rank(FindSurfaceRank())
  {
  }

  virtual void Update(const mfem::Coefficient & j) override;
  mfem::real_t SurfaceConcentration(const mfem::BlockVector & x);
  mfem::real_t AverageConcentration(const mfem::BlockVector & x);
  int FindSurfaceTrueDof();
  int FindSurfaceRank();
  int GetParticleRank() { return particle_rank; }
  int GetSurfaceRank() { return surface_rank; }
  bool IsParticleOwned() { return particle_owned; }
  bool IsSurfaceOwned() { return surface_owned; }
  Region GetParticleRegion() { return particle_region; }
  int GetParticleDof() { return particle_dof; }
};
