
#include "equations/ChargeBalance.hpp"

void
ChargeBalance::Update()
{
  const real_t frac = particle_region == PE ? -delta_p_scale/beta_p : -delta_n_scale/beta_n;

  mfem::ConstantCoefficient coeff(frac);

  if (!K)
  {
    K = new ParBilinearForm(&fespace);
    K->AddDomainIntegrator(new DiffusionIntegrator(coeff));
    K->Assemble(0);
    Kpmat = K->ParallelAssemble();
  }

  if (!C)
  {
    C = new ParMixedBilinearForm(&fespace, fespace2);
    C->AddDomainIntegrator(new MixedScalarMassIntegrator);
    C->Assemble(0);
    Cmat = C->ParallelAssemble();
  }

  if (!Q)
  {
    Q = new ParLinearForm(&fespace);
    Q->AddDomainIntegrator(new DomainLFIntegrator(coeff));
  }
  Q->Assemble();
  Q->ParallelAssemble(b);
}
