
#include "equations/CurrentConstraint.hpp"

void
CurrentConstraint::Update()
{

  if (!C)
  {
    C = new ParMixedBilinearForm(&fespace, fespace2);
    C->AddDomainIntegrator(new MixedScalarMassIntegrator);
    C->Assemble();
  }

  if (!Q)
  {
    Q = new ParLinearForm(&fespace);
    mfem::ConstantCoefficient one(1.0);
    Q->AddDomainIntegrator(new DomainLFIntegrator(one));
  }

  Q->Assemble();
  Q->ParallelAssemble(b);
  
}
