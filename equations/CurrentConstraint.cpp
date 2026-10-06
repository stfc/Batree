
#include "equations/CurrentConstraint.hpp"

void
CurrentConstraint::Update()
{

  if (!M)
  {
    M = new ParBilinearForm(fespace2);
    M->AddDomainIntegrator(new MassIntegrator);
    M->Assemble(0); // keep sparsity pattern of M and K the same
    Mpmat = M->ParallelAssemble();
  }

  HypreParVector one(*Mpmat);
  one = 1.0;

  Col = HypreParVector(fespace2->GetComm(), fespace2->GlobalTrueVSize(), fespace2->GetTrueDofOffsets());
  Mpmat->Mult(one, Col);


  if (!Q)
  {
    Q = new ParLinearForm(&fespace);
    mfem::ConstantCoefficient one(1.0);
    Q->AddDomainIntegrator(new DomainLFIntegrator(one));
  }

  Q->Assemble();
  Q->ParallelAssemble(b);
  
}
