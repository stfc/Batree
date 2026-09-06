
#include "equations/SolidPotential.hpp"

void
SolidPotential::Update(const mfem::Coefficient & j)
{
  // Source term.
  mfem::Vector source_vec({/* NE */ -AN /* length scaling */ * (LNE / NNE * NX),
                           /* SEP */ 0.,
                           /* PE */ -AP /* length scaling */ * (LPE / NPE * NX)});

  mfem::PWConstCoefficient source_part(source_vec);
  mfem::ProductCoefficient source(source_part, const_cast<mfem::Coefficient &>(j));

  // Effective conductivity (does not account for electrode filler).
  mfem::Vector sigma_vec({/* NE */ (1 - EPS_N) * SIGN /* length scaling */ / (LNE / NNE * NX),
                          /* SEP */ 0.,
                          /* PE */ (1 - EPS_P) * SIGP /* length scaling */ / (LPE / NPE * NX)});

  mfem::PWConstCoefficient sigma(sigma_vec);

  if (!K)
  {
    K = new mfem::ParBilinearForm(&fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(sigma));
    K->Assemble();
    K->FormSystemMatrix(ess_tdof_list, Kmat);
  }

  if (!Q)
  {
    Q = new mfem::ParLinearForm(&fespace);
    Q->AddDomainIntegrator(new mfem::DomainLFIntegrator(source));
  }
  Q->Assemble();
  Q->ParallelAssemble(b);
  b.SetSubVector(ess_tdof_list, 0.0);
}
