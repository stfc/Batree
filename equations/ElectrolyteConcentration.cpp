
#include "equations/ElectrolyteConcentration.hpp"

void
ElectrolyteConcentration::Update(const mfem::GridFunctionCoefficient & ec_gfc,
                                 const mfem::Coefficient & j)
{
  // Mass coefficient.
  mfem::Vector mass_vec({/* NE */ EPS_N /* length scaling */ * (LNE / NNE * NX),
                         /* SEP */ EPS_S /* length scaling */ * (LSEP / NSEP * NX),
                         /* PE */ EPS_P /* length scaling */ * (LPE / NPE * NX)});
  mfem::PWConstCoefficient mass_part(mass_vec);
  mfem::ProductCoefficient mass(te_scale, mass_part);

  // Source term.
  mfem::Vector source_vec({/* NE */ (1 - TPLUS) * AN /* length scaling */ * (LNE / NNE * NX),
                           /* SEP */ 0.,
                           /* PE */ (1 - TPLUS) * AP /* length scaling */ * (LPE / NPE * NX)});

  mfem::PWConstCoefficient source_part(source_vec);
  mfem::ProductCoefficient source(source_part, const_cast<mfem::Coefficient &>(j));

  // Diffusion coefficient.
  mfem::Vector D_scale_vec({/* NE */ BNE /* length scaling */ / (LNE / NNE * NX),
                            /* SEP */ BSEP /* length scaling */ / (LSEP / NSEP * NX),
                            /* PE */ BPE /* length scaling */ / (LPE / NPE * NX)});

  mfem::TransformedCoefficient D_coeff(&const_cast<mfem::GridFunctionCoefficient &>(ec_gfc), DE);
  mfem::PWConstCoefficient D_scale_coeff(D_scale_vec);
  mfem::ProductCoefficient D(D_scale_coeff, D_coeff);

  if (!M)
  {
    M = new mfem::ParBilinearForm(&fespace);
    M->AddDomainIntegrator(new mfem::MassIntegrator(mass));
    M->Assemble(0); // keep sparsity pattern of M and K the same
    M->FormSystemMatrix(ess_tdof_list, Mmat);
  }

  delete K;
  K = new mfem::ParBilinearForm(&fespace);
  K->AddDomainIntegrator(new mfem::DiffusionIntegrator(D));
  K->Assemble(0); // keep sparsity pattern of M and K the same
  K->FormSystemMatrix(ess_tdof_list, Kmat);

  if (!Q)
  {
    Q = new mfem::ParLinearForm(&fespace);
    Q->AddDomainIntegrator(new mfem::DomainLFIntegrator(source));
  }
  Q->Assemble();
  Q->ParallelAssemble(b);
}
