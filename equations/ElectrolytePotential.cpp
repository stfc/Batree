
#include "equations/ElectrolytePotential.hpp"

void
ElectrolytePotential::Update(const mfem::GridFunctionCoefficient & ec_gfc,
                             const mfem::Coefficient & j)
{
  // Source term.
  mfem::Vector source_vec({/* NE */ AN /* length scaling */ * (LNE / NNE * NX),
                           /* SEP */ 0.,
                           /* PE */ AP /* length scaling */ * (LPE / NPE * NX)});

  mfem::PWConstCoefficient source_part(source_vec);
  mfem::ProductCoefficient source(source_part, const_cast<mfem::Coefficient &>(j));

  mfem::Vector b_vec({/* NE */ BNE /* length scaling */ / (LNE / NNE * NX),
                      /* SEP */ BSEP /* length scaling */ / (LSEP / NSEP * NX),
                      /* PE */ BPE /* length scaling */ / (LPE / NPE * NX)});

  mfem::PWConstCoefficient b_part(b_vec);

  mfem::TransformedCoefficient kappa(&const_cast<mfem::GridFunctionCoefficient &>(ec_gfc), Kappa);
  mfem::ProductCoefficient kappa_eff(b_part, kappa);

  mfem::GradientGridFunctionCoefficient grad_ec(ec_gfc.GetGridFunction());
  mfem::RatioCoefficient ec_inv(2 * T * (1 - TPLUS),
                                const_cast<mfem::GridFunctionCoefficient &>(ec_gfc));
  mfem::ScalarVectorProductCoefficient grad_ln_ec(ec_inv, grad_ec);
  mfem::ScalarVectorProductCoefficient grad_ln_ec_kappad(kappa_eff, grad_ln_ec);

  if (!K)
  {
    K = new mfem::ParBilinearForm(&fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(kappa_eff));
    K->Assemble();
    K->FormSystemMatrix(ess_tdof_list, Kmat);
  }

  if (!Qc)
  {
    Qc = new mfem::ParLinearForm(&fespace);
    Qc->AddDomainIntegrator(new mfem::DomainLFGradIntegrator(grad_ln_ec_kappad));
    Qc->Assemble();
    Qc->ParallelAssemble(bc);
  }

  if (!Q)
  {
    Q = new mfem::ParLinearForm(&fespace);
    Q->AddDomainIntegrator(new mfem::DomainLFIntegrator(source));
  }
  Q->Assemble();
  Q->ParallelAssemble(b);
  b += bc;
  b.SetSubVector(ess_tdof_list, 0.0);
}
