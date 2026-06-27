
#include "equations/ElectrolytePotential.hpp"

void
ElectrolytePotential::Update(const GridFunctionCoefficient & ec_gfc, const Coefficient & j)
{
  // Source term.
  Vector source_vec({/* NE */ AN /* length scaling */ * (LNE / NNE * NX),
                     /* SEP */ 0.,
                     /* PE */ AP /* length scaling */ * (LPE / NPE * NX)});

  PWConstCoefficient source_part(source_vec);
  ProductCoefficient source(source_part, const_cast<Coefficient &>(j));

  Vector b_vec({/* NE */ BNE /* length scaling */ / (LNE / NNE * NX),
                /* SEP */ BSEP /* length scaling */ / (LSEP / NSEP * NX),
                /* PE */ BPE /* length scaling */ / (LPE / NPE * NX)});

  PWConstCoefficient b_part(b_vec);

  TransformedCoefficient kappa(&const_cast<GridFunctionCoefficient &>(ec_gfc), Kappa);
  ProductCoefficient kappa_eff(b_part, kappa);

  GradientGridFunctionCoefficient grad_ec(ec_gfc.GetGridFunction());
  RatioCoefficient ec_inv(1., const_cast<GridFunctionCoefficient &>(ec_gfc));
  ScalarVectorProductCoefficient grad_ln_ec(ec_inv, grad_ec);
  ScalarVectorProductCoefficient prod_part(kappa_eff, grad_ln_ec);
  ScalarVectorProductCoefficient grad_ln_ec_kappad(2 * T * (1 - TPLUS), prod_part);

  delete K;
  K = new ParBilinearForm(&fespace);
  K->AddDomainIntegrator(new DiffusionIntegrator(kappa_eff));
  K->Assemble();
  K->FormSystemMatrix(ess_tdof_list, Kmat);

  if (!Q)
  {
    Q = new ParLinearForm(&fespace);
    Q->AddDomainIntegrator(new DomainLFIntegrator(source));
    Q->AddDomainIntegrator(new DomainLFGradIntegrator(grad_ln_ec_kappad));
  }
  Q->Assemble();
  Q->ParallelAssemble(b);
  b.SetSubVector(ess_tdof_list, 0.0);
}
