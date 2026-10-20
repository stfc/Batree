
#include "equations/ChargeBalance.hpp"

void
ChargeBalance::Update()
{
  const real_t frac_one = particle_region == PE ? -delta_p_scale / beta_p : delta_n_scale / beta_n;
  const real_t frac_two =
      particle_region == PE ? -beta_p * (I / Acs_p_scale) : -beta_n * (I / Acs_n_scale);

  mfem::ConstantCoefficient coeff(frac_one);
  mfem::ConstantCoefficient proc(frac_one * frac_two);
  mfem::ConstantCoefficient sign(particle_region == PE ? 1.0 : -1.0);

  if (!K)
  {
    K = new ParBilinearForm(&fespace);
    K->AddDomainIntegrator(new DiffusionIntegrator(coeff));
    K->Assemble(0);
    K->FormSystemMatrix(ess_tdof_list, Kmat);
  }

  if (!C)
  {
    C = new ParMixedBilinearForm(&fespace, &l2_fespace);
    C->AddDomainIntegrator(new MixedScalarMassIntegrator(sign));
    C->Assemble(0);
    C->Finalize();
    Cmat = C->ParallelAssemble();
  }

  ParMesh *pmesh = fespace.GetParMesh();

  const int max_bdr_attr = pmesh->bdr_attributes.Max();
  Array<int> bdr_marker(max_bdr_attr);

 if (!Q)
  {
    Q = new ParLinearForm(&fespace);

    //bdr_marker.UseDevice(false);  // Here
    bdr_marker = 0;

    // Gmsh physical group "protrusion_top" has boundary attribute 3
    bdr_marker[3 - 1] = 1;

    Q->AddBoundaryIntegrator(
        new BoundaryLFIntegrator(proc), bdr_marker);
  }

Q->Assemble();
Q->ParallelAssemble(b);
}

ChargeBalance::~ChargeBalance()
{
  delete C;
  delete Cmat;
}
