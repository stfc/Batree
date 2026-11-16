#pragma once

#include "equations/Equation.hpp"

class ChargeBalance : public Equation
{
protected:
  mfem::ParFiniteElementSpace & l2_fespace;
  const Region particle_region;
  const mfem::real_t current; // scaled total current through the tabs
  mfem::ParMixedBilinearForm * C = nullptr;
  mfem::HypreParMatrix Cmat;

public:
  ChargeBalance(mfem::ParFiniteElementSpace & h1,
                mfem::ParFiniteElementSpace & l2,
                const Region & region,
                mfem::real_t current)
    : Equation(h1), l2_fespace(l2), particle_region(region), current(current)
  {

    if (particle_region == NE)
    {
      mfem::Array<int> ess_bdr(fespace.GetParMesh()->bdr_attributes.Max());
      ess_bdr = 0;
      ess_bdr[3 - 1] = 1;
      fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }
  }

  using Equation::Update; // keep the base-class overloads visible next to Update() below
  virtual void Update();
  const mfem::HypreParMatrix & GetC() const { return Cmat; };

  ~ChargeBalance() override;
};
