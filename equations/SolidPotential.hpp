#pragma once

#include "mfem.hpp"
#include "equations/Equation.hpp"

class SolidPotential : public Equation
{
public:
  SolidPotential(mfem::ParFiniteElementSpace & f) : Equation(f)
  {
    f.GetEssentialTrueDofs(mfem::Array<int>({1, 1}), ess_tdof_list);
  }
  virtual void Update(const mfem::Coefficient & j) override;
};
