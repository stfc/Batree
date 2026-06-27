#pragma once

#include "mfem.hpp"
#include "equations/Equation.hpp"

using namespace mfem;

class SolidPotential : public Equation
{
public:
  SolidPotential(ParFiniteElementSpace & f) : Equation(f)
  {
    f.GetEssentialTrueDofs(Array<int>({1, 1}), ess_tdof_list);
  }
  virtual void Update(const Coefficient & j) override;
};
