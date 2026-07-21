#pragma once

#include "mfem.hpp"
#include "equations/Equation.hpp"

using namespace mfem;

class SolidPotential : public Equation
{
protected:
  Coefficient &j;

public:
  SolidPotential(ParFiniteElementSpace & f, Coefficient & j) : Equation(f), j(j)
  {
    f.GetEssentialTrueDofs(Array<int>({1, 1}), ess_tdof_list);
  }
  virtual void Update();
};
