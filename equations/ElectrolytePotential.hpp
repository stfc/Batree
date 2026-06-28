#pragma once

#include "mfem.hpp"
#include "equations/Equation.hpp"

using namespace mfem;

class ElectrolytePotential : public Equation
{
protected:
  ParLinearForm * Qc = nullptr;
  Vector bc;

public:
  ElectrolytePotential(ParFiniteElementSpace & f) : Equation(f), bc(f.GetTrueVSize())
  {
    f.GetEssentialTrueDofs(Array<int>({1, 0}), ess_tdof_list);
  }
  virtual void Update(const GridFunctionCoefficient & ec_gfc, const Coefficient & j) override;
  void Reset()
  {
    delete K;
    delete Qc;
    K = nullptr;
    Qc = nullptr;
  }
};
