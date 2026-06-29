#pragma once

#include "mfem.hpp"
#include "equations/Equation.hpp"

using namespace mfem;

class ElectrolytePotential : public Equation
{
protected:
  Coefficient & j;
  GridFunctionCoefficient & ec_gfc;
  ParLinearForm * Qc = nullptr;
  Vector bc;

public:
  ElectrolytePotential(ParFiniteElementSpace & f, Coefficient & j, GridFunctionCoefficient & ec_gfc) : Equation(f), j(j), ec_gfc(ec_gfc), bc(f.GetTrueVSize())
  {
    f.GetEssentialTrueDofs(Array<int>({1, 0}), ess_tdof_list);
  }
  virtual void Update();
  void Reset()
  {
    delete K;
    delete Qc;
    K = nullptr;
    Qc = nullptr;
  }
};
