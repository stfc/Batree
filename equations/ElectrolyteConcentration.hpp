#pragma once

#include "mfem.hpp"
#include "equations/Equation.hpp"

using namespace mfem;

class ElectrolyteConcentration : public Equation
{
protected:
  Coefficient & j;
  GridFunctionCoefficient & ec_gfc;

public:
  ElectrolyteConcentration(ParFiniteElementSpace & f, Coefficient & j, GridFunctionCoefficient & ec_gfc) : Equation(f), j(j), ec_gfc(ec_gfc) {};
  virtual void Update();
};
