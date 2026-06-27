#pragma once

#include "mfem.hpp"
#include "equations/Equation.hpp"

using namespace mfem;

class ElectrolyteConcentration : public Equation
{
public:
  using Equation::Equation;
  virtual void Update(const GridFunctionCoefficient & ec_gfc, const Coefficient & j) override;
};
