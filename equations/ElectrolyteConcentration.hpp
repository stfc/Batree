#pragma once

#include "equations/Equation.hpp"

class ElectrolyteConcentration : public Equation
{
public:
  using Equation::Equation;
  virtual void Update(const mfem::GridFunctionCoefficient & ec_gfc,
                      const mfem::Coefficient & j) override;
};
