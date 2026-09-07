#pragma once

#include "equations/Equation.hpp"

class ChargeBalance : public Equation
{
public:
  using Equation::Equation;

  virtual void Update();
};
