#pragma once

#include "equations/Equation.hpp"

class CurrentConstraint : public Equation
{
public:
  using Equation::Equation;

  virtual void Update();
};
