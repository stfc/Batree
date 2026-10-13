#pragma once

#include "equations/Equation.hpp"

class ChargeBalance : public Equation
{
private:
  const Region particle_region;
public:
    ChargeBalance(const Region & region,
                  ParFiniteElementSpace & f,
                  ParFiniteElementSpace * f2 = nullptr)
        : Equation(f, f2), particle_region(region) {}

  virtual void Update();
};
