#pragma once

#include "equations/Equation.hpp"

class ChargeBalance : public Equation
{
protected:
  ParFiniteElementSpace & l2_fespace;
  const Region particle_region;
  ParMixedBilinearForm * C = nullptr;
  HypreParMatrix * Cmat = nullptr;
public:
    ChargeBalance(ParFiniteElementSpace & h1,
                  ParFiniteElementSpace & l2,
                  const Region & region)
        : Equation(h1), l2_fespace(l2), particle_region(region) {}

  virtual void Update();
  const HypreParMatrix & GetC() const { return *Cmat; };

  ~ChargeBalance() override;
};
