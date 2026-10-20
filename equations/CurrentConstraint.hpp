#pragma once

#include "equations/Equation.hpp"

class CurrentConstraint : public Equation
{
protected:
  HypreParMatrix * Product = nullptr;
  HypreParMatrix * ColMat = nullptr;
public:
  using Equation::Equation;

  virtual void Update(const Coefficient & j) override;
  const HypreParMatrix & GetProduct() const { return *Product; };

  virtual ~CurrentConstraint()
  {
    delete Product;
    delete ColMat;
  }
};
