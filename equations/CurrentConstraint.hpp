#pragma once

#include "equations/Equation.hpp"

class CurrentConstraint : public Equation
{
protected:
  mfem::HypreParMatrix * Product = nullptr;
  mfem::HypreParMatrix * ColMat = nullptr;

public:
  using Equation::Equation;

  virtual void Update(const mfem::Coefficient & j) override;
  const mfem::HypreParMatrix & GetProduct() const { return *Product; };

  virtual ~CurrentConstraint()
  {
    delete Product;
    delete ColMat;
  }
};
