#pragma once

#include "mfem.hpp"
#include "equations/Equation.hpp"

class ElectrolytePotential : public Equation
{
protected:
  mfem::ParLinearForm * Qc = nullptr;
  mfem::Vector bc;

public:
  ElectrolytePotential(mfem::ParFiniteElementSpace & f) : Equation(f), bc(f.GetTrueVSize())
  {
    f.GetEssentialTrueDofs(mfem::Array<int>({1, 0}), ess_tdof_list);
  }
  virtual void Update(const mfem::GridFunctionCoefficient & ec_gfc,
                      const mfem::Coefficient & j) override;
  void Reset()
  {
    delete K;
    delete Qc;
    K = nullptr;
    Qc = nullptr;
  }
};
