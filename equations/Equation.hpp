#pragma once

#include "parameters/settings.hpp"

using namespace settings;

class Equation
{
protected:
  mfem::ParFiniteElementSpace & fespace;

  mfem::Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

  mfem::ParBilinearForm * M = nullptr;
  mfem::ParBilinearForm * K = nullptr;
  mfem::ParLinearForm * Q = nullptr;

  mfem::HypreParMatrix Mmat;
  mfem::HypreParMatrix Kmat;

  mfem::Vector b; // auxiliary vector

public:
  Equation(mfem::ParFiniteElementSpace & f) : fespace(f), b(f.GetTrueVSize()) {};

  const mfem::HypreParMatrix & GetM() const { return Mmat; };
  const mfem::HypreParMatrix & GetK() const { return Kmat; };
  const mfem::Vector & GetZ() const { return b; };

  virtual void Update(const mfem::Coefficient & j) {}
  virtual void Update(const mfem::GridFunctionCoefficient & u, const mfem::Coefficient & j) {}

  virtual ~Equation()
  {
    delete M;
    delete K;
    delete Q;
  }
};
