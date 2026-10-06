#pragma once

#include "mfem.hpp"
#include "parameters/settings.hpp"

using namespace mfem;
using namespace settings;

class Equation
{
protected:
  ParFiniteElementSpace & fespace;
  ParFiniteElementSpace * fespace2 = nullptr;

  Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

  ParBilinearForm * M = nullptr;
  ParBilinearForm * K = nullptr;
  ParLinearForm * Q = nullptr;
  ParMixedBilinearForm * C = nullptr;

  HypreParMatrix Mmat;
  HypreParMatrix Kmat;
  HypreParMatrix * Mpmat = nullptr;
  HypreParMatrix * Kpmat = nullptr;
  HypreParMatrix * Cmat = nullptr;
  HypreParVector Col;
  HypreParMatrix * ColMat = nullptr;

  Vector b; // auxiliary vector


public:
  Equation(ParFiniteElementSpace & f,
           ParFiniteElementSpace * f2 = nullptr)
    : fespace(f),
      fespace2(f2),
      b(f.GetTrueVSize())
  {}

  const HypreParMatrix & GetM() const { return Mmat; };
  const HypreParMatrix & GetK() const { return Kmat; };
  const HypreParMatrix & GetMp() const { return *Mpmat; };
  const HypreParMatrix & GetKp() const { return *Kpmat; };
  const HypreParMatrix & GetC() const { return *Cmat; };
  const Vector & GetZ() const { return b; };
  const HypreParVector & GetCol() const { return Col; };

  virtual void Update(const Coefficient & j) {}
  virtual void Update(const GridFunctionCoefficient & u,
                      const Coefficient & j) {}

  virtual ~Equation()
  {
    delete M;
    delete K;
    delete Q;
    delete C;
    delete Mpmat;
    delete Kpmat;
    delete Cmat;
  }
};