#include "mfem.hpp"

using namespace std;
using namespace mfem;

#pragma once

class CurrentCollectorOperator : public Operator
{
protected:
   ParFiniteElementSpace &h1_fespace_current_collector_positive;
   ParFiniteElementSpace &h1_fespace_current_collector_negative;
   ParFiniteElementSpace &l2_fespace_currents;

   HypreParMatrix *T; // T = M + dt K
   real_t current_dt;

   CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
   HypreSmoother M_prec; // Preconditioner for the mass matrix M

   CGSolver T_solver;    // Implicit solver for T = M + dt K
   HypreSmoother T_prec; // Preconditioner for the implicit solver

   mutable Vector b; // auxiliary vector

   /// For solution true vector (4 MSMD equations)
   Array<int> _block_trueOffsets;

   Array2D<const HypreParMatrix *> _Bc{4, 4};

   /// Reference to solution true dof vector
   BlockVector & _x;

public:
   CurrentCollectorOperator(ParFiniteElementSpace * &f_pos, ParFiniteElementSpace * &f_neg, ParFiniteElementSpace * &f_l2, const unsigned & ndofs, BlockVector & x);

   virtual void Mult(const Vector & x, Vector &y) const;

   virtual ~CurrentCollectorOperator();
};