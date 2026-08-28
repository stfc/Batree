#include "mfem.hpp"

using namespace std;
using namespace mfem;

#pragma once

class CurrentCollectorOperator : public Operator
{
protected:
   ParFiniteElementSpace &fespace;
   Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

   HypreParMatrix *T; // T = M + dt K
   real_t current_dt;

   CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
   HypreSmoother M_prec; // Preconditioner for the mass matrix M

   CGSolver T_solver;    // Implicit solver for T = M + dt K
   HypreSmoother T_prec; // Preconditioner for the implicit solver

   mutable Vector b; // auxiliary vector

   const unsigned & np2ds;

   /// For solution true vector (4 MSMD equations)
   Array<int> _block_trueOffsets;

   /// Reference to solution true dof vector
   BlockVector & _x;

public:
   CurrentCollectorOperator(ParFiniteElementSpace * &f, const unsigned & ndofs, const unsigned & np2ds, BlockVector & x, const Array<int> &etl);

   virtual void Mult(const Vector & x, Vector &y) const;

   virtual ~CurrentCollectorOperator();
};