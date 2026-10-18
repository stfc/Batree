#pragma once

#include "mfem.hpp"

#include "equations/ChargeBalance.hpp"
#include "equations/CurrentConstraint.hpp"

using namespace std;
using namespace mfem;
//remove inheritence from Operator class
// possibly rename bcc by acc
class CurrentCollectorOperator 
{
protected:
   ParFiniteElementSpace &h1_fespace_current_collector_positive;
   ParFiniteElementSpace &h1_fespace_current_collector_negative;
   ParFiniteElementSpace &l2_fespace_currents;

   ChargeBalance * _cbp = nullptr; //charge balance for positive current collector
   ChargeBalance * _cbn = nullptr; //charge balance for negative current collector
   CurrentConstraint * _cc = nullptr;

   real_t current_dt;

   CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
   HypreSmoother M_prec; // Preconditioner for the mass matrix M

   CGSolver _Solver;    // Implicit solver for T = M + dt K
   HypreSmoother _I_prec; // Preconditioner for the implicit solver

   mutable Vector b; // auxiliary vector

   /// For solution true vector (4 MSMD equations)
   Array<int> _block_trueOffsets;

   /// System matrices for concentration and potential eqs
   HypreParMatrix *_Acc = nullptr;

   /// Auxiliary rhs vectors for MSMD equations
   mutable BlockVector _f;

   mutable Array2D<const HypreParMatrix *> _Bcc{4, 4};

   /// Reference to solution true dof vector
   BlockVector & _x;

public:
   CurrentCollectorOperator(ParFiniteElementSpace * &f_pos, ParFiniteElementSpace * &f_neg, ParFiniteElementSpace * &f_l2, const unsigned & ndofs, BlockVector & x);

   virtual void Solve();

   virtual ~CurrentCollectorOperator();
};