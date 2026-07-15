#include "CurrentOperator.hpp"

CurrentCollectorOperator::CurrentOperator(ParFiniteElementSpace &f, const Vector &u, const Array<int> &etl)
   : Operator(f.GetTrueVSize()), fespace(f),
     ess_tdof_list(etl), T(NULL), T_solver(f.GetComm())
{
   const real_t rel_tol = 1e-8;

   T_solver.iterative_mode = false;
   T_solver.SetRelTol(rel_tol);
   T_solver.SetAbsTol(0.0);
   T_solver.SetMaxIter(100);
   T_solver.SetPrintLevel(0);
   T_solver.SetPreconditioner(T_prec);
}

void CurrentCollectorOperator::Mult(const Vector &x, Vector &y)
{

}

CurrentCollectorOperator::~CurrentOperator()
{
   delete T;
}