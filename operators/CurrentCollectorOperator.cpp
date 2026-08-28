#include "CurrentCollectorOperator.hpp"

CurrentCollectorOperator::CurrentCollectorOperator(ParFiniteElementSpace *& f, const unsigned & ndofs, const unsigned & np2ds, BlockVector & x, const Array<int> &etl)
   : Operator(ndofs), fespace(*f), np2ds(np2ds), _x(x),
     ess_tdof_list(etl), T(NULL), T_solver(f->GetComm())
{
   const real_t rel_tol = 1e-8;

   T_solver.iterative_mode = false;
   T_solver.SetRelTol(rel_tol);
   T_solver.SetAbsTol(0.0);
   T_solver.SetMaxIter(100);
   T_solver.SetPrintLevel(0);
   T_solver.SetPreconditioner(T_prec);

  //_block_trueOffsets.SetSize(NMSMDEQS + 1);??
  _block_trueOffsets.SetSize(4 + 1);
  _block_trueOffsets[0] = 0;
  _block_trueOffsets[1] = fespace.GetTrueVSize();
  _block_trueOffsets[2] = fespace.GetTrueVSize();
  _block_trueOffsets[3] = np2ds;
  _block_trueOffsets[4] = 1;

  _block_trueOffsets.PartialSum();

  _x.Update(_block_trueOffsets);
}

void CurrentCollectorOperator::Mult(const Vector &x, Vector &y) const
{

}

CurrentCollectorOperator::~CurrentCollectorOperator()
{
   delete T;
}