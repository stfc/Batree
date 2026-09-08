#include "CurrentCollectorOperator.hpp"

CurrentCollectorOperator::CurrentCollectorOperator(ParFiniteElementSpace *& f_pos, ParFiniteElementSpace *& f_neg, ParFiniteElementSpace *& f_l2, const unsigned & ndofs, BlockVector & x)
   : Operator(ndofs), h1_fespace_current_collector_positive(*f_pos), h1_fespace_current_collector_negative(*f_neg), l2_fespace_currents(*f_l2), _x(x),
     T(NULL), T_solver(f_neg->GetComm())
{
   const real_t rel_tol = 1e-8;

   T_solver.iterative_mode = false;
   T_solver.SetRelTol(rel_tol);
   T_solver.SetAbsTol(0.0);
   T_solver.SetMaxIter(100);
   T_solver.SetPrintLevel(0);
   T_solver.SetPreconditioner(T_prec);

  _block_trueOffsets.SetSize(settings::NMSMDEQS + 1);
  _block_trueOffsets[0] = 0;
  _block_trueOffsets[1] = h1_fespace_current_collector_positive.GetTrueVSize();
  _block_trueOffsets[2] = h1_fespace_current_collector_negative.GetTrueVSize();
  _block_trueOffsets[3] = l2_fespace_currents.GetTrueVSize();
  _block_trueOffsets[4] = 1;

  _block_trueOffsets.PartialSum();

  _x.Update(_block_trueOffsets);
  _f.Update(_block_trueOffsets);
  //How to distinguish between the positive and negative current collector?
  _cb = new ChargeBalance(h1_fespace_current_collector_positive, &l2_fespace_currents);
  _cc = new CurrentConstraint(h1_fespace_current_collector_positive, &l2_fespace_currents);
}

void CurrentCollectorOperator::Mult(const Vector &x, Vector &y) const
{
   //Is the below correct place to populate the _Bcc matrix? Should it be done in the constructor?
   delete _Bcc(0, 0);
  _Bcc(0, 0) = &_cb->GetK();
   delete _Bcc(0, 2);
  _Bcc(0, 2) = &_cb->GetC();

  _f.GetBlock(0) = _cb->GetZ();
}

CurrentCollectorOperator::~CurrentCollectorOperator()
{
   delete T;
}