#include "CurrentCollectorOperator.hpp"

CurrentCollectorOperator::CurrentCollectorOperator(ParFiniteElementSpace *& f_pos,
                                                   ParFiniteElementSpace *& f_neg,
                                                   ParFiniteElementSpace *& f_l2,
                                                   const unsigned & ndofs,
                                                   BlockVector & x)
  : h1_fespace_current_collector_positive(*f_pos),
    h1_fespace_current_collector_negative(*f_neg),
    l2_fespace_currents(*f_l2),
    _x(x),
    _Solver(f_neg->GetComm())
{
  const real_t rel_tol = 1e-8;

  _Solver.iterative_mode = false;
  _Solver.SetRelTol(rel_tol);
  _Solver.SetAbsTol(0.0);
  _Solver.SetMaxIter(100);
  _Solver.SetPrintLevel(0);
  _Solver.SetPreconditioner(_I_prec);

  _block_trueOffsets.SetSize(settings::NMSMDEQS + 1);
  _block_trueOffsets[0] = 0;
  _block_trueOffsets[1] = h1_fespace_current_collector_positive.GetTrueVSize();
  _block_trueOffsets[2] = h1_fespace_current_collector_negative.GetTrueVSize();
  _block_trueOffsets[3] = l2_fespace_currents.GetTrueVSize();
  _block_trueOffsets[4] = 1;

  _block_trueOffsets.PartialSum();

  _x.Update(_block_trueOffsets);
  _f.Update(_block_trueOffsets);
  // How to distinguish between the positive and negative current collector?
  _cbp = new ChargeBalance(h1_fespace_current_collector_positive, l2_fespace_currents, PE);
  _cbn = new ChargeBalance(h1_fespace_current_collector_negative, l2_fespace_currents, NE);
  _cc = new CurrentConstraint(l2_fespace_currents);
}
// replace Mult with Solve, does not take any parameters!
void
CurrentCollectorOperator::Solve()
{
  // Is the below correct place to populate the _Bcc matrix? Should it be done in the constructor?
  delete _Bcc(0, 0);
  _Bcc(0, 0) = &_cbp->GetK();
  delete _Bcc(0, 2);
  _Bcc(0, 2) = _cbp->GetC().Transpose();
  _f.GetBlock(0) = _cbp->GetZ();

  delete _Bcc(1, 1);
  _Bcc(1, 1) = &_cbn->GetK();
  delete _Bcc(1, 2);
  _Bcc(1, 2) = _cbn->GetC().Transpose();
  _f.GetBlock(1) = _cbn->GetZ();

  delete _Bcc(2, 0);
  _Bcc(2, 0) = &_cbp->GetC();
  delete _Bcc(2, 1);
  _Bcc(2, 1) = &_cbn->GetC();
  delete _Bcc(2, 3);
  _Bcc(2, 3) = &_cc->GetProduct();
  _f.GetBlock(2) = _cc->GetZ();

  delete _Bcc(3, 2);
  _Bcc(3, 2) = _cc->GetProduct().Transpose();
  _f.GetBlock(3) = I;

  delete _Acc;
  _Acc = mfem::HypreParMatrixFromBlocks(_Bcc);

  // use 4.2 for voltage cofficient (j)
  _Solver.SetOperator(*_Acc);
  _Solver.Mult(_f, _x);
}

CurrentCollectorOperator::~CurrentCollectorOperator()
{
  for (int i = 0; i < _Bcc.NumRows(); i++)
    for (int j = 0; j < _Bcc.NumCols(); j++)
      delete _Bcc(i, j);
}