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
  _Solver.SetMaxIter(500);
  _Solver.SetPrintLevel(1);
  _Solver.SetKDim(50);
  _Solver.SetPreconditioner(_I_prec);

  _block_trueOffsets.SetSize(settings::NMSMDEQS + 1);
  _block_trueOffsets[0] = 0;
  _block_trueOffsets[1] = h1_fespace_current_collector_positive.GetTrueVSize();
  _block_trueOffsets[2] = h1_fespace_current_collector_negative.GetTrueVSize();
  _block_trueOffsets[3] = l2_fespace_currents.GetTrueVSize();
  _block_trueOffsets[4] = 1;

  std::cout << "size of _block_trueOffsets = " << _block_trueOffsets[0] << ", " << _block_trueOffsets[1] << ", " << _block_trueOffsets[2] << ", "
            << _block_trueOffsets[3] << ", " << _block_trueOffsets[4] << std::endl;

  _block_trueOffsets.PartialSum();

  _x.Update(_block_trueOffsets);
  _f.Update(_block_trueOffsets);
  // How to distinguish between the positive and negative current collector?
  _cbp = new ChargeBalance(h1_fespace_current_collector_positive, l2_fespace_currents, PE);
  _cbn = new ChargeBalance(h1_fespace_current_collector_negative, l2_fespace_currents, NE);
  _cc = new CurrentConstraint(l2_fespace_currents);

  _Bcc = nullptr;

   std::cout << "CurrentCollectorOperator constructed" << std::endl;
  _cbp->Update();
  _cbn->Update();
   mfem::ConstantCoefficient voltage(4.2/phi_scale); 
  _cc->Update(voltage);
  Solve(); 
  std::cout << "Solve completed" << std::endl;
}
// replace Mult with Solve, does not take any parameters!
void
CurrentCollectorOperator::Solve()
{   
  _Bcc(0, 0) = &_cbp->GetK();
  delete _Bcc(0, 2);
  _Bcc(0, 2) = _cbp->GetC().Transpose();

  _f.GetBlock(0) = _cbp->GetZ();

  _Bcc(1, 1) = &_cbn->GetK();
  delete _Bcc(1, 2);
  _Bcc(1, 2) = _cbn->GetC().Transpose();

  _f.GetBlock(1) = _cbn->GetZ();

  _Bcc(2, 0) = &_cbp->GetC();
  _Bcc(2, 1) = &_cbn->GetC();
  _Bcc(2, 3) = &_cc->GetProduct();

  _f.GetBlock(2) = _cc->GetZ();

   delete _Bcc(3, 2);
  _Bcc(3, 2) = _cc->GetProduct().Transpose();

  _f.GetBlock(3) = I;

  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
        std::cout << "_Bcc(" << i << "," << j << ") = "
                  << static_cast<const void *>(_Bcc(i, j))
                  << '\n';
    }
  }

  for (int i = 0; i < _Bcc.NumRows(); i++)
{
    for (int j = 0; j < _Bcc.NumCols(); j++)
    {
        const HypreParMatrix * block = _Bcc(i, j);

        if (block != nullptr)
        {
            std::cout << "_Bcc(" << i << "," << j << ") : "
                      << block->NumRows() << " x "
                      << block->NumCols() << '\n';
        }
        else
        {
            std::cout << "_Bcc(" << i << "," << j << ") : null\n";
        }
    }
}



  // ============================================================
  // DEBUG: Check block dimensions and matrix norms
  // ============================================================

  std::cout << "\n===== Current Collector Matrix Diagnostics =====\n";

  std::cout << "Kp = "
            << _cbp->GetK().GetGlobalNumRows()
            << " x "
            << _cbp->GetK().GetGlobalNumCols()
            << '\n';

  std::cout << "Kn = "
            << _cbn->GetK().GetGlobalNumRows()
            << " x "
            << _cbn->GetK().GetGlobalNumCols()
            << '\n';

  std::cout << "Cp = "
            << _cbp->GetC().GetGlobalNumRows()
            << " x "
            << _cbp->GetC().GetGlobalNumCols()
            << '\n';

  std::cout << "Cn = "
            << _cbn->GetC().GetGlobalNumRows()
            << " x "
            << _cbn->GetC().GetGlobalNumCols()
            << '\n';

  std::cout << "P = "
            << _cc->GetProduct().GetGlobalNumRows()
            << " x "
            << _cc->GetProduct().GetGlobalNumCols()
            << '\n';

  std::cout << "||Kp|| = " << _cbp->GetK().FNorm() << '\n';
  std::cout << "||Kn|| = " << _cbn->GetK().FNorm() << '\n';
  std::cout << "||Cp|| = " << _cbp->GetC().FNorm() << '\n';
  std::cout << "||Cn|| = " << _cbn->GetC().FNorm() << '\n';
  std::cout << "||P||  = " << _cc->GetProduct().FNorm() << '\n';

  std::cout << "===============================================\n";


  // ============================================================
  // Construct the global block matrix
  // ============================================================

  delete _Acc;
  _Acc = mfem::HypreParMatrixFromBlocks(_Bcc);

// Check symmetry: ||A - A^T|| / ||A||
HypreParMatrix * AT = _Acc->Transpose();

HypreParMatrix * A_minus_AT = Add(1.0, *_Acc, -1.0, *AT);

const real_t norm_A = _Acc->FNorm();
const real_t norm_diff = A_minus_AT->FNorm();

const real_t symmetry_error =
    norm_diff / std::max(norm_A, std::numeric_limits<real_t>::epsilon());


    std::cout << "||A||_F = " << norm_A << '\n';
    std::cout << "||A - A^T||_F = " << norm_diff << '\n';
    std::cout << "Relative symmetry error = " << symmetry_error << '\n';


delete A_minus_AT;
delete AT;

  // use 4.2 for voltage cofficient (j)
//_Solver.SetOperator(*_Acc);
//Solver.Mult(_f, _x);

std::cout << "\n===== RHS vector f =====\n";

for (int i = 0; i < 4; i++)
{
    std::cout << "\n===== RHS block " << i << " =====\n";
    _f.GetBlock(i).Print(std::cout);
}

MUMPSSolver direct_solver(_Acc->GetComm());
direct_solver.SetOperator(*_Acc);
direct_solver.Mult(_f, _x);

std::cout << "\n===== Solution =====\n";

for (int i = 0; i < 4; i++)
{
    std::cout << "\n===== Solution block " << i << " =====\n";
    _x.GetBlock(i).Print(std::cout);
}

}

CurrentCollectorOperator::~CurrentCollectorOperator()
{
  for (int i = 0; i < _Bcc.NumRows(); i++)
    for (int j = 0; j < _Bcc.NumCols(); j++)
      delete _Bcc(i, j);

    delete _Acc;

    delete _cbp;
    delete _cbn;
    delete _cc;
}