
#include "equations/CurrentConstraint.hpp"

void
CurrentConstraint::Update(const mfem::Coefficient & j)
{

  if (!M)
  {
    M = new mfem::ParBilinearForm(&fespace);
    M->AddDomainIntegrator(new mfem::MassIntegrator);
    M->Assemble(0); // keep sparsity pattern of M and K the same
    M->FormSystemMatrix(ess_tdof_list, Mmat);
  }

  // ------------------------------------------------------------
  // Construct ColMat = col as an n2 x 1 HypreParMatrix
  // ------------------------------------------------------------

  // 2. Local rows in CSR form: every local row has a single 1 in global column 0
  const int n = fespace.GetTrueVSize();
  mfem::Array<int> row_ptr(n + 1);
  mfem::Array<HYPRE_BigInt> col_idx(n);
  mfem::Vector data(n);
  for (int i = 0; i <= n; i++)
    row_ptr[i] = i;
  col_idx = 0;
  data = 1.0;

  // 3. Column partitioning: rank 0 owns global column 0, all other ranks own no columns
  const HYPRE_BigInt col_starts[2] = {fespace.GetMyRank() == 0 ? 0 : 1, 1};

  // 4. Construct the HypreParMatrix. This constructor splits each rank's rows into the columns
  // it owns and the ones other ranks own, so the column can live on rank 0 only.
  ColMat = new mfem::HypreParMatrix(fespace.GetComm(),
                                    n,                           // Local row count
                                    fespace.GlobalTrueVSize(),   // Global row count
                                    1,                           // Global col count
                                    row_ptr.GetData(),           // CSR row pointers
                                    col_idx.GetData(),           // CSR global column indices
                                    data.GetData(),              // CSR values
                                    fespace.GetTrueDofOffsets(), // Global row partitioning array
                                    col_starts);                 // Global col partitioning array

  Product = mfem::ParMult(&Mmat, ColMat);

  Product->Print("Product");

  if (!Q)
  {
    Q = new mfem::ParLinearForm(&fespace);
    Q->AddDomainIntegrator(new mfem::DomainLFIntegrator(const_cast<mfem::Coefficient &>(j)));
  }

  Q->Assemble();
  Q->ParallelAssemble(b);
}
