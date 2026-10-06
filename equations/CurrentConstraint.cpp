
#include "equations/CurrentConstraint.hpp"

void
CurrentConstraint::Update()
{

  if (!M)
  {
    M = new ParBilinearForm(fespace2);
    M->AddDomainIntegrator(new MassIntegrator);
    M->Assemble(0); // keep sparsity pattern of M and K the same
    Mpmat = M->ParallelAssemble();
  }

  HypreParVector one(*Mpmat);
  one = 1.0;

  Col = HypreParVector(fespace2->GetComm(), fespace2->GlobalTrueVSize(), fespace2->GetTrueDofOffsets());
  Mpmat->Mult(one, Col);

  // ------------------------------------------------------------
  // Construct ColMat = col as an n2 x 1 HypreParMatrix
  // ------------------------------------------------------------

  MPI_Comm comm = fespace2->GetComm();

  // 1. Get the uniform local row count for this processor
  int local_nrows = Col.Size(); 
  const double* col_data = Col.GetData(); // Raw vector values on this rank

  // 2. Local SparseMatrix setup: Local rows x 1 column
  // Every processor views this locally as a simple 1-column matrix
  mfem::SparseMatrix local_mat(local_nrows, 1);
  for (int i = 0; i < local_nrows; i++)
  {
      local_mat.Set(i, 0, col_data[i]);
  }
  local_mat.Finalize();

  // 3. Define global column partitioning (Must be identical on ALL ranks)
  // Rank 0 owns global column 0. All other ranks own 0 columns.
  HYPRE_BigInt col_starts[2];
  col_starts[0] = 0;
  col_starts[1] = 1;

  // 4. Construct the HypreParMatrix cleanly
  // This constructor auto-generates the correct diag, offd, and cmap layout for you.
  ColMat = new mfem::HypreParMatrix(
      comm,
      fespace2->GlobalTrueVSize(), // Global row count
      1,                           // Global col count
      fespace2->GetTrueDofOffsets(),// Global row partitioning array
      col_starts,                  // Global col partitioning array
      &local_mat
  );


  if (!Q)
  {
    Q = new ParLinearForm(&fespace);
    mfem::ConstantCoefficient one(1.0);
    Q->AddDomainIntegrator(new DomainLFIntegrator(one));
  }

  Q->Assemble();
  Q->ParallelAssemble(b);
  
}
