
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


  // ------------------------------------------------------------
  // Construct ColMat = col as an n2 x 1 HypreParMatrix
  // ------------------------------------------------------------

  // 2. Local SparseMatrix setup: Local rows x 1 column
  // Every processor views this locally as a simple 1-column matrix
  mfem::SparseMatrix local_mat(fespace2->GetTrueVSize(), 1);
  local_mat = 1;
  local_mat.Finalize();

  // 3. Define global column partitioning (Must be identical on ALL ranks)
  // Rank 0 owns global column 0. All other ranks own 0 columns.
  HYPRE_BigInt col_starts[2];
  col_starts[0] = 0;
  col_starts[1] = 1;

  // 4. Construct the HypreParMatrix cleanly
  // This constructor auto-generates the correct diag, offd, and cmap layout for you.
  HypreParMatrix * ColMat = new mfem::HypreParMatrix(
      fespace2->GetComm(),
      fespace2->GlobalTrueVSize(), // Global row count
      1,                           // Global col count
      fespace2->GetTrueDofOffsets(),// Global row partitioning array
      col_starts,                  // Global col partitioning array
      &local_mat
  );

  Product = mfem::ParMult(Mpmat, ColMat);

  if (!Q)
  {
    Q = new ParLinearForm(&fespace);
    mfem::ConstantCoefficient one(1.0);
    Q->AddDomainIntegrator(new DomainLFIntegrator(one));
  }

  Q->Assemble();
  Q->ParallelAssemble(b);
  
}
