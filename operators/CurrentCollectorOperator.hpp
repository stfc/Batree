#pragma once

#include "equations/ChargeBalance.hpp"
#include "equations/CurrentConstraint.hpp"
#include "operators/LocalModel.hpp"

// remove inheritence from Operator class
//  possibly rename bcc by acc
class CurrentCollectorOperator
{
protected:
  mfem::ParFiniteElementSpace & h1_fespace_current_collector_positive;
  mfem::ParFiniteElementSpace & h1_fespace_current_collector_negative;
  mfem::ParFiniteElementSpace & l2_fespace_currents;

  ChargeBalance * _cbp = nullptr; // charge balance for positive current collector
  ChargeBalance * _cbn = nullptr; // charge balance for negative current collector
  CurrentConstraint * _cc = nullptr;

  mfem::real_t current_dt;

  mfem::CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
  mfem::HypreSmoother M_prec; // Preconditioner for the mass matrix M

  // CGSolver _Solver;      // Implicit solver for T = M + dt K
  mfem::GMRESSolver _Solver;   // GMRES solver for the block system
  mfem::HypreSmoother _I_prec; // Preconditioner for the implicit solver

  mutable mfem::Vector b; // auxiliary vector

  /// For solution true vector (4 MSMD equations)
  mfem::Array<int> _block_trueOffsets;

  /// System matrices for concentration and potential eqs
  mfem::HypreParMatrix * _Acc = nullptr;

  /// Auxiliary rhs vectors for MSMD equations
  mutable mfem::BlockVector _f;

  // mutable Array2D<const HypreParMatrix *> _Bcc{4, 4};
  mutable mfem::Array2D<const mfem::HypreParMatrix *> _Bcc{3, 3}; // reference potential s removed

  /// (2,2) block D = M diag(R̄) of the coupled system. Owned here; _Bcc(2,2) points to it.
  mfem::HypreParMatrix * _D = nullptr;

  /// Element currents ī fed to the local cell models (true dofs of the L2 space).
  /// Kept between calls, so a later time step can start from the previous answer.
  mfem::Vector _ibar;

  /// Print the block, RHS and solution diagnostics in Solve()
  bool _verbose = true;

  /// Reference to solution true dof vector
  mfem::BlockVector & _x;

public:
  CurrentCollectorOperator(mfem::ParFiniteElementSpace *& f_pos,
                           mfem::ParFiniteElementSpace *& f_neg,
                           mfem::ParFiniteElementSpace *& f_l2,
                           const unsigned & ndofs,
                           mfem::BlockVector & x);

  /// One linear solve of the 3x3 system. _Bcc(2,2) and block 2 of the RHS must be set
  /// first (SolveCoupled does this).
  virtual void Solve();

  /// Couple the collectors to a local cell model V̄(ī), at one instant in time.
  /// Each pass solves the 3x3 system with D = M diag(R̄) in the (2,2) block and
  /// M (V̄ + R̄ ī) in block-row 2 of the RHS, then sets ī ← (1 − ω) ī + ω î.
  ///   slope = "exact"  : R̄ recomputed every pass (Newton)
  ///           "frozen" : R̄ from the first pass only (chord)
  ///           "none"   : R̄ = 0, the report's Option 1 (lagged Picard)
  void SolveCoupled(const LocalModel & model,
                    const std::string & slope,
                    mfem::real_t omega = 1.0,
                    mfem::real_t tol = 1e-10,
                    int max_it = 40);

  virtual ~CurrentCollectorOperator();
};
