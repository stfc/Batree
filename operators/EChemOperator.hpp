#pragma once

#include "equations/ElectrolytePotential.hpp"
#include "equations/ElectrolyteConcentration.hpp"
#include "equations/SolidPotential.hpp"
#include "equations/SolidConcentration.hpp"
#include "coefficients/ExchangeCurrentCoefficient.hpp"
#include "coefficients/OpenCircuitPotentialCoefficient.hpp"
#include "coefficients/OverPotentialCoefficient.hpp"
#include "coefficients/ReactionCurrentCoefficient.hpp"
#include "coefficients/ElectrodeReactionCurrentCoefficient.hpp"

class EChemOperator : public mfem::TimeDependentOperator
{
protected:
  mfem::ParFiniteElementSpace & _x_h1space;
  mfem::ParFiniteElementSpace & _r_h1space;

  ElectrolytePotential * _ep = nullptr;
  SolidPotential * _sp = nullptr;
  ElectrolyteConcentration * _ec = nullptr;
  mfem::Array<SolidConcentration *> _sc;

  /// Gridfunctions defined over _x_h1space (3 macro eqs plus _surface_ concentration)
  mfem::ParGridFunction _ep_gf;
  mfem::ParGridFunction _sp_gf;
  mfem::ParGridFunction _ec_gf;
  mfem::ParGridFunction _sc_gf;

  /// Coefficients for the gridfunction defined over _x_h1space (3 macro eqs plus _surface_ concentration)
  mfem::GridFunctionCoefficient _ep_gfc;
  mfem::GridFunctionCoefficient _sp_gfc;
  mfem::GridFunctionCoefficient _ec_gfc;
  mfem::GridFunctionCoefficient _sc_gfc;

  /// Big-enough array for the surface concentrations of the two SPM(e) particles
  mfem::Array<mfem::real_t> _sc_array{/* E */ 0., /* NE */ 0., /* SEP */ 0., /* PE */ 0.};

  /// Big-enough array for the reference potentials
  mfem::Array<mfem::real_t> _rp_array{/* E */ 0., /* NE */ 0., /* SEP */ 0., /* PE */ 0.};

  /// Coefficients for derived, i.e. not solved for, quantities
  ReactionCurrentCoefficient * _j;
  ExchangeCurrentCoefficient * _jex;
  OpenCircuitPotentialCoefficient * _ocp;
  OverPotentialCoefficient * _op;

  /// For solution true vector (3 macros eqs plus NPAR _radial_ concentrations)
  mfem::Array<int> _block_trueOffsets;
  /// For rhs true vectors (2 macro eqs)
  mfem::Array<int> _potential_trueOffsets;
  /// For rhs true vectors (1 macro eq plus NPAR _radial_ concentrations)
  mfem::Array<int> _concentration_trueOffsets;

  /// System matrices for concentration and potential eqs
  mfem::HypreParMatrix *_Ac = nullptr, *_Ap = nullptr;

  /// 2D array of pointers for each block in the system matrices
  mfem::Array2D<const mfem::HypreParMatrix *> _Bc{int(NPAR) + 1, int(NPAR) + 1}, _Bp{2, 2};

  /// Reference to solution true dof vector
  mfem::BlockVector & _x;

  /// Block vectors wrapping the concentration and potential solution true dof vectors
  mfem::BlockVector _xc, _xp;

  /// Implicit solver for T = M + dt K
  mfem::CGSolver _Solver;
  /// Preconditioner for the implicit solver
  mfem::HypreSmoother _Prec;

  /// Auxiliary rhs vectors for concentrations and potential eqs
  mfem::BlockVector _bc, _bp;

  /// Self-consistency loop "L2" error threshold (JuBat uses 1e-9)
  const mfem::real_t _scl_threshold = 1e-7;

  /// Self-consistency loop 4-point integration rule
  mfem::IntegrationRule _scl_ir = mfem::IntegrationRules().Get(mfem::Geometry::Type::SEGMENT, 7);

  /// Self-consistency loop quadrature space
  mfem::QuadratureSpace _scl_qspace = mfem::QuadratureSpace(*_x_h1space.GetParMesh(), _scl_ir);

  /// Self-consistency loop quadrature function for the reaction current
  mfem::QuadratureFunction _j_qfunction = mfem::QuadratureFunction(_scl_qspace);

  /// Self-consistency loop quadrature function vector for the reaction current
  mfem::Vector _j_vec = mfem::Vector(_j_qfunction.Size());

public:
  EChemOperator(mfem::ParFiniteElementSpace & x_h1space,
                mfem::ParFiniteElementSpace & r_h1space,
                const unsigned & ndofs,
                mfem::BlockVector & x);

  virtual void Mult(const mfem::Vector & x, mfem::Vector & dx_dt) const override {};

  /** Solve the Backward-Euler equation: k = f(x + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void
  ImplicitSolve(const mfem::real_t dt, const mfem::Vector & x, mfem::Vector & k) override;

  void SetPotentialGridFunctionsFromTrueVectors();
  void SetConcentrationGridFunctionsFromTrueVectors();
  void SetSurfaceConcentration();
  void SetReferencePotential();

  void UpdatePotentialEquations();
  void UpdateConcentrationEquations();

  mfem::real_t GetElectrodeReactionCurrent(const Region & r, const int & sign);
  mfem::Array<mfem::real_t> GetParticleReactionCurrent();

  /// Construct coefficients for derived quantities
  void ConstructReactionCurrent();
  void ConstructExchangeCurrent();
  void ConstructOpenCircuitPotential();
  void ConstructOverPotential();

  /// Helpers for quantities which are constant within a region
  const mfem::real_t & GetSurfaceConcentration(const Region & r);
  const mfem::real_t & GetReferencePotential(const Region & r);
  const mfem::real_t & GetReactionCurrent(const Region & r);
  const mfem::real_t & GetExchangeCurrent(const Region & r);
  const mfem::real_t & GetOpenCircuitPotential(const Region & r);
  const mfem::real_t & GetOverPotential(const Region & r);

  mfem::real_t GetVoltage();
  mfem::real_t GetVoltageMarquisCorrection();
  mfem::real_t GetSoC();

  virtual void GetParticleDofs(mfem::Array<int> & particle_dofs,
                               mfem::Array<Region> & particle_regions,
                               mfem::Array<int> & particle_offsets);

  virtual ~EChemOperator()
  {
    delete _ep;
    delete _sp;
    delete _ec;
    for (unsigned p = 0; p < NPAR; p++)
      delete _sc[p];

    delete _j;
    delete _jex;
    delete _ocp;
    delete _op;

    delete _Ap;
    delete _Ac;

    for (int i = 0; i < _Bc.NumRows(); i++)
      for (int j = 0; j < _Bc.NumCols(); j++)
        delete _Bc(i, j);
  }
};
