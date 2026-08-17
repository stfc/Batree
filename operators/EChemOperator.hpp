#pragma once

#include "mfem.hpp"
#include "equations/ElectrolytePotential.hpp"
#include "equations/ElectrolyteConcentration.hpp"
#include "equations/SolidPotential.hpp"
#include "equations/SolidConcentration.hpp"
#include "coefficients/ExchangeCurrentCoefficient.hpp"
#include "coefficients/OpenCircuitPotentialCoefficient.hpp"
#include "coefficients/OverPotentialCoefficient.hpp"
#include "coefficients/ReactionCurrentCoefficient.hpp"
#include "coefficients/ElectrodeReactionCurrentCoefficient.hpp"

using namespace mfem;

class EChemOperator : public TimeDependentOperator
{
protected:
  ParFiniteElementSpace & _x_h1space;
  ParFiniteElementSpace & _r_h1space;
  ParFiniteElementSpace * _x_l2space = nullptr;

  ElectrolytePotential * _ep = nullptr;
  SolidPotential * _sp = nullptr;
  ElectrolyteConcentration * _ec = nullptr;
  Array<SolidConcentration *> _sc;

  /// Gridfunctions defined over _x_h1space (3 macro eqs plus _surface_ concentration)
  ParGridFunction _ep_gf;
  ParGridFunction _sp_gf;
  ParGridFunction _ec_gf;
  ParGridFunction _sc_gf;

  /// Coefficients for the gridfunction defined over _x_h1space (3 macro eqs plus _surface_ concentration)
  GridFunctionCoefficient _ep_gfc;
  GridFunctionCoefficient _sp_gfc;
  GridFunctionCoefficient _ec_gfc;
  GridFunctionCoefficient _sc_gfc;

  /// Big-enough array for the surface concentrations of the two SPM(e) particles
  Array<real_t> _sc_array{/* E */ 0., /* NE */ 0., /* SEP */ 0., /* PE */ 0.};

  /// Big-enough array for the reference potentials
  Array<real_t> _rp_array{/* E */ 0., /* NE */ 0., /* SEP */ 0., /* PE */ 0.};

  /// Coefficients for derived, i.e. not solved for, quantities
  ReactionCurrentCoefficient * _j;
  ExchangeCurrentCoefficient * _jex;
  OpenCircuitPotentialCoefficient * _ocp;
  OverPotentialCoefficient * _op;

  /// For solution true vector (3 macros eqs plus NPAR _radial_ concentrations)
  Array<int> _block_trueOffsets;
  /// For rhs true vectors (2 macro eqs)
  Array<int> _potential_trueOffsets;
  /// For rhs true vectors (1 macro eq plus NPAR _radial_ concentrations)
  Array<int> _concentration_trueOffsets;

  /// System matrices for concentration and potential eqs
  HypreParMatrix *_Ac = nullptr, *_Ap = nullptr;

  /// 2D array of pointers for each block in the system matrices
  Array2D<const HypreParMatrix *> _Bc{int(NPAR) + 1, int(NPAR) + 1}, _Bp{2, 2};

  /// Reference to solution true dof vector
  BlockVector & _x;

  /// Block vectors wrapping the concentration and potential solution true dof vectors
  BlockVector _xc, _xp;

  /// Implicit solver for T = M + dt K
  CGSolver _Solver;
  /// Preconditioner for the implicit solver
  HypreSmoother _Prec;

  /// Auxiliary rhs vectors for concentrations and potential eqs
  BlockVector _bc, _bp;

  /// Self-consistency loop "L2" error threshold (JuBat uses 1e-9)
  const real_t _scl_threshold = 1e-7;

  /// Self-consistency loop 4-point integration rule
  IntegrationRule _scl_ir = IntegrationRules().Get(Geometry::Type::SEGMENT, 7);
  const IntegrationRule * _scl_irs[Geometry::Type::NUM_GEOMETRIES];

public:
  EChemOperator(ParFiniteElementSpace & x_h1space,
                ParFiniteElementSpace & r_h1space,
                const unsigned & ndofs,
                BlockVector & x);

  virtual void Mult(const Vector & x, Vector & dx_dt) const override {};

  /** Solve the Backward-Euler equation: k = f(x + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void ImplicitSolve(const real_t dt, const Vector & x, Vector & k) override;

  void SetPotentialGridFunctionsFromTrueVectors();
  void SetConcentrationGridFunctionsFromTrueVectors();
  void SetSurfaceConcentration();
  void SetReferencePotential();

  void UpdatePotentialEquations();
  void UpdateConcentrationEquations();

  real_t GetElectrodeReactionCurrent(const Region & r, const int & sign);
  Array<real_t> GetParticleReactionCurrent();

  /// Construct coefficients for derived quantities
  void ConstructReactionCurrent();
  void ConstructExchangeCurrent();
  void ConstructOpenCircuitPotential();
  void ConstructOverPotential();

  /// Helpers for quantities which are constant within a region
  const real_t & GetSurfaceConcentration(const Region & r);
  const real_t & GetReferencePotential(const Region & r);
  const real_t & GetReactionCurrent(const Region & r);
  const real_t & GetExchangeCurrent(const Region & r);
  const real_t & GetOpenCircuitPotential(const Region & r);
  const real_t & GetOverPotential(const Region & r);

  real_t GetVoltage();
  real_t GetVoltageMarquisCorrection();
  real_t GetSoC();

  virtual void GetParticleDofs(Array<int> & particle_dofs,
                               Array<Region> & particle_regions,
                               Array<int> & particle_offsets);

  virtual ~EChemOperator()
  {
    if (_x_l2space)
      delete _x_l2space->FEColl();
    delete _x_l2space;

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
