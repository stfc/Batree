#include "operators/EChemOperator.hpp"

EChemOperator::EChemOperator(ParFiniteElementSpace *& x_h1space,
                             Array<ParFiniteElementSpace *> & r_h1space,
                             const unsigned & ndofs,
                             BlockVector & x)
  : TimeDependentOperator(ndofs, (real_t)0.0),
    _x_h1space(x_h1space),
    _r_h1space(r_h1space),
    _ep_gf(_x_h1space),
    _sp_gf(_x_h1space),
    _ec_gf(_x_h1space),
    _sc_gf(_x_h1space),
    _ep_gfc(&_ep_gf),
    _sp_gfc(&_sp_gf),
    _ec_gfc(&_ec_gf),
    _sc_gfc(&_sc_gf),
    _x(x),
    _Solver(_x_h1space->GetComm())
{
  _Solver.iterative_mode = true;
  _Solver.SetRelTol(1e-11);
  _Solver.SetMaxIter(500);
  _Solver.SetPrintLevel(0);
  _Solver.SetPreconditioner(_Prec);

  _block_trueOffsets.SetSize(NEQS + 1);
  _potential_trueOffsets.SetSize(NMACROP + 1);
  _concentration_trueOffsets.SetSize(NMACROC + NPAR + 1);

  _block_trueOffsets[0] = 0;
  _block_trueOffsets[EP + 1] = _x_h1space->GetTrueVSize();
  _block_trueOffsets[SP + 1] = _x_h1space->GetTrueVSize();
  _block_trueOffsets[EC + 1] = _x_h1space->GetTrueVSize();

  _potential_trueOffsets[0] = 0;
  _potential_trueOffsets[EPP + 1] = _x_h1space->GetTrueVSize();
  _potential_trueOffsets[SPP + 1] = _x_h1space->GetTrueVSize();

  _concentration_trueOffsets[0] = 0;
  _concentration_trueOffsets[ECC + 1] = _x_h1space->GetTrueVSize();

  for (unsigned p = 0; p < NPAR; p++)
  {
    _block_trueOffsets[SC + 1 + p] = _r_h1space[p]->GetTrueVSize();
    _concentration_trueOffsets[SCC + 1 + p] = _r_h1space[p]->GetTrueVSize();
  }

  _block_trueOffsets.PartialSum();
  _potential_trueOffsets.PartialSum();
  _concentration_trueOffsets.PartialSum();

  if (!Mpi::WorldRank())
  {
    std::cout << "Variables: " << NEQS << std::endl;
    std::cout << "Unknowns (rank 0): " << _block_trueOffsets[NEQS] << std::endl;
  }

  // Set offsets for solution (complete, potential and concentration) true vectors
  _x.Update(_block_trueOffsets);
  _xp.Update(_x, _potential_trueOffsets);
  _xc.~BlockVector();
  new (&_xc) BlockVector(_x, _block_trueOffsets[C], _concentration_trueOffsets);

  // Set offsets for rhs (potential and concentration) true vectors
  _bp.Update(_potential_trueOffsets);
  _bc.Update(_concentration_trueOffsets);

  // Initialise 2D arrays of pointers for each block in the system matrices
  _Bp = nullptr;
  _Bc = nullptr;

  // Set initial conditions (in particular for electrolyte concentration)
  _x = 0.;
  _x.GetBlock(EC) = CE0;
  _ep_gf = 0.;
  _sp_gf = 0.;
  _ec_gf = CE0;
  _sc_gf = 0.;

  // Construct equation ojects, first the 3 macro equations, then the NPAR micro eqs
  if (P2D)
  {
    _ep = new ElectrolytePotential(*_x_h1space);
    _sp = new SolidPotential(*_x_h1space);
  }
  _ec = new ElectrolyteConcentration(*_x_h1space);

  if (SPM || SPMe)
  {
    _sc.Append(new SolidConcentration(*_r_h1space[0], 0, 0, -1, NE));
    _sc.Append(new SolidConcentration(*_r_h1space[1], 1, 0, -1, PE));
  }
  else
  {
    Array<int> particle_dofs, particle_offsets;
    Array<Region> particle_regions;
    GetParticleDofs(particle_dofs, particle_regions, particle_offsets);

    for (unsigned p = 0; p < NPAR; p++)
    {
      auto rank_iter = std::upper_bound(particle_offsets.begin(), particle_offsets.end(), p);
      int rank = std::distance(particle_offsets.begin(), rank_iter) - 1;
      bool owned = rank == Mpi::WorldRank();

      unsigned offset = particle_offsets[Mpi::WorldRank()];
      int dof = owned ? particle_dofs[p - offset] : -1;
      Region region = particle_regions[p];

      _sc.Append(new SolidConcentration(*_r_h1space[p], p, rank, dof, region));
    }
  }

  ConstructExchangeCurrent();
  ConstructOpenCircuitPotential();
  ConstructOverPotential();
  ConstructReactionCurrent();

  if (P2D)
  {
    // Construct space for discontinuous functions like the reaction current j
    _x_l2space = new ParFiniteElementSpace(_x_h1space->GetParMesh(),
                                           new L2_FECollection(_scl_ir.GetNPoints() - 1, 1));

    // Build special integration rule to be used only for self-consistency loop
    for (int i = 0; i < _scl_ir.GetNPoints(); i++)
      _scl_ir.IntPoint(i).weight = NX;
    _scl_irs[Geometry::Type::SEGMENT] = &_scl_ir;
  }
}

void
EChemOperator::ImplicitSolve(const real_t dt, const Vector & x, Vector & k)
{
  // Solve the equation:
  //   M x(t + dt) / dt = - K x(t + dt) + M x(t) / dt + b <=>
  //   (M / dt + K) x(t + dt) = M x(t) / dt + b
  // for k := x(t + dt), where K is linearized by using x from the previous timestep

  if (P2D)
  {
    ParGridFunction j_gf(_x_l2space);
    real_t j_norm;
    _ep->Reset();

    do
    {
      // save previous iteration reaction current and its l2 norm
      j_gf.ProjectCoefficient(*_j);
      j_norm = j_gf.Norml2();
      // get the global l2 norm if in parallel
      j_norm *= j_norm;
      MPI_Allreduce(MPI_IN_PLACE, &j_norm, 1, MFEM_MPI_REAL_T, MPI_SUM, MPI_COMM_WORLD);
      j_norm = sqrt(j_norm);

      // assemble each individual block of _Ap and _bp
      UpdatePotentialEquations();

      // put _Ap and _bp together
      _Bp(EPP, EPP) = &_ep->GetK();
      _Bp(SPP, SPP) = &_sp->GetK();
      delete _Ap;
      _Ap = mfem::HypreParMatrixFromBlocks(_Bp);
      _bp.GetBlock(EPP) = _ep->GetZ();
      _bp.GetBlock(SPP) = _sp->GetZ();

      // solve for xp(t + dt)
      _Solver.SetOperator(*_Ap);
      _Solver.Mult(_bp, _xp);

      SetGridFunctionsFromTrueVectors();
    } while (j_gf.ComputeL2Error(*_j, _scl_irs) > _scl_threshold * j_norm);
  }

  // assemble each individual block of _Ac and _bc
  UpdateConcentrationEquations();

  // put _Ac and _bc together
  delete _Bc(ECC, ECC);
  _Bc(ECC, ECC) = Add(1. / dt, _ec->GetM(), 1, _ec->GetK());
  _bc.GetBlock(ECC) = _ec->GetZ();
  _ec->GetM().AddMult(_xc.GetBlock(ECC), _bc.GetBlock(ECC), 1. / dt);
  for (unsigned p = 0; p < NPAR; p++)
  {
    delete _Bc(SCC + p, SCC + p);
    _Bc(SCC + p, SCC + p) = Add(1. / dt, _sc[p]->GetM(), 1, _sc[p]->GetK());
    _bc.GetBlock(SCC + p) = _sc[p]->GetZ();
    _sc[p]->GetM().AddMult(_xc.GetBlock(SCC + p), _bc.GetBlock(SCC + p), 1. / dt);
  }
  delete _Ac;
  _Ac = mfem::HypreParMatrixFromBlocks(_Bc);

  // solve for xc(t + dt)
  _Solver.SetOperator(*_Ac);
  _Solver.Mult(_bc, _xc);

  // k is the solution at t + dt, which we already have in _x
  k.MakeRef(_x, 0);

  SetGridFunctionsFromTrueVectors();
}

void
EChemOperator::SetGridFunctionsFromTrueVectors()
{
  _ep_gf.SetFromTrueDofs(_x.GetBlock(EP));
  _sp_gf.SetFromTrueDofs(_x.GetBlock(SP));
  _ec_gf.SetFromTrueDofs(_x.GetBlock(EC));
  SetSurfaceConcentration();
  SetReferencePotential();
}

void
EChemOperator::UpdatePotentialEquations()
{
  _ep->Update(_ec_gfc, *_j);
  _sp->Update(*_j);
}

void
EChemOperator::UpdateConcentrationEquations()
{
  _ec->Update(_ec_gfc, *_j);
  const Array<real_t> & j = GetParticleReactionCurrent();
  for (unsigned p = 0; p < NPAR; p++)
    _sc[p]->Update(ConstantCoefficient(j[p]));
}

//
// Surface Concentration
//

const real_t &
EChemOperator::GetSurfaceConcentration(const Region & r)
{
  MFEM_ASSERT(SPM || SPMe, "Cannot get surface concentration, only SPM and SPMe are supported.");
  MFEM_ASSERT(r == NE || r == PE,
              "Cannot get surface concentration, only negative (NE) and positive electrodes (PE) "
              "are supported.");
  return _sc_array[r];
}

void
EChemOperator::SetSurfaceConcentration()
{
  if (SPM || SPMe)
    for (unsigned p = 0; p < NPAR; p++)
    {
      Region r = _sc[p]->GetParticleRegion();
      real_t sc = (r == NE ? CN0 : CP0) + _sc[p]->SurfaceConcentration(_x);
      _sc_array[r] = sc;
      MPI_Bcast(&_sc_array[r], 1, MFEM_MPI_REAL_T, _sc[p]->GetParticleRank(), MPI_COMM_WORLD);
    }
  else if (P2D)
  {
    for (unsigned p = 0; p < NPAR; p++)
    {
      Region r = _sc[p]->GetParticleRegion();
      real_t sc = (r == NE ? CN0 : CP0) + _sc[p]->SurfaceConcentration(_x);
      if (_sc[p]->IsParticleOwned())
        _sc_gf(_sc[p]->GetParticleDof()) = sc;
    }
    // Apply prolongation after restriction. Might be unnecessary, but guarantees
    // all processors have the right information for all their local dofs.
    _sc_gf.SetTrueVector();
    _sc_gf.SetFromTrueVector();
  }
}

//
// Reference Potential
//

const real_t &
EChemOperator::GetReferencePotential(const Region & r)
{
  MFEM_ASSERT(P2D, "Cannot get reference potential, only P2D is supported.");
  MFEM_ASSERT(r == E || r == NE || r == PE,
              "Cannot get reference potential, only electrolyte (E), negative (NE) and positive "
              "electrodes (PE) are supported.");
  return _rp_array[r];
}

void
EChemOperator::SetReferencePotential()
{
  if (P2D)
  {
    _rp_array[E] = 0.;
    _rp_array[PE] = 0.;

    real_t Inp = GetElectrodeReactionCurrent(NE, 1.0);
    real_t Inn = GetElectrodeReactionCurrent(NE, -1.0);
    real_t Ipp = GetElectrodeReactionCurrent(PE, 1.0);
    real_t Ipn = GetElectrodeReactionCurrent(PE, -1.0);

    _rp_array[E] = -2.0 * T * log((I + sqrt(4.0 * Inp * Inn + I * I)) / (2.0 * Inp));
    _rp_array[PE] =
        2.0 * T * log((-I + sqrt(4.0 * Ipp * Ipn + I * I)) / (2.0 * Ipp)) + _rp_array[E];
  }
}

//
// Partial Reaction Current for each electrode
//

real_t
EChemOperator::GetElectrodeReactionCurrent(const Region & r, const int & sign)
{
  MFEM_ASSERT(r == NE || r == PE,
              "Cannot get partial electrode reaction current, only negative (NE) and positive "
              "electrodes (PE) are supported.");

  ElectrodeReactionCurrentCoefficient j(r, sign, *_jex, *_op);
  QuadratureSpace x_qspace(_x_h1space->GetParMesh(), 2 * _x_h1space->FEColl()->GetOrder());
  return x_qspace.Integrate(j);
}

//
// Reaction Current for each particle
//

Array<real_t>
EChemOperator::GetParticleReactionCurrent()
{
  Array<real_t> j(NPAR);

  if (SPM || SPMe)
    for (unsigned p = 0; p < NPAR; p++)
      j[p] = GetReactionCurrent(_sc[p]->GetParticleRegion());
  else if (P2D)
  {
    ParGridFunction j_gf(_x_h1space);

    { // ProjectDiscCoefficient uses the element w/ maximal attribute for shared dofs
      for (int elem = 0; elem < _x_h1space->GetParMesh()->GetNE(); elem++)
        if (_x_h1space->GetAttribute(elem) == SEP)
          _x_h1space->GetParMesh()->SetAttribute(elem, 0);

      j_gf.ProjectDiscCoefficient(*_j);

      for (int elem = 0; elem < _x_h1space->GetParMesh()->GetNE(); elem++)
        if (_x_h1space->GetAttribute(elem) == 0)
          _x_h1space->GetParMesh()->SetAttribute(elem, SEP);
    }

    for (unsigned p = 0; p < NPAR; p++)
    {
      j[p] = _sc[p]->IsParticleOwned() ? j_gf(_sc[p]->GetParticleDof()) : 0;

      if (_sc[p]->GetParticleRank() == _sc[p]->GetSurfaceRank())
        continue;

      if (_sc[p]->IsParticleOwned())
        MPI_Send(&j[p], 1, MFEM_MPI_REAL_T, _sc[p]->GetSurfaceRank(), p, MPI_COMM_WORLD);

      if (_sc[p]->IsSurfaceOwned())
        MPI_Recv(&j[p],
                 1,
                 MFEM_MPI_REAL_T,
                 _sc[p]->GetParticleRank(),
                 p,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }
  }

  return j;
}

//
// Reaction Current
//

const real_t &
EChemOperator::GetReactionCurrent(const Region & r)
{
  MFEM_ASSERT(SPM || SPMe,
              "Cannot get constant reaction current, only SPM and SPMe are supported.");
  MFEM_ASSERT(r == NE || r == PE,
              "Cannot get constant reaction current, only negative (NE) and positive electrodes "
              "(PE) are supported.");
  return _j->Eval()(r);
}

void
EChemOperator::ConstructReactionCurrent()
{
  if (SPM || SPMe)
    _j = new ReactionCurrentCoefficient();
  else if (P2D)
    _j = new ReactionCurrentCoefficient(T, *_jex, *_op);
}

//
// Exchange Current
//

const real_t &
EChemOperator::GetExchangeCurrent(const Region & r)
{
  MFEM_ASSERT(SPM || SPMe,
              "Cannot get constant exchange current, only SPM and SPMe are supported.");
  MFEM_ASSERT(r == NE || r == PE,
              "Cannot get constant exchange current, only negative (NE) and positive electrodes "
              "(PE) are supported.");
  return _jex->Eval()(r);
}

void
EChemOperator::ConstructExchangeCurrent()
{
  if (SPM)
    _jex = new ExchangeCurrentCoefficient(
        KN, KP, GetSurfaceConcentration(NE), GetSurfaceConcentration(PE), CE0);
  else if (SPMe)
    _jex = new ExchangeCurrentCoefficient(
        KN, KP, GetSurfaceConcentration(NE), GetSurfaceConcentration(PE), _ec_gfc);
  else if (P2D)
    _jex = new ExchangeCurrentCoefficient(KN, KP, _sc_gfc, _ec_gfc);
}

//
// Open Circuit Potential
//

const real_t &
EChemOperator::GetOpenCircuitPotential(const Region & r)
{
  MFEM_ASSERT(SPM || SPMe,
              "Cannot get constant open circuit potential, only SPM and SPMe are supported.");
  MFEM_ASSERT(r == NE || r == PE,
              "Cannot get constant open circuit potential, only negative (NE) and positive "
              "electrodes (PE) are supported.");
  return _ocp->Eval()(r);
}

void
EChemOperator::ConstructOpenCircuitPotential()
{
  if (SPM || SPMe)
    _ocp = new OpenCircuitPotentialCoefficient(
        UN, UP, GetSurfaceConcentration(NE), GetSurfaceConcentration(PE));
  else if (P2D)
    _ocp = new OpenCircuitPotentialCoefficient(UN, UP, _sc_gfc);
}

//
// OverPotential
//

const real_t &
EChemOperator::GetOverPotential(const Region & r)
{
  MFEM_ASSERT(SPM || SPMe, "Cannot get constant overpotential, only SPM and SPMe are supported.");
  MFEM_ASSERT(r == NE || r == PE,
              "Cannot get constant overpotential, only negative (NE) and positive electrodes (PE) "
              "are supported.");
  return _op->Eval()(r);
}

void
EChemOperator::ConstructOverPotential()
{
  if (SPM || SPMe)
    _op = new OverPotentialCoefficient(T, *_jex);
  else if (P2D)
    _op = new OverPotentialCoefficient(
        GetReferencePotential(E), GetReferencePotential(PE), _sp_gfc, _ep_gfc, *_ocp);
}

//
// Voltage
//

real_t
EChemOperator::GetVoltage()
{
  // Definition from JuBat: https://doi.org/10.1016/j.est.2023.107512
  if (SPM)
    return phi_scale * (GetOpenCircuitPotential(PE) - GetOpenCircuitPotential(NE) +
                        GetOverPotential(PE) - GetOverPotential(NE));
  else if (SPMe)
    return phi_scale *
           (GetOpenCircuitPotential(PE) - GetOpenCircuitPotential(NE) + GetOverPotential(PE) -
            GetOverPotential(NE) - GetVoltageMarquisCorrection());
  else if (P2D)
    return phi_scale * (GetReferencePotential(PE) - GetReferencePotential(NE));

  mfem_error("Unreachable as method must be one of SPM, SPMe or P2D.");
}

real_t
EChemOperator::GetVoltageMarquisCorrection()
{
  PWCoefficient ce_pwc;
  QuadratureSpace x_qspace(_x_h1space->GetParMesh(), _x_h1space->FEColl()->GetOrder());

  ce_pwc.UpdateCoefficient(NE, _ec_gfc);
  real_t ce_ne_int = x_qspace.Integrate(ce_pwc) * (NX / NNE);
  ce_pwc.ZeroCoefficient(NE);

  ce_pwc.UpdateCoefficient(PE, _ec_gfc);
  real_t ce_pe_int = x_qspace.Integrate(ce_pwc) * (NX / NPE);
  ce_pwc.ZeroCoefficient(PE);

  real_t eta_c = (2.0 * T / CE0) * (1 - TPLUS) * (ce_ne_int - ce_pe_int);
  real_t dphie = (I / KS) * (LNE / BNE / 3.0 + LSEP / BSEP + LPE / BPE / 3.0);
  real_t dphis = I / 3 * (LNE / SIGN + LPE / SIGP);

  return eta_c + dphie + dphis;
}

void
EChemOperator::GetParticleDofs(Array<int> & my_particle_dofs,
                               Array<Region> & particle_regions,
                               Array<int> & particle_offsets)
{
  Array<int> gtdofs;
  Array<Region> regions;
  for (int e = 0; e < _x_h1space->GetNE(); e++)
  {
    Array<int> dofs;
    _x_h1space->GetElementDofs(e, dofs);
    for (int d : dofs)
    {
      int gtdof = _x_h1space->GetGlobalTDofNumber(d);
      Region r = Region(_x_h1space->GetAttribute(e));
      gtdofs.Append(gtdof);
      regions.Append(r);
    }
  }

  const unsigned n_gtdofs = NX * (_x_h1space->FEColl()->GetOrder() + 1);

  gtdofs.SetSize(n_gtdofs, -1);
  Array<int> all_gtdofs(n_gtdofs * Mpi::WorldSize());
  MPI_Allgather(
      gtdofs.GetData(), n_gtdofs, MPI_INT, all_gtdofs.GetData(), n_gtdofs, MPI_INT, MPI_COMM_WORLD);

  regions.SetSize(n_gtdofs, UNKNOWN);
  Array<Region> all_regions(n_gtdofs * Mpi::WorldSize());
  MPI_Allgather(regions.GetData(),
                n_gtdofs,
                MPI_INT,
                all_regions.GetData(),
                n_gtdofs,
                MPI_INT,
                MPI_COMM_WORLD);

  Array<Region> my_particle_regions;
  for (int d = 0; d < _x_h1space->GetNDofs(); d++)
  {
    int ltdof = _x_h1space->GetLocalTDofNumber(d);
    int gtdof = _x_h1space->GetGlobalTDofNumber(d);
    Region r = UNKNOWN;
    if (ltdof != -1)
      for (unsigned i = 0; i < n_gtdofs * Mpi::WorldSize(); i++)
        if (gtdof == all_gtdofs[i] && (r = all_regions[i]) != SEP)
        {
          my_particle_dofs.Append(d);
          my_particle_regions.Append(r);
          break;
        }
  }
  my_particle_regions.SetSize(n_gtdofs, UNKNOWN);

  particle_regions.SetSize(n_gtdofs * Mpi::WorldSize(), UNKNOWN);
  MPI_Allgather(my_particle_regions.GetData(),
                n_gtdofs,
                MPI_INT,
                particle_regions.GetData(),
                n_gtdofs,
                MPI_INT,
                MPI_COMM_WORLD);
  while (particle_regions.Find(UNKNOWN) != -1)
    particle_regions.DeleteFirst(UNKNOWN);

  int my_particles = my_particle_dofs.Size();
  particle_offsets.SetSize(Mpi::WorldSize());
  MPI_Allgather(&my_particles, 1, MPI_INT, particle_offsets.GetData(), 1, MPI_INT, MPI_COMM_WORLD);

  particle_offsets.Prepend(0);
  particle_offsets.PartialSum();
  MFEM_ASSERT(unsigned(particle_offsets[Mpi::WorldSize()]) == NPAR &&
                  unsigned(particle_regions.Size()) == NPAR,
              "Failed to distribute particles across processors.");
}
