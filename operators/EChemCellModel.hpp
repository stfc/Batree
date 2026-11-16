#pragma once

#include "operators/EChemOperator.hpp"
#include "operators/LocalModel.hpp"

#include <algorithm>
#include <vector>

/// Couples every current-collector element to its own electrochemical cell (the SPM, SPMe or P2D
/// chosen with -m). One EChemOperator runs all of them: element e's cell is its saved state
/// vector. To step it, that state is copied into the operator's x, the current is set to iₑ and
/// one time step is taken. Since i0 is the cell's 1C current density, iₑ is element e's C-rate.
///
/// Within a time step, every Eval starts from the saved start-of-step states, so the coupling can
/// try as many currents as it needs. Commit() then accepts the end states of the last Eval.
///
/// Serial only: EChemOperator works on MPI_COMM_WORLD, so all ranks would have to step every
/// element together.
class EChemCellModel : public CellModel
{
public:
  /// t, dt: start time and time step. n: number of collector elements. h: finite-difference step
  /// in the scaled current for R̄ₑ = −(V̄ₑ(iₑ + h) − V̄ₑ(iₑ)) / h.
  EChemCellModel(EChemOperator & oper,
                 mfem::ODESolver & ode,
                 mfem::BlockVector & x,
                 mfem::real_t t,
                 mfem::real_t dt,
                 int n,
                 mfem::real_t h)
    : _oper(oper), _ode(ode), _x(x), _t(t), _dt(dt), _h(h), _start(n, x), _end(n, x), _Rbar(n)
  {
    _Rbar = 0.0;
  }

  void Eval(const mfem::Vector & i, mfem::Vector & Vbar, mfem::Vector & Rbar, bool slope) override
  {
    MFEM_VERIFY(i.Size() == int(_start.size()), "EChemCellModel: wrong number of elements.");
    Vbar.SetSize(i.Size());
    for (int e = 0; e < i.Size(); e++)
    {
      mfem::real_t Vh = 0.0;
      if (slope)
        Vh = Step(e, i(e) + _h);
      Vbar(e) = Step(e, i(e)); // taken last, so x holds the end state at iₑ
      _end[e] = _x;
      if (slope)
        _Rbar(e) = -(Vh - Vbar(e)) / _h;
    }
    Rbar.SetSize(i.Size());
    Rbar = _Rbar;
  }

  /// Accept the time step: the end states of the last Eval become the next start states.
  void Commit()
  {
    for (std::size_t e = 0; e < _start.size(); e++)
      _start[e] = _end[e];
    _t += _dt;
  }

  /// Area-weighted mean, minimum and maximum state of charge [%] of the committed states.
  void SoC(const mfem::Vector & area, mfem::real_t & mean, mfem::real_t & min, mfem::real_t & max)
  {
    mean = 0.0;
    min = mfem::infinity();
    max = -mfem::infinity();
    for (std::size_t e = 0; e < _start.size(); e++)
    {
      _x.Set(1.0, _start[e]);
      const mfem::real_t soc = _oper.GetSoC();
      mean += area(e) * soc;
      min = std::min(min, soc);
      max = std::max(max, soc);
    }
    mean /= area.Sum();
  }

  /// Cell time steps taken so far
  long Steps() const { return _steps; }

private:
  /// One time step of element e's cell from its start-of-step state at scaled current I.
  /// Returns the cell voltage in units of phi_scale.
  mfem::real_t Step(int e, mfem::real_t current)
  {
    _x.Set(1.0, _start[e]); // copy in place: the operator's block views of x stay valid
    _oper.SetConcentrationGridFunctionsFromTrueVectors();
    _oper.SetPotentialGridFunctionsFromTrueVectors();
    _oper.SetCurrent(current);
    mfem::real_t t = _t, dt = _dt;
    _ode.Step(_x, t, dt);
    _steps++;
    // Read the voltage now, before the next element's step changes the operator's state and
    // current (the SPM and SPMe voltages are computed from the current)
    return _oper.GetVoltage() / phi_scale;
  }

  EChemOperator & _oper;
  mfem::ODESolver & _ode;
  mfem::BlockVector & _x;
  mfem::real_t _t;
  const mfem::real_t _dt;
  const mfem::real_t _h;

  /// Start-of-step state of each element, and its end state from the last Eval
  std::vector<mfem::Vector> _start, _end;

  /// Last finite-difference slopes
  mfem::Vector _Rbar;

  long _steps = 0;
};
