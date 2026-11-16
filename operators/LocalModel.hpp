#pragma once

#include <cmath>
#include <string>

/// A cell model for every current-collector element e: given the element currents i, it returns
/// the cell voltage V̄ₑ(iₑ) and, when asked, the local resistance R̄ₑ = −dV̄ₑ/diₑ (scaled units).
class CellModel
{
public:
  /// Vbar(e) = V̄ₑ(i(e)) for every local element e. If slope is true, also Rbar(e) = −dV̄ₑ/di at
  /// i(e); otherwise Rbar holds the last slope the model computed (it is only reported then).
  virtual void
  Eval(const mfem::Vector & i, mfem::Vector & Vbar, mfem::Vector & Rbar, bool slope) = 0;

  virtual ~CellModel() = default;
};

/// Stand-in for the per-element P2D cell model, used to test the collector-cell coupling.
///
///   linear:     V = OCV − R0 i
///   nonlinear:  V = OCV − R0 i − 2 asinh(i / (2 j0))     (ohmic + Butler–Volmer-like)
///
/// These are the same models as the Python prototype (newton_proto.py). EChemCellModel couples
/// the real EChemOperator instead.
class LocalModel : public CellModel
{
public:
  LocalModel(const std::string & type, mfem::real_t ocv, mfem::real_t R0, mfem::real_t j0)
    : _nonlinear(type == "nonlinear"), _ocv(ocv), _R0(R0), _j0(j0)
  {
    if (type != "linear" && type != "nonlinear")
      mfem::mfem_error("LocalModel: type must be 'linear' or 'nonlinear'.");
  }

  /// The slope is exact and cheap, so it is always computed.
  void Eval(const mfem::Vector & i, mfem::Vector & Vbar, mfem::Vector & Rbar, bool) override
  {
    Vbar.SetSize(i.Size());
    Rbar.SetSize(i.Size());
    for (int e = 0; e < i.Size(); e++)
    {
      Vbar(e) = _ocv - _R0 * i(e);
      Rbar(e) = _R0;
      if (_nonlinear)
      {
        Vbar(e) -= 2.0 * std::asinh(i(e) / (2.0 * _j0));
        Rbar(e) += 2.0 / std::sqrt(4.0 * _j0 * _j0 + i(e) * i(e));
      }
    }
  }

private:
  const bool _nonlinear;
  const mfem::real_t _ocv;
  const mfem::real_t _R0;
  const mfem::real_t _j0;
};
