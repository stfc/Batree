#pragma once

#include <cmath>
#include <string>

/// Stand-in for the per-element P2D cell model, used to test the collector-cell coupling
/// before the real P2D is connected. For each current-collector element e it returns the
/// cell voltage V̄ₑ(īₑ) and the local resistance R̄ₑ = −dV̄ₑ/dīₑ (scaled units).
///
///   linear:     V = OCV − R0 i
///   nonlinear:  V = OCV − R0 i − 2 asinh(i / (2 j0))     (ohmic + Butler–Volmer-like)
///
/// These are the same models as the Python prototype (newton_proto.py). Later, Eval() is
/// replaced by one P2D step per element: V̄ₑ from the step, R̄ₑ from a finite difference.
class LocalModel
{
public:
  LocalModel(const std::string & type, mfem::real_t ocv, mfem::real_t R0, mfem::real_t j0)
    : _nonlinear(type == "nonlinear"), _ocv(ocv), _R0(R0), _j0(j0)
  {
    if (type != "linear" && type != "nonlinear")
      mfem::mfem_error("LocalModel: type must be 'linear' or 'nonlinear'.");
  }

  /// V(e) = V̄ₑ(i(e)) and R(e) = −dV̄ₑ/di at i(e), for every local element e.
  void Eval(const mfem::Vector & i, mfem::Vector & V, mfem::Vector & R) const
  {
    V.SetSize(i.Size());
    R.SetSize(i.Size());
    for (int e = 0; e < i.Size(); e++)
    {
      V(e) = _ocv - _R0 * i(e);
      R(e) = _R0;
      if (_nonlinear)
      {
        V(e) -= 2.0 * std::asinh(i(e) / (2.0 * _j0));
        R(e) += 2.0 / std::sqrt(4.0 * _j0 * _j0 + i(e) * i(e));
      }
    }
  }

private:
  const bool _nonlinear;
  const mfem::real_t _ocv;
  const mfem::real_t _R0;
  const mfem::real_t _j0;
};
