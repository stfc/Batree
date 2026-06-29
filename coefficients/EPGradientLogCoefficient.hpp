#include "mfem.hpp"
using namespace mfem;

class EPGradientLogCoefficient : public VectorCoefficient
{
private:
  GridFunctionCoefficient * _electrolyte_concentration_gfc = nullptr;
  GradientGridFunctionCoefficient _electrolyte_concentration_ggfc;

public:
  EPGradientLogCoefficient(GridFunctionCoefficient & ec)
    : VectorCoefficient(1),
      _electrolyte_concentration_gfc(&ec),
      _electrolyte_concentration_ggfc(ec.GetGridFunction())
  {
  }

  virtual void Eval(Vector &V, ElementTransformation & T, const IntegrationPoint & ip) override
  {
    _electrolyte_concentration_ggfc.Eval(V, T, ip);
    const real_t ec = _electrolyte_concentration_gfc->Eval(T, ip);
    V[0] *= 2 * constants::T * (1 - TPLUS) * Kappa(ec) / ec;

    switch (T.Attribute)
    {
      case NE:
        V[0] *= BNE / (LNE / NNE * NX);
        break;
      case SEP:
        V[0] *= BSEP / (LSEP / NSEP * NX);
        break;
      case PE:
        V[0] *= BPE / (LPE / NPE * NX);
        break;
    }
  }
};
