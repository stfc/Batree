#include "mfem.hpp"
using namespace mfem;

class OverPotentialCoefficient : public Coefficient
{
private:
  GridFunctionCoefficient * _solid_potential_gfc = nullptr;
  GridFunctionCoefficient * _electrolyte_potential_gfc = nullptr;

  ExchangeCurrentCoefficient * _jex = nullptr;
  OpenCircuitPotentialCoefficient * _ocp = nullptr;

  const real_t * _rpe = nullptr;
  const real_t * _rpp = nullptr;
  PWCoefficient _rp_pwc;

  PWConstCoefficient _op_pwcc;

public:
  /// SPM(e)
  OverPotentialCoefficient(const real_t & T, ExchangeCurrentCoefficient & jex)
    : _jex(&jex), _op_pwcc(3)
  {
  }

  /// P2D
  OverPotentialCoefficient(const real_t & rpe,
                           const real_t & rpp,
                           GridFunctionCoefficient & sp,
                           GridFunctionCoefficient & ep,
                           OpenCircuitPotentialCoefficient & ocp)
    : _solid_potential_gfc(&sp), _electrolyte_potential_gfc(&ep), _ocp(&ocp), _rpe(&rpe), _rpp(&rpp)
  {
  }

  /// SPM(e)
  virtual PWConstCoefficient & Eval()
  {
    _op_pwcc(NE) = 2 * T * asinh(+I / AN / LNE / 2.0 / _jex->Eval()(NE));
    _op_pwcc(PE) = 2 * T * asinh(-I / AP / LPE / 2.0 / _jex->Eval()(PE));
    return _op_pwcc;
  }

  /// P2D
  virtual real_t Eval(ElementTransformation & T, const IntegrationPoint & ip) override
  {
    switch (T.Attribute)
    {
      case NE:
      {
        const real_t sp = _solid_potential_gfc->Eval(T, ip);
        const real_t ep = _electrolyte_potential_gfc->Eval(T, ip);
        const real_t ocp = _ocp->Eval(T, ip);
        return sp - ep - ocp - *_rpe;
      }
      case PE:
      {
        const real_t sp = _solid_potential_gfc->Eval(T, ip);
        const real_t ep = _electrolyte_potential_gfc->Eval(T, ip);
        const real_t ocp = _ocp->Eval(T, ip);
        return sp - ep - ocp + *_rpp - *_rpe;
      }
      default:
        return 0;
    }
  }
};
