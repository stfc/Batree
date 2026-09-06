#include "mfem.hpp"

class OverPotentialCoefficient : public mfem::Coefficient
{
private:
  mfem::GridFunctionCoefficient * _solid_potential_gfc = nullptr;
  mfem::GridFunctionCoefficient * _electrolyte_potential_gfc = nullptr;

  ExchangeCurrentCoefficient * _jex = nullptr;
  OpenCircuitPotentialCoefficient * _ocp = nullptr;

  const mfem::real_t * _rpe = nullptr;
  const mfem::real_t * _rpp = nullptr;
  mfem::PWCoefficient _rp_pwc;

  mfem::PWConstCoefficient _op_pwcc;

public:
  /// SPM(e)
  OverPotentialCoefficient(ExchangeCurrentCoefficient & jex) : _jex(&jex), _op_pwcc(3) {}

  /// P2D
  OverPotentialCoefficient(const mfem::real_t & rpe,
                           const mfem::real_t & rpp,
                           mfem::GridFunctionCoefficient & sp,
                           mfem::GridFunctionCoefficient & ep,
                           OpenCircuitPotentialCoefficient & ocp)
    : _solid_potential_gfc(&sp), _electrolyte_potential_gfc(&ep), _ocp(&ocp), _rpe(&rpe), _rpp(&rpp)
  {
  }

  /// SPM(e)
  virtual mfem::PWConstCoefficient & Eval()
  {
    _op_pwcc(NE) = 2 * T * asinh(+I / AN / LNE / 2.0 / _jex->Eval()(NE));
    _op_pwcc(PE) = 2 * T * asinh(-I / AP / LPE / 2.0 / _jex->Eval()(PE));
    return _op_pwcc;
  }

  /// P2D
  virtual mfem::real_t Eval(mfem::ElementTransformation & Tr,
                            const mfem::IntegrationPoint & ip) override
  {
    switch (Tr.Attribute)
    {
      case NE:
      {
        const mfem::real_t sp = _solid_potential_gfc->Eval(Tr, ip);
        const mfem::real_t ep = _electrolyte_potential_gfc->Eval(Tr, ip);
        const mfem::real_t ocp = _ocp->Eval(Tr, ip);
        return sp - ep - ocp - *_rpe;
      }
      case PE:
      {
        const mfem::real_t sp = _solid_potential_gfc->Eval(Tr, ip);
        const mfem::real_t ep = _electrolyte_potential_gfc->Eval(Tr, ip);
        const mfem::real_t ocp = _ocp->Eval(Tr, ip);
        return sp - ep - ocp + *_rpp - *_rpe;
      }
      default:
        return 0;
    }
  }
};
