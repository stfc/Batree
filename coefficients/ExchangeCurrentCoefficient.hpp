#include "mfem.hpp"

class ExchangeCurrentCoefficient : public mfem::Coefficient
{
private:
  mfem::GridFunctionCoefficient * _surface_concentration_gfc = nullptr;
  mfem::GridFunctionCoefficient * _electrolyte_concentration_gfc = nullptr;

  const mfem::real_t * _scn = nullptr;
  const mfem::real_t * _scp = nullptr;

  mfem::TransformedCoefficient _jex_ne_tc;
  mfem::TransformedCoefficient _jex_pe_tc;

  const mfem::real_t * _kn = nullptr;
  const mfem::real_t * _kp = nullptr;

  mfem::Vector _jex_vec;
  mfem::PWConstCoefficient _jex_pwcc;
  mfem::PWCoefficient _jex_pwc;

public:
  /// SPM
  ExchangeCurrentCoefficient(const mfem::real_t & kn,
                             const mfem::real_t & kp,
                             const mfem::real_t & scn,
                             const mfem::real_t & scp,
                             const mfem::real_t & ec)
    : _scn(&scn),
      _scp(&scp),
      _jex_ne_tc(nullptr, [](mfem::real_t) { return 0; }),
      _jex_pe_tc(nullptr, [](mfem::real_t) { return 0; }),
      _jex_vec({kn * sqrt(ec), 0., kp * sqrt(ec)}),
      _jex_pwcc(3)
  {
  }

  /// SPMe
  ExchangeCurrentCoefficient(const mfem::real_t & kn,
                             const mfem::real_t & kp,
                             const mfem::real_t & scn,
                             const mfem::real_t & scp,
                             mfem::GridFunctionCoefficient & ec)
    : _electrolyte_concentration_gfc(&ec),
      _scn(&scn),
      _scp(&scp),
      _jex_ne_tc(_electrolyte_concentration_gfc, [=](mfem::real_t ec) { return kn * sqrt(ec); }),
      _jex_pe_tc(_electrolyte_concentration_gfc, [=](mfem::real_t ec) { return kp * sqrt(ec); }),
      _jex_pwcc(3)
  {
  }

  /// P2D
  ExchangeCurrentCoefficient(const mfem::real_t & kn,
                             const mfem::real_t & kp,
                             mfem::GridFunctionCoefficient & sc,
                             mfem::GridFunctionCoefficient & ec)
    : _surface_concentration_gfc(&sc),
      _electrolyte_concentration_gfc(&ec),
      _jex_ne_tc(nullptr, [](mfem::real_t) { return 0; }),
      _jex_pe_tc(nullptr, [](mfem::real_t) { return 0; }),
      _kn(&kn),
      _kp(&kp)
  {
  }

  /// SPM(e)
  virtual mfem::PWConstCoefficient & Eval()
  {
    /// SPMe
    if (_electrolyte_concentration_gfc)
    {
      mfem::ParFiniteElementSpace * x_h1space =
          static_cast<const mfem::ParGridFunction *>(
              _electrolyte_concentration_gfc->GetGridFunction())
              ->ParFESpace();
      mfem::QuadratureSpace x_qspace(x_h1space->GetParMesh(), 2 * x_h1space->FEColl()->GetOrder());

      /// NE
      _jex_pwc.UpdateCoefficient(NE, _jex_ne_tc);
      mfem::real_t integral_ne = x_qspace.Integrate(_jex_pwc);
      _jex_pwc.ZeroCoefficient(NE);

      /// PE
      _jex_pwc.UpdateCoefficient(PE, _jex_pe_tc);
      mfem::real_t integral_pe = x_qspace.Integrate(_jex_pwc);
      _jex_pwc.ZeroCoefficient(PE);

      _jex_pwcc(NE) = integral_ne * sqrt(*_scn * (1 - *_scn)) / NNE * NX;
      _jex_pwcc(PE) = integral_pe * sqrt(*_scp * (1 - *_scp)) / NPE * NX;
    }
    /// SPM
    else
    {
      _jex_pwcc(NE) = _jex_vec(NE - 1) * sqrt(*_scn * (1 - *_scn));
      _jex_pwcc(PE) = _jex_vec(PE - 1) * sqrt(*_scp * (1 - *_scp));
    }

    return _jex_pwcc;
  }

  /// P2D
  virtual mfem::real_t Eval(mfem::ElementTransformation & Tr,
                            const mfem::IntegrationPoint & ip) override
  {
    switch (Tr.Attribute)
    {
      case NE:
      {
        const mfem::real_t sc = _surface_concentration_gfc->Eval(Tr, ip);
        const mfem::real_t ec = _electrolyte_concentration_gfc->Eval(Tr, ip);
        return *_kn * sqrt(sc * ec * (1 - sc));
      }
      case PE:
      {
        const mfem::real_t sc = _surface_concentration_gfc->Eval(Tr, ip);
        const mfem::real_t ec = _electrolyte_concentration_gfc->Eval(Tr, ip);
        return *_kp * sqrt(sc * ec * (1 - sc));
      }
      default:
        return 0;
    }
  }
};
