#include "mfem.hpp"
using namespace mfem;

class OpenCircuitPotentialCoefficient : public Coefficient
{
private:
  GridFunctionCoefficient * _surface_concentration_gfc = nullptr;

  const std::function<real_t(real_t)> _un;
  const std::function<real_t(real_t)> _up;

  const real_t * _scn = nullptr;
  const real_t * _scp = nullptr;

  PWConstCoefficient _ocp_pwcc;

public:
  /// SPM(e)
  OpenCircuitPotentialCoefficient(const std::function<real_t(real_t)> & un,
                                  const std::function<real_t(real_t)> & up,
                                  const real_t & scn,
                                  const real_t & scp)
    : _un(un), _up(up), _scn(&scn), _scp(&scp), _ocp_pwcc(3)
  {
  }

  /// P2D
  OpenCircuitPotentialCoefficient(const std::function<real_t(real_t)> & un,
                                  const std::function<real_t(real_t)> & up,
                                  GridFunctionCoefficient & sc)
    : _surface_concentration_gfc(&sc), _un(un), _up(up)
  {
  }

  /// SPM(e)
  virtual PWConstCoefficient & Eval()
  {
    _ocp_pwcc(NE) = _un(*_scn);
    _ocp_pwcc(PE) = _up(*_scp);
    return _ocp_pwcc;
  }

  /// P2D
  virtual real_t Eval(ElementTransformation & Tr, const IntegrationPoint & ip) override
  {
    switch (Tr.Attribute)
    {
      case NE:
        return _un(_surface_concentration_gfc->Eval(Tr, ip));
      case PE:
        return _up(_surface_concentration_gfc->Eval(Tr, ip));
      default:
        return 0;
    }
  }
};
