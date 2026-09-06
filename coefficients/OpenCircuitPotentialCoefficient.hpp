class OpenCircuitPotentialCoefficient : public mfem::Coefficient
{
private:
  mfem::GridFunctionCoefficient * _surface_concentration_gfc = nullptr;

  const std::function<mfem::real_t(mfem::real_t)> _un;
  const std::function<mfem::real_t(mfem::real_t)> _up;

  const mfem::real_t * _scn = nullptr;
  const mfem::real_t * _scp = nullptr;

  mfem::PWConstCoefficient _ocp_pwcc;

public:
  /// SPM(e)
  OpenCircuitPotentialCoefficient(const std::function<mfem::real_t(mfem::real_t)> & un,
                                  const std::function<mfem::real_t(mfem::real_t)> & up,
                                  const mfem::real_t & scn,
                                  const mfem::real_t & scp)
    : _un(un), _up(up), _scn(&scn), _scp(&scp), _ocp_pwcc(3)
  {
  }

  /// P2D
  OpenCircuitPotentialCoefficient(const std::function<mfem::real_t(mfem::real_t)> & un,
                                  const std::function<mfem::real_t(mfem::real_t)> & up,
                                  mfem::GridFunctionCoefficient & sc)
    : _surface_concentration_gfc(&sc), _un(un), _up(up)
  {
  }

  /// SPM(e)
  virtual mfem::PWConstCoefficient & Eval()
  {
    _ocp_pwcc(NE) = _un(*_scn);
    _ocp_pwcc(PE) = _up(*_scp);
    return _ocp_pwcc;
  }

  /// P2D
  virtual mfem::real_t Eval(mfem::ElementTransformation & Tr,
                            const mfem::IntegrationPoint & ip) override
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
