class ReactionCurrentCoefficient : public mfem::Coefficient
{
private:
  ExchangeCurrentCoefficient * _jex = nullptr;
  OverPotentialCoefficient * _op = nullptr;

  /// SPM(e): the cell's scaled current, owned by its EChemOperator
  const mfem::real_t * _current = nullptr;

  mfem::PWConstCoefficient _j_pwcc;

public:
  /// SPM(e)
  ReactionCurrentCoefficient(const mfem::real_t & current) : _current(&current), _j_pwcc(3) {}

  /// P2D
  ReactionCurrentCoefficient(ExchangeCurrentCoefficient & jex, OverPotentialCoefficient & op)
    : _jex(&jex), _op(&op)
  {
  }

  /// SPM(e)
  virtual mfem::PWConstCoefficient & Eval()
  {
    _j_pwcc(NE) = +*_current / AN / LNE;
    _j_pwcc(PE) = -*_current / AP / LPE;
    return _j_pwcc;
  }

  /// P2D (and any integrators)
  virtual mfem::real_t Eval(mfem::ElementTransformation & Tr,
                            const mfem::IntegrationPoint & ip) override
  {
    if (_jex)
      return 2 * _jex->Eval(Tr, ip) * sinh(.5 * _op->Eval(Tr, ip) / T);
    else
      return Eval().Eval(Tr, ip);
  }
};
