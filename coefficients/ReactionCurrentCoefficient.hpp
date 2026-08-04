#include "mfem.hpp"
using namespace mfem;

class ReactionCurrentCoefficient : public Coefficient
{
private:
  ExchangeCurrentCoefficient * _jex = nullptr;
  OverPotentialCoefficient * _op = nullptr;

  PWConstCoefficient _j_pwcc;

public:
  /// SPM(e)
  ReactionCurrentCoefficient() : _j_pwcc(Vector{+I / AN / LNE, 0., -I / AP / LPE}) {}

  /// P2D
  ReactionCurrentCoefficient(ExchangeCurrentCoefficient & jex, OverPotentialCoefficient & op)
    : _jex(&jex), _op(&op)
  {
  }

  /// SPM(e)
  virtual PWConstCoefficient & Eval() { return _j_pwcc; }

  /// P2D (and any integrators)
  virtual real_t Eval(ElementTransformation & Tr, const IntegrationPoint & ip) override
  {
    if (_jex)
      return 2 * _jex->Eval(Tr, ip) * sinh(.5 * _op->Eval(Tr, ip) / T);
    else
      return _j_pwcc.Eval(Tr, ip);
  }
};
