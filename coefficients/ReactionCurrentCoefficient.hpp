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
  ReactionCurrentCoefficient(const real_t & T,
                             ExchangeCurrentCoefficient & jex,
                             OverPotentialCoefficient & op)
    : _jex(&jex), _op(&op)
  {
  }

  /// SPM(e)
  virtual PWConstCoefficient & Eval() { return _j_pwcc; }

  /// P2D (and any integrators)
  virtual real_t Eval(ElementTransformation & T, const IntegrationPoint & ip) override
  {
    if (_jex)
      return 2 * _jex->Eval(T, ip) * sinh(.5 * _op->Eval(T, ip) / constants::T);
    else
      return _j_pwcc.Eval(T, ip);
  }
};
