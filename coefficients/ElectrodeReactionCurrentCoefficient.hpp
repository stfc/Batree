#include "mfem.hpp"
using namespace mfem;

class ElectrodeReactionCurrentCoefficient : public Coefficient
{
private:
  ExchangeCurrentCoefficient * _jex = nullptr;
  OverPotentialCoefficient * _op = nullptr;

  const real_t _a;
  const Region _r;
  const int _sign;

public:
  /// P2D
  ElectrodeReactionCurrentCoefficient(const Region & r,
                                      const int & sign,
                                      ExchangeCurrentCoefficient & jex,
                                      OverPotentialCoefficient & op)
    : _jex(&jex),
      _op(&op),
      _a(r == NE   ? AN * LNE / NNE * NX
         : r == PE ? AP * LPE / NPE * NX
                   : 0),
      _r(r),
      _sign(sign)
  {
  }

  /// P2D
  virtual real_t Eval(ElementTransformation & Tr, const IntegrationPoint & ip) override
  {
    return Tr.Attribute == _r ? _a * _jex->Eval(Tr, ip) * exp(_sign * 0.5 * _op->Eval(Tr, ip)) : 0;
  }
};
