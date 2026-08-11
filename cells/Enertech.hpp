#pragma once

#include "Cell.hpp"
#include "mfem.hpp"

using namespace mfem;

class Enertech : public Cell
{
public:
  inline const real_t cn0() const override { return 24108.; }
  inline const real_t cnmax() const override { return 29700; } // PyBAMM uses 28700
  inline const real_t rn() const override { return 5e-6; }
  inline const real_t ln() const override { return 7.65e-5; }
  inline const real_t eps_n() const override { return 0.33; }
  inline const real_t eps_n_fi() const override { return 0.06; }
  inline const real_t kn_dim() const override { return 9.648533289e-7; }
  inline const real_t Dn() const override { return 3.9e-14; }
  inline const real_t sig_n() const override { return 100.; }
  inline const real_t brugg_n() const override { return 2.914; }
  inline const real_t Un(real_t cs) const override
  {
    return -2058.29865 * pow(cs, 9) + 10040.08960 * pow(cs, 8) - 20824.86740 * pow(cs, 7) +
           23911.86578 * pow(cs, 6) - 16576.36920 * pow(cs, 5) + 7098.09151 * pow(cs, 4) -
           1845.43634 * pow(cs, 3) + 275.31114 * pow(cs, 2) - 21.20097 * cs + 0.84498;
  }

  inline const real_t cp0() const override { return 21725.; }
  inline const real_t cpmax() const override { return 49943; }
  inline const real_t rp() const override { return 3e-6; }
  inline const real_t lp() const override { return 6.8e-5; }
  inline const real_t eps_p() const override { return 0.32; }
  inline const real_t eps_p_fi() const override { return 0.06; }
  inline const real_t kp_dim() const override { return 9.648533289e-7; }
  inline const real_t Dp() const override { return 5.387e-15; }
  inline const real_t sig_p() const override { return 10.; }
  inline const real_t brugg_p() const override { return 1.83; }
  inline const real_t Up(real_t cs) const override
  {
    return -107897.40 * pow(cs, 9) + 677406.28 * pow(cs, 8) - 1873803.91 * pow(cs, 7) +
           2996535.44 * pow(cs, 6) - 3052331.36 * pow(cs, 5) + 2053377.31 * pow(cs, 4) -
           912135.88 * pow(cs, 3) + 257964.35 * pow(cs, 2) - 42146.98 * cs + 3035.67;
  }

  inline const real_t ls() const override { return 2.5e-5; }
  inline const real_t eps_s() const override { return 0.5; }
  inline const real_t brugg_s() const override { return 1.5; }

  inline const real_t tplus() const override { return 0.38; }
  inline const real_t ce0() const override { return 1000.0; }
  inline const real_t De(real_t ce) const override
  {
    const real_t T = 298.15;
    return pow(10.0, (-4.43 - 54.0 / (T - 229.0 - 5e-3 * ce) - 0.22e-3 * ce));
  }
  inline const real_t kappa(real_t ce) const override
  {
    const real_t T = 298.15;
    real_t inner = (-10.5 + 0.668e-3 * ce + 0.494e-6 * ce * ce) +
                   (0.074 - 1.78e-5 * ce - 8.86e-10 * ce * ce) * T +
                   (-6.96e-5 + 2.8e-8 * ce) * T * T;
    return 1e-4 * ce * inner * inner;
  }

  inline const unsigned nlayer() const override { return 34; }
  inline const real_t llayer() const override { return 0.051; }
  inline const real_t wlayer() const override { return 0.047; }

  inline const real_t I_typ() const override { return 2.28; }
};
