#pragma once

#include "Cell.hpp"
#include "mfem.hpp"

class LGM50 : public Cell
{
public:
  inline const mfem::real_t cn0() const override { return 29866.; }
  inline const mfem::real_t cnmax() const override { return 33133; }
  inline const mfem::real_t rn() const override { return 5.86e-6; }
  inline const mfem::real_t ln() const override { return 85.2e-6; }
  inline const mfem::real_t eps_n() const override { return 0.25; }
  inline const mfem::real_t eps_n_fi() const override { return 0.0; }
  inline const mfem::real_t kn_dim() const override { return 6.48e-7; }
  inline const mfem::real_t Dn() const override { return 3.3e-14; }
  inline const mfem::real_t sig_n() const override { return 215.; }
  inline const mfem::real_t brugg_n() const override { return 1.5; }
  inline const mfem::real_t Un(mfem::real_t cs) const override
  {
    return 1.9793 * exp(-39.3631 * cs) + 0.2482 - 0.0909 * tanh(29.8538 * (cs - 0.1234)) -
           0.04478 * tanh(14.9159 * (cs - 0.2769)) - 0.0205 * tanh(30.4444 * (cs - 0.6103));
  }

  inline const mfem::real_t cp0() const override { return 17038.; }
  inline const mfem::real_t cpmax() const override { return 63104; }
  inline const mfem::real_t rp() const override { return 5.22e-6; }
  inline const mfem::real_t lp() const override { return 75.6e-6; }
  inline const mfem::real_t eps_p() const override { return 0.335; }
  inline const mfem::real_t eps_p_fi() const override { return 0.0; }
  inline const mfem::real_t kp_dim() const override { return 3.42e-6; }
  inline const mfem::real_t Dp() const override { return 4.0e-15; }
  inline const mfem::real_t sig_p() const override { return 0.18; }
  inline const mfem::real_t brugg_p() const override { return 1.5; }
  inline const mfem::real_t Up(mfem::real_t cs) const override
  {
    return -0.8090 * cs + 4.4875 - 0.0428 * tanh(18.5138 * (cs - 0.5542)) -
           17.7326 * tanh(15.7890 * (cs - 0.3117)) + 17.5842 * tanh(15.9308 * (cs - 0.3120));
  }

  inline const mfem::real_t ls() const override { return 12e-6; }
  inline const mfem::real_t eps_s() const override { return 0.47; }
  inline const mfem::real_t brugg_s() const override { return 1.5; }

  inline const mfem::real_t tplus() const override { return 0.2594; }
  inline const mfem::real_t ce0() const override { return 1000.0; }
  inline const mfem::real_t De(mfem::real_t ce) const override
  {
    return 8.794e-11 * pow(ce / 1000, 2) - 3.972e-10 * (ce / 1000) + 4.862e-10;
  }
  inline const mfem::real_t kappa(mfem::real_t ce) const override
  {
    return 0.1297 * pow(ce / 1000, 3) - 2.51 * pow(ce / 1000, 1.5) + 3.329 * (ce / 1000);
  }

  inline const unsigned nlayer() const override { return 1; }
  inline const mfem::real_t llayer() const override { return 1.58; }
  inline const mfem::real_t wlayer() const override { return 6.5e-2; }

  inline const mfem::real_t I_typ() const override { return 5.0; }

  inline const mfem::real_t lvoff() const override { return 2.5; }
  inline const mfem::real_t uvoff() const override { return 4.2; }
};
