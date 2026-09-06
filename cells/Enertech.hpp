#pragma once

#include "Cell.hpp"
#include "_EnertechUnSpline.hpp"
#include "_EnertechUpSpline.hpp"

class Enertech : public Cell
{
public:
  inline const mfem::real_t cn0() const override { return 24108.; }
  inline const mfem::real_t cnmax() const override { return 29700; } // PyBAMM uses 28700
  inline const mfem::real_t rn() const override { return 5e-6; }
  inline const mfem::real_t ln() const override { return 7.65e-5; }
  inline const mfem::real_t eps_n() const override { return 0.33; }
  inline const mfem::real_t eps_n_fi() const override { return 0.06; }
  inline const mfem::real_t kn_dim() const override { return 9.648533289e-7; }
  inline const mfem::real_t Dn() const override { return 3.9e-14; }
  inline const mfem::real_t sig_n() const override { return 100.; }
  inline const mfem::real_t brugg_n() const override { return 2.914; }
  inline const mfem::real_t Un(mfem::real_t cs) const override
  {
    const unsigned n = NBREAKS_GRAPHITE;
    const mfem::real_t * b = BREAKS_GRAPHITE;

    // locate interval: breaks[i] <= cs < breaks[i+1]
    const unsigned idx = std::distance(b, std::upper_bound(b, b + n, cs)) - 1;
    MFEM_ASSERT(idx >= 0 && idx < n - 1, "Enertech::Un: cs out of range");

    // compute spline value
    const mfem::real_t dx = cs - b[idx];
    const mfem::real_t * c = COEFS_GRAPHITE[idx];
    return ((c[0] * dx + c[1]) * dx + c[2]) * dx + c[3];
  }

  inline const mfem::real_t cp0() const override { return 21725.; }
  inline const mfem::real_t cpmax() const override { return 49943; }
  inline const mfem::real_t rp() const override { return 3e-6; }
  inline const mfem::real_t lp() const override { return 6.8e-5; }
  inline const mfem::real_t eps_p() const override { return 0.32; }
  inline const mfem::real_t eps_p_fi() const override { return 0.06; }
  inline const mfem::real_t kp_dim() const override { return 9.648533289e-7; }
  inline const mfem::real_t Dp() const override { return 5.387e-15; }
  inline const mfem::real_t sig_p() const override { return 10.; }
  inline const mfem::real_t brugg_p() const override { return 1.83; }
  inline const mfem::real_t Up(mfem::real_t cs) const override
  {
    const unsigned n = NBREAKS_LICO2;
    const mfem::real_t * b = BREAKS_LICO2;

    // locate interval: breaks[i] <= cs < breaks[i+1]
    const unsigned idx = std::distance(b, std::upper_bound(b, b + n, cs)) - 1;
    MFEM_ASSERT(idx >= 0 && idx < n - 1, "Enertech::Up: cs out of range");

    // compute spline value
    const mfem::real_t dx = cs - b[idx];
    const mfem::real_t * c = COEFS_LICO2[idx];
    return ((c[0] * dx + c[1]) * dx + c[2]) * dx + c[3];
  }

  inline const mfem::real_t ls() const override { return 2.5e-5; }
  inline const mfem::real_t eps_s() const override { return 0.5; }
  inline const mfem::real_t brugg_s() const override { return 1.5; }

  inline const mfem::real_t tplus() const override { return 0.38; }
  inline const mfem::real_t ce0() const override { return 1000.0; }
  inline const mfem::real_t De(mfem::real_t ce) const override
  {
    const mfem::real_t T = 298.15;
    return pow(10.0, (-4.43 - 54.0 / (T - 229.0 - 5e-3 * ce) - 0.22e-3 * ce));
  }
  inline const mfem::real_t kappa(mfem::real_t ce) const override
  {
    const mfem::real_t T = 298.15;
    const mfem::real_t inner = (-10.5 + 0.668e-3 * ce + 0.494e-6 * ce * ce) +
                               (0.074 - 1.78e-5 * ce - 8.86e-10 * ce * ce) * T +
                               (-6.96e-5 + 2.8e-8 * ce) * T * T;
    return 1e-4 * ce * inner * inner;
  }

  inline const unsigned nlayer() const override { return 34; }
  inline const mfem::real_t llayer() const override { return 0.051; }
  inline const mfem::real_t wlayer() const override { return 0.047; }

  inline const mfem::real_t I_typ() const override { return 2.28; }

  inline const mfem::real_t lvoff() const override { return 3.0; }
  inline const mfem::real_t uvoff() const override { return 4.2; }
};
