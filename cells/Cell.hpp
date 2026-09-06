#pragma once

#include "mfem.hpp"

class Cell
{
public:
  /// Negative electrode ///
  // Initial concentration, negative electrode [mol/m^3]
  virtual inline const mfem::real_t cn0() const = 0;
  // Maximum concentration, negative electrode [mol/m^3]
  virtual inline const mfem::real_t cnmax() const = 0;
  // Negative particle radius [m]
  virtual inline const mfem::real_t rn() const = 0;
  // Negative electrode thickness [m]
  virtual inline const mfem::real_t ln() const = 0;
  // Negative electrode electrolyte volume fraction
  virtual inline const mfem::real_t eps_n() const = 0;
  // Negative electrode filler volume fraction
  virtual inline const mfem::real_t eps_n_fi() const = 0;
  // Negative electrode solid volume fraction
  virtual inline const mfem::real_t eps_n_s() const { return 1 - eps_n() - eps_n_fi(); }
  // Negative electrode reaction rate [A m^{2.5}/mol^{1.5}]
  virtual inline const mfem::real_t kn_dim() const = 0;
  // Diffusion coefficient, negative electrode [m^2/s]
  virtual inline const mfem::real_t Dn() const = 0;
  // Negative electrode area [m^{-1}]
  virtual inline const mfem::real_t An() const { return 3 * eps_n_s() / rn(); }
  // Negative electrode conductivity [S/m]
  virtual inline const mfem::real_t sig_n() const = 0;
  // Negative electrode Bruggeman coefficient
  virtual inline const mfem::real_t brugg_n() const = 0;
  // Negative electrode open circuit potential [V]
  virtual inline const mfem::real_t Un(mfem::real_t) const = 0;

  /// Positive electrode ///
  // Initial concentration, positive electrode [mol/m^3]
  virtual inline const mfem::real_t cp0() const = 0;
  // Maximum concentration, positive electrode [mol/m^3]
  virtual inline const mfem::real_t cpmax() const = 0;
  // Positive particle radius [m]
  virtual inline const mfem::real_t rp() const = 0;
  // Positive electrode thickness [m]
  virtual inline const mfem::real_t lp() const = 0;
  // Positive electrode electrolyte volume fraction
  virtual inline const mfem::real_t eps_p() const = 0;
  // Positive electrode filler volume fraction
  virtual inline const mfem::real_t eps_p_fi() const = 0;
  // Positive electrode solid volume fraction
  virtual inline const mfem::real_t eps_p_s() const { return 1 - eps_p() - eps_p_fi(); }
  // Positive electrode reaction rate [A m^{2.5}/mol^{1.5}]
  virtual inline const mfem::real_t kp_dim() const = 0;
  // Diffusion coefficient, positive electrode [m^2/s]
  virtual inline const mfem::real_t Dp() const = 0;
  // Positive electrode area [m^2]
  virtual inline const mfem::real_t Ap() const { return 3 * eps_p_s() / rp(); }
  // Positive electrode conductivity [S/m]
  virtual inline const mfem::real_t sig_p() const = 0;
  // Positive electrode Bruggeman coefficient
  virtual inline const mfem::real_t brugg_p() const = 0;
  // Positive electrode open circuit potential [V]
  virtual inline const mfem::real_t Up(mfem::real_t) const = 0;

  /// Separator ///
  // Separator thickness [m]
  virtual inline const mfem::real_t ls() const = 0;
  // Separator electrolyte volume fraction
  virtual inline const mfem::real_t eps_s() const = 0;
  // Separator Bruggeman coefficient
  virtual inline const mfem::real_t brugg_s() const = 0;

  /// Electrolyte ///
  // Transference number of Li ions
  virtual inline const mfem::real_t tplus() const = 0;
  // Initial electrolyte concentration [mol/(m^3)]
  virtual inline const mfem::real_t ce0() const = 0;
  // Diffusivity of Li ions in the electrolyte [m^2/s]
  virtual inline const mfem::real_t De(mfem::real_t) const = 0;
  // Electrolyte conductivity [S/m]
  virtual inline const mfem::real_t kappa(mfem::real_t) const = 0;

  /// Geometric parameters ///
  // Number of layers (one layer for wound cells)
  virtual inline const unsigned nlayer() const = 0;
  // Layer length [m]
  virtual inline const mfem::real_t llayer() const = 0;
  // Layer width [m]
  virtual inline const mfem::real_t wlayer() const = 0;
  // Total area [m^2]
  virtual inline const mfem::real_t cell_area() const { return llayer() * wlayer() * nlayer(); }

  // Or I1C in Jubat. Reference total current [A]
  virtual inline const mfem::real_t I_typ() const = 0;

  // Lower voltage cut-off [V]
  virtual inline const mfem::real_t lvoff() const = 0;
  // Upper voltage cut-off [V]
  virtual inline const mfem::real_t uvoff() const = 0;
};
