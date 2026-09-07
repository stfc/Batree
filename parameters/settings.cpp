#include "mfem.hpp"
#include "cells/LGM50.hpp"
#include "cells/Enertech.hpp"

namespace settings
{
bool SPM = false;
bool SPMe = false;
bool P2D = false;

unsigned NNE = 0;
unsigned NSEP = 0;
unsigned NPE = 0;
unsigned NX = 0;
unsigned NR = 10;

unsigned NNEPAR = 0;
unsigned NPEPAR = 0;
unsigned NPAR = 0;

unsigned NMACROP = 2;
unsigned NMACROC = 1;
unsigned NMACRO = NMACROP + NMACROC;
unsigned NEQS = 0;

Cell * CELL = nullptr;
LGM50 LGM50_CELL = LGM50();
Enertech Enertech_CELL = Enertech();

mfem::real_t F = 96485.33289;
mfem::real_t R = 8.314;
mfem::real_t T_ref = 298.15;

mfem::real_t t0 = 1.0;
mfem::real_t r0 = 1e-6;
mfem::real_t L = 0;

mfem::real_t LNE = 0;
mfem::real_t LSEP = 0;
mfem::real_t LPE = 0;

mfem::real_t a0 = 0;

mfem::real_t tn = 0;
mfem::real_t tp = 0;

mfem::real_t te = 0;

mfem::real_t te_scale = 0;

mfem::real_t Dn_scale = 0;
mfem::real_t Dp_scale = 0;
mfem::real_t De_scale = 0;

mfem::real_t BPE = 0;
mfem::real_t BNE = 0;
mfem::real_t BSEP = 0;

mfem::real_t j_scale = 0;

mfem::real_t kn_scale = 0;
mfem::real_t kp_scale = 0;

mfem::real_t phi_scale = T_ref * R / F;

mfem::real_t tn_scale = 0;
mfem::real_t tp_scale = 0;

mfem::real_t ce_scale = 0;

mfem::real_t sig_scale = 0;
mfem::real_t kappa_scale = 0;

mfem::real_t DN = 0;
mfem::real_t DP = 0;

mfem::real_t AN = 0;
mfem::real_t AP = 0;

mfem::real_t KN = 0;
mfem::real_t KP = 0;

mfem::real_t CN0 = 0;
mfem::real_t CP0 = 0;

mfem::real_t RN = 0;
mfem::real_t RP = 0;

mfem::real_t SIGP = 0;
mfem::real_t SIGN = 0;

mfem::real_t CE0 = 0;
mfem::real_t I = 1.;
mfem::real_t T = 1.0;

mfem::real_t EPS_P = 0;
mfem::real_t EPS_N = 0;
mfem::real_t EPS_S = 0;

mfem::real_t TPLUS = 0;

void
init_settings(std::string m, std::string c, mfem::real_t c_rate, int order)
{
  std::transform(m.begin(), m.end(), m.begin(), [](unsigned char c) { return std::tolower(c); });

  if (m == "spm")
    SPM = true;
  else if (m == "spme")
    SPMe = true;
  else if (m == "p2d" || m == "dfn")
    P2D = true;
  else
    mfem::mfem_error("Unrecognised model.");

  std::transform(c.begin(), c.end(), c.begin(), [](unsigned char c) { return std::tolower(c); });

  if (c == "lgm50" || c == "chen2020")
    CELL = &LGM50_CELL;
  else if (c == "enertech" || c == "ai2020")
    CELL = &Enertech_CELL;
  else
    mfem::mfem_error("Unrecognised cell name.");

  if (SPM)
    NNE = NSEP = NPE = 0;
  else if (SPMe || P2D)
    NNE = NSEP = NPE = 10;
  NX = NNE + NSEP + NPE;

  if (SPM || SPMe)
    NNEPAR = NPEPAR = 1;
  else if (P2D)
  {
    NNEPAR = NNE * order + 1;
    NPEPAR = NPE * order + 1;
  }
  NPAR = NNEPAR + NPEPAR;
  NEQS = NMACRO + NPAR;

  L = CELL->ln() + CELL->ls() + CELL->lp();

  LNE = CELL->ln() / L;
  LSEP = CELL->ls() / L;
  LPE = CELL->lp() / L;

  a0 = 1.0 / r0;

  tn = F * CELL->cnmax() * CELL->cell_area() * L / CELL->I_typ();
  tp = F * CELL->cpmax() * CELL->cell_area() * L / CELL->I_typ();

  te = F * CELL->ce0() * CELL->cell_area() * L / CELL->I_typ();

  te_scale = te / t0;

  Dn_scale = r0 * r0 / t0;
  Dp_scale = r0 * r0 / t0;
  De_scale = L * L / te;

  BPE = pow(CELL->eps_p(), CELL->brugg_p());
  BNE = pow(CELL->eps_n(), CELL->brugg_n());
  BSEP = pow(CELL->eps_s(), CELL->brugg_s());

  j_scale = CELL->I_typ() / a0 / L / CELL->cell_area();

  kn_scale = j_scale / CELL->cnmax() / sqrt(CELL->ce0());
  kp_scale = j_scale / CELL->cpmax() / sqrt(CELL->ce0());

  tn_scale = tn / t0;
  tp_scale = tp / t0;

  ce_scale = CELL->ce0();

  sig_scale = L * CELL->I_typ() / (CELL->cell_area() * phi_scale);
  kappa_scale = L * CELL->I_typ() / (CELL->cell_area() * phi_scale);

  DN = CELL->Dn() / Dn_scale;
  DP = CELL->Dp() / Dp_scale;

  AN = CELL->An() / a0;
  AP = CELL->Ap() / a0;

  KN = CELL->kn_dim() / kn_scale;
  KP = CELL->kp_dim() / kp_scale;

  CN0 = CELL->cn0() / CELL->cnmax();
  CP0 = CELL->cp0() / CELL->cpmax();

  RN = CELL->rn() / r0;
  RP = CELL->rp() / r0;

  SIGP = CELL->sig_p() / sig_scale;
  SIGN = CELL->sig_n() / sig_scale;

  CE0 = CELL->ce0() / ce_scale;
  I = c_rate;

  EPS_P = CELL->eps_p();
  EPS_N = CELL->eps_n();
  EPS_S = CELL->eps_s();

  TPLUS = CELL->tplus();
}
}
