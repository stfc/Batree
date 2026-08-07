#include "cells/LGM50.hpp"
#include "mfem.hpp"

using namespace mfem;

namespace settings
{
enum PotentialBlock : int
{
  EPP, // Electrolyte Potential
  SPP, // Solid Potential (Postive and Negative Electrode region)
};

enum ConcentrationBlock : int
{
  ECC, // Electrolyte Concentration
  SCC  // Solid Concentration
};

enum Block : int
{
  P = 0, // Potential
  C = 2  // Concentration
};

enum XBlock : int
{
  EP = P + EPP, // Electrolyte Potential
  SP = P + SPP, // Solid Potential
  EC = C + ECC, // Electrolyte Concentration
  SC = C + SCC  // Solid Concentration
};

enum Region : int
{
  E,   // Electrolytea
  NE,  // Negative Electrode
  SEP, // Separator
  PE,  // Positive Electrode
  UNKNOWN
};

extern const bool SPM;  // Whether running SPM model
extern const bool SPMe; // Whether running SPMe model
extern const bool P2D;  // Whether running P2D model

extern const unsigned NNE;  // Number of elements in the Negative Electrode
extern const unsigned NSEP; // Number of elements in the Separator
extern const unsigned NPE;  // Number of elements in the Positive Electrode
extern const unsigned NX;   // Number of elements in the X-dimension (i.e Electrolye) (sum of above)
extern const unsigned NR;   // Number of elements in the R-dimension (i.e Particle)

extern const unsigned NNEPAR; // Number of Negative Electrode PARticle
extern const unsigned NPEPAR; // Number of Positive Electrode PARticle
extern const unsigned NPAR;   // Total number of  PARticles (Sum of the above two)

extern const unsigned NMACROP; // Number of Macro Potential equations //later change to NMACROPEQS!
extern const unsigned NMACROC; // Number of Macro Concentration equations
extern const unsigned NMACRO;  // Total Number of Macro Equations (Sum of the above two)
extern const unsigned NEQS;    // Total Number of Equations in the System

extern Cell * CELL;            // Chosen cell
extern const LGM50 LGM50_CELL; // LGM50 cell

// Dimensional constants
extern const real_t F;     // Faraday constant, C/mol
extern const real_t R;     // Universal gas constant, J/(mol*K)
extern const real_t T_ref; // Reference temperature, K

// Scalings
extern const real_t t0; // Time scale.
extern const real_t r0; // Length scale (particle)
extern const real_t L;  // Length scale (cell)

extern const real_t LNE;  // Length of Negative Electrode
extern const real_t LSEP; // Length of Separator
extern const real_t LPE;  // Length of Positive Electrode

extern const real_t a0;

extern const real_t tn; // Negative particle time scale.
extern const real_t tp; // Positive particle time scale.

extern const real_t te; // Electrolyte "particle" time scale.

extern const real_t te_scale;

extern const real_t Dn_scale; // Negative particle diffusion coefficient scale.  Units of m^2/s.
extern const real_t Dp_scale; // Positive particle diffusion coefficient scale.  Units of m^2/s.

extern const real_t De_scale;

// Transport efficiency (inverse MacMullin number). This is B(x) in Planella, and is absorbed into
// the definition of kappa_ne/kappa_pe/kappa_sp in JuBat.
extern const real_t BPE;
extern const real_t BNE;
extern const real_t BSEP;

extern const real_t j_scale;

extern const real_t kn_scale;
extern const real_t kp_scale;

extern const real_t phi_scale; // potential scale

// For scaling between particle time scale and cell time scale.
// Required for scaling the flux j in the SolidConcentration equation.
extern const real_t tn_scale; // time scale of negative electrode
extern const real_t tp_scale; // time scale of positive electrode

extern const real_t ce_scale;

extern const real_t sig_scale;   // Electrode conductivity scale.
extern const real_t kappa_scale; // Electrolyte conductivity scale.

// Scaled parameters
extern const real_t DN; // scaled diffusion coefficient of each Negative particle
extern const real_t DP; // scaled diffusion coefficient of each Positive particle

extern const real_t AN; // scaled surface Area of each Negative particle // later SAN!
extern const real_t AP; // scaled surface Area of each Positive particle

extern const real_t KN; // scaled reaction rate of each Negative particle
extern const real_t KP; // scaled reaction rate of each Positive particle

extern const real_t CN0; // scaled initial Concentration of Negative particle
extern const real_t CP0; // scaled initial Concentration of Positive particle

extern const real_t RN; // scaled Radius of Negative particle
extern const real_t RP; // scaled Radius of Positive particle

extern const real_t SIGP; // Scaled positive electrode conductivity.
extern const real_t SIGN; // Scaled negative electrode conductivity.

// Extras to be properly defined later.
extern const real_t CE0; // scaled initial Concentration of Electrolyte
extern const real_t I;   // scaled external current
extern const real_t T;   // scaled Temperature.

extern const real_t EPS_P;
extern const real_t EPS_N;
extern const real_t EPS_S;

extern const real_t TPLUS;

inline const real_t
UN(real_t ce)
{
  return CELL->Un(ce) / phi_scale;
}
inline const real_t
UP(real_t ce)
{
  return CELL->Up(ce) / phi_scale;
}
inline const real_t
DE(real_t ce)
{
  return CELL->De(ce * ce_scale) / De_scale;
}
inline const real_t
Kappa(real_t ce)
{
  return CELL->kappa(ce * ce_scale) / kappa_scale;
}

const real_t KS = Kappa(CE0); // / kappa_scale; // Scaled electrolyte conductivity.

void init_settings(std::string m, std::string c, int order);
}
