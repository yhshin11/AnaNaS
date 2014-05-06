#include "Analysis/utils/Constants.hh"

#include <cmath>

const float Constants::pi        = acos(-1);
const float Constants::twopi     = 2.*Constants::pi;
const float Constants::radToDeg  = 180./Constants::pi;

//
// From PDG, C. Amsler et al., Physics Letters B667, 1 (2008)
//      http://pdg.lbl.gov
//

const float Constants::Z0_m         = 91.188;     // GeV/c^2
const float Constants::Z0_Gamma     =  2.4952;    // GeV/c^2
const float Constants::Z0_br_had    = 69.910e-02;
const float Constants::Z0_br_bb     = 15.120e-02;
const float Constants::Z0_br_cc     = 12.030e-02;
const float Constants::Z0_br_inv    = 20.000e-02;
const float Constants::Z0_br_ll     =  3*3.366e-02;
const float Constants::Z0_br_ee     =  3.363e-02;
const float Constants::Z0_br_mumu   =  3.366e-02;
const float Constants::Z0_br_tautau =  3.370e-02;
const float Constants::Z0_br_nunu   =  Constants::Z0_br_inv;

const float Constants::W_m         = 80.398;     // GeV/c^2
const float Constants::W_Gamma     =  2.141;     // GeV/c^2
const float Constants::W_br_had    = 67.600e-02;
const float Constants::W_br_csbar  = 31.000e-02; // +/- 13!
const float Constants::W_br_lnu    = 3*10.800e-02;
const float Constants::W_br_enu    = 10.750e-02;
const float Constants::W_br_munu   = 10.570e-02;
const float Constants::W_br_taunu  = 11.350e-02;

const float Constants::tau_m         = 1.77684;    // GeV/c^2
const float Constants::tau_br_1p     = 85.360e-02;
const float Constants::tau_br_3p     = 14.560e-02;
const float Constants::tau_br_enunu  = 17.850e-02;
const float Constants::tau_br_mununu = 17.360e-02;
const float Constants::tau_br_hnu    = 11.600e-02;

