//======================================================================================
// Author Jens Chluba Sept/Oct 2010
// purpose: compute the 2s-1s two-photon profile and integral over it
//======================================================================================

#ifndef TWO_PHOTON_PROFILES_H
#define TWO_PHOTON_PROFILES_H

#include <string>

double phi_2s1s_em_spectrum(double y);

void phi_2s1s_function(double nu, double Tg, double nu0, double &phi_em, double &phi_abs);

void compute_stimulated_2s1s_integral(double Tgmin, double Tgmax, int nz, string filename);

#endif
