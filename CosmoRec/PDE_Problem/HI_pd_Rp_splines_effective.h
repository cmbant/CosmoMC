//===========================================================================================================
// Author Jens Chluba May 2009
// comment: HI populations have to be pre-loaded; Also the effective rates have to be initialized
//===========================================================================================================

#ifndef HI_PD_RP_SPLINES_EFFECTIVE_H
#define HI_PD_RP_SPLINES_EFFECTIVE_H

#include <string>
#include "Atom.h"
#include "Cosmos.h"

//===========================================================================================================
// For precalulated death probabilities and emission rates
//===========================================================================================================
void set_up_splines_for_HI_pd_Rp_effective(double zend, double zstart, 
                                           int nS_effective, Cosmos &cosmos, Gas_of_Atoms &HA);

void clear_HI_pd_Rp_effective_memory();
void Set_HI_pd_Rp_splines_effective_verbosity(int v);

//===========================================================================================================
// pd-functions
//===========================================================================================================
double calc_HI_pd_ns_splines_effective(double z, int n);
double calc_HI_pd_np_splines_effective(double z, int n);
double calc_HI_pd_nd_splines_effective(double z, int n);
double calc_HI_pd_nl_splines_effective(double z, int n, int l);
double calc_HI_pd_i_splines_effective(double z, int i);

//===========================================================================================================
// Rp-functions
//===========================================================================================================
double calc_HI_Rp_Rm_ns_splines_effective(double z, int n);
double calc_HI_Rp_Rm_np_splines_effective(double z, int n);
double calc_HI_Rp_Rm_nd_splines_effective(double z, int n);
double calc_HI_Rp_Rm_nl_splines_effective(double z, int n, int l);
double calc_HI_Rp_Rm_i_splines_effective(double z, int i);

//===========================================================================================================
// Dnem-functions
//===========================================================================================================
double calc_HI_Dnem_ns_splines_effective(double z, int n);
double calc_HI_Dnem_np_splines_effective(double z, int n);
double calc_HI_Dnem_nd_splines_effective(double z, int n);
double calc_HI_Dnem_nl_splines_effective(double z, int n, int l);
double calc_HI_Dnem_i_splines_effective(double z, int i);

#endif
