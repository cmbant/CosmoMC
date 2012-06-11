//======================================================================================
// Author Jens Chluba (July 2010)
//======================================================================================

#ifndef SOLVE_PDES_H
#define SOLVE_PDES_H

#include "Cosmos.h"
#include "Atom.h"

int compute_DPesc_with_diffusion_equation_effective(vector<double> &DF_vec_z, 
                                                    vector<vector<double> > &DF_2_gamma_vec, 
                                                    vector<vector<double> > &DF_Raman_vec, 
                                                    vector<double> &DI1_2s_vec, 
                                                    double zs, double ze, int nS_effective, 
                                                    int nmax_2g_corrs, int nmax_R_corrs,
                                                    Cosmos &cos, Gas_of_Atoms &HIA, 
                                                    vector<vector<double> > HI_Solution, 
                                                    int it_num);

#endif
