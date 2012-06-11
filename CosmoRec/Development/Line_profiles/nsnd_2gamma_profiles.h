//======================================================================================
// Author Jens Chluba Sept/Oct 2010
// purpose: compute the first few two-photon-profiles
// comment: a bit clunky but it works
//======================================================================================

#ifndef NSND_2GAMMA_PROFILES_H
#define NSND_2GAMMA_PROFILES_H

#include <string>

namespace nsnd_2gamma_profiles 
{
    void test_nsnd_2gamma_stuff();  
    void init_nsnd_2gamma_profiles(int nimax=5);

    double sigma_ns_1s_2gamma_ratio(int n, double y);
    double sigma_nd_1s_2gamma_ratio(int n, double y);

    //==================================================
    // 2s-1s two-photon profile functions
    //==================================================
    double sigma_2s_1s_non_res(double y);
    double sigma_2s_1s_2gamma(double y);
    void dump_2s_1s_2gamma_profile(string fname);
    
    //==================================================
    // 3s-1s two-photon profile functions
    //==================================================
    double sigma_3s_1s_non_res(double y);
    double sigma_3s_1s_res(double y);
    double sigma_3s_1s_res(double y, int choice);   
    double sigma_3s_1s_poles(double y);
    //
    double sigma_3s_1s_2gamma(double y);
    double sigma_3s_1s_2gamma_ratio(double y);
    void dump_3s_1s_2gamma_profile(string fname);

    //==================================================
    // 3d-1s two-photon profile functions
    //==================================================
    double sigma_3d_1s_non_res(double y);
    double sigma_3d_1s_res(double y);
    double sigma_3d_1s_res(double y, int choice);   
    double sigma_3d_1s_poles(double y);
    //
    double sigma_3d_1s_2gamma(double y);
    double sigma_3d_1s_2gamma_ratio(double y);
    void dump_3d_1s_2gamma_profile(string fname);

    //==================================================
    // 4s-1s and 4d-1s two-photon profile functions
    //==================================================
    double sigma_4s_1s_2gamma(double y);
    double sigma_4s_1s_2gamma_ratio(double y);
    void dump_4s_1s_2gamma_profile(string fname);
    //
    double sigma_4d_1s_2gamma(double y);
    double sigma_4d_1s_2gamma_ratio(double y);
    void dump_4d_1s_2gamma_profile(string fname);
    
    //==================================================
    // 5s-1s and 5d-1s two-photon profile functions
    //==================================================
    double sigma_5s_1s_2gamma(double y);
    double sigma_5s_1s_2gamma_ratio(double y);
    void dump_5s_1s_2gamma_profile(string fname);
    //
    double sigma_5d_1s_2gamma(double y);
    double sigma_5d_1s_2gamma_ratio(double y);
    void dump_5d_1s_2gamma_profile(string fname);

    //==================================================
    // 6s-1s and 6d-1s two-photon profile functions
    //==================================================
    double sigma_6s_1s_2gamma(double y);
    double sigma_6s_1s_2gamma_ratio(double y);
    void dump_6s_1s_2gamma_profile(string fname);
    //
    double sigma_6d_1s_2gamma(double y);
    double sigma_6d_1s_2gamma_ratio(double y);
    void dump_6d_1s_2gamma_profile(string fname);

    //==================================================
    // 7s-1s and 7d-1s two-photon profile functions
    //==================================================
    double sigma_7s_1s_2gamma(double y);
    double sigma_7s_1s_2gamma_ratio(double y);
    void dump_7s_1s_2gamma_profile(string fname);
    //
    double sigma_7d_1s_2gamma(double y);
    double sigma_7d_1s_2gamma_ratio(double y);
    void dump_7d_1s_2gamma_profile(string fname);

    //==================================================
    // 8s-1s and 8d-1s two-photon profile functions
    //==================================================
    double sigma_8s_1s_2gamma(double y);
    double sigma_8s_1s_2gamma_ratio(double y);
    void dump_8s_1s_2gamma_profile(string fname);
    //
    double sigma_8d_1s_2gamma(double y);
    double sigma_8d_1s_2gamma_ratio(double y);
    void dump_8d_1s_2gamma_profile(string fname);
}

#endif
