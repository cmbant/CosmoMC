//==========================================================================================
// Author Jens Chluba Aug/Sept 2010
// purpose: contains precomputed transitions rates for ns & nd-states of hydrogen
//==========================================================================================
#ifndef HI_TRANSITION_DATA_H
#define HI_TRANSITION_DATA_H

namespace HI_Transition_Data 
{
    //======================================================================================
    // Gamma_n=A_tot/(4 pi) for the first ten p-states (vacuum)
    //======================================================================================
    double Get_Gamma_np(int n);

    //======================================================================================
    // A_np1s for the first ten p-states (vacuum)
    //======================================================================================
    double Get_A_np1s(int n);

    //======================================================================================
    // A_npks for the first ten p-states (vacuum); n<k : Aksnp; n>k : Anpks*3;   k==n : 0
    //======================================================================================
    double Get_A_np2s(int n);
    double Get_A_np3s(int n);
    double Get_A_np4s(int n);
    double Get_A_np5s(int n);
    double Get_A_np6s(int n);
    double Get_A_np7s(int n);
    double Get_A_np8s(int n);
    double Get_A_npks(int k, int n);

    //======================================================================================
    // A_npkd for the first ten p-states (vacuum); n<k : Akdnp; n>k : Anpkd*3/5; k==n : 0
    //======================================================================================
    double Get_A_np3d(int n);
    double Get_A_np4d(int n);
    double Get_A_np5d(int n);
    double Get_A_np6d(int n);
    double Get_A_np7d(int n);
    double Get_A_np8d(int n);
    double Get_A_npkd(int k, int n);
}   

#endif
