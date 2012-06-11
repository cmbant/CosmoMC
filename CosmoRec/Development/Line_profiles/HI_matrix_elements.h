//==========================================================================================
// Author Jens Chluba Aug/Sept 2010
// purpose: explicit expression for matrix elements of the first few ns & nd-states in HI
//==========================================================================================

#ifndef HI_MATRIX_ELEMENTS_H
#define HI_MATRIX_ELEMENTS_H

namespace HI_matrix_elements 
{
    double Rnsnp(int n);
    double Rndnp(int n);
    
    //======================================================================================
    // bound-bound matrix elements
    //======================================================================================
    double R1snp(int n);
    double R2snp(int n);
    //
    double R3snp(int n);
    double R3dnp(int n);
    //
    double R4snp(int n);
    double R4dnp(int n);
    //
    double R5snp(int n);
    double R5dnp(int n);
    //
    double R6snp(int n);
    double R6dnp(int n);
    //
    double R7snp(int n);
    double R7dnp(int n);
    //
    double R8snp(int n);
    double R8dnp(int n);
    
    //======================================================================================
    // bound-free matrix elements
    //======================================================================================
    double C1s(double x);
    double C2s(double x);
    //
    double C3s(double x);
    double C3d(double x);
    //
    double C4s(double x);
    double C4d(double x);
    //
    double C5s(double x);
    double C5d(double x);   
    //
    double C6s(double x);
    double C6d(double x);   
    //
    double C7s(double x);
    double C7d(double x);   
    //
    double C8s(double x);
    double C8d(double x);   
}

#endif
