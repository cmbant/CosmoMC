//======================================================================================
// Author Jens Chluba Sept/Oct 2010
// purpose: compute the first few two-photon-profiles
// comment: a bit clunky but it works
//======================================================================================

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h>

#include "physical_consts.h"
#include "routines.h"
#include "Patterson.h"

#include "HI_matrix_elements.h"
#include "HI_Transition_Data.h"
#include "nsnd_2gamma_profiles.h"

using namespace std;
using namespace HI_matrix_elements;
using namespace HI_Transition_Data;

//======================================================================================
// Atomic Data for HI atom
//======================================================================================
const double nsnd_2gamma_profiles_C_sd   =9.0/1024.0*pow(const_alpha, 6)*const_cl
                                          *const_Ry_inf_icm/(1.0+const_me_mp);

const double nsnd_2gamma_profiles_nu_1sc =const_EH_inf_Hz/(1.0+const_me_mp);

//======================================================================================
// variables
//======================================================================================
const int nsnd_2gamma_profiles_nmax=500;
const int nsnd_2gamma_profiles_xpts=2*256;
const double nsnd_2gamma_profiles_xmin=1.0e-4;

struct nsnd_2gamma_profiles_Data
{
    int xnpt;
    int memindex;
    double xmax;
    vector<double> kappa_res;
    //
    vector<double> f_re;
    vector<double> f_im;
    vector<double> yr;
};

nsnd_2gamma_profiles_Data nsnd_2gamma_profiles_Data_2s_1s;
//
nsnd_2gamma_profiles_Data nsnd_2gamma_profiles_Data_3s_1s;
nsnd_2gamma_profiles_Data nsnd_2gamma_profiles_Data_3d_1s;
//
nsnd_2gamma_profiles_Data nsnd_2gamma_profiles_Data_4s_1s;
nsnd_2gamma_profiles_Data nsnd_2gamma_profiles_Data_4d_1s;
//
nsnd_2gamma_profiles_Data nsnd_2gamma_profiles_Data_5s_1s;
nsnd_2gamma_profiles_Data nsnd_2gamma_profiles_Data_5d_1s;
//
nsnd_2gamma_profiles_Data nsnd_2gamma_profiles_Data_6s_1s;
nsnd_2gamma_profiles_Data nsnd_2gamma_profiles_Data_6d_1s;
//
nsnd_2gamma_profiles_Data nsnd_2gamma_profiles_Data_7s_1s;
nsnd_2gamma_profiles_Data nsnd_2gamma_profiles_Data_7d_1s;
//
nsnd_2gamma_profiles_Data nsnd_2gamma_profiles_Data_8s_1s;
nsnd_2gamma_profiles_Data nsnd_2gamma_profiles_Data_8d_1s;

//======================================================================================
// local functions
//======================================================================================
namespace nsnd_2gamma_profiles_local 
{
    double nuij(int n, int np)
    { return nsnd_2gamma_profiles_nu_1sc*(pow(1.0*np,-2) - pow(1.0*n,-2)); } 
    
    //======================================================================================
    // resonance frequencies (n<ni)
    //======================================================================================
    double y_res(int ni, int n){ return -(1.0/ni/ni-1.0/n/n)/(1.0-1.0/ni/ni); }
    
    //======================================================================================
    // energy factor
    //======================================================================================
    double fn_nsnd_2gamma(int ni, int n, int nf, double y)
    { 
        return 1.0/((pow(1.0*ni,-2) -pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y) 
             + 1.0/((pow(1.0*nf,-2) -pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y);
    }

    double fn_nsnd_2gamma_cont(int ni, double x, int nf, double y)
    { 
        return 1.0/(( x*x + pow(1.0*ni,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y) 
             + 1.0/(( x*x + pow(1.0*nf,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y);
    }
    
    //======================================================================================
    // energy factors for the resonances
    //======================================================================================
    double Lorentzian(double a, double b){ return a/(a*a+b*b); }
    
    double fn_nsnd_2gamma_r(int ni, int n, int nf, double y)
    {
        double dy1=(pow(1.0*ni,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y;
        double dy2=(pow(1.0*nf,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y;
        double nuscale=nuij(ni, nf);
        
        //==================================================================================
        // comment: the '+' sign is due to change of phase when using 
        // f= 1/(yp+y-i*d) + 1/(ym-y-i*d) instead of f= 1/(yp+y-i*d) - 1/(y-ym-i*d) as in 
        // CS 2009 paper
        //==================================================================================
        return Lorentzian(dy1, Get_Gamma_np(n)/nuscale) + Lorentzian(dy2, Get_Gamma_np(n)/nuscale);
    }

    double fn_nsnd_2gamma_i(int ni, int n, int nf, double y)
    {
        double dy1=(pow(1.0*ni,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y;
        double dy2=(pow(1.0*nf,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y;
        double nuscale=nuij(ni, nf);
        
        return Lorentzian(Get_Gamma_np(n)/nuscale, dy1) - Lorentzian(Get_Gamma_np(n)/nuscale, dy2);
    }
        
    //======================================================================================
    // normalization factor
    //======================================================================================
    double G_ns(int n){ return nsnd_2gamma_profiles_C_sd    *pow(4.0/3.0*(1.0-1.0/n/n), 5); }
    double G_nd(int n){ return nsnd_2gamma_profiles_C_sd/2.5*pow(4.0/3.0*(1.0-1.0/n/n), 5); }
    
    //======================================================================================
    //======================================================================================

    //======================================================================================
    // ns1s, and nd1s integrals over free states
    //======================================================================================
    double (*Cj_x_ptr__)(double);
    
    double dIntC1sCj_dx(double logx, void *p)
    {
        double *y=(double *) p;
        double x=exp(logx);
        return x*C1s(x)*Cj_x_ptr__(x)*fn_nsnd_2gamma_cont(y[1], x, 1, y[0]);
    }

    double IntC1sC_nsd(double y, int ni, double (*Cj_x)(double))
    {
        if(y>=1.0 || y<=0.0) return 0.0;
        
        double r=0.0;
        double a=log(1.0e-8), b=log(0.5), epsrel=1.0e-8, epsabs=1.0e-10;
        double pars[]={y, ni};
        void *p=&pars;
        Cj_x_ptr__=Cj_x;
        
        r =Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dIntC1sCj_dx, p);
        //
        a=b; b=log(1.0); epsabs=max(epsabs, r*epsrel/2.0);
        r+=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dIntC1sCj_dx, p);
        //
        a=b; b=log(10.0); epsabs=max(epsabs, r*epsrel/2.0);
        r+=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dIntC1sCj_dx, p);
        //
        a=b; b=log(100.0); epsabs=max(epsabs, r*epsrel/2.0);
        r+=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dIntC1sCj_dx, p);
        //
        a=b; b=log(500.0); epsabs=max(epsabs, r*epsrel/2.0);
        r+=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dIntC1sCj_dx, p);
        
        return r;
    }

    //======================================================================================
    //======================================================================================

    //======================================================================================
    // Matrix element ns-->np
    //======================================================================================
    double Matrix_element_ns_np_1s(double y, int n)
    { return R1snp(n)*Rnsnp(n)*fn_nsnd_2gamma(n, n, 1, y); }

    //======================================================================================
    // Matrix element nd-->np
    //======================================================================================
    double Matrix_element_nd_np_1s(double y, int n)
    { return R1snp(n)*Rndnp(n)*fn_nsnd_2gamma(n, n, 1, y); }

    //======================================================================================
    // general Matrix element ns/d-->np
    //======================================================================================
    double Matrix_element_j_np_1s(int n, double y, int ni, double (*Rj_np)(int))
    { return R1snp(n)*Rj_np(n)*fn_nsnd_2gamma(ni, n, 1, y); }
    
    double Matrix_element_j_1s_total_non_res(double y, int ni, double (*Rj_np)(int), 
                                             double (*Cj)(double))
    {
        double r=IntC1sC_nsd(y, ni, Cj) + Matrix_element_ns_np_1s(y, ni);

        for(int nint=nsnd_2gamma_profiles_nmax; nint>ni; nint--) 
            r+=Matrix_element_j_np_1s(nint, y, ni, Rj_np);
        
        return r;
    }
    //======================================================================================
    //======================================================================================

    //======================================================================================
    // Matrix element Mn_2s_1s
    //======================================================================================
    double Matrix_element_2s_1s_total_non_res(double y)
    { return Matrix_element_j_1s_total_non_res(y, 2, R2snp, C2s); }
    
    //======================================================================================
    // Matrix element Mn_3s_1s && Mn_3d_1s
    //======================================================================================
    double Matrix_element_3s_1s_total_non_res(double y)
    { return Matrix_element_j_1s_total_non_res(y, 3, R3snp, C3s); }

    double Matrix_element_3d_1s_total_non_res(double y)
    { return Matrix_element_j_1s_total_non_res(y, 3, R3dnp, C3d); }
    
    //======================================================================================
    // Matrix element Mn_4s_1s && M_4d_1s
    //======================================================================================
    double Matrix_element_4s_1s_total_non_res(double y)
    { return Matrix_element_j_1s_total_non_res(y, 4, R4snp, C4s); }
    
    double Matrix_element_4d_1s_total_non_res(double y)
    { return Matrix_element_j_1s_total_non_res(y, 4, R4dnp, C4d); }
    
    //======================================================================================
    // Matrix element Mn_5s_1s && M_5d_1s
    //======================================================================================
    double Matrix_element_5s_1s_total_non_res(double y)
    { return Matrix_element_j_1s_total_non_res(y, 5, R5snp, C5s); }
    
    double Matrix_element_5d_1s_total_non_res(double y)
    { return Matrix_element_j_1s_total_non_res(y, 5, R5dnp, C5d); }
    
    //======================================================================================
    // Matrix element Mn_6s_1s && M_6d_1s
    //======================================================================================
    double Matrix_element_6s_1s_total_non_res(double y)
    { return Matrix_element_j_1s_total_non_res(y, 6, R6snp, C6s); }
    
    double Matrix_element_6d_1s_total_non_res(double y)
    { return Matrix_element_j_1s_total_non_res(y, 6, R6dnp, C6d); }
    
    //======================================================================================
    // Matrix element Mn_7s_1s && M_7d_1s
    //======================================================================================
    double Matrix_element_7s_1s_total_non_res(double y)
    { return Matrix_element_j_1s_total_non_res(y, 7, R7snp, C7s); }
    
    double Matrix_element_7d_1s_total_non_res(double y)
    { return Matrix_element_j_1s_total_non_res(y, 7, R7dnp, C7d); }
    
    //======================================================================================
    // Matrix element Mn_8s_1s && M_8d_1s
    //======================================================================================
    double Matrix_element_8s_1s_total_non_res(double y)
    { return Matrix_element_j_1s_total_non_res(y, 8, R8snp, C8s); }
    
    double Matrix_element_8d_1s_total_non_res(double y)
    { return Matrix_element_j_1s_total_non_res(y, 8, R8dnp, C8d); }
    
    //======================================================================================
    //======================================================================================
    
    //======================================================================================
    // Matrix element using the spline functions after setup
    //======================================================================================
    double compute_Mnr_2s_1s(double y)
    { return calc_spline_JC(min(nsnd_2gamma_profiles_Data_2s_1s.xmax, max(nsnd_2gamma_profiles_xmin, y)), 
                            nsnd_2gamma_profiles_Data_2s_1s.memindex)/y/(1.0-y); }
    
    double compute_Mnr_3s_1s(double y)
    { return calc_spline_JC(min(nsnd_2gamma_profiles_Data_3s_1s.xmax, max(nsnd_2gamma_profiles_xmin, y)), 
                            nsnd_2gamma_profiles_Data_3s_1s.memindex)/y/(1.0-y); }
    
    double compute_Mnr_3d_1s(double y)
    { return calc_spline_JC(min(nsnd_2gamma_profiles_Data_3d_1s.xmax, max(nsnd_2gamma_profiles_xmin, y)), 
                            nsnd_2gamma_profiles_Data_3d_1s.memindex)/y/(1.0-y); }

    double compute_Mnr_4s_1s(double y)
    { return calc_spline_JC(min(nsnd_2gamma_profiles_Data_4s_1s.xmax, max(nsnd_2gamma_profiles_xmin, y)), 
                            nsnd_2gamma_profiles_Data_4s_1s.memindex)/y/(1.0-y); }
    
    double compute_Mnr_4d_1s(double y)
    { return calc_spline_JC(min(nsnd_2gamma_profiles_Data_4d_1s.xmax, max(nsnd_2gamma_profiles_xmin, y)), 
                            nsnd_2gamma_profiles_Data_4d_1s.memindex)/y/(1.0-y); }

    double compute_Mnr_5s_1s(double y)
    { return calc_spline_JC(min(nsnd_2gamma_profiles_Data_5s_1s.xmax, max(nsnd_2gamma_profiles_xmin, y)), 
                            nsnd_2gamma_profiles_Data_5s_1s.memindex)/y/(1.0-y); }
    
    double compute_Mnr_5d_1s(double y)
    { return calc_spline_JC(min(nsnd_2gamma_profiles_Data_5d_1s.xmax, max(nsnd_2gamma_profiles_xmin, y)), 
                            nsnd_2gamma_profiles_Data_5d_1s.memindex)/y/(1.0-y); }

    double compute_Mnr_6s_1s(double y)
    { return calc_spline_JC(min(nsnd_2gamma_profiles_Data_6s_1s.xmax, max(nsnd_2gamma_profiles_xmin, y)), 
                            nsnd_2gamma_profiles_Data_6s_1s.memindex)/y/(1.0-y); }
    
    double compute_Mnr_6d_1s(double y)
    { return calc_spline_JC(min(nsnd_2gamma_profiles_Data_6d_1s.xmax, max(nsnd_2gamma_profiles_xmin, y)), 
                            nsnd_2gamma_profiles_Data_6d_1s.memindex)/y/(1.0-y); }
    
    double compute_Mnr_7s_1s(double y)
    { return calc_spline_JC(min(nsnd_2gamma_profiles_Data_7s_1s.xmax, max(nsnd_2gamma_profiles_xmin, y)), 
                            nsnd_2gamma_profiles_Data_7s_1s.memindex)/y/(1.0-y); }
    
    double compute_Mnr_7d_1s(double y)
    { return calc_spline_JC(min(nsnd_2gamma_profiles_Data_7d_1s.xmax, max(nsnd_2gamma_profiles_xmin, y)), 
                            nsnd_2gamma_profiles_Data_7d_1s.memindex)/y/(1.0-y); }
    
    double compute_Mnr_8s_1s(double y)
    { return calc_spline_JC(min(nsnd_2gamma_profiles_Data_8s_1s.xmax, max(nsnd_2gamma_profiles_xmin, y)), 
                            nsnd_2gamma_profiles_Data_8s_1s.memindex)/y/(1.0-y); }
    
    double compute_Mnr_8d_1s(double y)
    { return calc_spline_JC(min(nsnd_2gamma_profiles_Data_8d_1s.xmax, max(nsnd_2gamma_profiles_xmin, y)), 
                            nsnd_2gamma_profiles_Data_8d_1s.memindex)/y/(1.0-y); }
    
    //======================================================================================
    //======================================================================================
}

using namespace nsnd_2gamma_profiles_local;

namespace nsnd_2gamma_profiles 
{

    //======================================================================================
    // testing 
    //======================================================================================
    void test_nsnd_2gamma_stuff()
    {
        wait_f_r();
    }

    //======================================================================================
    // setup routines
    //======================================================================================
    void alloc_memory(int ni, nsnd_2gamma_profiles_Data &nsnd_2gamma_profiles_Data_i)
    {
        nsnd_2gamma_profiles_Data_i.kappa_res.clear();
        nsnd_2gamma_profiles_Data_i.kappa_res.resize(ni, 0.0);
        nsnd_2gamma_profiles_Data_i.yr.clear();
        nsnd_2gamma_profiles_Data_i.yr.       resize(ni, 0.0);

        nsnd_2gamma_profiles_Data_i.f_re.resize(ni);
        nsnd_2gamma_profiles_Data_i.f_im.resize(ni);
        
        return;
    }
    
    void init_splines(vector<double> &xarr, vector<double> &yarr, 
                      nsnd_2gamma_profiles_Data &nsnd_2gamma_profiles_Data_i, double (*M_nr_i)(double))
    {
        for(int k=0; k<nsnd_2gamma_profiles_xpts; k++) yarr[k]=xarr[k]*(1.0-xarr[k])*M_nr_i(xarr[k]);
        
        nsnd_2gamma_profiles_Data_i.xnpt=nsnd_2gamma_profiles_xpts;
        nsnd_2gamma_profiles_Data_i.xmax=xarr[nsnd_2gamma_profiles_xpts-1];
        nsnd_2gamma_profiles_Data_i.memindex=calc_spline_coeffies_JC(nsnd_2gamma_profiles_Data_i.xnpt, 
                                                                     &xarr[0], &yarr[0], 
                                                                     "nsnd_2gamma_profiles");
        
        return;
    }
    
    void init_mem_and_data(int ni, nsnd_2gamma_profiles_Data &nsnd_2gamma_profiles_Data_i, 
                           double (*Rj_np)(int), vector<double> &xarr, vector<double> &yarr, 
                           double (*M_nr_i)(double))
    {
        alloc_memory(ni, nsnd_2gamma_profiles_Data_i);
        
        for(int n=2; n<ni; n++) 
        {
            nsnd_2gamma_profiles_Data_i.kappa_res[n]=R1snp(n)*Rj_np(n);
            nsnd_2gamma_profiles_Data_i.yr[n]=y_res(ni, n);
        }
        
        init_splines(xarr, yarr, nsnd_2gamma_profiles_Data_i, M_nr_i);
        
        return;
    }
    
    void dump_mess(string mess)
    {
        cout << " dumping 2gamma profile for " << mess << endl;
        return;
    }
    
    //======================================================================================
    // setup non-resonant part of 2gamma-profile; this depends on chosen xmax
    //======================================================================================
    void init_nsnd_2gamma_profiles(int nimax)
    {       
        if(nimax<2) return;
        if(nimax>8){ cerr << "\n init_nsnd_2gamma_profiles:: two-photon profiles are only available up to nmax=8 " << endl; exit(0);}
        
        cout << "\n init_nsnd_2gamma_profiles:: Initializing ns-1s and nd-1s 2gamma-profiles up to nmax= " << nimax << endl;

        //======================================================================
        // setup splines for non-resonant part
        //======================================================================
        vector<double> xarr(nsnd_2gamma_profiles_xpts);
        vector<double> yarr(nsnd_2gamma_profiles_xpts);
        double xmax=1.0-nsnd_2gamma_profiles_xmin;
        //
        init_xarr(nsnd_2gamma_profiles_xmin, xmax, &xarr[0], nsnd_2gamma_profiles_xpts, 0, 0);

        //======================================================================
        // 2s
        //======================================================================
        init_splines(xarr, yarr, nsnd_2gamma_profiles_Data_2s_1s, Matrix_element_2s_1s_total_non_res);
        
        //======================================================================
        // compute all profiles splines and kappa_n
        //======================================================================
        if(nimax>2)
        {
            init_mem_and_data(3, nsnd_2gamma_profiles_Data_3s_1s, R3snp, xarr, yarr, Matrix_element_3s_1s_total_non_res);
            init_mem_and_data(3, nsnd_2gamma_profiles_Data_3d_1s, R3dnp, xarr, yarr, Matrix_element_3s_1s_total_non_res);
            
            if(nimax>3)
            {
                init_mem_and_data(4, nsnd_2gamma_profiles_Data_4s_1s, R4snp, xarr, yarr, Matrix_element_4s_1s_total_non_res);
                init_mem_and_data(4, nsnd_2gamma_profiles_Data_4d_1s, R4dnp, xarr, yarr, Matrix_element_4d_1s_total_non_res);
                
                if(nimax>4)
                {
                    init_mem_and_data(5, nsnd_2gamma_profiles_Data_5s_1s, R5snp, xarr, yarr, Matrix_element_5s_1s_total_non_res);
                    init_mem_and_data(5, nsnd_2gamma_profiles_Data_5d_1s, R5dnp, xarr, yarr, Matrix_element_5d_1s_total_non_res);

                    if(nimax>5)
                    {
                        init_mem_and_data(6, nsnd_2gamma_profiles_Data_6s_1s, R6snp, xarr, yarr, Matrix_element_6s_1s_total_non_res);
                        init_mem_and_data(6, nsnd_2gamma_profiles_Data_6d_1s, R6dnp, xarr, yarr, Matrix_element_6d_1s_total_non_res);

                        if(nimax>6)
                        {
                            init_mem_and_data(7, nsnd_2gamma_profiles_Data_7s_1s, R7snp, xarr, yarr, Matrix_element_7s_1s_total_non_res);
                            init_mem_and_data(7, nsnd_2gamma_profiles_Data_7d_1s, R7dnp, xarr, yarr, Matrix_element_7d_1s_total_non_res);

                            if(nimax>7)
                            {
                                init_mem_and_data(8, nsnd_2gamma_profiles_Data_8s_1s, R8snp, xarr, yarr, Matrix_element_8s_1s_total_non_res);
                                init_mem_and_data(8, nsnd_2gamma_profiles_Data_8d_1s, R8dnp, xarr, yarr, Matrix_element_8d_1s_total_non_res);
                            }
                        }
                    }
                }
            }
        }

        cout << " init_nsnd_2gamma_profiles:: done " << endl;

        return;
    }
    
    //======================================================================================
    //
    // access to different profiles
    //
    //======================================================================================
    double sigma_ns_1s_2gamma_ratio(int n, double y)
    {
        if(n==3) return sigma_3s_1s_2gamma_ratio(y);
        if(n==4) return sigma_4s_1s_2gamma_ratio(y);
        if(n==5) return sigma_5s_1s_2gamma_ratio(y);
        if(n==6) return sigma_6s_1s_2gamma_ratio(y);
        if(n==7) return sigma_7s_1s_2gamma_ratio(y);
        if(n==8) return sigma_8s_1s_2gamma_ratio(y);
        
        return 1.0;
    }
    
    double sigma_nd_1s_2gamma_ratio(int n, double y)
    {
        if(n==3) return sigma_3d_1s_2gamma_ratio(y);
        if(n==4) return sigma_4d_1s_2gamma_ratio(y);
        if(n==5) return sigma_5d_1s_2gamma_ratio(y);
        if(n==6) return sigma_6d_1s_2gamma_ratio(y);
        if(n==7) return sigma_7d_1s_2gamma_ratio(y);
        if(n==8) return sigma_8d_1s_2gamma_ratio(y);
        
        return 1.0;
    }   

    //======================================================================================
    //
    // different cross sections ns/d divided by Gnl(!!!)
    //
    //======================================================================================
    double sigma_nsd_1s_res(double y, int ni, int choice, double Mnr, nsnd_2gamma_profiles_Data &nsnd_2gamma_profiles_Data_nsd_1s)
    { 
        double r=0.0;

        // prepare fn tables
        for(int n=2; n<ni; n++) 
        {
            nsnd_2gamma_profiles_Data_nsd_1s.f_re[n]=fn_nsnd_2gamma_r(ni, n, 1, y);
            nsnd_2gamma_profiles_Data_nsd_1s.f_im[n]=fn_nsnd_2gamma_i(ni, n, 1, y);
        }
        
        // resonance part
        if(choice==0 || choice==1)
            for(int n=2; n<ni; n++) 
            {
                r+=pow(nsnd_2gamma_profiles_Data_nsd_1s.kappa_res[n], 2)*( pow(nsnd_2gamma_profiles_Data_nsd_1s.f_re[n], 2) 
                                                                         + pow(nsnd_2gamma_profiles_Data_nsd_1s.f_im[n], 2) );
                
                for(int m=2; m<n; m++) 
                    r+=2.0*nsnd_2gamma_profiles_Data_nsd_1s.kappa_res[n]*nsnd_2gamma_profiles_Data_nsd_1s.kappa_res[m]
                        *( nsnd_2gamma_profiles_Data_nsd_1s.f_re[n]*nsnd_2gamma_profiles_Data_nsd_1s.f_re[m] 
                         + nsnd_2gamma_profiles_Data_nsd_1s.f_im[n]*nsnd_2gamma_profiles_Data_nsd_1s.f_im[m] );
            }
        
        // interference part
        if(choice==0 || choice==2)
        {
            double r_int=0.0;
            for(int n=2; n<ni; n++) r_int+=nsnd_2gamma_profiles_Data_nsd_1s.kappa_res[n]*nsnd_2gamma_profiles_Data_nsd_1s.f_re[n];
            
            r+=2.0*Mnr*r_int;
        }
        
        return r*pow(y*(1.0-y), 3);
    }

    //======================================================================================
    double sigma_nsd_1s_poles(double y, int ni, nsnd_2gamma_profiles_Data &nsnd_2gamma_profiles_Data_nsd_1s)
    { 
        double r=0.0;
        
        // prepare fn tables
        for(int n=2; n<ni; n++) 
        {
            nsnd_2gamma_profiles_Data_nsd_1s.f_re[n]=fn_nsnd_2gamma_r(ni, n, 1, y);
            nsnd_2gamma_profiles_Data_nsd_1s.f_im[n]=fn_nsnd_2gamma_i(ni, n, 1, y);
        }
        
        // sum of resonances
        for(int n=2; n<ni; n++) 
            r+=pow(nsnd_2gamma_profiles_Data_nsd_1s.kappa_res[n], 2)*( pow(nsnd_2gamma_profiles_Data_nsd_1s.f_re[n], 2) 
                                                                     + pow(nsnd_2gamma_profiles_Data_nsd_1s.f_im[n], 2) );
        
        return r*pow(y*(1.0-y), 3);
    }

    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_nsd_1s_nsnd_2gamma(double y, int ni, double Mnr, nsnd_2gamma_profiles_Data &nsnd_2gamma_profiles_Data_nsd_1s)
    { 
        double r=Mnr*Mnr;
        
        // prepare fn tables
        for(int n=2; n<ni; n++) 
        {
            nsnd_2gamma_profiles_Data_nsd_1s.f_re[n]=fn_nsnd_2gamma_r(ni, n, 1, y);
            nsnd_2gamma_profiles_Data_nsd_1s.f_im[n]=fn_nsnd_2gamma_i(ni, n, 1, y);
        }
        
        // resonance part
        for(int n=2; n<ni; n++) 
        {
            r+=pow(nsnd_2gamma_profiles_Data_nsd_1s.kappa_res[n], 2)*( pow(nsnd_2gamma_profiles_Data_nsd_1s.f_re[n], 2) 
                                                                     + pow(nsnd_2gamma_profiles_Data_nsd_1s.f_im[n], 2) );
            
            for(int m=2; m<n; m++) 
                r+=2.0*nsnd_2gamma_profiles_Data_nsd_1s.kappa_res[n]*nsnd_2gamma_profiles_Data_nsd_1s.kappa_res[m]
                     *(nsnd_2gamma_profiles_Data_nsd_1s.f_re[n]*nsnd_2gamma_profiles_Data_nsd_1s.f_re[m] 
                      +nsnd_2gamma_profiles_Data_nsd_1s.f_im[n]*nsnd_2gamma_profiles_Data_nsd_1s.f_im[m] );
        }
        
        // interference part
        double r_int=0.0;
        for(int n=2; n<ni; n++) r_int+=nsnd_2gamma_profiles_Data_nsd_1s.kappa_res[n]*nsnd_2gamma_profiles_Data_nsd_1s.f_re[n];
        r+=2.0*Mnr*r_int;
        
        return r*pow(y*(1.0-y), 3);
    }

    //======================================================================================
    //
    // different cross sections 2s
    //
    //======================================================================================
    double sigma_2s_1s_non_res(double y)
    { 
        double Mnr=compute_Mnr_2s_1s(y);
        return G_ns(2)*Mnr*Mnr*pow(y*(1.0-y), 3);
    }

    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_2s_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return sigma_2s_1s_non_res(y);
    }
    
    //======================================================================================
    // plot profile for 2s
    //======================================================================================
    void dump_2s_1s_2gamma_profile(string fname)
    {
        dump_mess("2s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int npy=1000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, nsnd_2gamma_profiles_Data_2s_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << sigma_2s_1s_2gamma(xarr[k]) << " " << sigma_2s_1s_non_res(xarr[k]) << endl;

        ofile.close();
        
        return;
    }

    //======================================================================================
    //
    // different cross sections 3s
    //
    //======================================================================================
    double sigma_3s_1s_non_res(double y)
    { 
        double Mnr=compute_Mnr_3s_1s(y);
        return G_ns(3)*Mnr*Mnr*pow(y*(1.0-y), 3);
    }
    
    //======================================================================================
    double sigma_3s_1s_res(double y, int choice)
    { return G_ns(3)*sigma_nsd_1s_res(y, 3, choice, compute_Mnr_3s_1s(y), nsnd_2gamma_profiles_Data_3s_1s); }
    
    //======================================================================================
    double sigma_3s_1s_poles(double y)
    { return G_ns(3)*sigma_nsd_1s_poles(y, 3, nsnd_2gamma_profiles_Data_3s_1s); }
    
    //======================================================================================
    double sigma_3s_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<3; n++) 
            r+=Get_A_np3s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*( Lorentzian(Get_Gamma_np(n)/nuij(3, 1), y-nsnd_2gamma_profiles_Data_3s_1s.yr[n])
                                                        +Lorentzian(Get_Gamma_np(n)/nuij(3, 1), y-(1.0-nsnd_2gamma_profiles_Data_3s_1s.yr[n])) );
        
        return r/PI;
    }
    
    //======================================================================================
    double sigma_3s_1s_res(double y){ return sigma_3s_1s_res(y, 0); }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_3s_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return G_ns(3)*sigma_nsd_1s_nsnd_2gamma(y, 3, compute_Mnr_3s_1s(y), nsnd_2gamma_profiles_Data_3s_1s);
    }
    
    double sigma_3s_1s_2gamma_ratio(double y){ return sigma_3s_1s_2gamma(y)/sigma_3s_1s_sum_of_Lorentzians(y); }
    
    //======================================================================================
    // plot profile for 3s
    //======================================================================================
    void dump_3s_1s_2gamma_profile(string fname)
    {
        dump_mess("3s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=3;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, nsnd_2gamma_profiles_Data_3s_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << xarr[k] /( nsnd_2gamma_profiles_Data_3s_1s.yr[which_res] ) << " " << sigma_3s_1s_2gamma(xarr[k]) << " " << sigma_3s_1s_non_res(xarr[k]) << " " 
                  << sigma_3s_1s_res(xarr[k]) << " " << sigma_3s_1s_res(xarr[k], 1) << " " << sigma_3s_1s_res(xarr[k], 2) << " " 
                  << sigma_3s_1s_2gamma(xarr[k])/sigma_3s_1s_poles(xarr[k]) << " " << sigma_3s_1s_poles(xarr[k]) << " "
                  << sigma_3s_1s_2gamma(xarr[k])/sigma_3s_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_3s_1s_sum_of_Lorentzians(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }

    //======================================================================================
    //
    // different cross sections 3d
    //
    //======================================================================================
    double sigma_3d_1s_non_res(double y)
    { 
        double Mnr=compute_Mnr_3d_1s(y);
        return G_nd(3)*Mnr*Mnr*pow(y*(1.0-y), 3);
    }
    
    //======================================================================================
    double sigma_3d_1s_res(double y, int choice)
    { return G_nd(3)*sigma_nsd_1s_res(y, 3, choice, compute_Mnr_3d_1s(y), nsnd_2gamma_profiles_Data_3d_1s); }
    
    //======================================================================================
    double sigma_3d_1s_poles(double y)
    { return G_nd(3)*sigma_nsd_1s_poles(y, 3, nsnd_2gamma_profiles_Data_3d_1s); }
    
    //======================================================================================
    double sigma_3d_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<3; n++) 
            r+=Get_A_np3d(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*( Lorentzian(Get_Gamma_np(n)/nuij(3, 1), y-nsnd_2gamma_profiles_Data_3d_1s.yr[n])
                                                        +Lorentzian(Get_Gamma_np(n)/nuij(3, 1), y-(1.0-nsnd_2gamma_profiles_Data_3d_1s.yr[n])) );
        
        return r/PI;
    }
    
    //======================================================================================
    double sigma_3d_1s_res(double y){ return sigma_3d_1s_res(y, 0); }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_3d_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return G_nd(3)*sigma_nsd_1s_nsnd_2gamma(y, 3, compute_Mnr_3d_1s(y), nsnd_2gamma_profiles_Data_3d_1s);
    }
    
    double sigma_3d_1s_2gamma_ratio(double y){ return sigma_3d_1s_2gamma(y)/sigma_3d_1s_sum_of_Lorentzians(y); }
    
    //======================================================================================
    // plot profile for 3d
    //======================================================================================
    void dump_3d_1s_2gamma_profile(string fname)
    {
        dump_mess("3d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=3;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, nsnd_2gamma_profiles_Data_3d_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << xarr[k] /( nsnd_2gamma_profiles_Data_3d_1s.yr[which_res] ) << " " << sigma_3d_1s_2gamma(xarr[k]) << " " << sigma_3d_1s_non_res(xarr[k]) << " " 
                  << sigma_3d_1s_res(xarr[k]) << " " << sigma_3d_1s_res(xarr[k], 1) << " " << sigma_3d_1s_res(xarr[k], 2) << " " 
                  << sigma_3d_1s_2gamma(xarr[k])/sigma_3d_1s_poles(xarr[k]) << " " << sigma_3d_1s_poles(xarr[k]) << " "
                  << sigma_3d_1s_2gamma(xarr[k])/sigma_3d_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_3d_1s_sum_of_Lorentzians(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }

    //======================================================================================
    //
    // different cross sections 4s && 4d
    //
    //======================================================================================
    double sigma_4s_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<4; n++) 
            r+=Get_A_np4s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*( Lorentzian(Get_Gamma_np(n)/nuij(4, 1), y-nsnd_2gamma_profiles_Data_4s_1s.yr[n])
                                                        +Lorentzian(Get_Gamma_np(n)/nuij(4, 1), y-(1.0-nsnd_2gamma_profiles_Data_4s_1s.yr[n])) );
        
        return r/PI;
    }
    
    double sigma_4d_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<4; n++) 
            r+=Get_A_np4d(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*( Lorentzian(Get_Gamma_np(n)/nuij(4, 1), y-nsnd_2gamma_profiles_Data_4d_1s.yr[n])
                                                        +Lorentzian(Get_Gamma_np(n)/nuij(4, 1), y-(1.0-nsnd_2gamma_profiles_Data_4d_1s.yr[n])) );
        
        return r/PI;
    }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_4s_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return G_ns(4)*sigma_nsd_1s_nsnd_2gamma(y, 4, compute_Mnr_4s_1s(y), nsnd_2gamma_profiles_Data_4s_1s);
    }

    double sigma_4d_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return G_nd(4)*sigma_nsd_1s_nsnd_2gamma(y, 4, compute_Mnr_4d_1s(y), nsnd_2gamma_profiles_Data_4d_1s);
    }
    
    double sigma_4s_1s_2gamma_ratio(double y){ return sigma_4s_1s_2gamma(y)/sigma_4s_1s_sum_of_Lorentzians(y); }
    double sigma_4d_1s_2gamma_ratio(double y){ return sigma_4d_1s_2gamma(y)/sigma_4d_1s_sum_of_Lorentzians(y); }
    
    //======================================================================================
    // plot profile for 4s && 4d
    //======================================================================================
    void dump_4s_1s_2gamma_profile(string fname)
    {
        dump_mess("4s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, nsnd_2gamma_profiles_Data_4s_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << sigma_4s_1s_2gamma(xarr[k]) << " " << sigma_4s_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_4s_1s_2gamma_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }

    void dump_4d_1s_2gamma_profile(string fname)
    {
        dump_mess("4d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, nsnd_2gamma_profiles_Data_4d_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << sigma_4d_1s_2gamma(xarr[k]) << " " << sigma_4d_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_4d_1s_2gamma_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }
    
    //======================================================================================
    //
    // different cross sections 5s && 5d
    //
    //======================================================================================
    double sigma_5s_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<5; n++) 
            r+=Get_A_np5s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*( Lorentzian(Get_Gamma_np(n)/nuij(5, 1), y-nsnd_2gamma_profiles_Data_5s_1s.yr[n])
                                                        +Lorentzian(Get_Gamma_np(n)/nuij(5, 1), y-(1.0-nsnd_2gamma_profiles_Data_5s_1s.yr[n])) );
        
        return r/PI;
    }
    
    double sigma_5d_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<5; n++) 
            r+=Get_A_np5d(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*( Lorentzian(Get_Gamma_np(n)/nuij(5, 1), y-nsnd_2gamma_profiles_Data_5d_1s.yr[n])
                                                        +Lorentzian(Get_Gamma_np(n)/nuij(5, 1), y-(1.0-nsnd_2gamma_profiles_Data_5d_1s.yr[n])) );
        
        return r/PI;
    }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_5s_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return G_ns(5)*sigma_nsd_1s_nsnd_2gamma(y, 5, compute_Mnr_5s_1s(y), nsnd_2gamma_profiles_Data_5s_1s);
    }
    
    double sigma_5d_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return G_nd(5)*sigma_nsd_1s_nsnd_2gamma(y, 5, compute_Mnr_5d_1s(y), nsnd_2gamma_profiles_Data_5d_1s);
    }
    
    double sigma_5s_1s_2gamma_ratio(double y){ return sigma_5s_1s_2gamma(y)/sigma_5s_1s_sum_of_Lorentzians(y); }
    double sigma_5d_1s_2gamma_ratio(double y){ return sigma_5d_1s_2gamma(y)/sigma_5d_1s_sum_of_Lorentzians(y); }
    
    //======================================================================================
    // plot profile for 5s && 5d
    //======================================================================================
    void dump_5s_1s_2gamma_profile(string fname)
    {
        dump_mess("5s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, nsnd_2gamma_profiles_Data_5s_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << sigma_5s_1s_2gamma(xarr[k]) << " " << sigma_5s_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_5s_1s_2gamma_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }
    
    void dump_5d_1s_2gamma_profile(string fname)
    {
        dump_mess("5d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, nsnd_2gamma_profiles_Data_5d_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << sigma_5d_1s_2gamma(xarr[k]) << " " << sigma_5d_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_5d_1s_2gamma_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }

    //======================================================================================
    //
    // different cross sections 6s && 6d
    //
    //======================================================================================
    double sigma_6s_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<6; n++) 
            r+=Get_A_np6s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*( Lorentzian(Get_Gamma_np(n)/nuij(6, 1), y-nsnd_2gamma_profiles_Data_6s_1s.yr[n])
                                                        +Lorentzian(Get_Gamma_np(n)/nuij(6, 1), y-(1.0-nsnd_2gamma_profiles_Data_6s_1s.yr[n])) );
        
        return r/PI;
    }
    
    double sigma_6d_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<6; n++) 
            r+=Get_A_np6d(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*( Lorentzian(Get_Gamma_np(n)/nuij(6, 1), y-nsnd_2gamma_profiles_Data_6d_1s.yr[n])
                                                        +Lorentzian(Get_Gamma_np(n)/nuij(6, 1), y-(1.0-nsnd_2gamma_profiles_Data_6d_1s.yr[n])) );
        
        return r/PI;
    }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_6s_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return G_ns(6)*sigma_nsd_1s_nsnd_2gamma(y, 6, compute_Mnr_6s_1s(y), nsnd_2gamma_profiles_Data_6s_1s);
    }
    
    double sigma_6d_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return G_nd(6)*sigma_nsd_1s_nsnd_2gamma(y, 6, compute_Mnr_6d_1s(y), nsnd_2gamma_profiles_Data_6d_1s);
    }
    
    double sigma_6s_1s_2gamma_ratio(double y){ return sigma_6s_1s_2gamma(y)/sigma_6s_1s_sum_of_Lorentzians(y); }
    double sigma_6d_1s_2gamma_ratio(double y){ return sigma_6d_1s_2gamma(y)/sigma_6d_1s_sum_of_Lorentzians(y); }
    
    //======================================================================================
    // plot profile for 6s && 6d
    //======================================================================================
    void dump_6s_1s_2gamma_profile(string fname)
    {
        dump_mess("6s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, nsnd_2gamma_profiles_Data_6s_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << sigma_6s_1s_2gamma(xarr[k]) << " " << sigma_6s_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_6s_1s_2gamma_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }
    
    void dump_6d_1s_2gamma_profile(string fname)
    {
        dump_mess("6d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, nsnd_2gamma_profiles_Data_6d_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << sigma_6d_1s_2gamma(xarr[k]) << " " << sigma_6d_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_6d_1s_2gamma_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }

    //======================================================================================
    //
    // different cross sections 7s && 7d
    //
    //======================================================================================
    double sigma_7s_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<7; n++) 
            r+=Get_A_np7s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*( Lorentzian(Get_Gamma_np(n)/nuij(7, 1), y-nsnd_2gamma_profiles_Data_7s_1s.yr[n])
                                                        +Lorentzian(Get_Gamma_np(n)/nuij(7, 1), y-(1.0-nsnd_2gamma_profiles_Data_7s_1s.yr[n])) );
        
        return r/PI;
    }
    
    double sigma_7d_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<7; n++) 
            r+=Get_A_np7d(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*( Lorentzian(Get_Gamma_np(n)/nuij(7, 1), y-nsnd_2gamma_profiles_Data_7d_1s.yr[n])
                                                        +Lorentzian(Get_Gamma_np(n)/nuij(7, 1), y-(1.0-nsnd_2gamma_profiles_Data_7d_1s.yr[n])) );
        
        return r/PI;
    }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_7s_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return G_ns(7)*sigma_nsd_1s_nsnd_2gamma(y, 7, compute_Mnr_7s_1s(y), nsnd_2gamma_profiles_Data_7s_1s);
    }
    
    double sigma_7d_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return G_nd(7)*sigma_nsd_1s_nsnd_2gamma(y, 7, compute_Mnr_7d_1s(y), nsnd_2gamma_profiles_Data_7d_1s);
    }
    
    double sigma_7s_1s_2gamma_ratio(double y){ return sigma_7s_1s_2gamma(y)/sigma_7s_1s_sum_of_Lorentzians(y); }
    double sigma_7d_1s_2gamma_ratio(double y){ return sigma_7d_1s_2gamma(y)/sigma_7d_1s_sum_of_Lorentzians(y); }
    
    //======================================================================================
    // plot profile for 7s && 7d
    //======================================================================================
    void dump_7s_1s_2gamma_profile(string fname)
    {
        dump_mess("7s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int npy=70000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, nsnd_2gamma_profiles_Data_7s_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << sigma_7s_1s_2gamma(xarr[k]) << " " << sigma_7s_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_7s_1s_2gamma_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }
    
    void dump_7d_1s_2gamma_profile(string fname)
    {
        dump_mess("7d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int npy=70000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, nsnd_2gamma_profiles_Data_7d_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << sigma_7d_1s_2gamma(xarr[k]) << " " << sigma_7d_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_7d_1s_2gamma_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }

    //======================================================================================
    //
    // different cross sections 8s && 8d
    //
    //======================================================================================
    double sigma_8s_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<8; n++) 
            r+=Get_A_np8s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*( Lorentzian(Get_Gamma_np(n)/nuij(8, 1), y-nsnd_2gamma_profiles_Data_8s_1s.yr[n])
                                                        +Lorentzian(Get_Gamma_np(n)/nuij(8, 1), y-(1.0-nsnd_2gamma_profiles_Data_8s_1s.yr[n])) );
        
        return r/PI;
    }
    
    double sigma_8d_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<8; n++) 
            r+=Get_A_np8d(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*( Lorentzian(Get_Gamma_np(n)/nuij(8, 1), y-nsnd_2gamma_profiles_Data_8d_1s.yr[n])
                                                        +Lorentzian(Get_Gamma_np(n)/nuij(8, 1), y-(1.0-nsnd_2gamma_profiles_Data_8d_1s.yr[n])) );
        
        return r/PI;
    }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_8s_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return G_ns(8)*sigma_nsd_1s_nsnd_2gamma(y, 8, compute_Mnr_8s_1s(y), nsnd_2gamma_profiles_Data_8s_1s);
    }
    
    double sigma_8d_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return G_nd(8)*sigma_nsd_1s_nsnd_2gamma(y, 8, compute_Mnr_8d_1s(y), nsnd_2gamma_profiles_Data_8d_1s);
    }
    
    double sigma_8s_1s_2gamma_ratio(double y){ return sigma_8s_1s_2gamma(y)/sigma_8s_1s_sum_of_Lorentzians(y); }
    double sigma_8d_1s_2gamma_ratio(double y){ return sigma_8d_1s_2gamma(y)/sigma_8d_1s_sum_of_Lorentzians(y); }
    
    //======================================================================================
    // plot profile for 8s && 8d
    //======================================================================================
    void dump_8s_1s_2gamma_profile(string fname)
    {
        dump_mess("8s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int npy=80000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, nsnd_2gamma_profiles_Data_8s_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << sigma_8s_1s_2gamma(xarr[k]) << " " << sigma_8s_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_8s_1s_2gamma_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }
    
    void dump_8d_1s_2gamma_profile(string fname)
    {
        dump_mess("8d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int npy=80000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, nsnd_2gamma_profiles_Data_8d_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << sigma_8d_1s_2gamma(xarr[k]) << " " << sigma_8d_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_8d_1s_2gamma_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }
}






