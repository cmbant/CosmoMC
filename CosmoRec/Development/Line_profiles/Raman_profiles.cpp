//======================================================================================
// Author Jens Chluba Aug/Sept 2010
// purpose: compute the first few Raman-profiles
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
#include "Raman_profiles.h"

using namespace std;
using namespace HI_matrix_elements;
using namespace HI_Transition_Data;

//======================================================================================
// Atomic Data for HI atom
//======================================================================================
const double Raman_profiles_C_sd   =9.0/1024.0*pow(const_alpha, 6)*const_cl*const_Ry_inf_icm/(1.0+const_me_mp);
const double Raman_profiles_nu_1sc =const_EH_inf_Hz/(1.0+const_me_mp);

//======================================================================================
// variables
//======================================================================================
const int Raman_profiles_nmax=500;
const int Raman_profiles_xpts=2*256;
const double Raman_profiles_xmin=1.0e-4;

struct Raman_profiles_Data
{
    int res_up;
    int xnpt;
    int memindex;
    double xmax;
    vector<double> kappa_res;
    //
    vector<double> f_re;
    vector<double> f_im;
    vector<double> yr;
};

Raman_profiles_Data Raman_profiles_Data_2s_1s;
Raman_profiles_Data Raman_profiles_Data_3s_1s;
Raman_profiles_Data Raman_profiles_Data_3d_1s;
Raman_profiles_Data Raman_profiles_Data_4s_1s;
Raman_profiles_Data Raman_profiles_Data_4d_1s;
Raman_profiles_Data Raman_profiles_Data_5s_1s;
Raman_profiles_Data Raman_profiles_Data_5d_1s;
Raman_profiles_Data Raman_profiles_Data_6s_1s;
Raman_profiles_Data Raman_profiles_Data_6d_1s;
Raman_profiles_Data Raman_profiles_Data_7s_1s;
Raman_profiles_Data Raman_profiles_Data_7d_1s;

//======================================================================================
// local functions
//======================================================================================
namespace Raman_profiles_local 
{
    double nuij(int n, int np){ return Raman_profiles_nu_1sc*(pow(1.0*np,-2) - pow(1.0*n,-2)); } 
    
    //======================================================================================
    // resonance frequencies
    //======================================================================================
    double y_res(int ni, int n){ return (1.0/ni/ni-1.0/n/n)/(1.0-1.0/ni/ni); }
    
    //======================================================================================
    // energy factor
    //======================================================================================
    double fn_Raman(int ni, int n, int nf, double y)
    { 
        return 1.0/((pow(1.0*ni,-2) -pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y) 
             + 1.0/((pow(1.0*nf,-2) -pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y);
    }

    double fn_Raman_cont(int ni, double x, int nf, double y)
    { 
        return 1.0/(( x*x + pow(1.0*ni,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y) 
             + 1.0/(( x*x + pow(1.0*nf,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y);
    }
    
    //======================================================================================
    // energy factors for the resonances
    //======================================================================================
    double Lorentzian(double a, double b){ return a/(a*a+b*b); }
    
    double fn_Raman_r(int ni, int n, int nf, double y)
    {
        double dy1=(pow(1.0*ni,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y;
        double dy2=(pow(1.0*nf,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y;
        double nuscale=nuij(ni, nf);
        
        //==================================================================================
        // comment: the '+' sign is due to change of phase when using 
        // f= 1/(yp+y-i*d) + 1/(ym-y-i*d) instead of f= 1/(yp+y-i*d) - 1/(y-ym-i*d) as in 
        // CS 2009 paper
        //==================================================================================
        return Lorentzian(dy1, Get_Gamma_np(n)/nuscale) + Lorentzian(dy2, Get_Gamma_np(n)/nuscale);
    }

    double fn_Raman_i(int ni, int n, int nf, double y)
    {
        double dy1=(pow(1.0*ni,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y;
        double dy2=(pow(1.0*nf,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y;
        double nuscale=nuij(ni, nf);
        
        return Lorentzian(Get_Gamma_np(n)/nuscale, dy1) - Lorentzian(Get_Gamma_np(n)/nuscale, dy2);
    }
        
    //======================================================================================
    // normalization factor
    //======================================================================================
    double G_ns(int n){ return Raman_profiles_C_sd    *pow(4.0/3.0*(1.0-1.0/n/n), 5); }
    double G_nd(int n){ return Raman_profiles_C_sd/2.5*pow(4.0/3.0*(1.0-1.0/n/n), 5); }
    
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
        return x*C1s(x)*Cj_x_ptr__(x)*fn_Raman_cont(y[1], x, 1, y[0]);
    }
    
    double IntC1sC_nsd(double y, int ni, double (*Cj_x)(double))
    {
        // y=nu/nu_n1 --> pole at y=yc
        if(y>=1.0/(ni*ni-1.0)) return 0.0;
        
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
    { return R1snp(n)*Rnsnp(n)*fn_Raman(n, n, 1, y); }

    //======================================================================================
    // Matrix element nd-->np
    //======================================================================================
    double Matrix_element_nd_np_1s(double y, int n)
    { return R1snp(n)*Rndnp(n)*fn_Raman(n, n, 1, y); }
    
    //======================================================================================
    // general Matrix element ns/d-->np
    //======================================================================================
    double Matrix_element_j_np_1s(int n, double y, int ni, double (*Rj_np)(int))
    { return R1snp(n)*Rj_np(n)*fn_Raman(ni, n, 1, y); }
    
    double Matrix_element_j_1s_total_non_res(int nmin, double y, int ni, double (*Rj_np)(int), double (*Cj)(double))
    {
        double r=IntC1sC_nsd(y, ni, Cj) + Matrix_element_ns_np_1s(y, ni);
        
        nmin=(int)max(ni+1, nmin);
        
        for(int nint=Raman_profiles_nmax; nint>nmin; nint--) r+=Matrix_element_j_np_1s(nint, y, ni, Rj_np);
        for(int nint=ni-1; nint>=2; nint--)                  r+=Matrix_element_j_np_1s(nint, y, ni, Rj_np);
        
        return r;
    }

    //======================================================================================
    //======================================================================================
    
    //======================================================================================
    // Matrix element Mn_2s_1s
    //======================================================================================
    double Matrix_element_2s_1s_total_non_res(int nmin, double y)
    { return Matrix_element_j_1s_total_non_res(nmin, y, 2, R2snp, C2s); }
    
    //======================================================================================
    // Matrix element Mn_3s_1s && M_3d_1s
    //======================================================================================
    double Matrix_element_3s_1s_total_non_res(int nmin, double y)
    { return Matrix_element_j_1s_total_non_res(nmin, y, 3, R3snp, C3s); }
    
    double Matrix_element_3d_1s_total_non_res(int nmin, double y)
    { return Matrix_element_j_1s_total_non_res(nmin, y, 3, R3dnp, C3d); }
    
    //======================================================================================
    // Matrix element Mn_4s_1s && M_4d_1s
    //======================================================================================
    double Matrix_element_4s_1s_total_non_res(int nmin, double y)
    { return Matrix_element_j_1s_total_non_res(nmin, y, 4, R4snp, C4s); }
    
    double Matrix_element_4d_1s_total_non_res(int nmin, double y)
    { return Matrix_element_j_1s_total_non_res(nmin, y, 4, R4dnp, C4d); }

    //======================================================================================
    // Matrix element Mn_5s_1s && M_5d_1s
    //======================================================================================
    double Matrix_element_5s_1s_total_non_res(int nmin, double y)
    { return Matrix_element_j_1s_total_non_res(nmin, y, 5, R5snp, C5s); }

    double Matrix_element_5d_1s_total_non_res(int nmin, double y)
    { return Matrix_element_j_1s_total_non_res(nmin, y, 5, R5dnp, C5d); }
    
    //======================================================================================
    // Matrix element Mn_6s_1s && M_6d_1s
    //======================================================================================
    double Matrix_element_6s_1s_total_non_res(int nmin, double y)
    { return Matrix_element_j_1s_total_non_res(nmin, y, 6, R6snp, C6s); }
    
    double Matrix_element_6d_1s_total_non_res(int nmin, double y)
    { return Matrix_element_j_1s_total_non_res(nmin, y, 6, R6dnp, C6d); }
    
    //======================================================================================
    // Matrix element Mn_7s_1s && M_7d_1s
    //======================================================================================
    double Matrix_element_7s_1s_total_non_res(int nmin, double y)
    { return Matrix_element_j_1s_total_non_res(nmin, y, 7, R7snp, C7s); }
    
    double Matrix_element_7d_1s_total_non_res(int nmin, double y)
    { return Matrix_element_j_1s_total_non_res(nmin, y, 7, R7dnp, C7d); }
    
    //======================================================================================
    //======================================================================================
    
    //======================================================================================
    // Matrix element using the spline functions after setup
    //======================================================================================
    double compute_Mnr_2s_1s(double y)
    { return calc_spline_JC(max(Raman_profiles_xmin, y), Raman_profiles_Data_2s_1s.memindex)/y/(1.0/3.0-y); }
    
    double compute_Mnr_3s_1s(double y)
    { return calc_spline_JC(max(Raman_profiles_xmin, y), Raman_profiles_Data_3s_1s.memindex)/y/(1.0/8.0-y); }
    
    double compute_Mnr_3d_1s(double y)
    { return calc_spline_JC(max(Raman_profiles_xmin, y), Raman_profiles_Data_3d_1s.memindex)/y/(1.0/8.0-y); }

    double compute_Mnr_4s_1s(double y)
    { return calc_spline_JC(max(Raman_profiles_xmin, y), Raman_profiles_Data_4s_1s.memindex)/y/(1.0/15.0-y); }
    
    double compute_Mnr_4d_1s(double y)
    { return calc_spline_JC(max(Raman_profiles_xmin, y), Raman_profiles_Data_4d_1s.memindex)/y/(1.0/15.0-y); }

    double compute_Mnr_5s_1s(double y)
    { return calc_spline_JC(max(Raman_profiles_xmin, y), Raman_profiles_Data_5s_1s.memindex)/y/(1.0/24.0-y); }
    
    double compute_Mnr_5d_1s(double y)
    { return calc_spline_JC(max(Raman_profiles_xmin, y), Raman_profiles_Data_5d_1s.memindex)/y/(1.0/24.0-y); }

    double compute_Mnr_6s_1s(double y)
    { return calc_spline_JC(max(Raman_profiles_xmin, y), Raman_profiles_Data_6s_1s.memindex)/y/(1.0/35.0-y); }
    
    double compute_Mnr_6d_1s(double y)
    { return calc_spline_JC(max(Raman_profiles_xmin, y), Raman_profiles_Data_6d_1s.memindex)/y/(1.0/35.0-y); }
    
    double compute_Mnr_7s_1s(double y)
    { return calc_spline_JC(max(Raman_profiles_xmin, y), Raman_profiles_Data_7s_1s.memindex)/y/(1.0/48.0-y); }
    
    double compute_Mnr_7d_1s(double y)
    { return calc_spline_JC(max(Raman_profiles_xmin, y), Raman_profiles_Data_7d_1s.memindex)/y/(1.0/48.0-y); }
    
    //======================================================================================
    //======================================================================================
}

using namespace Raman_profiles_local;

namespace Raman_profiles 
{

    //======================================================================================
    // testing 
    //======================================================================================
    void test_Raman_stuff()
    {
        cout << C1s(0.000001) << " " << C2s(0.000001) << endl;
        cout << C1s(10.0) << " " << C2s(30) << endl;
        cout << R1snp(10) << " " << R2snp(130) << endl;

        cout << endl;
        cout << R2snp(10) << " " << R2snp(130) << endl;
        cout << R3snp(10) << " " << R3snp(130) << endl;
        cout << R4snp(10) << " " << R4snp(130) << endl;
        cout << R5snp(10) << " " << R5snp(130) << endl;

        cout << endl;
        cout << R3dnp(10) << " " << R3dnp(130) << endl;
        cout << R4dnp(10) << " " << R4dnp(130) << endl;
        cout << R5dnp(10) << " " << R5dnp(130) << endl;
        
        cout << endl;
        cout << C2s(0.00001) << " " << C2s(30) << endl;
        cout << C3s(0.00001) << " " << C3s(30) << endl;
        cout << C4s(0.00001) << " " << C4s(30) << endl;
        cout << C5s(0.00001) << " " << C5s(30) << endl;

        cout << endl;
        cout << C3d(0.00001) << " " << C3d(30) << endl;
        cout << C4d(0.00001) << " " << C4d(30) << endl;
        cout << C5d(0.00001) << " " << C5d(30) << endl;
        
        wait_f_r();
    }

    //======================================================================================
    // setup routines
    //======================================================================================
    void alloc_memory(int nres_up, Raman_profiles_Data &Raman_profiles_Data_i)
    {
        Raman_profiles_Data_i.res_up=nres_up;
        Raman_profiles_Data_i.kappa_res.clear();
        Raman_profiles_Data_i.kappa_res.resize(Raman_profiles_Data_i.res_up+1, 0.0);
        Raman_profiles_Data_i.yr.clear();
        Raman_profiles_Data_i.yr.       resize(Raman_profiles_Data_i.res_up+1, 0.0);

        Raman_profiles_Data_i.f_re.resize(Raman_profiles_Data_i.res_up+1);
        Raman_profiles_Data_i.f_im.resize(Raman_profiles_Data_i.res_up+1);
        
        return;
    }
    
    void init_splines(int ni, vector<double> &xarr, vector<double> &yarr, Raman_profiles_Data &Raman_profiles_Data_i, double (*M_nr_i)(int, double))
    {
        for(int k=0; k<Raman_profiles_xpts; k++) yarr[k]=xarr[k]*(1.0/(ni*ni-1.0)-xarr[k])*M_nr_i(Raman_profiles_Data_i.res_up, xarr[k]);
        
        Raman_profiles_Data_i.xnpt=Raman_profiles_xpts;
        Raman_profiles_Data_i.xmax=xarr[Raman_profiles_xpts-1];
        Raman_profiles_Data_i.memindex=calc_spline_coeffies_JC(Raman_profiles_Data_i.xnpt, 
                                                               &xarr[0], &yarr[0],
                                                               "Raman_profiles");
        
        return;
    }
    
    void init_mem_and_data(int ni, int nres_up, double xmax, Raman_profiles_Data &Raman_profiles_Data_i, double (*Rj_np)(int),
                           vector<double> &xarr, vector<double> &yarr, double (*M_nr_i)(int, double))
    {
        alloc_memory(nres_up, Raman_profiles_Data_i);
        
        for(int n=ni+1; n<=nres_up; n++) 
        {
            Raman_profiles_Data_i.kappa_res[n]=R1snp(n)*Rj_np(n);
            Raman_profiles_Data_i.yr[n]=y_res(ni, n);
        }
        
        double loc_xmax=(xmax+1.0)*nuij(2, 1)/nuij(ni, 1)-1.0;
        
        init_xarr(Raman_profiles_xmin, loc_xmax, &xarr[0], Raman_profiles_xpts, 0, 0);
        init_splines(ni, xarr, yarr, Raman_profiles_Data_i, M_nr_i);
        
        return;
    }
    
    void dump_mess(string mess)
    {
        cout << " dumping Raman profile for " << mess << endl;
        return;
    }
    
    //======================================================================================
    // setup non-resonant part of Raman-profile; this depends on chosen xmax
    //======================================================================================
    void init_Raman_profiles(double xmax, int nimax)
    {       
        if(nimax<2) return;
        if(nimax>7){ cout << "\n init_Raman_profiles:: Raman profiles are only available up to nmax=7 " << endl; exit(0);}
        
        cout << "\n init_Raman_profiles:: Initializing ns-1s and nd-1s Raman-profiles up to nmax= " << nimax << " xmax= " << xmax<< endl;

        int nres_up=floor( 2.0/sqrt(1.0-3.0*xmax) ); // upper limit on resonances
        nres_up+=1;
        if(nres_up>10){ cout << " More than 10 resonances in computational domain. " << endl; exit(0); }
        
        //======================================================================
        // setup splines for non-resonant part
        //======================================================================
        vector<double> xarr(Raman_profiles_xpts);
        vector<double> yarr(Raman_profiles_xpts);
        
        //======================================================================
        // ns && nd
        //======================================================================
        init_mem_and_data(2, nres_up, xmax, Raman_profiles_Data_2s_1s, R2snp, xarr, yarr, Matrix_element_2s_1s_total_non_res);
        if(nimax>2)
        {
            init_mem_and_data(3, nres_up, xmax, Raman_profiles_Data_3s_1s, R3snp, xarr, yarr, Matrix_element_3s_1s_total_non_res);
            init_mem_and_data(3, nres_up, xmax, Raman_profiles_Data_3d_1s, R3dnp, xarr, yarr, Matrix_element_3d_1s_total_non_res);
            
            if(nimax>3)
            {
                init_mem_and_data(4, nres_up, xmax, Raman_profiles_Data_4s_1s, R4snp, xarr, yarr, Matrix_element_4s_1s_total_non_res);
                init_mem_and_data(4, nres_up, xmax, Raman_profiles_Data_4d_1s, R4dnp, xarr, yarr, Matrix_element_4d_1s_total_non_res);
                
                if(nimax>4)
                {
                    init_mem_and_data(5, nres_up, xmax, Raman_profiles_Data_5s_1s, R5snp, xarr, yarr, Matrix_element_5s_1s_total_non_res);
                    init_mem_and_data(5, nres_up, xmax, Raman_profiles_Data_5d_1s, R5dnp, xarr, yarr, Matrix_element_5d_1s_total_non_res);

                    if(nimax>5)
                    {
                        init_mem_and_data(6, nres_up, xmax, Raman_profiles_Data_6s_1s, R6snp, xarr, yarr, Matrix_element_6s_1s_total_non_res);
                        init_mem_and_data(6, nres_up, xmax, Raman_profiles_Data_6d_1s, R6dnp, xarr, yarr, Matrix_element_6d_1s_total_non_res);

                        if(nimax>6)
                        {
                            init_mem_and_data(7, nres_up, xmax, Raman_profiles_Data_7s_1s, R7snp, xarr, yarr, Matrix_element_7s_1s_total_non_res);
                            init_mem_and_data(7, nres_up, xmax, Raman_profiles_Data_7d_1s, R7dnp, xarr, yarr, Matrix_element_7d_1s_total_non_res);
                        }
                    }
                }
            }
        }
        
        cout << " init_Raman_profiles:: done " << endl;

        return;
    }
    
    //======================================================================================
    //
    // access to different profiles
    //
    //======================================================================================
    double sigma_ns_1s_Raman_ratio(int n, double y)
    {
        if(n==2) return sigma_2s_1s_Raman_ratio(y);
        if(n==3) return sigma_3s_1s_Raman_ratio(y);
        if(n==4) return sigma_4s_1s_Raman_ratio(y);
        if(n==5) return sigma_5s_1s_Raman_ratio(y);
        if(n==6) return sigma_6s_1s_Raman_ratio(y);
        if(n==7) return sigma_7s_1s_Raman_ratio(y);
        
        return 1.0;
    }
    
    double sigma_nd_1s_Raman_ratio(int n, double y)
    {
        if(n==3) return sigma_3d_1s_Raman_ratio(y);
        if(n==4) return sigma_4d_1s_Raman_ratio(y);
        if(n==5) return sigma_5d_1s_Raman_ratio(y);
        if(n==6) return sigma_6d_1s_Raman_ratio(y);
        if(n==7) return sigma_7d_1s_Raman_ratio(y);
        
        return 1.0;
    }   

    //======================================================================================
    //
    // different cross sections ns/d divided by Gnl(!!!)
    //
    //======================================================================================
    double sigma_nsd_1s_res(double y, int ni, int choice, double Mnr, Raman_profiles_Data &Raman_profiles_Data_nsd_1s)
    { 
        double r=0.0;
        int nmax=(int)min(10, Raman_profiles_Data_nsd_1s.res_up);
        
        // prepare fn tables
        for(int n=ni+1; n<=nmax; n++) 
        {
            Raman_profiles_Data_nsd_1s.f_re[n]=fn_Raman_r(ni, n, 1, y);
            Raman_profiles_Data_nsd_1s.f_im[n]=fn_Raman_i(ni, n, 1, y);
        }
        
        // resonance part
        if(choice==0 || choice==1)
            for(int n=ni+1; n<=nmax; n++) 
            {
                r+=pow(Raman_profiles_Data_nsd_1s.kappa_res[n], 2)*( pow(Raman_profiles_Data_nsd_1s.f_re[n], 2) + pow(Raman_profiles_Data_nsd_1s.f_im[n], 2) );
                
                for(int m=ni+1; m<n; m++) 
                    r+=2.0*Raman_profiles_Data_nsd_1s.kappa_res[n]*Raman_profiles_Data_nsd_1s.kappa_res[m]
                        *( Raman_profiles_Data_nsd_1s.f_re[n]*Raman_profiles_Data_nsd_1s.f_re[m] + Raman_profiles_Data_nsd_1s.f_im[n]*Raman_profiles_Data_nsd_1s.f_im[m] );
            }
        
        // interference part
        if(choice==0 || choice==2)
        {
            double r_int=0.0;
            for(int n=ni+1; n<=nmax; n++) r_int+=Raman_profiles_Data_nsd_1s.kappa_res[n]*Raman_profiles_Data_nsd_1s.f_re[n];
            
            r+=2.0*Mnr*r_int;
        }
        
        return r*pow(y*(1.0+y), 3);
    }

    //======================================================================================
    double sigma_nsd_1s_poles(double y, int ni, Raman_profiles_Data &Raman_profiles_Data_nsd_1s)
    { 
        double r=0.0;
        int nmax=(int)min(10, Raman_profiles_Data_nsd_1s.res_up);
        
        // prepare fn tables
        for(int n=ni+1; n<=nmax; n++) 
        {
            Raman_profiles_Data_nsd_1s.f_re[n]=fn_Raman_r(ni, n, 1, y);
            Raman_profiles_Data_nsd_1s.f_im[n]=fn_Raman_i(ni, n, 1, y);
        }
        
        // sum of resonances
        for(int n=ni+1; n<=nmax; n++) 
            r+=pow(Raman_profiles_Data_nsd_1s.kappa_res[n], 2)*( pow(Raman_profiles_Data_nsd_1s.f_re[n], 2) + pow(Raman_profiles_Data_nsd_1s.f_im[n], 2) );
        
        return r*pow(y*(1.0+y), 3);
    }

    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_nsd_1s_Raman(double y, int ni, double Mnr, Raman_profiles_Data &Raman_profiles_Data_nsd_1s)
    { 
        double r=Mnr*Mnr;
        int nmax=(int)min(10, Raman_profiles_Data_nsd_1s.res_up);
        
        // prepare fn tables
        for(int n=ni+1; n<=nmax; n++) 
        {
            Raman_profiles_Data_nsd_1s.f_re[n]=fn_Raman_r(ni, n, 1, y);
            Raman_profiles_Data_nsd_1s.f_im[n]=fn_Raman_i(ni, n, 1, y);
        }
        
        // resonance part
        for(int n=ni+1; n<=nmax; n++) 
        {
            r+=pow(Raman_profiles_Data_nsd_1s.kappa_res[n], 2)*( pow(Raman_profiles_Data_nsd_1s.f_re[n], 2) + pow(Raman_profiles_Data_nsd_1s.f_im[n], 2) );
            
            for(int m=ni+1; m<n; m++) 
                r+=2.0*Raman_profiles_Data_nsd_1s.kappa_res[n]*Raman_profiles_Data_nsd_1s.kappa_res[m]
                     *(Raman_profiles_Data_nsd_1s.f_re[n]*Raman_profiles_Data_nsd_1s.f_re[m] 
                      +Raman_profiles_Data_nsd_1s.f_im[n]*Raman_profiles_Data_nsd_1s.f_im[m] );
        }
        
        // interference part
        double r_int=0.0;
        for(int n=ni+1; n<=nmax; n++) r_int+=Raman_profiles_Data_nsd_1s.kappa_res[n]*Raman_profiles_Data_nsd_1s.f_re[n];
        r+=2.0*Mnr*r_int;
        
        return r*pow(y*(1.0+y), 3);
    }

    //======================================================================================
    //
    // different cross sections 2s
    //
    //======================================================================================
    double sigma_2s_1s_non_res(double y)
    { 
        double Mnr=compute_Mnr_2s_1s(y);
        return G_ns(2)*Mnr*Mnr*pow(y*(1.0+y), 3);
    }

    //======================================================================================
    double sigma_2s_1s_res(double y, int choice)
    { return G_ns(2)*sigma_nsd_1s_res(y, 2, choice, compute_Mnr_2s_1s(y), Raman_profiles_Data_2s_1s); }
    
    //======================================================================================
    double sigma_2s_1s_poles(double y)
    { return G_ns(2)*sigma_nsd_1s_poles(y, 2, Raman_profiles_Data_2s_1s); }

    //======================================================================================
    double sigma_2s_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=3; n<=(int)min(10, Raman_profiles_Data_2s_1s.res_up); n++) 
            r+=Get_A_np2s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*Lorentzian(Get_Gamma_np(n)/nuij(2, 1), y-Raman_profiles_Data_2s_1s.yr[n]);
        
        return r/PI;
    }

    //======================================================================================
    double sigma_2s_1s_res(double y){ return sigma_2s_1s_res(y, 0); }

    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_2s_1s_Raman(double y)
    { 
        if(y<0.0) return 0.0;
        return G_ns(2)*sigma_nsd_1s_Raman(y, 2, compute_Mnr_2s_1s(y), Raman_profiles_Data_2s_1s);
    }
    
    //======================================================================================
    // total cross section with motion
    //======================================================================================
    double sigma_2s_1s_Raman_motion(double y, vector<double> phi_i_y)
    { 
        double phi_V_tot=0.0;
        for(int n=3; n<=(int)min(Raman_profiles_Data_2s_1s.res_up, min(10, phi_i_y.size()+3)); n++) 
            phi_V_tot+=Get_A_np2s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*phi_i_y[n-3];
            
        return phi_V_tot*sigma_2s_1s_Raman(y)/sigma_2s_1s_sum_of_Lorentzians(y);
    }
    
    double sigma_2s_1s_Raman_ratio(double y){ return sigma_2s_1s_Raman(y)/sigma_2s_1s_sum_of_Lorentzians(y); }

    //======================================================================================
    // plot profile for 2s
    //======================================================================================
    void dump_2s_1s_Raman_profile(string fname)
    {
        dump_mess("2s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=3;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_2s_1s.xmax, &xarr[0], npy, 1, 0);
//      init_xarr(Raman_profiles_Data_2s_1s.yr[which_res]-1.0e-3, Raman_profiles_Data_2s_1s.yr[which_res]+1.0e-3, &xarr[0], npy, 0, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << xarr[k] /( Raman_profiles_Data_2s_1s.yr[which_res] ) << " " << sigma_2s_1s_Raman(xarr[k]) << " " << sigma_2s_1s_non_res(xarr[k]) << " " 
                                    << sigma_2s_1s_res(xarr[k]) << " " << sigma_2s_1s_res(xarr[k], 1) << " " << sigma_2s_1s_res(xarr[k], 2) << " " 
                                    << sigma_2s_1s_Raman(xarr[k])/sigma_2s_1s_poles(xarr[k]) << " " << sigma_2s_1s_poles(xarr[k]) << " "
                                    << sigma_2s_1s_Raman(xarr[k])/sigma_2s_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_2s_1s_sum_of_Lorentzians(xarr[k]) << endl;

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
        return G_ns(3)*Mnr*Mnr*pow(y*(1.0+y), 3);
    }
    
    //======================================================================================
    double sigma_3s_1s_res(double y, int choice)
    { return G_ns(3)*sigma_nsd_1s_res(y, 3, choice, compute_Mnr_3s_1s(y), Raman_profiles_Data_3s_1s); }
    
    //======================================================================================
    double sigma_3s_1s_poles(double y)
    { return G_ns(3)*sigma_nsd_1s_poles(y, 3, Raman_profiles_Data_3s_1s); }
    
    //======================================================================================
    double sigma_3s_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=4; n<=(int)min(10, Raman_profiles_Data_3s_1s.res_up); n++) 
            r+=Get_A_np3s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*Lorentzian(Get_Gamma_np(n)/nuij(3, 1), y-Raman_profiles_Data_3s_1s.yr[n]);
        
        return r/PI;
    }
    
    //======================================================================================
    double sigma_3s_1s_res(double y){ return sigma_3s_1s_res(y, 0); }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_3s_1s_Raman(double y)
    { 
        if(y<0.0) return 0.0;
        return G_ns(3)*sigma_nsd_1s_Raman(y, 3, compute_Mnr_3s_1s(y), Raman_profiles_Data_3s_1s);
    }
    
    //======================================================================================
    // total cross section with motion
    //======================================================================================
    double sigma_3s_1s_Raman_motion(double y, vector<double> phi_i_y)
    { 
        double phi_V_tot=0.0;
        for(int n=4; n<=(int)min(Raman_profiles_Data_3s_1s.res_up, min(10, phi_i_y.size()+4)); n++) 
            phi_V_tot+=Get_A_np3s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*phi_i_y[n-4];
        
        return phi_V_tot*sigma_3s_1s_Raman(y)/sigma_3s_1s_sum_of_Lorentzians(y);
    }
    
    double sigma_3s_1s_Raman_ratio(double y){ return sigma_3s_1s_Raman(y)/sigma_3s_1s_sum_of_Lorentzians(y); }
    
    //======================================================================================
    // plot profile for 3s
    //======================================================================================
    void dump_3s_1s_Raman_profile(string fname)
    {
        dump_mess("3s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=4;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_3s_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << xarr[k] /( Raman_profiles_Data_3s_1s.yr[which_res] ) << " " << sigma_3s_1s_Raman(xarr[k]) << " " << sigma_3s_1s_non_res(xarr[k]) << " " 
                  << sigma_3s_1s_res(xarr[k]) << " " << sigma_3s_1s_res(xarr[k], 1) << " " << sigma_3s_1s_res(xarr[k], 2) << " " 
                  << sigma_3s_1s_Raman(xarr[k])/sigma_3s_1s_poles(xarr[k]) << " " << sigma_3s_1s_poles(xarr[k]) << " "
                  << sigma_3s_1s_Raman(xarr[k])/sigma_3s_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_3s_1s_sum_of_Lorentzians(xarr[k]) << endl;
        
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
        return G_nd(3)*Mnr*Mnr*pow(y*(1.0+y), 3);
    }
    
    //======================================================================================
    double sigma_3d_1s_res(double y, int choice)
    { return G_nd(3)*sigma_nsd_1s_res(y, 3, choice, compute_Mnr_3d_1s(y), Raman_profiles_Data_3d_1s); }
    
    //======================================================================================
    double sigma_3d_1s_poles(double y)
    { return G_nd(3)*sigma_nsd_1s_poles(y, 3, Raman_profiles_Data_3d_1s); }
    
    //======================================================================================
    double sigma_3d_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=4; n<=(int)min(10, Raman_profiles_Data_3d_1s.res_up); n++) 
            r+=Get_A_np3d(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*Lorentzian(Get_Gamma_np(n)/nuij(3, 1), y-Raman_profiles_Data_3d_1s.yr[n]);
        
        return r/PI;
    }
    
    //======================================================================================
    double sigma_3d_1s_res(double y){ return sigma_3d_1s_res(y, 0); }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_3d_1s_Raman(double y)
    { 
        if(y<0.0) return 0.0;
        return G_nd(3)*sigma_nsd_1s_Raman(y, 3, compute_Mnr_3d_1s(y), Raman_profiles_Data_3d_1s);
    }
    
    //======================================================================================
    // total cross section with motion
    //======================================================================================
    double sigma_3d_1s_Raman_motion(double y, vector<double> phi_i_y)
    { 
        double phi_V_tot=0.0;
        for(int n=4; n<=(int)min(Raman_profiles_Data_3d_1s.res_up, min(10, phi_i_y.size()+4)); n++) 
            phi_V_tot+=Get_A_np3d(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*phi_i_y[n-4];
        
        return phi_V_tot*sigma_3d_1s_Raman(y)/sigma_3d_1s_sum_of_Lorentzians(y);
    }
    
    double sigma_3d_1s_Raman_ratio(double y){ return sigma_3d_1s_Raman(y)/sigma_3d_1s_sum_of_Lorentzians(y); }
    
    //======================================================================================
    // plot profile for 3d
    //======================================================================================
    void dump_3d_1s_Raman_profile(string fname)
    {
        dump_mess("3d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=4;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_3d_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << xarr[k] /( Raman_profiles_Data_3d_1s.yr[which_res] ) << " " << sigma_3d_1s_Raman(xarr[k]) << " " << sigma_3d_1s_non_res(xarr[k]) << " " 
                  << sigma_3d_1s_res(xarr[k]) << " " << sigma_3d_1s_res(xarr[k], 1) << " " << sigma_3d_1s_res(xarr[k], 2) << " " 
                  << sigma_3d_1s_Raman(xarr[k])/sigma_3d_1s_poles(xarr[k]) << " " << sigma_3d_1s_poles(xarr[k]) << " "
                  << sigma_3d_1s_Raman(xarr[k])/sigma_3d_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_3d_1s_sum_of_Lorentzians(xarr[k]) << endl;
        
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
        for(int n=5; n<=(int)min(10, Raman_profiles_Data_4s_1s.res_up); n++) 
            r+=Get_A_np4s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*Lorentzian(Get_Gamma_np(n)/nuij(4, 1), y-Raman_profiles_Data_4s_1s.yr[n]);
        
        return r/PI;
    }

    double sigma_4d_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=5; n<=(int)min(10, Raman_profiles_Data_4d_1s.res_up); n++) 
            r+=Get_A_np4d(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*Lorentzian(Get_Gamma_np(n)/nuij(4, 1), y-Raman_profiles_Data_4d_1s.yr[n]);
        
        return r/PI;
    }

    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_4s_1s_Raman(double y)
    { 
        if(y<0.0) return 0.0;
        return G_ns(4)*sigma_nsd_1s_Raman(y, 4, compute_Mnr_4s_1s(y), Raman_profiles_Data_4s_1s);
    }
    double sigma_4d_1s_Raman(double y)
    { 
        if(y<0.0) return 0.0;
        return G_nd(4)*sigma_nsd_1s_Raman(y, 4, compute_Mnr_4d_1s(y), Raman_profiles_Data_4d_1s);
    }

    double sigma_4s_1s_Raman_ratio(double y){ return sigma_4s_1s_Raman(y)/sigma_4s_1s_sum_of_Lorentzians(y); }
    double sigma_4d_1s_Raman_ratio(double y){ return sigma_4d_1s_Raman(y)/sigma_4d_1s_sum_of_Lorentzians(y); }

    //======================================================================================
    // plot profile for 4s && 4d
    //======================================================================================
    void dump_4s_1s_Raman_profile(string fname)
    {
        dump_mess("4s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=5;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_4s_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << xarr[k] /( Raman_profiles_Data_4s_1s.yr[which_res] ) << " " 
                  << sigma_4s_1s_Raman(xarr[k]) << " " << sigma_4s_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_4s_1s_Raman_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }

    void dump_4d_1s_Raman_profile(string fname)
    {
        dump_mess("4d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=5;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_4d_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << xarr[k] /( Raman_profiles_Data_4d_1s.yr[which_res] ) << " " 
                             << sigma_4d_1s_Raman(xarr[k]) << " " << sigma_4d_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_4d_1s_Raman_ratio(xarr[k]) << endl;
        
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
        for(int n=6; n<=(int)min(10, Raman_profiles_Data_5s_1s.res_up); n++) 
            r+=Get_A_np5s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*Lorentzian(Get_Gamma_np(n)/nuij(5, 1), y-Raman_profiles_Data_5s_1s.yr[n]);
        
        return r/PI;
    }
    
    double sigma_5d_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=6; n<=(int)min(10, Raman_profiles_Data_5d_1s.res_up); n++) 
            r+=Get_A_np5d(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*Lorentzian(Get_Gamma_np(n)/nuij(5, 1), y-Raman_profiles_Data_5d_1s.yr[n]);
        
        return r/PI;
    }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_5s_1s_Raman(double y)
    { 
        if(y<0.0) return 0.0;
        return G_ns(5)*sigma_nsd_1s_Raman(y, 5, compute_Mnr_5s_1s(y), Raman_profiles_Data_5s_1s);
    }
    double sigma_5d_1s_Raman(double y)
    { 
        if(y<0.0) return 0.0;
        return G_nd(5)*sigma_nsd_1s_Raman(y, 5, compute_Mnr_5d_1s(y), Raman_profiles_Data_5d_1s);
    }
    
    double sigma_5s_1s_Raman_ratio(double y){ return sigma_5s_1s_Raman(y)/sigma_5s_1s_sum_of_Lorentzians(y); }
    double sigma_5d_1s_Raman_ratio(double y){ return sigma_5d_1s_Raman(y)/sigma_5d_1s_sum_of_Lorentzians(y); }
    
    //======================================================================================
    // plot profile for 5s && 5d
    //======================================================================================
    void dump_5s_1s_Raman_profile(string fname)
    {
        dump_mess("5s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=6;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_5s_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << xarr[k] /( Raman_profiles_Data_5s_1s.yr[which_res] ) << " " 
                  << sigma_5s_1s_Raman(xarr[k]) << " " << sigma_5s_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_5s_1s_Raman_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }
    
    void dump_5d_1s_Raman_profile(string fname)
    {
        dump_mess("5d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=6;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_5d_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << xarr[k] /( Raman_profiles_Data_5d_1s.yr[which_res] ) << " " 
                  << sigma_5d_1s_Raman(xarr[k]) << " " << sigma_5d_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_5d_1s_Raman_ratio(xarr[k]) << endl;
        
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
        for(int n=6; n<=(int)min(10, Raman_profiles_Data_6s_1s.res_up); n++) 
            r+=Get_A_np6s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*Lorentzian(Get_Gamma_np(n)/nuij(6, 1), y-Raman_profiles_Data_6s_1s.yr[n]);
        
        return r/PI;
    }
    
    double sigma_6d_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=6; n<=(int)min(10, Raman_profiles_Data_6d_1s.res_up); n++) 
            r+=Get_A_np6d(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*Lorentzian(Get_Gamma_np(n)/nuij(6, 1), y-Raman_profiles_Data_6d_1s.yr[n]);
        
        return r/PI;
    }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_6s_1s_Raman(double y)
    { 
        if(y<0.0) return 0.0;
        return G_ns(6)*sigma_nsd_1s_Raman(y, 6, compute_Mnr_6s_1s(y), Raman_profiles_Data_6s_1s);
    }
    double sigma_6d_1s_Raman(double y)
    { 
        if(y<0.0) return 0.0;
        return G_nd(6)*sigma_nsd_1s_Raman(y, 6, compute_Mnr_6d_1s(y), Raman_profiles_Data_6d_1s);
    }
    
    double sigma_6s_1s_Raman_ratio(double y){ return sigma_6s_1s_Raman(y)/sigma_6s_1s_sum_of_Lorentzians(y); }
    double sigma_6d_1s_Raman_ratio(double y){ return sigma_6d_1s_Raman(y)/sigma_6d_1s_sum_of_Lorentzians(y); }
    
    //======================================================================================
    // plot profile for 6s && 6d
    //======================================================================================
    void dump_6s_1s_Raman_profile(string fname)
    {
        dump_mess("6s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=7;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_6s_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << xarr[k] /( Raman_profiles_Data_6s_1s.yr[which_res] ) << " " 
            << sigma_6s_1s_Raman(xarr[k]) << " " << sigma_6s_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_6s_1s_Raman_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }
    
    void dump_6d_1s_Raman_profile(string fname)
    {
        dump_mess("6d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=7;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_6d_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << xarr[k] /( Raman_profiles_Data_6d_1s.yr[which_res] ) << " " 
            << sigma_6d_1s_Raman(xarr[k]) << " " << sigma_6d_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_6d_1s_Raman_ratio(xarr[k]) << endl;
        
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
        for(int n=6; n<=(int)min(10, Raman_profiles_Data_7s_1s.res_up); n++) 
            r+=Get_A_np7s(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*Lorentzian(Get_Gamma_np(n)/nuij(7, 1), y-Raman_profiles_Data_7s_1s.yr[n]);
        
        return r/PI;
    }
    
    double sigma_7d_1s_sum_of_Lorentzians(double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=6; n<=(int)min(10, Raman_profiles_Data_7d_1s.res_up); n++) 
            r+=Get_A_np7d(n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*Lorentzian(Get_Gamma_np(n)/nuij(7, 1), y-Raman_profiles_Data_7d_1s.yr[n]);
        
        return r/PI;
    }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_7s_1s_Raman(double y)
    { 
        if(y<0.0) return 0.0;
        return G_ns(7)*sigma_nsd_1s_Raman(y, 7, compute_Mnr_7s_1s(y), Raman_profiles_Data_7s_1s);
    }
    double sigma_7d_1s_Raman(double y)
    { 
        if(y<0.0) return 0.0;
        return G_nd(7)*sigma_nsd_1s_Raman(y, 7, compute_Mnr_7d_1s(y), Raman_profiles_Data_7d_1s);
    }
    
    double sigma_7s_1s_Raman_ratio(double y){ return sigma_7s_1s_Raman(y)/sigma_7s_1s_sum_of_Lorentzians(y); }
    double sigma_7d_1s_Raman_ratio(double y){ return sigma_7d_1s_Raman(y)/sigma_7d_1s_sum_of_Lorentzians(y); }
    
    //======================================================================================
    // plot profile for 7s && 7d
    //======================================================================================
    void dump_7s_1s_Raman_profile(string fname)
    {
        dump_mess("7s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=8;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_7s_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << xarr[k] /( Raman_profiles_Data_7s_1s.yr[which_res] ) << " " 
            << sigma_7s_1s_Raman(xarr[k]) << " " << sigma_7s_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_7s_1s_Raman_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }
    
    void dump_7d_1s_Raman_profile(string fname)
    {
        dump_mess("7d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=8;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_7d_1s.xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " << xarr[k] /( Raman_profiles_Data_7d_1s.yr[which_res] ) << " " 
            << sigma_7d_1s_Raman(xarr[k]) << " " << sigma_7d_1s_sum_of_Lorentzians(xarr[k]) << " " << sigma_7d_1s_Raman_ratio(xarr[k]) << endl;
        
        ofile.close();
        
        return;
    }
}






