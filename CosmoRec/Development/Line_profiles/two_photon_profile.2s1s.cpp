//======================================================================================
// Author Jens Chluba Sept/Oct 2010
// purpose: compute the 2s-1s two-photon profile and integral over it
//======================================================================================

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

#include "physical_consts.h"
#include "routines.h"
#include "Patterson.h"

#include "two_photon_profile.2s1s.h"

using namespace std;

//======================================================================================
// Atomic data
//======================================================================================
const double two_photon_profiles_nu_2s1s=const_EH_inf_Hz/(1.0+const_me_mp)*0.75; // ~ 2.466038424e+15;

//======================================================================================
// planck function
//======================================================================================
double nPlanck_function(double nu, double Tg)
{
    double dum=exp(-const_h_kb*nu/Tg);
    return dum/(1.0-dum);
}

//======================================================================================
// 2s two photon decay profile function, normalised to 2
//======================================================================================
double phi_2s1s_em_spectrum(double y)
{
    if(y<0.0000001 || y> 0.9999999) return 0.0;
    double w=y*(1.0-y);
    
    double C=24.5561, a=0.88, b=1.53, c=0.8;
    return C*( w*(1-pow(4.0*w, c))+a*pow(w, b)*pow(4.0*w, c) );
}

void phi_2s1s_function(double nu, double Tg, double nu0, double &phi_em, double &phi_abs)
{
    double y=nu/nu0;
    phi_abs=phi_2s1s_em_spectrum(y);
    
    double np1=nPlanck_function(nu, Tg);
    double np2=nPlanck_function(nu0-nu, Tg);
    phi_em=phi_abs*(1.0+np1)*(1.0+np2);
    
    return;
}
    
//======================================================================================
// integral over 2s-1s emission spectrum with stimulated term
//======================================================================================
double dI2_2s1s_em(double y, void *p)
{
    double *var=(double *) p;
    double Tg=var[0];
    double nu21=var[1];
    double nu=y*nu21;
    double phi_abs=phi_2s1s_em_spectrum(y); 
    double np1=nPlanck_function(nu, Tg);
    double np2=nPlanck_function(nu21-nu, Tg);
    
    return phi_abs*( np1 + np2 + np1*np1 );
}

double Dphi_2s1s_em_int(double Tg)
{
    double a=0.5, b=0.9999999, epsrel=1.0e-8, epsabs=1.0e-60;
    double pars[]={Tg, two_photon_profiles_nu_2s1s};
    void *p=&pars;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dI2_2s1s_em, p);
}

//======================================================================================
// dump I2_planck to disk
//======================================================================================
void compute_stimulated_2s1s_integral(double Tgmin, double Tgmax, int np, string filename)
{
    cout << " compute_stimulated_2s1s_integral :: Tabulating 2s-1s stimulated rate. "  << endl;
    
    double *Tgarr=new double[np];
    double dum;
    init_xarr(Tgmin, Tgmax, Tgarr, np, 1, 0);
    
    ofstream ofile;
    ofile.open(filename.c_str());
    ofile.precision(10);
    
    for(int k=0; k<np; k++)
    {
        dum=Dphi_2s1s_em_int(Tgarr[k]);
        ofile << Tgarr[k] << " " << dum << " " << log(dum) << endl;
    }

    ofile.close();
    
    return;
}


