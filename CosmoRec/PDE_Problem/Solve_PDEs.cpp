//======================================================================================
// Author Jens Chluba (July 2010)
//======================================================================================
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <sstream>

//======================================================================================
// my stuff
//======================================================================================
#include "physical_consts.h"
#include "routines.h"
#include "Cosmos.h"

#include "Sobolev.h"
#include "Patterson.h"

#include "Load.Populations.HI.h"
#include "HI_pd_Rp_splines_effective.h"

#include "PDE_solver.h"
#include "define_PDE.h"
#include "Solve_PDEs.h"

#include "two_photon_profile.2s1s.h"
#include "HI_Transition_Data.h"
#include "Raman_profiles.h"
#include "nsnd_2gamma_profiles.h"


//====================================================================================================================
using namespace std;
using namespace HI_Transition_Data;
using namespace Raman_profiles;
using namespace nsnd_2gamma_profiles;

//====================================================================================================================
//
// flags for message output and aVoigt setting
//
//====================================================================================================================
bool use_full_width=1;
int show_messages=0;


//====================================================================================================================
//
// important functions
//
//====================================================================================================================
double HI_Dnem_for_nl_effective(double z, int n, int l){ return calc_HI_Dnem_nl_splines_effective(z, n, l); }
double HI_pd_for_nl_effective(double z, int n, int l){ return calc_HI_pd_nl_splines_effective(z, n, l); }

//====================================================================================================================
double Dn_nu_ref_effective(int Lyn, double z, double Dn_ref, double pd, double nu, Cosmos &cos, Gas_of_Atoms &HIA)
{
    double Tm=calc_HI_rho(z)*cos.TCMB(z);
    double tau_d=pd*3.0*HIA.HI_Lyn_profile(Lyn).Get_A21()*pow(HIA.HI_Lyn_profile(Lyn).Get_lambda21(), 3)
                   *calc_HI_X1s(z)*cos.NH(z)/8.0/PI/cos.H(z);
    
    double xi_I=HIA.HI_Lyn_profile(Lyn).xi_Int(HIA.HI_Lyn_profile(Lyn).nu2x(nu, Tm), 1.0e+5, HIA.HI_Lyn_profile(Lyn).aVoigt(Tm));
    return Dn_ref*(1.0-exp(-tau_d*xi_I));
}

double Dn_nu_ref_effective_sum(int Lyn_max, double z, double nu, double nu21, Cosmos &cos, Gas_of_Atoms &HIA)
{
    double r=0.0;
    
    for(int ni=2; ni<=Lyn_max; ni++) r+=Dn_nu_ref_effective(ni, z, HI_Dnem_for_nl_effective(z, ni, 1)*nu21, HI_pd_for_nl_effective(z, ni, 1), nu, cos, HIA);
    
    return r;
}



//====================================================================================================================
//
// different Integrals over the spectral distortion
//
//====================================================================================================================
#include "./Solve_PDEs_integrals.cpp"

//====================================================================================================================
// different Grids
//====================================================================================================================
#include "./Solve_PDEs_grid.cpp"


//====================================================================================================================
//
// functions to setup PDESolver; this should only happen the first time
//
//====================================================================================================================
PDE_solver_functions PDE_funcs;
PDE_Stepper_Data PDE_D;
bool PDE_funcs_are_set=0;

//====================================================================================================================
void arm_PDE_solver(Cosmos &cos, Gas_of_Atoms &HIA, vector<double> &xarr, vector<double> &resonances)
{
    if(!PDE_funcs_are_set)
    {
        //==============================================================
        // here time is wasted a bit, but this is only done once...
        //==============================================================
        PDE_funcs.HI_Xe=calc_HI_Xe;
        PDE_funcs.HI_rho=calc_HI_rho;
        PDE_funcs.HI_X1s=calc_HI_X1s;
        //
        PDE_funcs.HI_Xnl=calc_HI_Xnl;
        PDE_funcs.HI_Dnem_nl=HI_Dnem_for_nl_effective;
        PDE_funcs.HI_pd_nl=HI_pd_for_nl_effective;
        
        Setup_PDE_Data(cos, HIA, PDE_funcs);
        
        //==============================================================
        int npts=xarr.size();
        int nresmax=resonances.size()+1;
        
        //==============================================================
        // save ratio of Raman profiles once grid is set
        //==============================================================
        init_Raman_profiles(xarr[npts-1]-1.0, Get_nmax_Raman_correction());
        //==============================================================
        //  dump_2s_1s_Raman_profile(COSMORECDIR+"./temp/test.2s.R.dat");
        //  dump_3s_1s_Raman_profile(COSMORECDIR+"./temp/test.3s.R.dat");
        //  dump_3d_1s_Raman_profile(COSMORECDIR+"./temp/test.3d.R.dat");
        //  dump_4s_1s_Raman_profile(COSMORECDIR+"./temp/test.4s.R.dat");
        //  dump_4d_1s_Raman_profile(COSMORECDIR+"./temp/test.4d.R.dat");
        //  dump_5s_1s_Raman_profile(COSMORECDIR+"./temp/test.5s.R.dat");
        //  dump_5d_1s_Raman_profile(COSMORECDIR+"./temp/test.5d.R.dat");
        
        init_nsnd_2gamma_profiles(Get_nmax_two_g_correction());
        //==============================================================
        //  dump_2s_1s_2gamma_profile(COSMORECDIR+"./temp/test.2s.2g.dat");
        //  dump_3s_1s_2gamma_profile(COSMORECDIR+"./temp/test.3s.2g.dat");
        //  dump_3d_1s_2gamma_profile(COSMORECDIR+"./temp/test.3d.2g.dat");
        //  dump_4s_1s_2gamma_profile(COSMORECDIR+"./temp/test.4s.2g.dat");
        //  dump_4d_1s_2gamma_profile(COSMORECDIR+"./temp/test.4d.2g.dat");
        //  dump_5s_1s_2gamma_profile(COSMORECDIR+"./temp/test.5s.2g.dat");
        //  dump_5d_1s_2gamma_profile(COSMORECDIR+"./temp/test.5d.2g.dat");
        //  dump_6s_1s_2gamma_profile(COSMORECDIR+"./temp/test.6s.2g.dat");
        //  dump_6d_1s_2gamma_profile(COSMORECDIR+"./temp/test.6d.2g.dat");
        
        init_profile_memory(npts, nresmax, PDE_funcs);
        
        //==============================================================
        // setting profile ratios
        //==============================================================
        double xres_loc;

        //==============================================================
        // Raman profiles
        //==============================================================
        if(Get_nmax_Raman_correction()>=2) 
            for(int k=0; k<npts; k++) PDE_funcs.Raman_ns_1s[0].ratio[k]=sigma_2s_1s_Raman_ratio(xarr[k]-1.0);

        for(int ni=3; ni<=Get_nmax_Raman_correction(); ni++)
        {
            xres_loc=(1.0-1.0/ni/ni)/0.75;
            
            for(int k=0; k<npts; k++)
            {
                PDE_funcs.Raman_ns_1s[ni-2].ratio[k]=sigma_ns_1s_Raman_ratio(ni, xarr[k]/xres_loc-1.0);
                PDE_funcs.Raman_nd_1s[ni-2].ratio[k]=sigma_nd_1s_Raman_ratio(ni, xarr[k]/xres_loc-1.0);
            }
        }

        //==============================================================
        for(int ni=(int)max(2, Get_nmax_Raman_correction()+1); ni<=nresmax; ni++)
        {
            for(int k=0; k<npts; k++)
            {
                PDE_funcs.Raman_ns_1s[ni-2].ratio[k]=1.0;
                PDE_funcs.Raman_nd_1s[ni-2].ratio[k]=1.0;
            }
        }
        
        //==============================================================
        // two-photon profiles
        //==============================================================
        if(Get_nmax_two_g_correction()>=2) 
            for(int k=0; k<npts; k++) PDE_funcs.two_g_ns_1s[0].ratio[k]=sigma_2s_1s_2gamma(xarr[k]);
        
        for(int ni=3; ni<=Get_nmax_two_g_correction(); ni++)
        {
            xres_loc=(1.0-1.0/ni/ni)/0.75;
        
            for(int k=0; k<npts; k++)
            {
                PDE_funcs.two_g_ns_1s[ni-2].ratio[k]=sigma_ns_1s_2gamma_ratio(ni, xarr[k]/xres_loc);
                PDE_funcs.two_g_nd_1s[ni-2].ratio[k]=sigma_nd_1s_2gamma_ratio(ni, xarr[k]/xres_loc);
            }
        }

        //==============================================================
        for(int ni=(int)max(2, Get_nmax_two_g_correction()+1); ni<=nresmax; ni++)
        {
            for(int k=0; k<npts; k++)
            {
                PDE_funcs.two_g_ns_1s[ni-2].ratio[k]=1.0;
                PDE_funcs.two_g_nd_1s[ni-2].ratio[k]=1.0;
            }
        }
        
        //==============================================================
        // A-coefficients and splitting points
        //==============================================================
        int index=0;
        PDE_funcs.index_emission.clear();
        for(; index<npts; index++) if(xarr[index]>=0.5){ PDE_funcs.index_2=index; break; }
        //
        index=0;
        for(int k=0; k<(int)resonances.size(); k++)
        {
            for(int m=0; m<(int)resonances.size(); m++)
            {
                PDE_funcs.Voigt_profiles_A_npns[m][k]=Get_A_npks(m+2, k+2);
                PDE_funcs.Voigt_profiles_A_npnd[m][k]=Get_A_npkd(m+2, k+2);
            }
            //
            for(; index<npts; index++) if(xarr[index]>=resonances[k]){ PDE_funcs.index_emission.push_back(index); break; }
        }   
        
        //==============================================================
        // setup for main run
        //==============================================================
        init_PDE_Stepper_Data(PDE_D, npts);
        setup_Lagrange_interpolation_coefficients_O2(PDE_D, xarr);
    
        PDE_funcs_are_set=1;
    }
    
    reset_PDE_solver_variables();
    
    return;
}

//====================================================================================================================
struct Solve_PDE_Data
{
    int npts;
    int nresmax;
    vector<double> resonances;
    vector<double> Fwork;
    vector<double> xarr, yarr;
};

Solve_PDE_Data SPDE_D;
bool SPDE_D_is_set=0;

//====================================================================================================================
void init_Solve_PDE_Data(int nresmax, int npts_res)
{
    //==================================================================
    if(SPDE_D_is_set && nresmax!=SPDE_D.nresmax)
    { 
        cout << " init_Solve_PDE_Data :: the number of resonances has changed."
             << " Please check your code... Resetting. " << endl; 
        
        SPDE_D_is_set=0;
        PDE_funcs_are_set=0;
    } 
    
    //==================================================================
    if(!SPDE_D_is_set)
    {
        SPDE_D.resonances.clear();
        SPDE_D.Fwork.clear();
        SPDE_D.xarr.clear(); SPDE_D.yarr.clear();
        
        //==============================================================
        // set resonances
        //==============================================================
        SPDE_D.nresmax=nresmax;
        for(int n=2; n<=nresmax; n++) SPDE_D.resonances.push_back((1.0-1.0/n/n)/0.75);
        SPDE_D.npts=npts_res*SPDE_D.resonances.size()+2200;
        //
        SPDE_D.Fwork.resize(SPDE_D.npts);
        SPDE_D.xarr.resize(SPDE_D.npts);
        SPDE_D.yarr.resize(SPDE_D.npts);
        
        //==============================================================
        // prepare grid
        //==============================================================
        init_PDE_xarr_cores(SPDE_D.npts, npts_res, SPDE_D.xarr, 0.001, ( 1.0-1.0/(nresmax+1)/(nresmax+1) )/0.75, SPDE_D.resonances);
        
        SPDE_D_is_set=1;
    }   
    
    //==================================================================
    // initial solution
    //==================================================================
    for(int k=0; k<(int)SPDE_D.yarr.size(); k++) SPDE_D.yarr[k]=0.0;
    
    return;
}

//====================================================================================================================
//
// to output solution for spectrum
//
//====================================================================================================================
void output_solution(int output_count, string endname, double nu21, Cosmos &cos, Gas_of_Atoms &HIA, 
                     double zout, int nresmax, bool boundary_up)
{
    string fname=COSMORECDIR+"./temp/DPesc/sol."+int_to_string(output_count, 5)+endname;
    ofstream ofile(fname.c_str());
    ofile.precision(8);
    
    fname=COSMORECDIR+"./temp/DPesc/sol.approx."+int_to_string(output_count, 5)+endname;
    ofstream ofileappr(fname.c_str());
    ofileappr.precision(8);
    
    //==============================================================
    // reference spectrum (at zout)
    //==============================================================
    double x_c=const_h_kb*nu21/cos.TCMB(zout);
    double Dnem=HI_Dnem_for_nl_effective(zout, 2, 1);
    double DnL_ref=Dnem*nu21;
    double rho=calc_HI_rho(zout);
    double Tm=rho*cos.TCMB(zout);
    int nmax=(boundary_up ? nresmax+1 : nresmax);
    //
    double exp_xc=exp(x_c);
    double exp_xc_rho=exp(x_c/rho);
    double exp_xc_2=exp(x_c*32.0/27.0);
    //
    for(int i=0; i<SPDE_D.npts; i++)
    { 
        double exp_fac_Tg=PDE_funcs.exp_x[i]*exp_xc;
        double exp_fac_Tg_2=PDE_funcs.exp_x[i]*exp_xc_2;
        double exp_fac_Te=PDE_funcs.exp_x[i]*exp_xc_rho;
        
        ofile << SPDE_D.xarr[i] << " " << HIA.HI_Lyn_profile(2).nu2x(SPDE_D.xarr[i]*nu21, Tm) << " " << HIA.HI_Lyn_profile(3).nu2x(SPDE_D.xarr[i]*nu21, Tm) 
              << " " << SPDE_D.yarr[i] << " " << SPDE_D.yarr[i]*pow(SPDE_D.xarr[i], 3) << " " << DnL_ref << " " << DnL_ref*exp_fac_Te << " ";
        
        double Dnu_2=Dn_nu_ref_effective(2, zout, HI_Dnem_for_nl_effective(zout, 2, 1)*nu21, HI_pd_for_nl_effective(zout, 2, 1), SPDE_D.xarr[i]*nu21, cos, HIA);
        double Dnu_3=Dn_nu_ref_effective(3, zout, HI_Dnem_for_nl_effective(zout, 3, 1)*nu21, HI_pd_for_nl_effective(zout, 3, 1), SPDE_D.xarr[i]*nu21, cos, HIA);
        double Dnu_tot=Dn_nu_ref_effective_sum(nmax, zout, SPDE_D.xarr[i]*nu21, nu21, cos, HIA);
        
        ofileappr << SPDE_D.xarr[i] << " " << HIA.HI_Lyn_profile(2).nu2x(SPDE_D.xarr[i]*nu21, Tm) << " " << HIA.HI_Lyn_profile(3).nu2x(SPDE_D.xarr[i]*nu21, Tm) << " " 
        //
        << Dnu_2 << " " << Dnu_3 << " " << Dnu_tot << " "
        //
        << Dnu_2*pow(SPDE_D.xarr[i], 3) << " " << Dnu_3*pow(SPDE_D.xarr[i], 3) << " " << Dnu_tot*pow(SPDE_D.xarr[i], 3) << " "
        //
        << Dnu_2*exp_fac_Tg << " " << Dnu_3*exp_fac_Tg_2 << endl;
        
        ofile << SPDE_D.yarr[i]/exp_fac_Te - DnL_ref << " " << SPDE_D.yarr[i] << " " << SPDE_D.yarr[i]/exp_fac_Tg << " " << SPDE_D.yarr[i]/exp_fac_Te << " "; 
        
        ofile << endl;
    }           
    
    ofile.close();  
    ofileappr.close();
    
    return;
}

//==============================================================
//
// code with effective rates
//
//==============================================================
int compute_DPesc_with_diffusion_equation_effective(vector<double> &DF_vec_z, 
                                                    vector<vector<double> > &DF_2_gamma_vec, 
                                                    vector<vector<double> > &DF_Raman_vec, 
                                                    vector<double> &DI1_2s_vec, 
                                                    double zs, double ze, int nS_effective, 
                                                    int nmax_2g_corrs, int nmax_R_corrs,
                                                    Cosmos &cos, Gas_of_Atoms &HIA, 
                                                    vector<vector<double> > HI_Solution, 
                                                    int it_num)
{
    Set_Load_Populations_HI_verbosity(show_messages);
    Set_HI_pd_Rp_splines_effective_verbosity(show_messages);
    
    //==============================================================
    // extension for filename
    //==============================================================
    string exten=".test";

    //==============================================================
    // general flags
    //==============================================================
    bool do_output=0;
    bool output_DI1=0;
    bool output_2gR=0;
    bool write_solution_z0=0;

    //==============================================================
    int nShells=HIA.Get_nShells();
    int npts_res=1000;            // points per resonance
    bool boundary_up=0;           // change upper boundary condition
    
    //==============================================================
    // step size
    //==============================================================
    double dz_out=20.0;
    int output_count=0;
    int output_step=10;
    int DF_step=1;
    
    //==============================================================
    // switching on/off Raman/two-gamma profile corrections
    //==============================================================
    // always make sure that Raman-processes are not switched on 
    // above two-gamma processes!
    //--------------------------------------------------------------
    nmax_R_corrs=(nmax_R_corrs<nmax_2g_corrs ? nmax_R_corrs : nmax_2g_corrs);
    //
    if(nmax_2g_corrs<3) switch_off_two_g_corrections();
    else switch_on_two_g_corrections(nmax_2g_corrs, nShells);
    //
    if(nmax_R_corrs<2) switch_off_Raman_corrections();
    else switch_on_Raman_corrections(nmax_R_corrs, nShells);
    //
    int nresmax=(int)min(8, nShells);
    if(nmax_2g_corrs>=3) nresmax=(int)min(nmax_2g_corrs, nresmax); 

    //==============================================================
    // for more that 4 shells the PDE stepper is not precise enough 
    // to 'see' the difference. Higher accuracy will be set. 
    // The precision of the out put should then be like ~<0.01%.
    //==============================================================
    if(nresmax>4){ dz_out=2.0; DF_step=10; }
        
    //==============================================================
    // prepare memory for back-communication
    //==============================================================
    DF_vec_z.clear();
    DI1_2s_vec.clear();
    DF_2_gamma_vec.clear();
    DF_Raman_vec.clear();
    //
    vector<double> DPesc_dum_z;

    //==============================================================
    // two-gamma integrals can always be computed to take into accout 
    // the correction due to time-dependence etc
    //==============================================================
    for(int k=0; k<nresmax-2; k++)
    { 
        DF_2_gamma_vec.push_back(DPesc_dum_z); 
        DF_2_gamma_vec.push_back(DPesc_dum_z); 
    }

    //==============================================================
    // only include the Raman integrals when the profile correction 
    // is switched on
    //==============================================================
    if(Get_nmax_Raman_correction()>=2) DF_Raman_vec.push_back(DPesc_dum_z);
    for(int k=0; k<Get_nmax_Raman_correction()-2; k++)
    { 
        DF_Raman_vec.push_back(DPesc_dum_z); 
        DF_Raman_vec.push_back(DPesc_dum_z); 
    }
    
    //==============================================================
    // load precomputed solutions (load could be avoided...)
    //==============================================================
    compute_Xi_HI_splines(zs, ze, HIA, HI_Solution);

    set_up_splines_for_HI_pd_Rp_effective(ze, zs, nS_effective, cos, HIA);
    
    
    //==============================================================
    // memory for Solve_PDE
    //==============================================================
    init_Solve_PDE_Data(nresmax, npts_res);
    

    //==============================================================
    // memory for PDE Solver
    //==============================================================
    arm_PDE_solver(cos, HIA, SPDE_D.xarr, SPDE_D.resonances);
            

    //==============================================================
    // outputs/names/etc
    //==============================================================
    string fname, endname=exten+".it_"+int_to_string(it_num)+".dat";
    ofstream two_g_file, Raman_file, Ifile;

    if(output_2gR)
    {
        fname=COSMORECDIR+"./temp/DPesc/DF.2_gamma"+endname;
        two_g_file.open(fname.c_str());
        two_g_file.precision(8);
        
        fname=COSMORECDIR+"./temp/DPesc/DF.Raman"+endname;
        Raman_file.open(fname.c_str());
        Raman_file.precision(8);
    }
    
    if(output_DI1)
    {
        fname=COSMORECDIR+"./temp/DPesc/DI1.2s1s"+endname;
        Ifile.open(fname.c_str());
        Ifile.precision(8);
    }
    
    //==============================================================
    // main run
    //==============================================================
    if(show_messages>=1) cout << "\n entering diffusion part " << endl << endl;

    double nu21=HIA.Level(2, 1).Get_Dnu_1s();
    double zin=zs, zout, dz;
    do
    {
        //==============================================================
        // set step-size
        //==============================================================
        if(zin>zs-10.0) dz=-0.5;
        if(zin>zs-100.0) dz=-2.0;
        else dz=-dz_out;
        
        zout=max(ze, zin+dz);
        
        //==============================================================
        // define lower & upper boundary
        //==============================================================
        double xeval=SPDE_D.xarr[0]*(1.0+zin)/(1.0+zout); // (zin > zout)
        int i;
        for(i=0; i<SPDE_D.npts; i++) if(SPDE_D.xarr[i]>=xeval){ i--; break; }
        double y_lower, dydum;
        polint_JC(&SPDE_D.xarr[0], &SPDE_D.yarr[0], SPDE_D.npts, xeval, i, 6, &y_lower, &dydum);
        double y_upper=( boundary_up ? PDE_funcs.HI_Dnem_nl(zout, nresmax+1, 1)*nu21 : 0.0);
        
        //==============================================================
        // solve PDE
        //==============================================================
        Step_PDE_O2t(0.55, zin, zout, SPDE_D.xarr, SPDE_D.yarr, y_lower, y_upper, 
                     PDE_D, def_PDE_Lyn_and_2s1s);       
        
        if(show_messages>=3) 
            cout << " " << zs << " --> " << zout << " # " << output_count 
                 << " y_low= " << y_lower << " i= " << i << endl;

        //==============================================================
        // write results only after some initial evolution
        //==============================================================
        if(zout<zs-200.0)
        {
            //==============================================================
            // output solution
            //==============================================================
            if(!(output_count%output_step) && do_output) 
               output_solution(output_count, endname, nu21, cos, HIA, zout, nresmax, boundary_up);
 
            //==============================================================
            // compute correction
            //==============================================================
            if(!(output_count%DF_step))
            {
                //==============================================================
                DF_vec_z.push_back(zout);
                
                //==============================================================
                // SPDE_D.yarr is Dn_x
                //==============================================================
                compute_DI1_2s_and_dump_it(zout, SPDE_D.xarr, SPDE_D.yarr, HIA, cos, Ifile, DI1_2s_vec, output_DI1);
                
                //==============================================================
                // 2-gamma integrals
                //==============================================================
                double DF;
                int count=0;
                for(int ni=3; ni<=nresmax; ni++)
                {
                    DF=compute_DF_2gamma(ni, 0, zout, cos, HIA, two_g_file, SPDE_D.xarr, SPDE_D.yarr, SPDE_D.Fwork, PDE_funcs.exp_x, phi_ns1s_2g, PDE_funcs);
                    DF_2_gamma_vec[count++].push_back(DF);
                    DF=compute_DF_2gamma(ni, 2, zout, cos, HIA, two_g_file, SPDE_D.xarr, SPDE_D.yarr, SPDE_D.Fwork, PDE_funcs.exp_x, phi_nd1s_2g, PDE_funcs);
                    DF_2_gamma_vec[count++].push_back(DF);
                }
                
                //==============================================================
                // Raman integrals
                //==============================================================
                if(Get_nmax_Raman_correction()>=2)
                {
                    count=0;
                    
                    DF=compute_DF_Raman(2, 0, nresmax, zout, cos, HIA, Raman_file, SPDE_D.xarr, SPDE_D.yarr, SPDE_D.Fwork, PDE_funcs.exp_x, phi_ns1s_Raman, PDE_funcs);
                    DF_Raman_vec[count++].push_back(DF);
                    
                    for(int ni=3; ni<=Get_nmax_Raman_correction(); ni++)
                    {
                        DF=compute_DF_Raman(ni, 0, nresmax, zout, cos, HIA, Raman_file, SPDE_D.xarr, SPDE_D.yarr, SPDE_D.Fwork, PDE_funcs.exp_x, phi_ns1s_Raman, PDE_funcs);
                        DF_Raman_vec[count++].push_back(DF);
                        DF=compute_DF_Raman(ni, 2, nresmax, zout, cos, HIA, Raman_file, SPDE_D.xarr, SPDE_D.yarr, SPDE_D.Fwork, PDE_funcs.exp_x, phi_nd1s_Raman, PDE_funcs);
                        DF_Raman_vec[count++].push_back(DF);
                    }
                }
                
                if(output_2gR)
                {
                    two_g_file << zout << " ";
                    for(int m=0; m<(int)DF_2_gamma_vec.size(); m++) two_g_file << DF_2_gamma_vec[m].back() << " ";
                    two_g_file << endl;
                    //
                    Raman_file << zout << " ";
                    for(int m=0; m<(int)DF_Raman_vec.size(); m++) Raman_file << DF_Raman_vec[m].back() << " ";
                    Raman_file << endl;
                }
            }
            
            output_count++;
            
            if(show_messages>=2) cout << endl;
        }
        
        zin=zout;
    }
    while(zout>ze);
    
    //=====================================================================================
    // write solution for spectrum at z=0
    //=====================================================================================
    if(write_solution_z0)
    {
        fname=COSMORECDIR+"./temp/DPesc/sol.z_0"+endname;
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        for(int i=0; i<SPDE_D.npts; i++) 
            ofile << SPDE_D.xarr[i] << " " << SPDE_D.yarr[i] << " " 
                  << SPDE_D.yarr[i]*pow(SPDE_D.xarr[i], 3) << " " 
                  // spectral distortion at z=0
                  << SPDE_D.xarr[i]*nu21*1.0e-9/(1.0+ze) << " " 
                  << SPDE_D.yarr[i]/nu21*2.0*const_h*const_cl*pow(SPDE_D.xarr[i]*nu21/const_cl/(1.0+ze), 3) 
                  << endl;
        
        ofile.close();      
    }
    
    //=====================================================================================
    // clean up
    //=====================================================================================
    if(output_DI1) Ifile.close();
    if(output_2gR)
    {
        two_g_file.close();
        Raman_file.close();
    }
    
    return 0;
}


