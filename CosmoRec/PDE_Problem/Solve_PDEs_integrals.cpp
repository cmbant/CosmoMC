//====================================================================================================================
// to compute the integrals over the solution for the photon field
//====================================================================================================================
struct HI_nx_spline_Data
{
    int Ly_n_res;
    int n_nu;
    int memindex;
    double xmin, xmax;
    double Tm, aV, nu21, x_c;
    Gas_of_Atoms *HIA;
    Cosmos *cos;
    
    //==============================================
    // for 2gamma and Raman integrals
    //==============================================
    double tau_S[11];
    double exp[11]; 
};

bool DF_integrals_splines_are_set=0;
HI_nx_spline_Data Spline_for_nx_integral;

//====================================================================================================================
// integrals for resonances
//====================================================================================================================
void define_intervals_Diffusion(double xmin, double xmax, double eps, 
                                vector<double > &intervals, double *crit_points, int np)
{
    // no interval with proper ordering
    if(xmin>=xmax) return; 
    // points outside the crit_point array
    if( (xmin>=crit_points[np-1] && xmax>=crit_points[np-1]) || (xmin<=crit_points[0] && xmax<=crit_points[0]))
    {
        intervals.push_back(xmin);
        intervals.push_back(xmax);
        return;
    }
    
    double xmax0=xmax;
    if(xmax>crit_points[np-1]) xmax=crit_points[np-1]; // Only the intervals with upper boundary strickly 
                                                       // inside the array crit_points should be used
    
    intervals.clear();
    intervals.push_back(xmin);
    int i=0, j=0;
    
    for(i=0; i<np; i++)
    {
        if(crit_points[i]>intervals[j] && crit_points[i]<xmax){ intervals.push_back(crit_points[i]); j++; }
        else if(crit_points[i]>=xmax){ intervals.push_back(xmax); break; }
    }
    // check size of first and last interval
    if(intervals.size()>2)
    {
        double i1=intervals[1]-intervals[0];
        double i2=intervals[2]-intervals[1];
        if(i1/i2<=eps) intervals.erase(intervals.begin()+1);
    }
    if(intervals.size()>2)
    {
        double i1=intervals[intervals.size()-1]-intervals[intervals.size()-2];
        double i2=intervals[intervals.size()-2]-intervals[intervals.size()-3];
        if(i1/i2<=eps) intervals.erase(intervals.end()-2);
    }
    // add last interval if needed. Here the exact distance to the last point is neglected
    if(xmax0>crit_points[np-1]) intervals.push_back(xmax0); 
    
    if(intervals[0]!=xmin || intervals[intervals.size()-1]!=xmax0)
    { 
        cout << " ALERT!!! " << endl; 
        wait_f_r(); 
    }
    
    return;
}


//====================================================================================================================
// simple functions for Sobolev approximation
//====================================================================================================================
double tau_S_function(int Lyn, double z, Cosmos &cos, Gas_of_Atoms &HIA)
{ 
    return 3.0*HIA.HI_Lyn_profile(Lyn).Get_A21()*pow(HIA.HI_Lyn_profile(Lyn).Get_lambda21(), 3)
              *calc_HI_X1s(z)*cos.NH(z)/8.0/PI/cos.H(z); 
}

//====================================================================================================================
// two-gamma & Raman process
//====================================================================================================================
double dnbar_dx_2g_R(double x)
{ return calc_spline_JC(x, Spline_for_nx_integral.memindex); }

double calc_DF(double xmin, double xmax)
{
    //=========================================================================
    // Integrals
    //=========================================================================
    double epsrel=1.0e-6, epsabs=1.0e-50;
    double r=Integrate_using_Patterson_adaptive(xmin, xmax, epsrel, epsabs, dnbar_dx_2g_R);
    
    return r;   
}

//====================================================================================================================
double compute_DF_2gamma(int ni, int li, double z, Cosmos &cos, Gas_of_Atoms &HIA, ofstream &Pfile, 
                         vector<double> &xarr, vector<double> &yarr, vector<double> &Fwork, vector<double> &exp_x,
                         double (*phi_2g)(int n, int k), PDE_solver_functions &PDE_funcs)
{
    if(ni>11){ cerr << " compute_DF_2gamma:: too many intermediate resonances..." << endl; exit(0); }
    
    int nnu_pts=PDE_funcs.index_emission[ni-2];
    double nuik, nui1=HIA.Level(ni, li).Get_Dnu_1s();
    double nu21=HIA.Level(2, 1).Get_Dnu_1s();
    double Tg=cos.TCMB(z);
    //
    double w=(2.0*li+1);
    double exp_i1=exp(-const_h_kb*nui1/Tg);
    double Xnl=calc_HI_Xnl(z, ni, li);
    double X1s=calc_HI_X1s(z);
    double DnLj=Xnl/X1s/w-exp_i1;
    //
    double Dnem=HI_Dnem_for_nl_effective(z, ni, li);
    double DRtot=0.0;
    
    //================================================================================
    // create data for part from Sobolev approximation
    //================================================================================
    for(int n=2; n<ni; n++) 
    {
        nuik=HIA.Level(ni, li).Get_nu21(n, 1);
        double Aik=HIA.Level(ni, li).Get_A21(n, 1);
        //
        Spline_for_nx_integral.exp[n]=exp(-const_h_kb*HIA.HI_Lyn_profile(n).Get_nu21()/Tg);
        Spline_for_nx_integral.tau_S[n]=tau_S_function(n, z, cos, HIA);
        //
        double pem=1.0-calc_HI_pd_nl_splines_effective(z, n, 1);
        double exp_ik=exp(-const_h_kb*nuik/Tg);
        double PS=p_ij(Spline_for_nx_integral.tau_S[n]);
        double Xkp=calc_HI_Xnl(z, n, 1);
        double Atilde=Aik*pem/(1.0-exp_ik);
        //
        double DnLkp=Xkp/X1s/3.0-Spline_for_nx_integral.exp[n];
        double Del=DnLkp-DnLj/exp_ik;
        //
        DRtot+=Atilde*exp_ik*(Del-PS*DnLkp);
    }
        
    //==========================================================================
    // tabulate Dn(nu) for splines
    //==========================================================================
    for(int k=0; k<nnu_pts; k++)
    { 
        double npl=exp_x[k]/(1.0-exp_x[k]);
        double f_x=exp_i1/npl; 

        Fwork[k]=(yarr[k]*f_x-nu21*Dnem)*phi_2g(ni, k); 
    }       
    
    //==========================================================================
    // setup splines
    //==========================================================================
    if(!DF_integrals_splines_are_set) 
    {
        Spline_for_nx_integral.memindex=calc_spline_coeffies_JC(nnu_pts, &xarr[0], &Fwork[0], "DF-integrals");
        DF_integrals_splines_are_set=1;
    }
    else 
        update_spline_coeffies_JC(Spline_for_nx_integral.memindex, nnu_pts, &xarr[0], &Fwork[0]);
    
    //==========================================================================
    //==========================================================================
    double xmin=max(0.5*nui1/nu21, xarr[0]);
    double xmax=min(nui1/nu21, xarr[nnu_pts-2]);
    double a=xmin, b;
    
    double dx_x=1.0e-4;
    double r=0.0, r1, xres;
    
    for(int n=2; n<ni; n++)
    {
        xres=HIA.HI_Lyn_profile(n).Get_nu21()/nu21;
        //----------------------------------------------------------------------
        // below Lyn-resonance
        //----------------------------------------------------------------------
        b=xres*(1.0-dx_x);
        r1=calc_DF(a, b);
        r+=r1;  
        a=b;
        
        //----------------------------------------------------------------------
        // across Lyn-resonance
        //----------------------------------------------------------------------
        b=xres*(1.0+dx_x);
        r1=calc_DF(a, b);
        r+=r1;  
        a=b;
    }
    
    b=xmax*(1.0-dx_x);
    r1=calc_DF(a, b);
    r+=r1;
    
    r1=calc_DF(b, xmax);
    r+=r1;
    
    if(show_messages>=2) cout << " difference:: (ni, li)= " << ni << " " << li << " Dr-2g = " << r      
                              << " DRtot= " << DRtot << " " << w*(r-DRtot)/DnLj 
                              << " ( error-check " << DnLj/Dnem-1.0 << " )" << endl;
    
    return w*(r-DRtot)/DnLj;
}

//====================================================================================================================
// Raman-process
//====================================================================================================================
double compute_DF_Raman(int ni, int li, int nmax, double z, Cosmos &cos, Gas_of_Atoms &HIA, ofstream &Pfile, 
                        vector<double> &xarr, vector<double> &yarr, vector<double> &Fwork, vector<double> &exp_x,
                        double (*phi_R)(int n, int k), PDE_solver_functions &PDE_funcs)
{
    if(nmax>10){ cerr << " compute_DF_Raman:: too many intermediate resonances..." << endl; exit(0); }
    
    int k0=PDE_funcs.index_emission[ni-2];
    int nnu_pts=xarr.size()-k0;
    double nuik, nui1=HIA.Level(ni, li).Get_Dnu_1s();
    double nu21=HIA.Level(2, 1).Get_Dnu_1s();
    double Tg=cos.TCMB(z);
    //
    double w=(2.0*li+1.0);
    double exp_i1=exp(-const_h_kb*nui1/Tg);
    double Xnl=calc_HI_Xnl(z, ni, li);
    double X1s=calc_HI_X1s(z);
    double DnLj=Xnl/X1s/w-exp_i1;
    //
    double Dnem=HI_Dnem_for_nl_effective(z, ni, li);
    double DRtot=0.0;
    
    //================================================================================
    // create data for part from Sobolev approximation
    //================================================================================
    for(int n=ni+1; n<=nmax; n++) 
    {
        nuik=HIA.Level(n, 1).Get_nu21(ni, li);
        double Aik=HIA.Level(n, 1).Get_A21(ni, li);
        //
        Spline_for_nx_integral.exp[n]=exp(-const_h_kb*HIA.HI_Lyn_profile(n).Get_nu21()/Tg);
        Spline_for_nx_integral.tau_S[n]=tau_S_function(n, z, cos, HIA);
        //
        double pem=1.0-calc_HI_pd_nl_splines_effective(z, n, 1);
        double exp_ik=exp(-const_h_kb*nuik/Tg);
        double PS=p_ij(Spline_for_nx_integral.tau_S[n]);
        double Xkp=calc_HI_Xnl(z, n, 1);
        double Atilde=3.0/w*Aik*pem/(1.0-exp_ik);
        //
        double DnLkp=Xkp/X1s/3.0-Spline_for_nx_integral.exp[n];
        double Del=DnLkp-DnLj*exp_ik;
        //
        DRtot+=Atilde*(Del-PS*DnLkp);
    }
    
    //==========================================================================
    // tabulate Dn(nu) for splines
    //==========================================================================
    for(int k=0; k<nnu_pts; k++)
    { 
        double npl=exp_x[k+k0]/(1.0-exp_x[k+k0]);
        double f_x=exp_i1/npl; 

        Fwork[k+k0]=(yarr[k+k0]*f_x-nu21*Dnem)*phi_R(ni, k+k0); 
    }       
    
    //==========================================================================
    // setup splines
    //==========================================================================
    if(!DF_integrals_splines_are_set) 
    {
        Spline_for_nx_integral.memindex=calc_spline_coeffies_JC(nnu_pts, &xarr[0+k0], &Fwork[0+k0], "DF-integrals");
        DF_integrals_splines_are_set=1;
    }
    else update_spline_coeffies_JC(Spline_for_nx_integral.memindex, nnu_pts, &xarr[0+k0], &Fwork[0+k0]);
    
    //==========================================================================
    //==========================================================================
    double xmin=max(nui1/nu21, xarr[0+k0]);
    double xmax=xarr[nnu_pts+k0-2];
    double a=xmin, b;
    
    double dx_x=1.0e-4;
    double r=0.0, r1, xres;

    //----------------------------------------------------------------------
    // lower boundary contains large phase space density (better to separate it)
    //----------------------------------------------------------------------
    b=a*(1.0+0.01);
    r1=calc_DF(a, b);
    r+=r1;  
    a=b;
    
    for(int n=ni+1; n<=nmax; n++)
    {
        xres=HIA.HI_Lyn_profile(n).Get_nu21()/nu21;
        //----------------------------------------------------------------------
        // below Lyn-resonance (term from Sobolev xmin-->0)
        //----------------------------------------------------------------------
        b=xres*(1.0-dx_x);
        r1=calc_DF(a, b);
        r+=r1;  
        a=b;
        
        //----------------------------------------------------------------------
        // across Lyn-resonance
        //----------------------------------------------------------------------
        b=xres*(1.0+dx_x);
        r1=calc_DF(a, b);
        r+=r1;  
        a=b;
    }
    
    b=xmax*(1.0-dx_x);
    r1=calc_DF(a, b);
    r+=r1;
    
    r1=calc_DF(b, xmax);
    r+=r1;
    
    if(show_messages>=2) cout << " difference:: (ni, li)= " << ni << " " << li << " Dr-Raman = " << r 
                              << " DRtot= " << DRtot << " " << w*(r-DRtot)/DnLj 
                              << " ( error-check " << DnLj/Dnem-1.0 << " )" << endl;
    
    return w*(r-DRtot)/DnLj;
}

//====================================================================================================================
// DI1 output function
//====================================================================================================================
double dI1_2s1s_abs_Patt(double y, void *p)
{
    double *var=(double *) p;
    double xc=var[1];
    //
    double npl_prim=1.0/( exp(xc*(1.0-y)) - 1.0 );
    //
    double Dn=calc_spline_JC(y, Spline_for_nx_integral.memindex);
    //
    double phi_em=phi_2s1s_em_spectrum(y);  
    return phi_em*npl_prim*Dn;
}

double calc_DI1_2s_Patt(HI_nx_spline_Data &SD)
{
    //=========================================================================
    // 1s-->2s Integral; only high nu part (nu/nu21>=0.5) matters
    //=========================================================================
    double r=0.0;
    double epsrel=5.0e-7, epsabs=1.0e-60;
    
    double pars[]={SD.nu21, SD.x_c, SD.xmin, SD.xmax};
    void *p=&pars;
    
    r=Integrate_using_Patterson_adaptive(max(SD.xmin, 0.5), 1.0, epsrel, epsabs, dI1_2s1s_abs_Patt, p);
    
    return r/SD.nu21;   
}

void compute_DI1_2s_and_dump_it(double z, vector<double> &xarr, vector<double> &nx_arr, Gas_of_Atoms &HIA, Cosmos &cos, 
                                ofstream &Ifile, vector<double> &DI1_2s_vec, bool output=0)
{
    int nnu_pts=xarr.size();
    
    Spline_for_nx_integral.nu21=HIA.Level(2, 0).Get_Dnu_1s();
    Spline_for_nx_integral.x_c=const_h_kb*Spline_for_nx_integral.nu21/cos.TCMB(z);
    Spline_for_nx_integral.xmin=xarr[0];
    Spline_for_nx_integral.xmax=xarr[nnu_pts-1];
    
    if(!DF_integrals_splines_are_set) 
    {
        Spline_for_nx_integral.memindex=calc_spline_coeffies_JC(nnu_pts, &xarr[0], &nx_arr[0], "DF-integrals");
        DF_integrals_splines_are_set=1;
    }
    else update_spline_coeffies_JC(Spline_for_nx_integral.memindex, nnu_pts, &xarr[0], &nx_arr[0]);
    
    double DI1=calc_DI1_2s_Patt(Spline_for_nx_integral);
    
    double Xnl=calc_HI_Xnl(z, 2, 0);
    double X1s=calc_HI_X1s(z);
    double Del=Xnl-exp(-Spline_for_nx_integral.x_c)*X1s;    
    DI1/=Del;
    
    if(show_messages>=2) cout << " 2s-1s-2g " << z << " DI1= " << DI1 << endl;
    if(output) Ifile << z << " " << DI1 << endl;
    
    DI1_2s_vec.push_back(DI1);
    
    return;
}

