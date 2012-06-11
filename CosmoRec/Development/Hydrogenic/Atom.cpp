//===================================================================================
// Author: Jens Chluba 
// last change: Oct 2010
//===================================================================================

#include <iostream>
#include <string>
#include <cmath>
#include <stdlib.h>

#include "Atom.h"
#include "Rec_Phot_BB.h"
#include "Voigtprofiles.h"
#include "Oscillator_strength.h"
#include "routines.h"
#include "physical_consts.h"

using namespace std;

//===================================================================================
// Konstructors and Destructors for Electron_Level
//===================================================================================
void Electron_Level::init(int n, int l, int nm, int ZZ, double NNp, int Rec_flag, int mflag)
{
    nn=n; ll=l; nmax=nm;
    gw=2.0*(2.0*l+1.0);
    Z=ZZ;
    Np=NNp;
    mu_red=1.0/(1.0+const_me_mp/Np);
    
    Recombination_flag=Rec_flag;
    mess_flag=mflag;
    
    DE=DE_ul(nn, 1);
    Dnu=nu_ul(nn, 1);
    Eion=E_ion(nn);
    Eion_ergs=E_ion_ergs(nn);
    nuion=nu_ion(nn);

    if(mess_flag>=1)
    {
        if(Recombination_flag==0) cout << " Electron_Level::init: no Recombination routines are set " << endl;
        else if(Recombination_flag==1) Interaction_with_Photons_SH_QSP.init(nn, ll, Z, Np, mess_flag); 
        else if(Recombination_flag==2) Interaction_with_Photons_SH.init(nn, ll, Z, Np, mess_flag); 
        else{ cerr << " Electron_Level::init: This option does not exist" << endl; exit(0); }
    }
    
    calculate_transition_Data();
    
    A=calculate_A_tot();
    Gamma=A/FOURPI;
    
    return;
}

Electron_Level::Electron_Level(int n, int l, int nm, int Z, double Np, int Rec_flag, int mflag)
{ init(n, l, nm, Z, Np, Rec_flag, mflag); }

Electron_Level::~Electron_Level()
{ 
    A_values_down_lm1.clear();
    A_values_down_lp1.clear();
    B_values.clear(); 
}

void Electron_Level::calculate_transition_Data()
{
    //---------------------------------------------
    // transition data stored according to n
    //---------------------------------------------
    A_values_down_lm1.clear();
    A_values_down_lp1.clear();
    ZERO_Data.np=ZERO_Data.lp=0; ZERO_Data.Dnu=ZERO_Data.A21=ZERO_Data.lambda21=0.0;

    B_values.clear();
    for(int n=1; n<=nmax; n++)
    {
        calculate_transition_to_level(n, ll-1);
        calculate_transition_to_level(n, ll+1);
    }
    
    if(mess_flag>1 && nn>1) display_all_Trans_Data();
    
    return;
}

//-----------------------------------------------------------------------------------------
void Electron_Level::calculate_transition_to_level(int n, int l)
{
    if(n==nn || abs(ll-l)!=1 || l>=n || l < 0) return;
    
    if(n<nn)
    {
        //---------------------------------------------
        // transition data stored according to n
        //---------------------------------------------
        Transition_Data_A vA;
        vA.np=n;
        vA.lp=l;
        vA.jp=l+n*(n-1)/2;
        vA.Dnu=nu_ul(nn, n);
        vA.A21=A_SH(Z, Np, nn, ll, n, l);
        vA.lambda21=const_cl/vA.Dnu;
        
        if(l<ll) A_values_down_lm1.push_back(vA);
        if(l>ll) A_values_down_lp1.push_back(vA);
    }
    
    if(n>nn)
    {
        Transition_Data_B v;
        v.np=n; v.lp=l;
        // level index directly
        v.jp=l+n*(n-1)/2;
        B_values.push_back(v);
    }
    
    return;
}

void Electron_Level::display_Trans_Data(Transition_Data_A &v)
{
    cout << "\n %==============================================================%" << endl;
    cout << " % Electron_Level::display_Trans_Data:\n % (" << nn << ", " << ll << ") --> (" << v.np << ", " << v.lp << ") | gw= " << gw << " gwp= " << 2.0*(2.0*v.lp+1.0) << endl;
    cout << " %==============================================================%" << endl;
    cout << " Dnu: " << v.Dnu << " Hz, lambda21: " << v.lambda21 << " cm, A21: " << v.A21  << " sec^-1 )" 
    << endl << endl; 
    
    return;
}

void Electron_Level::display_Trans_Data(Transition_Data_B &v)
{
    cout << "\n %==============================================================%" << endl;
    cout << " % Electron_Level::display_Trans_Data:\n % (" << nn << ", " << ll 
         << ") --> (" << v.np << ", " << v.lp << ") | gw= " << gw << " gwp= " << 2.0*(2.0*v.lp+1.0) << endl;
    cout << " %==============================================================%" << endl;
    
    return;
}

void Electron_Level::display_all_Trans_Data()
{
    if(A_values_down_lm1.size()>0)
    {
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Electron_Level::display_all_Trans_Data: A_values_down_lm1:\n % " << endl;
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%" << endl;
        for(unsigned int i=0; i<A_values_down_lm1.size(); i++) display_Trans_Data(A_values_down_lm1[i]); 
    }
    
    if(A_values_down_lp1.size()>0)
    {
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Electron_Level::display_all_Trans_Data: A_values_down_lp1:\n % " << endl;
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%" << endl;
        for(unsigned int i=0; i<A_values_down_lp1.size(); i++) display_Trans_Data(A_values_down_lp1[i]); 
    }
    
    return;
}

void Electron_Level::display_all_upward_Trans_Data()
{
    if(B_values.size()>0)
    {
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Electron_Level::: display_all_upward_Trans_Data: ntrans= " 
             << B_values.size() << " nmax= " << nmax << "\n % " << endl;
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%" << endl;
        for(unsigned int i=0; i<B_values.size(); i++) display_Trans_Data(B_values[i]); 
    }
    
    return;
}

//=====================================================================================================
// access to photoionization and recombination rates
// Author: Jens Chluba
// changed: 12.04.2006
//=====================================================================================================
void Electron_Level::clear_Interaction_w_photons()
{
    if(Recombination_flag==0)
    { cout << " Electron_Level::clear_Interaction_w_photons: Recombination rates were not setup!!! " << endl; exit(0); }
    else if(Recombination_flag==1) Interaction_with_Photons_SH_QSP.clear();
    else if(Recombination_flag==2) Interaction_with_Photons_SH.clear();
    
    return;
}

double Electron_Level::Get_nu_ionization()
{
    double r=0.0;
    if(Recombination_flag==0)
    { cout << " Electron_Level::Get_nu_ionization: Recombination rates were not setup!!! " << endl; exit(0); }
    else if(Recombination_flag==1) r=Interaction_with_Photons_SH_QSP.Get_nu_ionization();
    else if(Recombination_flag==2) r=Interaction_with_Photons_SH.Get_nu_ionization();
    
    return r;
}        

double Electron_Level::R_nl_c(double T_g)
{
    double r=0.0;
    if(Recombination_flag==0)
    { cout << " Electron_Level::R_nl_c: Recombination rates were not setup!!! " << endl; exit(0); } 
    else if(Recombination_flag==1) r=Interaction_with_Photons_SH_QSP.R_nl_c_Int(T_g);
    else if(Recombination_flag==2) r=Interaction_with_Photons_SH.R_nl_c_Int(T_g);
    return r;
}

double Electron_Level::R_c_nl(double T_g, double rho)
{
    double r=0.0;
    if(Recombination_flag==0)
    { cout << " Electron_Level::R_c_nl: Recombination rates were not setup!!! " << endl; exit(0); } 
    else if(Recombination_flag==1) r=Interaction_with_Photons_SH_QSP.R_c_nl_Int(T_g, rho);
    else if(Recombination_flag==2) r=Interaction_with_Photons_SH.R_c_nl_Int(T_g, rho);
    
    return r;
}

double Electron_Level::dR_c_nl_dTe(double T_g, double rho)
{
    double r=0.0;
    if(Recombination_flag==0)
    { cout << " Electron_Level::R_c_nl: Recombination rates were not setup!!! " << endl; exit(0); } 
    else if(Recombination_flag==1) r=Interaction_with_Photons_SH_QSP.dR_c_nl_dTe_Int(T_g, rho);
    else if(Recombination_flag==2) r=Interaction_with_Photons_SH.dR_c_nl_dTe_Int(T_g, rho);
    
    return r;
}

double Electron_Level::sig_phot_ion_nuc()
{ 
    double r=0.0;
    if(Recombination_flag==0)
    { cout << " Electron_Level::g_phot_ion_E: Recombination rates were not setup!!! " << endl; exit(0); } 
    else if(Recombination_flag==1) r=Interaction_with_Photons_SH_QSP.sig_phot_ion_nuc();
    else if(Recombination_flag==2) r=Interaction_with_Photons_SH.sig_phot_ion_nuc();
    
    return r;
}

double Electron_Level::sig_phot_ion(double nu)
{ 
    double r=0.0;
    if(Recombination_flag==0)
    { cout << " Electron_Level::sig_phot_ion: Recombination rates were not setup!!! " << endl; exit(0); } 
    else if(Recombination_flag==1) r=Interaction_with_Photons_SH_QSP.sig_phot_ion(nu);
    else if(Recombination_flag==2) r=Interaction_with_Photons_SH.sig_phot_ion(nu);
    
    return r;
}

double Electron_Level::g_phot_ion(double nu)
{ 
    double r=0.0;
    if(Recombination_flag==0)
    { cout << " Electron_Level::g_phot_ion: Recombination rates were not setup!!! " << endl; exit(0); } 
    else if(Recombination_flag==1) r=Interaction_with_Photons_SH_QSP.g_phot_ion(nu);
    else if(Recombination_flag==2) r=Interaction_with_Photons_SH.g_phot_ion(nu);
    
    return r;
}

//=====================================================================================================
double Electron_Level::calculate_A_tot()
{
    double r=0.0;
    for(unsigned int i=0; i<A_values_down_lm1.size(); i++) r+=A_values_down_lm1[i].A21;
    for(unsigned int i=0; i<A_values_down_lp1.size(); i++) r+=A_values_down_lp1[i].A21;
    return r;
}

double Electron_Level::A21(double f_ul, int nu, int lu, int nl, int ll)                        // in 1/sec
{ 
    if(nu<nl)
    {
        cout << " Upper and lower levels are confused: nu: " << nu << " nl: " << nl << endl;
        exit(0);
    }
    
    return -2.0*FOURPI*const_PIe2_mec*f_ul/pow(lambda_ul(nu, nl), 2); 
}

double Electron_Level::B12(double f_lu, int nl, int nu)
{
    if(nl>nu)
    {
        cout << " Upper and lower levels are confused: nup: " << nu << " nlow: " << nl << endl;
        exit(0);
    }
    return FOURPI/(nu_ul(nu, nl)*const_h)*const_PIe2_mec*f_lu; 
}  

double Electron_Level::B21(double f_ul, int nu, int nl)
{ 
    if(nu<nl)
    {
        cout << " Upper and lower levels are confused: nup: " << nu << " nlow: " << nl << endl;
        exit(0);
    }
    
    return FOURPI/(nu_ul(nl, nu)*const_h)*const_PIe2_mec*f_ul; 
}

//======================================================================
// functions with fast access to transition data
//======================================================================
const Transition_Data_A& Electron_Level::Get_Data_A(const int &np, const int &lp) const
{
    if(np>=nn) return ZERO_Data;            // no transition to level itself
    if(lp<0 || lp>np-1) return ZERO_Data;   // level (np, lp) does not exist
    
    if(lp==ll-1) return A_values_down_lm1[np-ll];
    if(lp==ll+1) return A_values_down_lp1[np-(ll+2)];
    
    return ZERO_Data;
}

double Electron_Level::Get_A21(const int &np, const int &lp) const
{   
    if(np>=nn) return 0.0;            // no transition to level itself
    if(lp<0 || lp>np-1) return 0.0;   // level (np, lp) does not exist

    if(lp==ll-1) return A_values_down_lm1[np-ll].A21;
    if(lp==ll+1) return A_values_down_lp1[np-(ll+2)].A21;
        
    return 0.0;
}

double Electron_Level::Get_nu21(const int &np, const int &lp) const
{   
    if(np>=nn) return 0.0;            // no transition to level itself
    if(lp<0 || lp>np-1) return 0.0;   // level (np, lp) does not exist
    
    if(lp==ll-1) return A_values_down_lm1[np-ll].Dnu;
    if(lp==ll+1) return A_values_down_lp1[np-(ll+2)].Dnu;
    
    return 0.0;
}

double Electron_Level::Get_lambda21(const int &np, const int &lp) const
{   
    if(np>=nn) return 0.0;            // no transition to level itself
    if(lp<0 || lp>np-1) return 0.0;   // level (np, lp) does not exist
    
    if(lp==ll-1) return A_values_down_lm1[np-ll].lambda21;
    if(lp==ll+1) return A_values_down_lp1[np-(ll+2)].lambda21;
    
    return 0.0;
}

//===================================================================================
// Konstructors and Destructors for class: Atomic_Shell
//===================================================================================
void Atomic_Shell::init(int n, int nm, int ZZ, double NNp, int Rec_flag, int mflag)
{
    nn=n; nmax=nm;
    Z=ZZ;
    Np=NNp;
    mess_flag=mflag;
    
    create_Electron_Levels(Rec_flag);
}

Atomic_Shell::Atomic_Shell(int n, int nm, int ZZ, double NNp, int Rec_flag, int mflag)
{ init(n, nm, ZZ, NNp, Rec_flag, mflag); }

Atomic_Shell::~Atomic_Shell(){ Angular_Momentum_Level.clear(); }

//-----------------------------------------------------------------------------------
void Atomic_Shell::create_Electron_Levels(int Rec_flag)
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atomic_Shell::create_Electron_Levels: filling shell: " << nn <<"\n % " << endl;
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    // free memory (if necessary)
    Angular_Momentum_Level.clear();
    
    Electron_Level v;
    int l;
    // fill with empty electron-states
    for(l=0; l<nn; l++) Angular_Momentum_Level.push_back(v);
    // now initialize each state
    for(l=0; l<nn; l++) Angular_Momentum_Level[l].init(nn, l, nmax, Z, Np, Rec_flag, mess_flag);

    return;
}

void Atomic_Shell::display_general_data_of_level(unsigned int i)
{
    if(i>=(unsigned int)(nn) || Angular_Momentum_Level.size()==0)
    {
        cout << " This level does not exist inside shell " << nn << endl;
        return;
    }
    
    cout << "\n %==============================================================%" << endl;
    cout << " % Atomic_Shell::_general_data_of_level:\n % (" << nn << ", " << i << ")" << endl;
    cout << " %==============================================================%" << endl;
    cout << " Z: " << Z << " Np: " << Np << endl; 
    cout << " Dnu: " << Angular_Momentum_Level[i].Get_Dnu_1s() 
    << " Hz, DE: " << Angular_Momentum_Level[i].Get_DE_1s() << " eV"<< endl;
    cout << " A_tot: " << Angular_Momentum_Level[i].Get_A() 
    << " sec^-1, Gamma: " << Angular_Momentum_Level[i].Get_Gamma() << " sec^-1" << endl << endl; 
    
    return;
}

//===================================================================================
// Konstructors and Destructors for Atom
//===================================================================================
void Atom::init(int nS, int ZZ, double NNp, int Rec_flag, int mflag)
{
    nShells=nS;
    Z=ZZ;
    Np=NNp;
    mess_flag=mflag;
    
    fill_Level_Map();
    if(mess_flag>2) display_Level_Map();
    
    create_Shells(Rec_flag);
}

Atom:: Atom(){}

Atom:: Atom(int nS, int ZZ, double NNp, int Rec_flag, int mflag)
{ init(nS, ZZ, NNp, Rec_flag, mflag); }

Atom::~Atom()
{ 
    Level_Map.clear(); 
    Shell.clear();
}

//-----------------------------------------------------------------------------------
void Atom::create_Shells(int Rec_flag)
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atom::create_Shells:creating all the shells up to n=" << nShells << endl;
        cout << " %                     The hydrogenic atom has charge Z=" << Z << endl;        
        cout << " %                     Effective number of protons in the nucleus Np=" << Np << endl;        
        if(Rec_flag==0) 
            cout << " %                     No recombination rates will be setup" << "\n % " << endl;
        else if(Rec_flag==1) 
            cout << " %                     Storey & Hummer Recombination rates with interpolation will be setup" 
            << "\n % " << endl;
        else if(Rec_flag==2) 
            cout << " %                     Storey & Hummer Recombination rates will be setup" << "\n % " << endl;
        else if(Rec_flag==3) 
            cout << " %                     Karzas & Latter Recombination rates will be setup" << "\n % " << endl;
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    int m=mess_flag-1;
    if(mess_flag>=1) m++;
    if(mess_flag>=2) m++;
    
    // free memory (if necessary)
    Shell.clear();
    
    Atomic_Shell v;
    // fill with empty shells
    for(unsigned int n=0; n<=nShells; n++) Shell.push_back(v);
    // create each shells
    for(unsigned int n=1; n<=nShells; n++) Shell[n].init(n, nShells, Z, Np, Rec_flag, m); 
    
    return;
}

void Atom::display_general_data_of_Shell(unsigned int i)
{
    if(i>=(unsigned int)(nShells) || i==0 || Shell.size()==0)
    {
        cout << " This Shell does not exist inside Atom " << endl;
        return;
    }
    
    cout << " Function not specified yet..." << endl; 
    
    return;
}

//------------------------------------------------------------------------------------------------------
void Atom::fill_Level_Map()
{
    // free memory (if necessary)
    Level_Map.clear(); 
    
    // this function creates the map of indicies i-> (n,l)
    Level_nl v;
    
    for(unsigned int n=1; n<=nShells; n++) 
    {
        for(unsigned int l=0; l<n; l++)
        {
            v.n=n; v.l=l;       
            Level_Map.push_back(v);
        }
    }
    
    return;
}

void Atom::display_Level_Map()
{
    cout << "\n %==============================================================%" << endl;
    cout << " % Atom::display_Level_Map:" << endl; 
    cout << " %==============================================================%" << endl;
    
    for(unsigned int i=0; i< Level_Map.size(); i++)
        cout << " " << i << " --> (" << Level_Map[i].n << ", " << Level_Map[i].l <<")" << endl; 
    cout << endl;
    
    return;
}

const Electron_Level& Atom::Level(const unsigned int &i) const
{ 
    if(i> Get_total_number_of_Levels())
    {
        cout << " Atom::Level: you are trying to access a non-existing level: " << i 
        << " total number of levels: " << Get_total_number_of_Levels() << endl;
        exit(0);
    }

    return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l]; 
}

const Electron_Level& Atom::Level(const unsigned int &n, const unsigned int &l) const
{ 
    if(n==0 || n>nShells || l>=n)
    {
        cout << " Atom::Level: you are trying to access a non-existing level: " << n 
        << ", " << l << endl;
        exit(0);
    }
    
    return Shell[n].Angular_Momentum_Level[l]; 
}

Electron_Level& Atom::Level(const unsigned int &i) 
{ 
    if(i> Get_total_number_of_Levels())
    {
        cout << " Atom::Level: you are trying to access a non-existing level: " << i 
        << " total number of levels: " << Get_total_number_of_Levels() << endl;
        exit(0);
    }
    
    return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l]; 
}

Electron_Level& Atom::Level(const unsigned int &n, const unsigned int &l)
{ 
    if(n==0 || n>nShells || l>=n)
    {
        cout << " Atom::Level: you are trying to access a non-existing level: " << n 
        << ", " << l << endl;
        exit(0);
    }
    
    return Shell[n].Angular_Momentum_Level[l]; 
}

//===================================================================================
// Konstructors and Destructors for class: Gas_of_Atoms
//===================================================================================
void Gas_of_Atoms::init(int nS, int Z, double Np, int Rec_flag, int mflag)
{ 
    Atom::init(nS, Z, Np, Rec_flag, mflag);
    create_vectors_X(nS);
    
    //===========================================================
    // create the Voigt profiles for Ly-n 
    //===========================================================
    double A21, f=0.0, nu21, lam21, Gamma=0;
    //
    for(int k=2; k<=(int)min(Get_nShells(), 10); k++)
    {
        A21=Level(k, 1).Get_A21(1, 0);
        nu21=Level(k, 1).Get_nu21(1, 0);
        lam21=Level(k, 1).Get_lambda21(1, 0);
        Gamma=Level(k, 1).Get_A(); // in Voigt-profile Gamma == A_tot not A_tot/4 pi !!!
        //
        if(mflag>=1) cout << " Initializing profile for " << k << "P - 1S transition : " << nu21 << " " << lam21 << " " << A21 << " " << Gamma << " " << f << endl;
        // initialize Voigt-profiles for HI 1s-np series
        phi_HI_Lyn[k].Set_atomic_data_Gamma(nu21, lam21, A21, f, Gamma, 1.0);
    }
    
    return;
}

Gas_of_Atoms::Gas_of_Atoms():Atom(){}

Gas_of_Atoms::Gas_of_Atoms(int nS, int Z, double Np, int Rec_flag, int mflag):Atom()
{ init(nS, Z, Np, Rec_flag, mflag); }

Gas_of_Atoms::~Gas_of_Atoms()
{ 
    Xvec.clear();
    Ric_BB_vec.clear();
    Rci_BB_vec.clear();
    dRci_dTe_BB_vec.clear();
}

const Voigtprofile_Dawson& Gas_of_Atoms::HI_Lyn_profile(int n) const 
{ 
    if(n>1 && n<=(int)min(Get_nShells(), 10)) return phi_HI_Lyn[n]; 
    
    return phi_HI_Lyn[0]; 
}

Voigtprofile_Dawson& Gas_of_Atoms::HI_Lyn_profile(int n)
{ 
    if(n>1 && n<=(int)min(Get_nShells(), 10)) return phi_HI_Lyn[n]; 
    
    return phi_HI_Lyn[0]; 
}

void Gas_of_Atoms::create_vectors_X(int imax)
{
    vector<double> v;
    
    Xvec.clear();
    Ric_BB_vec.clear();
    Rci_BB_vec.clear();
    dRci_dTe_BB_vec.clear();
    
    for(int i=0; i<=imax; i++)
    {
        v.assign(i, 0.0);
        Xvec.push_back(v);
        Ric_BB_vec.push_back(v);
        Rci_BB_vec.push_back(v);
        dRci_dTe_BB_vec.push_back(v);
        v.clear();
    }
    
    return;
}
//-----------------------------------------------------------------------------------
double Gas_of_Atoms::Get_A(unsigned int n, unsigned int l, unsigned int np, unsigned int lp) const
{ return Level(n, l).Get_A21(np, lp); }

double Gas_of_Atoms::Get_population_of_level(const unsigned int &i) const
{ 
    if(i>= Get_total_number_of_Levels())
    {
        cout << " Gas_of_Atoms::Get_population_of_level: you are trying to access a non-existing level: " << i 
        << " total number of levels: " << Get_total_number_of_Levels() << endl;
        exit(0);
    }
    
    return Xvec[Get_n_of_Level(i)][Get_l_of_Level(i)]; 
}

void Gas_of_Atoms::Set_population_of_level(const unsigned int &i, const double &XX)
{ 
    if(i>= Get_total_number_of_Levels())
    {
        cout << " Gas_of_Atoms::Set_population_of_level: you are trying to access a non-existing level: " << i 
        << " total number of levels: " << Get_total_number_of_Levels() << endl;
        exit(0);
    }
    
    Xvec[Get_n_of_Level(i)][Get_l_of_Level(i)]=XX;
    return;
}

double Gas_of_Atoms::Get_population_of_level(const unsigned int &n, const unsigned int &l) const
{ 
    if(n==0 || n>Xvec.size() || l>=Xvec[n].size())
    {
        cout << " Gas_of_Atoms::Get_population_of_level: you are trying to access a non-existing level: " 
        << "( " << n << ", " << l << " )" << " nmax= " << Get_nShells() << endl;
        exit(0);
    }
    
    return Xvec[n][l];
}

void Gas_of_Atoms::Set_population_of_level(const unsigned int &n, const unsigned int &l, const double &XX)
{ 
    if(n==0 || n>Xvec.size() || l>=Xvec[n].size())
    {
        cout << " Gas_of_Atoms::Set_population_of_level: you are trying to access a non-existing level: " 
        << "( " << n << ", " << l << " )" << " nmax= " << Get_nShells() << endl;
        exit(0);
    }
    
    Xvec[n][l]=XX;
    return;
}

//-------------------------------------------------------------------------------------------------------------------
double Gas_of_Atoms::X_tot() const
{
    double r=0.0;
    // smallest numbers first
    for(int n=Get_nShells(); n>=1; n--)
        for(int l=0; l<n; l++) r+=Xvec[n][l];
    
    return r; 
}

double Gas_of_Atoms::X_tot(int nlow) const
{
    if(nlow<1) return X_tot();
    
    double r=0.0;
    // smallest numbers first
    for(int n=Get_nShells(); n>=nlow; n--)
        for(int l=0; l<n; l++) r+=Xvec[n][l];
    
    return r; 
}

double Gas_of_Atoms::R_rec_tot_BB(int nlow, double T, double rho)       // total recombination rate from nlow upwards
{
    if(nlow<1) nlow=1;
    double r=0.0;
    
    // smallest numbers first
    for(int n=Get_nShells(); n>=nlow; n--)
        for(int l=n-1; l>=0; l--) r+=Shell[n].Angular_Momentum_Level[l].R_c_nl(T, rho);
    
    return r;
}

double Gas_of_Atoms::R_phot_tot_BB(int nlow, double T)                  // total photoionization rate from nlow upwards
{
    if(nlow<1) nlow=1;
    double r=0.0;
    
    // smallest numbers first
    for(int n=Get_nShells(); n>=nlow; n--)
        for(int l=n-1; l>=0; l--) r+=Shell[n].Angular_Momentum_Level[l].R_nl_c(T)*Xvec[n][l];
    
    return r;
}

//=================================================================================
// routines for numerical applications
//=================================================================================
void Gas_of_Atoms::clear_Interaction_w_photons()
{ 
    for(int n=1; n<=Get_nShells(); n++) 
        for(int l=0; l<n; l++) Shell[n].Angular_Momentum_Level[l].clear_Interaction_w_photons();
    
    return;
}

void Gas_of_Atoms::update_Ric_BB(double TT)
{
    update_Ric_BB(1, TT);
    return;
}

void Gas_of_Atoms::update_Rci_BB(double TT, double rho)
{
    update_Rci_BB(1, TT, rho);
    return;
}

void Gas_of_Atoms::update_dRci_dTe_BB(double TT, double rho)
{
    update_dRci_dTe_BB(1, TT, rho);
    return;
}

void Gas_of_Atoms::update_Ric_BB(int nm, double TT)
{
    T=TT;
    for(int n=nm; n<=Get_nShells(); n++) 
        for(int l=0; l<n; l++) Ric_BB_vec[n][l]=Shell[n].Angular_Momentum_Level[l].R_nl_c(T);

    return;
}

void Gas_of_Atoms::update_Rci_BB(int nm, double TT, double rho)
{
    T=TT;
    for(int n=nm; n<=Get_nShells(); n++) 
        for(int l=0; l<n; l++) Rci_BB_vec[n][l]=Shell[n].Angular_Momentum_Level[l].R_c_nl(T, rho);

    return;
}

void Gas_of_Atoms::update_dRci_dTe_BB(int nm, double TT, double rho)
{
    T=TT;
    for(int n=nm; n<=Get_nShells(); n++) 
            for(int l=0; l<n; l++) dRci_dTe_BB_vec[n][l]=Shell[n].Angular_Momentum_Level[l].dR_c_nl_dTe(T, rho);

    return;
}

void Gas_of_Atoms::Set_Ric_rate(const unsigned int &n, const unsigned int &l, const double &Ric)
{ 
    if(n==0 || n>Ric_BB_vec.size() || l>=Ric_BB_vec[n].size())
    {
        cout << " Gas_of_Atoms::Set_Ric_rate: you are trying to access a non-existing level: " 
        << "( " << n << ", " << l << " )" << endl;
        exit(0);
    }
    
    Ric_BB_vec[n][l]=Ric;
    return;
}

void Gas_of_Atoms::Set_Rci_rate(const unsigned int &n, const unsigned int &l, const double &Rci) 
{ 
    if(n==0 || n>Rci_BB_vec.size() || l>=Rci_BB_vec[n].size())
    {
        cout << " Gas_of_Atoms::Set_Rci_rate: you are trying to access a non-existing level: " 
        << "( " << n << ", " << l << " )" << endl;
        exit(0);
    }
    
    Rci_BB_vec[n][l]=Rci;
    return;
}

void Gas_of_Atoms::Set_dRci_dTe_rate(const unsigned int &n, const unsigned int &l, const double &dRci_dTe) 
{ 
    if(n==0 || n>dRci_dTe_BB_vec.size() || l>=dRci_dTe_BB_vec[n].size())
    {
        cout << " Gas_of_Atoms::Set_dRci_dTe_rate: you are trying to access a non-existing level: " 
        << "( " << n << ", " << l << " )" << endl;
        exit(0);
    }
    
    dRci_dTe_BB_vec[n][l]=dRci_dTe;
    return;
}

//----------------------------------------------------------------
double Gas_of_Atoms::Ric_BB(unsigned int i) 
{
    if(i>= Get_total_number_of_Levels())
    {
        cout << " Gas_of_Atoms::Ric_BB: you are trying to access a non-existing level: " << i 
        << " total number of levels: " << Get_total_number_of_Levels() << endl;
        exit(0);
    }
    
    return Ric_BB_vec[Get_n_of_Level(i)][Get_l_of_Level(i)]; 
}

double Gas_of_Atoms::Rci_BB(unsigned int i)
{
    if(i>= Get_total_number_of_Levels())
    {
        cout << " Gas_of_Atoms::Rci_BB: you are trying to access a non-existing level: " << i 
        << " total number of levels: " << Get_total_number_of_Levels() << endl;
        exit(0);
    }
    
    return Rci_BB_vec[Get_n_of_Level(i)][Get_l_of_Level(i)]; 
}

double Gas_of_Atoms::dRci_dTe_BB(unsigned int i)
{
    if(i>= Get_total_number_of_Levels())
    {
        cout << " Gas_of_Atoms::dRci_dTe_BB: you are trying to access a non-existing level: " << i 
        << " total number of levels: " << Get_total_number_of_Levels() << endl;
        exit(0);
    }
    
    return dRci_dTe_BB_vec[Get_n_of_Level(i)][Get_l_of_Level(i)]; 
}

double Gas_of_Atoms::Ric_BB(unsigned int n, unsigned int l)
{
    if(n==0 || n>Ric_BB_vec.size() || l>=Ric_BB_vec[n].size())
    {
        cout << " Gas_of_Atoms::Ric_BB: you are trying to access a non-existing level: " 
        << "( " << n << ", " << l << " )" << endl;
        exit(0);
    }
    
    return Ric_BB_vec[n][l]; 
} 

double Gas_of_Atoms::Rci_BB(unsigned int n, unsigned int l) 
{
    if(n==0 || n>Rci_BB_vec.size() || l>=Rci_BB_vec[n].size())
    {
        cout << " Gas_of_Atoms::Rci_BB: you are trying to access a non-existing level: " 
        << "( " << n << ", " << l << " )" << endl;
        exit(0);
    }
    
    return Rci_BB_vec[n][l]; 
} 

double Gas_of_Atoms::dRci_dTe_BB(unsigned int n, unsigned int l) 
{
    if(n==0 || n>dRci_dTe_BB_vec.size() || l>=dRci_dTe_BB_vec[n].size())
    {
        cout << " Gas_of_Atoms::dRci_dTe_BB: you are trying to access a non-existing level: " 
        << "( " << n << ", " << l << " )" << endl;
        exit(0);
    }
    
    return dRci_dTe_BB_vec[n][l]; 
} 

double Gas_of_Atoms::R_rec_tot_BB(int nlow)                   // total recombination rate from nlow upwards
{
    if(nlow<1) nlow=1;
    double r=0.0;
    
    for(int n=nlow; n<=Get_nShells(); n++) 
        for(int l=0; l<n; l++) r+=Rci_BB_vec[n][l];
    
    return r;
}

double Gas_of_Atoms::R_phot_tot_BB(int nlow)                  // total photoionization rate from nlow upwards
{
    if(nlow<1) nlow=1;
    double r=0.0;
    
    for(int n=nlow; n<=Get_nShells(); n++)
        for(int l=0; l<n; l++) r+=Ric_BB_vec[n][l]*Xvec[n][l];
    
    return r;
}

//==============================================================================
// Saha-relations with continuum
//==============================================================================
// f(T) == (Ni/[Ne Nc])_LTE
//
double Gas_of_Atoms::Ni_NeNc_LTE(unsigned int i, double TM)
{ 
    double gi=2.0*(2.0*Get_l_of_Level(i)+1.0), gc=1.0;
    return gi/2.0/gc*pow(const_lambdac, 3)*pow(2.0*PI*const_kb_mec2*TM*Level(i).Get_mu_red(), -1.5)*exp(Level(i).Get_E_ion_ergs()/const_kB/TM );
}

double Gas_of_Atoms::Ni_NeNc_LTE(unsigned int n, unsigned int l, double TM)
{ return Ni_NeNc_LTE(Get_Level_index(n, l), TM); }

double Gas_of_Atoms::Xi_Saha(unsigned int i, double Xe, double Xc, double NH, double TM, double z)
{ return Xe*Xc*NH*Ni_NeNc_LTE(i, TM); }

double Gas_of_Atoms::Xi_Saha(unsigned int n, unsigned int l, double Xe, double Xc, double NH, double TM, double z)
{ return Xi_Saha(Get_Level_index(n, l), Xe, Xc, NH, T, z); }

double Gas_of_Atoms::Ni_Saha(unsigned int i, double Ne, double Nc, double TM)
{ return Ne*Nc*Ni_NeNc_LTE(i, TM); }

double Gas_of_Atoms::Ni_Saha(unsigned int n, unsigned int l, double Ne, double Nc, double TM)
{ return Ni_Saha(Get_Level_index(n, l), Ne, Nc, TM); }

