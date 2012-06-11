//===========================================================================================================
// CosmoRec version
//===========================================================================================================
string CosmoRec_version="v1.4.1 beta";

//===========================================================================================================
// setting for two-gamma and Raman corrections
//===========================================================================================================
int ntg_max=0;
int nR_max=0;
bool hydrogen_and_helium_atom_are_set=0;

//===========================================================================================================
// for data output
//===========================================================================================================
string path, detailedname;
string addname;
ofstream *Xfiles=new ofstream[6];
vector<vector<double> > output_CosmoRec;
vector<vector<double> > pass_on_the_Solution_CosmoRec;

//===========================================================================================================
// level of verbosity
//===========================================================================================================
int show_CosmoRec_mess=0;

//===========================================================================================================
// output files
//===========================================================================================================
bool write_populations=0;
bool write_Xe_Te_sol=1;

//===========================================================================================================
// path to the effective rate database
//===========================================================================================================
string Rec_database_path=COSMORECDIR+COSMOREC_BASE_PATH;

//===========================================================================================================
// to include HeI Lya absorption by hydrogen
//===========================================================================================================
Photoionization_cross_section_SH H1s(1, 0, 1, 1, 1);
Photoionization_Lyc Lyc1s(1);

//===========================================================================================================
// below z=200 extend the solution with Recfast system
//===========================================================================================================
int flag_compute_Recfast_part=0;
double ze_Recfast=0.0;
int nz_Recfast_part=100;

//==========================================================================================================
// Variables for module that accounts for the diffusion corrections to escape probabilities
//==========================================================================================================
bool Diffusion_correction=0;
bool Diffusion_correction_is_on=0;
int Diffusion_flag;
//
int Diff_iteration_max=1, Diff_iteration_min=0;
int iteration_count=Diff_iteration_min;

//==========================================================================================================
// Variables for dark matter annihilation module
//==========================================================================================================
bool DM_annihilation=0;
double fDM_CosmoRec=0.0;

//===========================================================================================================
// memory for CosmoRec-ODE solver
//===========================================================================================================
ODE_solver_Rec::ODE_solver_Solution Sz;
ODE_solver_Rec::ODE_solver_accuracies tols;
ODE_solver_Rec::ODE_Solver_data ODE_Solver_info_CR;

//===========================================================================================================
// for passing on the data about cosmology
//===========================================================================================================
struct Parameters_form_initfile {
    int nz;
    double YP;
    double fHe;
    double To;
    double Nnu;
    double Omega0;    /* Omega matter    */
    double OmegaB;    /* Omega Baryons   */
    double OmegaK;    /* Omega Curvature */
    double OmegaL;    /* Omega Lambda    */
    double zstart;
    double zend;
    double F;
    double H0;
    double h100;
};

Parameters_form_initfile parameters;

//===========================================================================================================
// Level I: Xi=Ni/NHtot 
//===========================================================================================================
struct Data_Level_I {
    int neq;
    int nHIeq;
    int nHeIeq;
    int index_HI;
    int index_HeI;
    double *X;
    double *g;
};

//===========================================================================================================
// Level III: 
//===========================================================================================================
struct Data_Level_III {
    int neq;
    int nHIeq;
    int nHeIeq;
    int index_HI;
    int index_HeI;
    double *y;
    double *g;
};

Data_Level_I Level_I; 
Data_Level_III Level_III; 

Data_Level_I Level_I_temp; 
Data_Level_III Level_III_temp; 
Data_Level_III Level_III_temp_JAC; 

//===========================================================================================================
double *zarr;               // redshift array
double *f0;                 // rescaling factors Yi=Xi/fi

//===========================================================================================================
// for Cosmology object
//===========================================================================================================
Cosmos cosmos;

