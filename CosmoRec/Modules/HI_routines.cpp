//========================================================================
// data for hydrogen atom
//========================================================================

Gas_of_Atoms *HAtoms=NULL;
#define Hydrogen_Atoms (*HAtoms)

//========================================================================
// setting for hydrogen effective rates
//========================================================================
int nS_effective=500;                    // can be changed with start-file
int Lyn_max_parameter=40;

string effective_rate_path_arr[3]=
{
    Rec_database_path+"/Effective_Rates.HI/Effective_Rate_Tables.nS_3/",
    Rec_database_path+"/Effective_Rates.HI/Effective_Rate_Tables.nS_4/",
    Rec_database_path+"/Effective_Rates.HI/Effective_Rate_Tables.nS_8.np_10/"
};

string effective_rate_path;

//========================================================================
// vectors for effective rates
//========================================================================
vector<double> Ai_df_dx, Bi_df_dx;
vector<vector<double> > RijVec_df_dx;


//========================================================================
// different parameters
//========================================================================
const double zref=1100.0;
const double zcrit=1000.0;

int nShells=5;

//========================================================================
// setup routines
//========================================================================
void setup_Hydrogen_atom()
{
    print_message(" entering the setup for neutral hydrogen ");
    
    HAtoms=new Gas_of_Atoms(nShells, 1, 1.0, 0, -2);

    
    //=====================================================================
    // allocate memory to pass on the solution to the PDE solver
    //=====================================================================
    pass_on_the_Solution_CosmoRec.clear();
    vector<double> dum(2+Hydrogen_Atoms.Get_total_number_of_Levels()+3);
    for(int k=0; k<parameters.nz; k++) 
        pass_on_the_Solution_CosmoRec.push_back(dum);
    
    //=====================================================================
    // Load effective rates
    //=====================================================================
    if(nShells==3) effective_rate_path=effective_rate_path_arr[0];
    else if(nShells==4) effective_rate_path=effective_rate_path_arr[1];
    else if(nShells==10) effective_rate_path=effective_rate_path_arr[2];
    else
    { 
        cerr << " setup_Hydrogen_atom :: please choose nSHI = 3, 4 or 10! Exiting... \n" << endl; 
        exit(0); 
    }
    load_rates(effective_rate_path, nS_effective, Hydrogen_Atoms);
    
    //=====================================================================
    // memory for effective rates
    //=====================================================================
    Ai_df_dx.resize(get_number_of_resolved_levels());
    Bi_df_dx.resize(get_number_of_resolved_levels());
    RijVec_df_dx.clear();
    for(int k=0; k<(int)Ai_df_dx.size(); k++) RijVec_df_dx.push_back(Ai_df_dx);

    return;
}

//========================================================================
void Set_Hydrogen_Levels_to_Saha(double z)
{
    Level_I.X[0]=min(cosmos.Xe_Seager(z), 1.0+parameters.fHe);// electrons
    double Xe=Level_I.X[0];
    double Xp=min(cosmos.Xp(z), 1.0);
    
    if(show_CosmoRec_mess>=1) cout << " Xe= " << Xe << " Xp= " << Xp << endl;
    
    // follow all sub-level 
    for(long int i=0; i<Level_I.nHIeq; i++)
    {
        Level_I.X[i+1]=Hydrogen_Atoms.Xi_Saha(i, Xe, Xp, cosmos.NH(z), cosmos.Te(z), z);        
        
        Hydrogen_Atoms.Set_population_of_level(i, Level_I.X[i+1]);
        
        if(show_CosmoRec_mess>=2) cout << " H-Level: " << i << " " << Level_I.X[i+1] << endl;
    }
    Level_I.X[Level_I.neq-1]=1.0;

    if(show_CosmoRec_mess>=2)
        cout << " total population: " << Hydrogen_Atoms.X_tot() << endl << endl;
    
    return;
}

