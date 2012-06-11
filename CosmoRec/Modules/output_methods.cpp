//--------------------------------------------------------------------------------------------
void set_detailed_name()
{
    detailedname="CosmoRec.";
    detailedname+="HI.nS_eff_"+int_to_string(nS_effective)+"."; 
    //
    detailedname+="nH_"+int_to_string(nShells)+".";
    detailedname+="nHe_"+int_to_string(nShellsHeI)+".";
    //
    if(flag_HI_absorption==1 && _HI_abs_appr_flag==0) detailedname+="H_abs.";
    if(flag_HI_absorption==1 && _HI_abs_appr_flag==1) detailedname+="H_abs.fcorr.";
    //
    if(flag_spin_forbidden==0) detailedname+="no_spin_forbidden.";
    if(set_HeI_Quadrupole_lines) detailedname+="Q_lines.";
    //
    if(HeISTfeedback==1) detailedname+="HeI_feedback_nf_"+int_to_string(nHeFeedback)+".";
    if(_HeI_feedback_w_HI_abs==1) detailedname+="w_HI_abs.";
    //
    if(DM_annihilation) detailedname+="DM_annihil.";
        
    if(show_CosmoRec_mess>=0) 
        cout << " set_detailed_name: detailedname: " << detailedname << endl;
    
    return; 
}

//--------------------------------------------------------------------------------------------
void open_files(ofstream *Xfiles, string path, int fresh=1)
{
    ios::openmode wrmode;
    if(fresh==1) wrmode=ios::out;
    else wrmode=ios::app;
    
    string prename=path+detailedname;
    if(write_Xe_Te_sol)
    {
        if(iteration_count==0) 
        {
            Xfiles[0].open((prename+"X_Recfast"+addname).c_str(), wrmode);
            Xfiles[0] << "# Recfast++ output. The columns are z, Xe=Ne/NH, Te\n#" << endl;
        }
        
        Xfiles[1].open((prename+"Xe"+addname).c_str(), wrmode);
        Xfiles[1] << "# CosmoRec "+CosmoRec_version+" output. The columns are z, Xe=Ne/NH, Te\n#" 
                  << endl;        
    }
    
    //  Xfiles[2].open((prename+"rho"+addname).c_str(), wrmode); 
    //
    if(write_populations) Xfiles[3].open((prename+"HI.X"+addname).c_str(), wrmode); 
//  Xfiles[4].open((prename+"HeI.X"+addname).c_str(), wrmode); 
//  Xfiles[5].open((prename+"HI.HeI.1s"+addname).c_str(), wrmode); 

    if(write_Xe_Te_sol && iteration_count==0) Xfiles[0].precision(8);
    if(write_Xe_Te_sol) Xfiles[1].precision(8);
    if(write_populations) Xfiles[3].precision(12);
    
    return;
}

void close_files(ofstream *Xfiles)
{
    if(write_Xe_Te_sol && iteration_count==0) Xfiles[0].close();
    if(write_Xe_Te_sol) Xfiles[1].close();
    if(write_populations) Xfiles[3].close();
    
    return;
}

//--------------------------------------------------------------------------------------------
void write_HI_HeI_1s_array_to_file(ofstream &ofile, double zval)
{
    ofile << zval;
    ofile << " " << Hydrogen_Atoms.X(1, 0);
    ofile << " " << HeI_Atoms.Xi(0) << endl;
    
    return;
}

void write_HI_X_array_to_file(ofstream &ofile, double zval, double *X, int nHIeq, int neq)
{
    ofile << zval;
    // Xe & HI Xi
    for(int l=0; l<Level_I.nHIeq+1; l++) ofile << " " << X[l];
    // rho
    ofile << " " << X[neq-1];
    // Xp, XeHeII
    ofile << " " << 1.0-Hydrogen_Atoms.X_tot();
    ofile << " " << parameters.fHe-HeI_Atoms.X_tot() << endl;
    
    return;
}

//--------------------------------------------------------------------------------------------
void write_array_to_file(ofstream &ofile, double zval, double *X, int neq)
{
    ofile << zval;
    // write all the variables
    for(int l=0; l<neq; l++) ofile << " " << X[l];
    ofile << endl;
    
    return;
}

//--------------------------------------------------------------------------------------------
void write_X_Recfast_z_to_file(ofstream &ofile, double zval, Cosmos &cosmos)
{
    ofile << zval << " " 
          << cosmos.Xe_Seager(zval) << " " 
          << cosmos.Te(zval) << " "
//        << cosmos.dXe_dz_Seager(zval) << " " << cosmos.X1s(zval) << " " 
//        << cosmos.dX1s_dz(zval) << " " << cosmos.Te_Tg(zval) << " " << cosmos.H(zval)
          << endl;
    
    return;
}

//--------------------------------------------------------------------------------------------
void write_Xe_Te_to_file(ofstream &ofile, double zval, double Xe, double Te)
{
    ofile << zval << " " << Xe << " " << Te << endl;
    
    return;
}

//--------------------------------------------------------------------------------------------
void write_files(double z, ofstream *Xfiles)
{
    if(write_Xe_Te_sol)
    {
        //==========================================================
        // only Xe and Te
        //==========================================================
        if(iteration_count==0) write_X_Recfast_z_to_file(Xfiles[0], z, cosmos);
        write_Xe_Te_to_file(Xfiles[1], z, Level_I.X[0], Level_I.X[Level_I.neq-1]*cosmos.TCMB(z));
    }
    //==============================================================
    // only rho
    //==============================================================
//  write_array_to_file(Xfiles[2], z, &Level_I.X[Level_I.neq-1], 1);
 
    //==============================================================
    // Xe & hydrogen populations
    //==============================================================
    if(write_populations) write_HI_X_array_to_file(Xfiles[3], z, Level_I.X, Level_I.nHIeq, Level_I.neq);

    //==============================================================
    // Helium populations
    //==============================================================
//  write_array_to_file(Xfiles[4], z, &Level_I.X[Level_I.nHIeq+1], Level_I.neq-(Level_I.nHIeq+1)-1);
    //
//  write_HI_HeI_1s_array_to_file(Xfiles[5], z);    
    
    return;
}


