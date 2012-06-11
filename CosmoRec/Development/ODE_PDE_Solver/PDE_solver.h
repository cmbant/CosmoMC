//============================================================================================================
//
// Author Jens Chluba July 2010
// last modification: Feb 2011; added possibility to take integral part over +/-2 off-diags into account
//
//============================================================================================================

#ifndef PDE_SOLVER_H
#define PDE_SOLVER_H

#include <vector>

using namespace std;

struct Lagrange_interpolation_coefficients
{
    vector<vector<double> > dli_dx;
    vector<vector<double> > d2li_d2x;
};

struct PDE_Stepper_Data
{
    vector<double> *ptr;    
    //
    vector<double> Ai;
    vector<double> Bi;
    vector<double> Ci;
    vector<double> Di;
    // used as workspace
    vector<double> Ui;
    vector<double> Vi;
    vector<double> bi;
    vector<double> zi; 
    //
    vector<double> Aip;
    vector<double> Bip;
    vector<double> Cip;
    vector<double> Dip;
    //
    vector<double> *Ai_ptr; 
    vector<double> *Bi_ptr; 
    vector<double> *Ci_ptr; 
    vector<double> *Di_ptr; 
    //
    vector<double> *Aip_ptr;    
    vector<double> *Bip_ptr;    
    vector<double> *Cip_ptr;    
    vector<double> *Dip_ptr;    
    //
    vector<vector<double> > lambda; // temporary data for banded solve
    vector<vector<double> > lambdap;

    Lagrange_interpolation_coefficients LG;
};

void init_PDE_Stepper_Data(PDE_Stepper_Data &PDE_D, int npts);
void reset_PDE_solver_variables();

//============================================================================================================
void setup_Lagrange_interpolation_coefficients_O1(PDE_Stepper_Data &PDE_D, vector<double> &xarr);
void setup_Lagrange_interpolation_coefficients_O2(PDE_Stepper_Data &PDE_D, vector<double> &xarr);

void Step_PDE_O1(double zs, double ze, const vector<double> &xi, vector<double> &yi, 
                 double yi_low, double yi_up, PDE_Stepper_Data &PDE_D,
                 void (*func)(double z, const vector<double> &xi, vector<double> &Ai, 
                              vector<double> &Bi, vector<double> &Ci, vector<double> &Di));

void Step_PDE_O2t(double theta, double zs, double ze, const vector<double> &xi, vector<double> &yi, 
                  double yi_low, double yi_up, PDE_Stepper_Data &PDE_D,
                  void (*func)(double z, const vector<double> &xi, vector<double> &Ai, 
                               vector<double> &Bi, vector<double> &Ci, vector<double> &Di));

#endif
//============================================================================================================
//============================================================================================================

