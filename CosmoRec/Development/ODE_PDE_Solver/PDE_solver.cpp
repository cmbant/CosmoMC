//=====================================================================================
//
// Author Jens Chluba July 2010
//
//=====================================================================================

//=====================================================================================
//
// Purpose: Solve linear (!) parabolic differential equation of the form
//
// du/dt = A(x,t) d^2u/dx^2 + B(x,t) du/dx + C(x,t) u + D(x,t)
//
// The mesh of grid-points given by xi[...] is not assumed to be uniform.
//=====================================================================================

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

#include "PDE_solver.h"
#include "routines.h"

using namespace std;

//=====================================================================================
//
// global structure
//
//=====================================================================================
struct previous_solution_info
{
    vector<double> Fsol;
    bool is_set_up;
    double zsol;
};

bool PDE_solver_initial_call_of_solver=0;

//=====================================================================================
void init_PDE_Stepper_Data(PDE_Stepper_Data &PDE_D, int npts)
{
    PDE_D.ptr=NULL;
    //
    PDE_D.Ai.resize(npts);
    PDE_D.Bi.resize(npts);
    PDE_D.Ci.resize(npts);
    PDE_D.Di.resize(npts);
    //
    PDE_D.Ai_ptr=&PDE_D.Ai;
    PDE_D.Bi_ptr=&PDE_D.Bi;
    PDE_D.Ci_ptr=&PDE_D.Ci;
    PDE_D.Di_ptr=&PDE_D.Di;
    //
    PDE_D.Ui.resize(npts);
    PDE_D.Vi.resize(npts);
    PDE_D.bi.resize(npts);
    PDE_D.zi.resize(npts);
    //
    PDE_D.Aip.resize(npts);
    PDE_D.Bip.resize(npts);
    PDE_D.Cip.resize(npts);
    PDE_D.Dip.resize(npts);
    //
    PDE_D.Aip_ptr=&PDE_D.Aip;
    PDE_D.Bip_ptr=&PDE_D.Bip;
    PDE_D.Cip_ptr=&PDE_D.Cip;
    PDE_D.Dip_ptr=&PDE_D.Dip;
        
    PDE_D.lambda. resize(npts, vector<double> (5));
    PDE_D.lambdap.resize(npts, vector<double> (5));

    return;
}

//=====================================================================================
//
// reset solver
//
//=====================================================================================
void reset_PDE_solver_variables()
{
    PDE_solver_initial_call_of_solver=0;
    return;
}

//-------------------------------------------------------------------------------------
// 1. derivative
//-------------------------------------------------------------------------------------
double dli_dx_xj(int i, int j, double *xarr, int np)
{
    double dl_dx=0.0;
    
    if(i==j) 
    {
        for(int k=0; k<i; k++)    dl_dx+=1.0/(xarr[i]-xarr[k]);
        for(int k=i+1; k<np; k++) dl_dx+=1.0/(xarr[i]-xarr[k]);
    }
    else if(i<j) 
    {
        dl_dx=1.0/(xarr[i]-xarr[j]);
        for(int k=0; k<i; k++)    dl_dx*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
        for(int k=i+1; k<j; k++)  dl_dx*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
        for(int k=j+1; k<np; k++) dl_dx*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
    }
    else 
    {
        dl_dx=1.0/(xarr[i]-xarr[j]);
        for(int k=0; k<j; k++)    dl_dx*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
        for(int k=j+1; k<i; k++)  dl_dx*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
        for(int k=i+1; k<np; k++) dl_dx*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
    }
    
    return dl_dx;
}

//-------------------------------------------------------------------------------------
// 2. derivative
//-------------------------------------------------------------------------------------
double d2li_d2x_xj(int i, int j, double *xarr, int np)
{
    double d2l_d2x=0.0, dum;
    
    if(i==j) 
    {
        for(int k=0; k<i; k++)
        {
            dum=0.0;
            
            for(int m=0; m<k; m++) dum+=1.0/(xarr[i]-xarr[m]);
            for(int m=k+1; m<i; m++) dum+=1.0/(xarr[i]-xarr[m]);
            for(int m=i+1; m<np; m++) dum+=1.0/(xarr[i]-xarr[m]);
            
            d2l_d2x+=1.0/(xarr[i]-xarr[k])*dum;
        }
        
        for(int k=i+1; k<np; k++)
        {
            dum=0.0;
            
            for(int m=0; m<i; m++) dum+=1.0/(xarr[i]-xarr[m]);
            for(int m=i+1; m<k; m++) dum+=1.0/(xarr[i]-xarr[m]);
            for(int m=k+1; m<np; m++) dum+=1.0/(xarr[i]-xarr[m]);
            
            d2l_d2x+=1.0/(xarr[i]-xarr[k])*dum;
        }
    }
    else 
    {
        for(int m=0; m<np; m++) 
        {
            if(m!=i && m!=j)
            {
                dum=1.0;
                for(int k=0; k<np; k++) 
                    if(k!=i && k!=j && k!=m) dum*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
                d2l_d2x+=dum/(xarr[i]-xarr[m]);
            }
        }
        
        d2l_d2x*=2.0/(xarr[i]-xarr[j]);
    }
    
    return d2l_d2x;
}

//=====================================================================================
//
// 1. order scheme
//
//=====================================================================================

//-------------------------------------------------------------------------------------
// precompute coefficients
//-------------------------------------------------------------------------------------
void setup_Lagrange_interpolation_coefficients_O1(PDE_Stepper_Data &PDE_D, vector<double> &xarr)
{
    int np=xarr.size();
    
    cout << "\n setup_Lagrange_interpolation_coefficients_O1:: setting up grid " << endl;
    
    PDE_D.LG.dli_dx.clear();
    PDE_D.LG.d2li_d2x.clear();
    
    //====================================================================
    // allocate memory
    //====================================================================
    vector<double> dum(3, 0.0);
    for(int k=0; k<np; k++)
    {
        PDE_D.LG.dli_dx  .push_back(dum);
        PDE_D.LG.d2li_d2x.push_back(dum);
    }
    
    //====================================================================
    // all intereor points
    //====================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int i=1; i<np-1; i++)
        for(int j=0; j<3; j++)
        {
            PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 1, &xarr[i-1], 3);
            PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 1, &xarr[i-1], 3);
        }
    
    cout << " setup_Lagrange_interpolation_coefficients_O1:: done " << endl;

    return;
}


//-------------------------------------------------------------------------------------
// compute lambda_i
//-------------------------------------------------------------------------------------
void compute_lambda_i_01(double dz, int i, vector<double> &Ai, 
                         vector<double> &Bi, vector<double> &Ci, 
                         PDE_Stepper_Data &PDE_D, double *lambda)
{       
    for(int j=0; j<3; j++) lambda[j]=PDE_D.LG.d2li_d2x[i][j]*Ai[i] +PDE_D.LG.dli_dx[i][j]*Bi[i];
    lambda[1]+=Ci[i];
    for(int j=0; j<3; j++) lambda[j]*=-dz;
    lambda[1]+=1.0;
}

//-------------------------------------------------------------------------------------
// compute solution
//-------------------------------------------------------------------------------------
void Step_PDE_O1(double zs, double ze, const vector<double> &xi, vector<double> &yi, 
                 double yi_low, double yi_up, PDE_Stepper_Data &PDE_D,
                 void (*func)(double z, const vector<double> &xi, vector<double> &Ai, 
                              vector<double> &Bi, vector<double> &Ci, vector<double> &Di))
{
    double dz=ze-zs;
    int npxi=xi.size();
    
    //=============================================
    // evaluate coefficients of PDE at future time
    //=============================================
    func(ze, xi, PDE_D.Ai, PDE_D.Bi, PDE_D.Ci, PDE_D.Di);
    
    //=============================================
    // lower & upper boundary condition
    //=============================================
    yi[0]=yi_low;
    yi[npxi-1]=yi_up;
    
    //=============================================
    // compute Ui* and bi*
    //=============================================
    double lambda[3], dum=0.0; lambda[2]=0.0;
    //=============================================
    PDE_D.Ui[0]=0.0;
    PDE_D.bi[0]=yi[0];  
    //=============================================
    for(int i=1; i<npxi-1; i++)
    {
        compute_lambda_i_01(dz, i, PDE_D.Ai, PDE_D.Bi, PDE_D.Ci, PDE_D, lambda);
        //
        dum=(lambda[1]-lambda[0]*PDE_D.Ui[i-1]);
        PDE_D.Ui[i]=lambda[2]/dum;
        PDE_D.bi[i]=(yi[i]+dz*PDE_D.Di[i]-lambda[0]*PDE_D.bi[i-1])/dum;
    }
    //=============================================
    PDE_D.bi[npxi-2]-=yi[npxi-1]*lambda[2]/dum;
    
    //=============================================
    // compute solution
    //=============================================
    yi[npxi-2]=PDE_D.bi[npxi-2];
    for(int i=npxi-3; i>0; i--) yi[i]=PDE_D.bi[i]-PDE_D.Ui[i]*yi[i+1];  
    
    return;
}


//=====================================================================================
//
// 2. order scheme (in spacial derivative)
//
//=====================================================================================

//-------------------------------------------------------------------------------------
// precompute coefficients
//-------------------------------------------------------------------------------------
void setup_Lagrange_interpolation_coefficients_O2(PDE_Stepper_Data &PDE_D, vector<double> &xarr)
{
    int np=xarr.size();

    cout << "\n setup_Lagrange_interpolation_coefficients_O2:: setting up grid " << endl;
    
    PDE_D.LG.dli_dx.clear();
    PDE_D.LG.d2li_d2x.clear();

    //====================================================================
    // allocate memory
    //====================================================================
    vector<double> dum(5, 0.0);
    for(int k=0; k<np; k++)
    {
        PDE_D.LG.dli_dx  .push_back(dum);
        PDE_D.LG.d2li_d2x.push_back(dum);
    }
    
    //====================================================================
    // 1. interior point (lower boundary not needed)
    //====================================================================
    int i=1;
    for(int j=0; j<5; j++)
    {
        PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 1, &xarr[i-1], 5);
        PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 1, &xarr[i-1], 5);
    }
    
    //====================================================================
    // last interior point (upper boundary not needed)
    //====================================================================
    i=np-2;
    for(int j=0; j<5; j++)
    {
        PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 3, &xarr[i-3], 5);
        PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 3, &xarr[i-3], 5);
    }

    //====================================================================
    // the rest
    //====================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(i=2; i<np-2; i++)
        for(int j=0; j<5; j++)
        {
            PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 2, &xarr[i-2], 5);
            PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 2, &xarr[i-2], 5);
        }

    cout << " setup_Lagrange_interpolation_coefficients_O2:: done " << endl;

    return;
}
        

//-------------------------------------------------------------------------------------
// compute lambda_i
//-------------------------------------------------------------------------------------
void compute_lambda_i_02(const double dz, const int i, const int center, 
                         const vector<double> &Ai, const vector<double> &Bi, 
                         const vector<double> &Ci, 
                         PDE_Stepper_Data &PDE_D, vector<double> &lambda)
{       
    for(int j=0; j<5; j++) lambda[j]=PDE_D.LG.d2li_d2x[i][j]*Ai[i]+PDE_D.LG.dli_dx[i][j]*Bi[i];
    lambda[center]+=Ci[i];
    lambda[center]-=1.0/(dz+1.0e-100);
}

//-------------------------------------------------------------------------------------
// compute solution
//
// 27.04.2011: rearranged some of the evaluation to make things more compact
// 22.04.2011: Added simple openmp support for matrix element evaluations
// 28.02.2011: Found small bug in the matrix solving part. Affected the boundary.
// 24.02.2011: This routine was optimized to reduce the number of operations per call
//
//-------------------------------------------------------------------------------------
void Step_PDE_O2t(double theta, double zs, double ze, const vector<double> &xi, vector<double> &yi, 
                  double yi_low, double yi_up, PDE_Stepper_Data &PDE_D,
                  void (*func)(double z, const vector<double> &xi, vector<double> &Ai, 
                               vector<double> &Bi, vector<double> &Ci, vector<double> &Di))
{
    double kap=theta-1.0, eta=kap/theta;
    double dz=ze-zs;
    int npxi=xi.size();
    
    //=================================================================================
    // evaluate coefficients of PDE at future time
    //================================================================================= 
    if(!PDE_solver_initial_call_of_solver)
    {
        func(zs, xi, *PDE_D.Aip_ptr, *PDE_D.Bip_ptr, *PDE_D.Cip_ptr, *PDE_D.Dip_ptr);
        PDE_solver_initial_call_of_solver=1;
    }
    
    func(ze, xi, *PDE_D.Ai_ptr, *PDE_D.Bi_ptr, *PDE_D.Ci_ptr, *PDE_D.Di_ptr);
    
    //=================================================================================
    // compute all matrix elements
    //=================================================================================
    int i=1;
    compute_lambda_i_02(dz*theta, i, 1, *PDE_D.Ai_ptr, *PDE_D.Bi_ptr, *PDE_D.Ci_ptr , PDE_D, PDE_D.lambda [i]);
    compute_lambda_i_02(dz*kap  , i, 1, *PDE_D.Aip_ptr,*PDE_D.Bip_ptr,*PDE_D.Cip_ptr, PDE_D, PDE_D.lambdap[i]); 
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(i=2; i<npxi-2; i++)
    {
        compute_lambda_i_02(dz*theta, i, 2, *PDE_D.Ai_ptr, *PDE_D.Bi_ptr, *PDE_D.Ci_ptr , PDE_D, PDE_D.lambda [i]);
        compute_lambda_i_02(dz*kap  , i, 2, *PDE_D.Aip_ptr,*PDE_D.Bip_ptr,*PDE_D.Cip_ptr, PDE_D, PDE_D.lambdap[i]); 
    }
    
    i=npxi-2;
    compute_lambda_i_02(dz*theta, i, 3, *PDE_D.Ai_ptr, *PDE_D.Bi_ptr, *PDE_D.Ci_ptr , PDE_D, PDE_D.lambda [i]);
    compute_lambda_i_02(dz*kap  , i, 3, *PDE_D.Aip_ptr,*PDE_D.Bip_ptr,*PDE_D.Cip_ptr, PDE_D, PDE_D.lambdap[i]); 
    
    //================================================================================= 
    // initialize all b-vectors with common parts
    //================================================================================= 
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(i=1; i<npxi-1; i++) PDE_D.bi[i]=eta*(*PDE_D.Dip_ptr)[i]-(*PDE_D.Di_ptr)[i];
    
    //================================================================================= 
    i=1;
    for(int m=0; m<5; m++) PDE_D.bi[i]+=PDE_D.lambdap[i][m]*eta*yi[i+m-1];
    PDE_D.bi[i]-=PDE_D.lambda[i][0]*yi_low; // JC changed sign 28.02.2011
    
    i=2;
    PDE_D.bi[i]-=PDE_D.lambda[i][0]*yi_low; // JC changed sign 28.02.2011
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(i=2; i<npxi-2; i++) 
    {
        for(int m=0; m<5; m++) PDE_D.bi[i]+=PDE_D.lambdap[i][m]*eta*yi[i+m-2];
    }
    
    i=npxi-3;
    PDE_D.bi[i]-=PDE_D.lambda[i][4]*yi_up;
    
    i=npxi-2;
    for(int m=0; m<5; m++) PDE_D.bi[i]+=PDE_D.lambdap[i][m]*eta*yi[i+m-3];
    PDE_D.bi[i]-=PDE_D.lambda[i][4]*yi_up;
    //================================================================================= 
    
    //================================================================================= 
    // compute Ui* and bi*
    //================================================================================= 
    double dum1, dum2, dum3;
    double Zlow;
    
    //================================================================================= 
    // 1. intereor point
    //=================================================================================  
    i=1;
    //
    PDE_D.bi[i]/=PDE_D.lambda[i][1];
    PDE_D.Ui[i] =PDE_D.lambda[i][2]/PDE_D.lambda[i][1];
    PDE_D.Vi[i] =PDE_D.lambda[i][3]/PDE_D.lambda[i][1];
    Zlow        =PDE_D.lambda[i][4]/PDE_D.lambda[i][1];
    
    //================================================================================= 
    // 2. intereor point
    //=================================================================================  
    i=2;
    //
    dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][1]*PDE_D.Ui[i-1];
    //
    PDE_D.bi[i]-=PDE_D.lambda[i][1]*PDE_D.bi[i-1];
    PDE_D.bi[i]/=dum2;
    //
    PDE_D.Ui[i] =(PDE_D.lambda[i][3]-PDE_D.lambda[i][1]*PDE_D.Vi[i-1])/dum2;
    PDE_D.Vi[i] =(PDE_D.lambda[i][4]-PDE_D.lambda[i][1]*Zlow)/dum2;
    
    //================================================================================= 
    // 3. intereor point
    //=================================================================================  
    i=3;
    //
    dum1=PDE_D.lambda[i][1]-PDE_D.lambda[i][0]*PDE_D.Ui[i-2];
    dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][0]*PDE_D.Vi[i-2];
    //
    PDE_D.bi[i]-=PDE_D.lambda[i][0]*PDE_D.bi[i-2];
    PDE_D.Ui[i] =PDE_D.lambda[i][3]-PDE_D.lambda[i][0]*Zlow;
    //---------------------------------------------
    dum2-=dum1*PDE_D.Ui[i-1];
    //
    PDE_D.bi[i]-=dum1*PDE_D.bi[i-1];
    PDE_D.bi[i]/=dum2;
    //
    PDE_D.Ui[i]-=dum1*PDE_D.Vi[i-1];
    PDE_D.Ui[i]/=dum2;
    PDE_D.Vi[i] =PDE_D.lambda[i][4]/dum2;
    
    //================================================================================= 
    // intereor point 4...n-4
    //=================================================================================  
    for(i=4; i<npxi-3; i++)
    {
        dum1=PDE_D.lambda[i][1]-PDE_D.lambda[i][0]*PDE_D.Ui[i-2];
        dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][0]*PDE_D.Vi[i-2];
        //
        PDE_D.bi[i]-=PDE_D.lambda[i][0]*PDE_D.bi[i-2];
        PDE_D.Ui[i] =PDE_D.lambda[i][3];
        //---------------------------------------------
        dum2-=dum1*PDE_D.Ui[i-1];
        //
        PDE_D.bi[i]-=dum1*PDE_D.bi[i-1];
        PDE_D.bi[i]/=dum2;
        //
        PDE_D.Ui[i]-=dum1*PDE_D.Vi[i-1];
        PDE_D.Ui[i]/=dum2;
        PDE_D.Vi[i] =PDE_D.lambda[i][4]/dum2;
    }
    
    //================================================================================= 
    // intereor point n-3
    //=================================================================================  
    i=npxi-3;
    //
    dum1=PDE_D.lambda[i][1]-PDE_D.lambda[i][0]*PDE_D.Ui[i-2];
    dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][0]*PDE_D.Vi[i-2];
    //
    PDE_D.bi[i]-=PDE_D.lambda[i][0]*PDE_D.bi[i-2];
    PDE_D.Ui[i] =PDE_D.lambda[i][3];
    //---------------------------------------------
    dum2-=dum1*PDE_D.Ui[i-1];
    //
    PDE_D.bi[i]-=dum1*PDE_D.bi[i-1];
    PDE_D.bi[i]/=dum2;
    //
    PDE_D.Ui[i]-=dum1*PDE_D.Vi[i-1];
    PDE_D.Ui[i]/=dum2;
    PDE_D.Vi[i] =0.0;
    
    //================================================================================= 
    // intereor point n-2
    //=================================================================================  
    i=npxi-2;
    //
    dum1=PDE_D.lambda[i][1]-PDE_D.lambda[i][0]*PDE_D.Ui[i-3];
    dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][0]*PDE_D.Vi[i-3];
    dum3=PDE_D.lambda[i][3];
    //
    PDE_D.bi[i]-=PDE_D.lambda[i][0]*PDE_D.bi[i-3];
    //---------------------------------------------
    dum2-=dum1*PDE_D.Ui[i-2];
    dum3-=dum1*PDE_D.Vi[i-2];
    //
    PDE_D.bi[i]-=dum1*PDE_D.bi[i-2];
    //---------------------------------------------
    dum3-=dum2*PDE_D.Ui[i-1];
    //
    PDE_D.bi[i]-=dum2*PDE_D.bi[i-1];
    PDE_D.bi[i]/=dum3;
    PDE_D.Ui[i]=PDE_D.Vi[i]=0.0;
    
    //================================================================================= 
    // lower & upper boundary condition
    //================================================================================= 
    yi[0]=yi_low;
    yi[npxi-1]=yi_up;
    
    //================================================================================= 
    // compute solution
    //================================================================================= 
    yi[npxi-2]=PDE_D.bi[npxi-2];
    //
    for(i=npxi-3; i>0; i--) yi[i]=PDE_D.bi[i]-PDE_D.Ui[i]*yi[i+1]-PDE_D.Vi[i]*yi[i+2];
    //
    yi[1]-=Zlow*yi[4];  
    
    //================================================================================= 
    // copy old coefficients
    //================================================================================= 
    PDE_D.ptr=PDE_D.Aip_ptr;
    PDE_D.Aip_ptr=PDE_D.Ai_ptr;
    PDE_D.Ai_ptr=PDE_D.ptr;
    //
    PDE_D.ptr=PDE_D.Bip_ptr;
    PDE_D.Bip_ptr=PDE_D.Bi_ptr;
    PDE_D.Bi_ptr=PDE_D.ptr;
    //
    PDE_D.ptr=PDE_D.Cip_ptr;
    PDE_D.Cip_ptr=PDE_D.Ci_ptr;
    PDE_D.Ci_ptr=PDE_D.ptr;
    //  
    PDE_D.ptr=PDE_D.Dip_ptr;
    PDE_D.Dip_ptr=PDE_D.Di_ptr;
    PDE_D.Di_ptr=PDE_D.ptr;
    //
    PDE_D.ptr=NULL;
    
    return;
}
