//==================================================================================================
// ODE_solver routines 
//
// Purpose: solve a system of stiff ODEs with moderate dimension (n<~100). 
// 
// Basic Aspects: The solver is based on a 6th order Gears method (implicit && stiffly-stable) with 
// variable time-step. The linear algebra system is solved iteratively using a biconjugate gradient 
// method. A guess for the next solution is extrapolated from the previous solutions; 
//
// Author: Jens Chluba with contributions from Geoff Vasil at CITA
// First implementation: 12/01/2010
// Last modification   : 15/06/2011
//==================================================================================================
// 15/06/2011: set max order to 5; This was making the code a bit faster; 
//             tidied the code a bit;
//==================================================================================================

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "ODE_solver_Rec.h"

using namespace std;

int verbosity_ODE_solver=0;

//=======================================================================================
namespace ODE_solver_Rec
{
    //===================================================================================
    void wait_f_r_loc(string mess)
    {
        char c;
        if(mess!="") cout << " " << mess << endl;
        cout << " To go on press return! " << endl;
        cin.get(c);
        return;
    }
    
    void wait_f_r_loc()
    {
        wait_f_r_loc("");
        return;
    }
    
    bool isnan_loc(double a)
    {
#if defined(__APPLE__) 
        return !(a==a);
#else 
        return isnan(a);
#endif
    }
    //===================================================================================

    
    //===================================================================================
    // pointers to gobal variables and memory
    //===================================================================================
    static ODE_Solver_data *global_ODE_Solver_info_ptr_for_function=NULL;
    
    
    //===================================================================================
    //
    // Linear-Algebra part
    //
    //===================================================================================
    
    //===================================================================================
    // All these routines are for simple linear algebra operations assuming that the
    // dimension of the problem is small (i.e. comparable to a few like in RECFAST)
    //===================================================================================
    
    //===================================================================================
    // scalar product c=a*b
    //===================================================================================
    double dot(const vector<double> &a, const vector<double> &b)
    {
        double scalar=0.0;
        for(unsigned int i=0; i<a.size(); i++) scalar+=a[i]*b[i];
        return scalar;
    }
    
    //===================================================================================
    // norm of vector a
    //===================================================================================
    double norm(const vector<double> &a){ return sqrt(dot(a, a)); }
    
    //===================================================================================
    // compute c=A*b (i.e. matrix times vector)
    //===================================================================================
    void A_times_b_is_c(const ODE_solver_matrix &Jac, 
                        const vector<double> &b, 
                        vector<double> &c)
    {
        //------------------------------------------------------
        // here it is assumed that A is symmetric with dimension 
        // equal to b & c; c is overwritten
        //------------------------------------------------------
        for(unsigned int j=0; j<c.size(); j++) c[j]=0.0;

        for(unsigned int i=0; i<Jac.A.size(); i++) 
            c[Jac.row[i]]+= Jac.A[i] * b[Jac.col[i]];
        
        return;
    }
    
    //===================================================================================
    // get inverse diagonal elements of matrix A
    //===================================================================================
    void Get_inverse_diags(int dim, const ODE_solver_matrix &Jac, vector<double> &d)
    {
        //-----------------------------------------------------------
        // here it is assumed that A is symmetric with dimension dim
        //-----------------------------------------------------------
        for(int i=0; i<dim; i++)
        {
            d[i]=Jac.A[Jac.diags[i]];
            
            if(d[i]==0)
            { 
                cerr << " error in preconditioner for element: " << i << endl; 
                exit(0); 
            }
            else d[i]=1.0/d[i];
        }
        
        return;
    }
    
    //===================================================================================
    // Iterative biconjugate gradiant routine -- BiCGSTAB
    //
    // BiCGSTAB solves the unsymmetric linear system Ax = b 
    // using the Preconditioned BiConjugate Gradient Stabilized method
    //
    // BiCGSTAB follows the algorithm described on p. 27 of the 
    // SIAM Templates book.
    //
    // The return value indicates convergence within max_iter (input)
    // iterations (0), or no convergence within max_iter iterations (1).
    //
    // Upon successful return, output arguments have the following values:
    //  
    //        x  --  approximate solution to Ax = b
    // max_iter  --  the number of iterations performed before the
    //               tolerance was reached
    //      tol  --  the residual after the final iteration
    //  
    // This routine was adapted from the IML++ http://math.nist.gov/iml++/
    //===================================================================================
    int BiCGSTAB_JC(const ODE_solver_matrix &Jac, vector<double> &x, 
                    const vector<double> &b, int &max_iter, double &tol)
    {
        int neq=b.size();
        double resid;
        double rho_1=1.0, rho_2=1.0, alpha=1.0, beta=1.0, omega=1.0;
        //
        vector<double> p(neq), phat(neq), s(neq), shat(neq);
        vector<double> t(neq), v(neq), r(neq), rtilde(neq);
        vector<double> invdiag(neq);
        
        //------------------------------------------------------
        // this is to precondition the Matrix (i.e. A=M*Atilde)
        // here simply the inverse diagonal elements are used
        //------------------------------------------------------
        Get_inverse_diags(neq, Jac, invdiag);
        
        double normb = norm(b);
        //------------------------------------------------------
        // r = b - A*x
        //------------------------------------------------------
        A_times_b_is_c(Jac, x, r);
        for(int i=0; i<neq; i++) r[i]=b[i]-r[i]; 
        //
        rtilde = r;
        
        if (normb == 0.0) normb = 1;
        
        if ((resid = norm(r) / normb) <= tol) 
        {
            tol = resid;
            max_iter = 0;
            return 0;
        }
        
        for (int k = 1; k <= max_iter; k++) 
        {
            rho_1 = dot(rtilde, r);
            if (rho_1 == 0) 
            {
                tol = norm(r) / normb;
                return 2;
            }
            
            if (k == 1) p = r;
            else 
            {
                beta = (rho_1/rho_2) * (alpha/omega);
                for(int i=0; i<neq; i++) p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
            //------------------------------------------------------
            // preconditioning phat
            //------------------------------------------------------
            for(int i=0; i<neq; i++) phat[i]=invdiag[i]*p[i];
            
            //------------------------------------------------------
            // v = A * phat;
            //------------------------------------------------------
            A_times_b_is_c(Jac, phat, v);
            
            alpha = rho_1 / dot(rtilde, v);
            for(int i=0; i<neq; i++) s[i] = r[i] - alpha * v[i];
            
            if ((resid = norm(s)/normb) < tol) 
            {
                for(int i=0; i<neq; i++) x[i] += alpha * phat[i];
                tol = resid;
                return 0;
            }
            //------------------------------------------------------
            // preconditioning shat
            //------------------------------------------------------
            for(int i=0; i<neq; i++) shat[i]=invdiag[i]*s[i];
            
            //------------------------------------------------------
            // t = A * shat;
            //------------------------------------------------------
            A_times_b_is_c(Jac, shat, t);
            
            omega = dot(t,s) / dot(t,t);
            //
            for(int i=0; i<neq; i++) x[i] += alpha * phat[i] + omega * shat[i];
            for(int i=0; i<neq; i++) r[i] = s[i] - omega * t[i];
            
            rho_2 = rho_1;
            if ((resid = norm(r) / normb) < tol) 
            {
                tol = resid;
                max_iter = k;
                return 0;
            }
            
            if (omega == 0) 
            {
                tol = norm(r) / normb;
                return 3;
            }
        }
        
        tol = resid;
        return 1;
    }
    
    //===================================================================================
    // Solving the Matrix Equation A*x == b for x
    //===================================================================================
    int ODE_Solver_Solve_LA_system(ODE_solver_matrix &Jac, vector<double > &b, 
                                   vector<double > &x, double tol, int verbose=0)
    {
        int ifail=0, maxit = 10000;                      // Maximum iterations
        
        for(unsigned int i=0; i<x.size(); i++) x[i]=0.0; // reset x-vector
        
        ifail = BiCGSTAB_JC(Jac, x, b, maxit, tol);      // Solve linear system
        
        if(verbose>0 && ifail!=0)
        {
            cout << endl << endl;
            cout << " BiCG flag = " << ifail << endl;
            cout << " iterations performed: " << maxit << endl;
            cout << " tolerance achieved  : " << tol << endl;
            cout << " Something bad happened to the ODE system. Restarting solver. " 
                 << endl;
            
            if(isnan_loc(tol)) ifail=10;
        }
        
        return ifail;
    }
    
        
    //===================================================================================
    //
    // Jacobian and ODE function part
    //
    //===================================================================================
    
    //===================================================================================
    // right hand side of the differential equation system dX/dz = g(z, X)
    //===================================================================================
    void ODE_solver_f(int neq, double z, double *y, double *f)
    {
        //==================================================================
        // evaluate the rhs of system  dy/dz == g(t,y)
        //==================================================================
        global_ODE_Solver_info_ptr_for_function->fcn_col_ptr(&neq, &z, &y[0], &f[0], -1); 
        return;
    }
    
    void ODE_solver_f_col(int neq, double z, double *y, double *f, int col)
    {
        //==================================================================
        // evaluate the rhs of system  dy/dz == g(t,y)
        //==================================================================
//      global_ODE_Solver_info_ptr_for_function->fcn_col_ptr(&neq, &z, &y[0], &f[0], col); 
        global_ODE_Solver_info_ptr_for_function->fcn_col_ptr(&neq, &z, &y[0], &f[0], -1); 
        return;
    }
    
    //===================================================================================
    //
    // to get the numerical Jacobian on the system
    //
    //===================================================================================
    inline double ODE_solver_df_dx_2point(const double &fp1, const double &fm1, 
                                          const double &h, double eps)
    { 
        if(fm1==0.0) return (fp1-fm1)/(2.0*h);        
        double dum=fp1/fm1-1.0;
        return fabs(dum)<=eps ? 0.0 : dum*fm1/(2.0*h);
    }
    
    //===================================================================================
    void ODE_solver_jac_2p(int neq, int col, double z, double Dz, double *y, double *r)
    { 
        //===============================================================================
        // r[i] is reset here
        // y[i] should contain the solution at z
        //===============================================================================
        double y0=y[col], Dyj=y[col]*1.0e-8, eps=1.0e-12;
        if(y0==0.0) Dyj=1.0e-8;
        
        vector<double> fp1(neq), fm1(neq);
        
        //===============================================================================
        // derivatives with respect to Xj.
        // A two-point formula is used. It has accuracy O(h^2)
        //===============================================================================
        // get f(yj+Dyj)
        y[col]=y0+Dyj;
        ODE_solver_f_col(neq, z, &y[0], &fp1[0], col);
        // get f(yj-Dyj)
        y[col]=y0-Dyj;
        ODE_solver_f_col(neq, z, &y[0], &fm1[0], col);
        // restore y again
        y[col]=y0;
        
        //===============================================================================
        // define numerical derivative
        //===============================================================================
        for(int k=0; k<neq; k++) 
            r[k]=-Dz*ODE_solver_df_dx_2point(fp1[k], fm1[k], Dyj, eps);
        
        //===============================================================================
        // add unity to the i==j element
        //===============================================================================
        r[col]+=1.0;
        
        return; 
    }
    
    //===================================================================================
    void ODE_solver_jac_user(int neq, int col, double z, double Dz, 
                             double *y, double *r)
    { 
        global_ODE_Solver_info_ptr_for_function->fcn_col_ptr(&neq, &z, y, r, col);
        
        //===============================================================================
        // define derivative
        //===============================================================================
        for(int k=0; k<neq; k++) r[k]=-Dz*r[k];
        
        //===============================================================================
        // add unity to the i==j element
        //===============================================================================
        r[col]+=1.0;
        
        return; 
    }
    
    //===================================================================================
    void ODE_solver_compute_Jacobian_Matrix(double z, double Dz, int neq, double *y, 
                                            ODE_solver_matrix &Jac, int mess=0)
    {
        //===============================================================================
        // z is the current redshift
        // y[] should contain the solution for which the jacobian is needed
        // Jac[] should have dimension neq*neq
        //===============================================================================
        if(mess==1) cout << " entering full Jacobinan update. " << endl;
        
        Jac.A.clear();
        Jac.row.clear();
        Jac.col.clear();
        Jac.diags.clear();
        
        //===============================================================================
        // this vector will be used to store the jacobian after calls of f_i(X)
        //===============================================================================
        vector<double> Jij(neq);
        
        //===============================================================================
        // fill Matrix with values
        //===============================================================================
        if(global_ODE_Solver_info_ptr_for_function->num_jac) 
        {
            for(int J=0; J<neq; J++)
            {   
                ODE_solver_jac_2p(neq, J, z, Dz, y, &Jij[0]);
                
                for(int row=0; row<neq; row++)
                    if(Jij[row]!=0.0)
                    {
                        Jac.A.push_back(Jij[row]);
                        Jac.col.push_back(J);
                        Jac.row.push_back(row);
                        
                        if(J==row) Jac.diags.push_back((int)(Jac.A.size()-1));
                    }
            }
        }
        else 
        {
            for(int J=0; J<neq; J++)
            {   
                ODE_solver_jac_user(neq, J, z, Dz, y, &Jij[0]);
                
                for(int row=0; row<neq; row++)
                    if(Jij[row]!=0.0)
                    {
                        Jac.A.push_back(Jij[row]);
                        Jac.col.push_back(J);
                        Jac.row.push_back(row);
                        
                        if(J==row) Jac.diags.push_back((int)(Jac.A.size()-1));
                    }
            }
        }
        
        if(mess==1) cout << " Number of non-zero elements = " << Jac.A.size() 
                         << ". Full matrix has " << neq*neq << " elements " << endl;
        
        if((int)Jac.diags.size()!=neq) 
            cout << " Hmmm. That should not happen... " << endl;
        
        return;
    }
    
    
    
    //===================================================================================
    //
    // this function has to be called to set/reset errors
    //
    //===================================================================================
    void ODE_Solver_set_errors(ODE_solver_accuracies &tols, 
                               ODE_Solver_data &ODE_Solver_info)
    {
        int neq=tols.rel.size();
        
        ODE_Solver_info.abs_vector.resize(neq);
        ODE_Solver_info.rel_vector.resize(neq);
        
        //===============================================================================
        // Copy absolute & relative accuracies
        //===============================================================================
        ODE_Solver_info.abs_vector=tols.abs;
        ODE_Solver_info.rel_vector=tols.rel;
        
        // find minimal relative tolerance
        double minrtol=1.0;
        for(int k=0; k<neq; k++)
        {
            if(minrtol>tols.rel[k] && tols.rel[k]>0.0) minrtol=tols.rel[k];
        }
        
        if(verbosity_ODE_solver>=1) 
            cout << " setting tolerance for Jacobian iteration to " << minrtol << endl;
        ODE_Solver_info.tolSol=minrtol;
        
        return;
    }
    
    //===================================================================================
    //
    // memory setup for ODE_solver_Solution and ODE_solver_accuracies
    //
    //===================================================================================
    void allocate_memory(ODE_solver_Solution &Sz, int neq)
    {
        Sz.y.resize(neq);
        Sz.dy.resize(neq);
        return;
    }
        
    void allocate_memory(ODE_solver_accuracies &tols, int neq)
    {
        tols.rel.resize(neq);
        tols.abs.resize(neq);
        return;
    }

    void allocate_memory(ODE_solver_Solution &Sz, ODE_solver_accuracies &tols, int neq)
    {
        allocate_memory(Sz, neq);
        allocate_memory(tols, neq);
        return;
    }
    
    //===================================================================================
    void clear_memory(ODE_solver_Solution &Sz)
    {
        Sz.y.clear();
        Sz.dy.clear();
        return;
    }

    void clear_memory(ODE_solver_accuracies &tols)
    {
        tols.rel.clear();
        tols.abs.clear();
        return;
    }

    void clear_memory(ODE_solver_Solution &Sz, ODE_solver_accuracies &tols)
    {
        clear_memory(Sz);
        clear_memory(tols);
        return;
    }
    
    //===================================================================================
    void clear_SolverData(ODE_Solver_data &SD)
    {
        clear_memory(SD.Stemp);
        clear_memory(*SD.Snewptr);
        clear_memory(*SD.Sguessptr);
        clear_memory(*SD.Snptr);
        clear_memory(*SD.Snm1ptr);
        clear_memory(*SD.Snm2ptr);
        clear_memory(*SD.Snm3ptr);
        clear_memory(*SD.Snm4ptr);
        clear_memory(*SD.Snm5ptr);
        
        SD.F.clear();
        SD.dY.clear();
        SD.abs_vector.clear();
        SD.rel_vector.clear();
        
        SD.Jacobian_ptr->A.clear();
        SD.Jacobian_ptr->row.clear();
        SD.Jacobian_ptr->col.clear();
        SD.Jacobian_ptr->diags.clear();
        
        delete SD.Snewptr;
        delete SD.Sguessptr;
        delete SD.Snptr;
        delete SD.Snm1ptr;
        delete SD.Snm2ptr;
        delete SD.Snm3ptr;
        delete SD.Snm4ptr;
        delete SD.Snm5ptr;
        delete SD.Jacobian_ptr;

        return;
    }

    //===================================================================================
    void clear_all_memory(ODE_solver_Solution &Sz, 
                          ODE_solver_accuracies &tols, 
                          ODE_Solver_data &SD)
    {
        clear_memory(Sz, tols);
        clear_SolverData(SD);
        return;
    }

    //===================================================================================
    //
    // this function has to be called to set up the memory and the initial solution
    //
    //===================================================================================
    void ODE_Solver_setup(ODE_solver_Solution &Sz, ODE_solver_accuracies &tols, 
                          ODE_Solver_data &ODE_Solver_info)
    {
        int neq=Sz.y.size();
        
        if(verbosity_ODE_solver>=1)
        {
            if(ODE_Solver_info.Snewptr==NULL) cout << "\n Setting up memory " << endl;
            else  cout << "\n Resetting memory " << endl;
            cout << " Total number of equations = " << neq << endl;
        }
        
        //-------------------------------------------------------------------------------
        // setting this pointer is only important for the ODE-function routine. All the 
        // other parts of the code do not use global variables, but pass the memory along.
        //-------------------------------------------------------------------------------
        global_ODE_Solver_info_ptr_for_function=&ODE_Solver_info;
        
        //-------------------------------------------------------------
        // create vectors & allocate memory
        //------------------------------------------------------------- 
        ODE_Solver_info.F.resize(neq);
        ODE_Solver_info.dY.resize(neq);
        ODE_Solver_info.Sptr=NULL;
        
        if(ODE_Solver_info.Snewptr==NULL)
        {
            //wait_f_r_loc("Memory");
            ODE_Solver_info.Snewptr=new ODE_solver_Solution;
            ODE_Solver_info.Sguessptr=new ODE_solver_Solution;
            ODE_Solver_info.Snptr=new ODE_solver_Solution;
            ODE_Solver_info.Snm1ptr=new ODE_solver_Solution;
            ODE_Solver_info.Snm2ptr=new ODE_solver_Solution;
            ODE_Solver_info.Snm3ptr=new ODE_solver_Solution;
            ODE_Solver_info.Snm4ptr=new ODE_solver_Solution;
            ODE_Solver_info.Snm5ptr=new ODE_solver_Solution;
            ODE_Solver_info.Jacobian_ptr= new ODE_solver_matrix;
        }
        
        ODE_Solver_info.Stemp.y.resize(neq); ODE_Solver_info.Stemp.dy.resize(neq); 
        //
        ODE_Solver_info.Snewptr->y.resize(neq); ODE_Solver_info.Snewptr->dy.resize(neq); 
        ODE_Solver_info.Sguessptr->y.resize(neq); ODE_Solver_info.Sguessptr->dy.resize(neq);
        // 
        ODE_Solver_info.Snptr->y.resize(neq); ODE_Solver_info.Snptr->dy.resize(neq); 
        ODE_Solver_info.Snm1ptr->y.resize(neq); ODE_Solver_info.Snm1ptr->dy.resize(neq); 
        ODE_Solver_info.Snm2ptr->y.resize(neq); ODE_Solver_info.Snm2ptr->dy.resize(neq); 
        ODE_Solver_info.Snm3ptr->y.resize(neq); ODE_Solver_info.Snm3ptr->dy.resize(neq); 
        ODE_Solver_info.Snm4ptr->y.resize(neq); ODE_Solver_info.Snm4ptr->dy.resize(neq);
        ODE_Solver_info.Snm5ptr->y.resize(neq); ODE_Solver_info.Snm5ptr->dy.resize(neq);
        
        //===============================================================================
        // copy solution
        //===============================================================================
        ODE_Solver_info.Snptr->z=Sz.z;
        ODE_Solver_info.Snptr->y=Sz.y;
        ODE_Solver_info.Snptr->dy=Sz.dy;
        
        ODE_solver_f(neq, ODE_Solver_info.Snptr->z, 
                     &ODE_Solver_info.Snptr->y[0], 
                     &ODE_Solver_info.Snptr->dy[0]);
        
        //===============================================================================
        // other data
        //===============================================================================
        ODE_Solver_info.Dz_z_last=0;    
        ODE_Solver_info.zstart_solver=ODE_Solver_info.Snptr->z; 
        ODE_Solver_info.order=1;    
        ODE_Solver_info.count=1;
        ODE_Solver_info.Jac_is_set=0;
        //  
        ODE_Solver_info.n_up=0;
        ODE_Solver_info.n_down=0;
        
        //===============================================================================
        // Copy absolute & relative accuracies
        //===============================================================================
        ODE_Solver_set_errors(tols, ODE_Solver_info);
        
        global_ODE_Solver_info_ptr_for_function=NULL;
        return;
    }
    
    //===================================================================================
    void ODE_Solver_set_up_solution_and_memory(ODE_solver_Solution &Sz, 
                                               ODE_solver_accuracies &tols, 
                                               ODE_Solver_data &ODE_Solver_info, 
                                               void (*fcn_ptr)(int *neq, double *z, 
                                                               double *y, double *f, 
                                                               int col), 
                                               bool num_jac)
    {
        ODE_Solver_info.fcn_col_ptr=fcn_ptr;
        ODE_Solver_info.num_jac=num_jac;
        ODE_Solver_setup(Sz, tols, ODE_Solver_info);
        
        return;
    }
    
    //===================================================================================
    //
    // time-step using Gears-method
    // order can be ==1; 2; 3; 4; 5; 6;
    //
    //===================================================================================
    
    //===================================================================================
    // aux-functions alpha_i coefficients
    //===================================================================================
    double ODE_Solver_Gears_fk1(double rk)
    { return 2.0+rk; } 
    
    double ODE_Solver_Gears_fk2(double r1, double rk)
    {
        double r=ODE_Solver_Gears_fk1(r1)*ODE_Solver_Gears_fk1(rk);
        return r-1.0;
    } 
    
    double ODE_Solver_Gears_fk3(double r1, double r2, double rk)
    {
        double r=ODE_Solver_Gears_fk2(r1, r2)*ODE_Solver_Gears_fk1(rk);
        return r-ODE_Solver_Gears_fk1(r1+r2);
    } 
    
    double ODE_Solver_Gears_fk4(double r1, double r2, double r3, double rk)
    {
        double r=ODE_Solver_Gears_fk3(r1, r2, r3)*ODE_Solver_Gears_fk1(rk);
        r-=ODE_Solver_Gears_fk1(r1+r2)*r3;
        return r-ODE_Solver_Gears_fk2(r1, r2);
    } 
    
    inline double ODE_Solver_Gears_fk5(double r1, double r2, double r3, 
                                       double r4, double rk)
    {
        double r=ODE_Solver_Gears_fk4(r1, r2, r3, r4)*ODE_Solver_Gears_fk1(rk);
        r-=(ODE_Solver_Gears_fk2(r1, r2+r3)+r2*r3)*r4;
        return r-ODE_Solver_Gears_fk3(r1, r2, r3);
    } 
    
    //===================================================================================
    // alpha_i coefficients
    //===================================================================================
    double ODE_Solver_delta0(double r1, double r2, double r3, double r4, double r5, 
                             double alp1, double alp2, double alp3, double alp4, double alp5)
    { return 1.0 + alp1*r1 + alp2*r2 + alp3*r3 + alp4*r4 + alp5*r5; }
    
    double ODE_Solver_alp0(double alp1, double alp2, double alp3, double alp4, double alp5)
    { return 1.0 - alp1 - alp2 - alp3 - alp4 - alp5; }
    
    double ODE_Solver_alp1(double r1, double r2, double r3, double r4, double r5, 
                           double alp2, double alp3, double alp4, double alp5)
    { return -((1.0 + alp2*r2*(2.0+r2) + alp3*r3*(2.0+r3) + alp4*r4*(2.0+r4) 
                    + alp5*r5*(2.0+r5))/(r1*(2.0+r1))); }
    
    double ODE_Solver_alp2(double r1, double r2, double r3, double r4, double r5, 
                           double alp3, double alp4, double alp5)
    {   
        double t=(1.0+r1);
        return -( t*t 
                 +alp3*r3*(r1-r3)*ODE_Solver_Gears_fk2(r1, r3)
                 +alp4*r4*(r1-r4)*ODE_Solver_Gears_fk2(r1, r4)
                 +alp5*r5*(r1-r5)*ODE_Solver_Gears_fk2(r1, r5) 
                 )
                /(     r2*(r1-r2)*ODE_Solver_Gears_fk2(r1, r2)  ); 
    }
    
    double ODE_Solver_alp3(double r1, double r2, double r3, double r4, double r5, 
                           double alp4, double alp5)
    { 
        double t=(1.0+r1)*(1.0+r2);
        return -( t*t
                 +alp4*r4*(r1-r4)*(r2-r4)*ODE_Solver_Gears_fk3(r1, r2, r4) 
                 +alp5*r5*(r1-r5)*(r2-r5)*ODE_Solver_Gears_fk3(r1, r2, r5) 
                 ) 
                /(     r3*(r1-r3)*(r2-r3)*ODE_Solver_Gears_fk3(r1, r2, r3)  ) ; 
    }
    
    double ODE_Solver_alp4(double r1, double r2, double r3, double r4, double r5, double alp5)
    { 
        double t=(1.0+r1)*(1.0+r2)*(1.0+r3);
        return -( t*t 
                 +alp5*r5*(r1-r5)*(r2-r5)*(r3-r5)*ODE_Solver_Gears_fk4(r1, r2, r3, r5)
                 )
                /(     r4*(r1-r4)*(r2-r4)*(r3-r4)*ODE_Solver_Gears_fk4(r1, r2, r3, r4)  ); 
    }
    
    double ODE_Solver_alp5(double r1, double r2, double r3, double r4, double r5)
    { 
        double t=(1.0+r1)*(1.0+r2)*(1.0+r3)*(1.0+r4);
        return -t*t/( r5*(r1-r5)*(r2-r5)*(r3-r5)*(r4-r5)
                    *ODE_Solver_Gears_fk5(r1, r2, r3, r4, r5) ); 
    }
    
    //===================================================================================
    // Gear's corrector (implicit)
    //===================================================================================
    double ODE_Solver_compute_ynp1(int order, ODE_solver_Solution &Snp1_new, 
                                   const ODE_solver_Solution &Snp1_guess, 
                                   const ODE_solver_Solution &Sn, 
                                   const ODE_solver_Solution &Snm1, 
                                   const ODE_solver_Solution &Snm2, 
                                   const ODE_solver_Solution &Snm3, 
                                   const ODE_solver_Solution &Snm4, 
                                   const ODE_solver_Solution &Snm5)
    {
        if(order<1 || order>6)
        { 
            cerr << " check order for ODE_Solver_compute_ynp1 " << endl; 
            exit(0); 
        }
        
        //===============================================================================
        // the structures yn, ynm1, ynm2 contain the information from the previous time 
        // steps the structure ynp1 contains the current version of ynp1 and dynm1
        //===============================================================================
        double delta=1.0, a0=0.0, a1=0.0, a2=0.0, a3=0.0, a4=0.0, a5=0.0;
        double Dznp1=Snp1_guess.z-Sn.z, Dzn=0.0, Dznm1=0.0, Dznm2=0.0, Dznm3=0.0, Dznm4=0.0;
        double rho[5]={0.0, 0.0, 0.0, 0.0, 0.0};
        int neq=Snp1_new.y.size();
        
        if(order>=2){ Dzn  =Sn.z-Snm1.z; rho[0]=Dzn  /Dznp1; }
        if(order>=3){ Dznm1=Sn.z-Snm2.z; rho[1]=Dznm1/Dznp1; }
        if(order>=4){ Dznm2=Sn.z-Snm3.z; rho[2]=Dznm2/Dznp1; }
        if(order>=5){ Dznm3=Sn.z-Snm4.z; rho[3]=Dznm3/Dznp1; }
        if(order>=6){ Dznm4=Sn.z-Snm5.z; rho[4]=Dznm4/Dznp1; 
        //    
                     a5=ODE_Solver_alp5(rho[0], rho[1], rho[2], rho[3], rho[4]); }
        if(order>=5) a4=ODE_Solver_alp4(rho[0], rho[1], rho[2], rho[3], rho[4], a5);
        if(order>=4) a3=ODE_Solver_alp3(rho[0], rho[1], rho[2], rho[3], rho[4], a4, a5);
        if(order>=3) a2=ODE_Solver_alp2(rho[0], rho[1], rho[2], rho[3], rho[4], a3, a4, a5);
        if(order>=2) a1=ODE_Solver_alp1(rho[0], rho[1], rho[2], rho[3], rho[4], a2, a3, a4, a5);
        //
        a0=ODE_Solver_alp0(a1, a2, a3, a4, a5); 
        delta=ODE_Solver_delta0(rho[0], rho[1], rho[2], rho[3], rho[4], a1, a2, a3, a4, a5);
        //
        for(int i=0; i<neq; i++) Snp1_new.y[i]=Dznp1*delta*Snp1_guess.dy[i]+a0*Sn.y[i]+a1*Snm1.y[i]
                                               +a2*Snm2.y[i]+a3*Snm3.y[i]+a4*Snm4.y[i]+a5*Snm5.y[i];    
        
        return delta;   
    }
    
    //===================================================================================
    // beta_i coefficients for extrapolation with avariable step-size
    // These expression where derived by Geoff Vasil (CITA)
    //===================================================================================
    double ODE_Solver_beta0(double beta1, double beta2, double beta3, double beta4, double beta5)
    { return 1.0 - beta1 - beta2 - beta3 - beta4 - beta5; }
    
    double ODE_Solver_beta1(double r1, double r2, double r3, double r4, double r5, 
                            double beta2, double beta3, double beta4, double beta5)
    { return -(1.0 + beta2*r2 + beta3*r3 + beta4*r4 + beta5*r5)/r1; }
    
    double ODE_Solver_beta2(double r1, double r2, double r3, double r4, double r5, 
                            double beta3, double beta4, double beta5)
    { return -( (1.0+r1) + beta3*r3*(r1-r3) + beta4*r4*(r1-r4) + beta5*r5*(r1-r5) )/((r1-r2)*r2); }
    
    double ODE_Solver_beta3(double r1, double r2, double r3, double r4, double r5, 
                            double beta4, double beta5)
    { return ( (1.0+r1)*(1.0+r2) + beta4*r4*(r1-r4)*(r2-r4) 
                                 + beta5*r5*(r1-r5)*(r2-r5) )/((r1-r3)*r3*(r3-r2)); }
    
    double ODE_Solver_beta4(double r1, double r2, double r3, double r4, double r5, double beta5)
    { return -( (1.0+r1)*(1.0+r2)*(1.0+r3) + beta5*r5*(r1-r5)*(r2-r5)*(r3-r5) )
              /((r3-r4)*r4*(r4-r1)*(r4-r2)); }
    
    double ODE_Solver_beta5(double r1, double r2, double r3, double r4, double r5)
    { return (1.0+r1)*(1.0+r2)*(1.0+r3)*(1.0+r4)/((r2-r5)*r5*(r5-r1)*(r5-r3)*(r5-r4)); }
    
    //===================================================================================
    //
    // extrapolation using old function values
    //
    //===================================================================================
    void ODE_Solver_extrapolate_ynp1(int order, ODE_solver_Solution &Snp1, 
                                     const ODE_solver_Solution &Sn, 
                                     const ODE_solver_Solution &Snm1, 
                                     const ODE_solver_Solution &Snm2, 
                                     const ODE_solver_Solution &Snm3,
                                     const ODE_solver_Solution &Snm4, 
                                     const ODE_solver_Solution &Snm5)
    {
        if(order<1 || order>6)
        { 
            cerr << " check order for ODE_Solver_extrapolate_ynp1 " << endl; 
            exit(0); 
        }
        
        //================================================================================
        // the structures yn, ynm1,... contain the information from the previous time 
        // steps the structure ynp1 contains the current version of ynp1 and dynm1
        //================================================================================
        double b0=1.0,  b1=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.0;
        double Dznp1=Snp1.z-Sn.z, Dzn=0.0, Dznm1=0.0, Dznm2=0.0, Dznm3=0.0, Dznm4=0.0;
        double rho[5]={0.0, 0.0, 0.0, 0.0, 0.0};
        int neq=Sn.y.size();
        
        if(order>=2){ Dzn  =Sn.z-Snm1.z; rho[0]=Dzn  /Dznp1; }
        if(order>=3){ Dznm1=Sn.z-Snm2.z; rho[1]=Dznm1/Dznp1; }
        if(order>=4){ Dznm2=Sn.z-Snm3.z; rho[2]=Dznm2/Dznp1; }
        if(order>=5){ Dznm3=Sn.z-Snm4.z; rho[3]=Dznm3/Dznp1; }
        if(order>=6){ Dznm4=Sn.z-Snm5.z; rho[4]=Dznm4/Dznp1; 
        //    
                     b5=ODE_Solver_beta5(rho[0], rho[1], rho[2], rho[3], rho[4]); }
        if(order>=5) b4=ODE_Solver_beta4(rho[0], rho[1], rho[2], rho[3], rho[4], b5);
        if(order>=4) b3=ODE_Solver_beta3(rho[0], rho[1], rho[2], rho[3], rho[4], b4, b5);
        if(order>=3) b2=ODE_Solver_beta2(rho[0], rho[1], rho[2], rho[3], rho[4], b3, b4, b5);
        if(order>=2) b1=ODE_Solver_beta1(rho[0], rho[1], rho[2], rho[3], rho[4], b2, b3, b4, b5);
        //
        b0=ODE_Solver_beta0(b1, b2, b3, b4, b5);
        //
        for(int i=0; i<neq; i++) Snp1.y[i]=b0*Sn.y[i]+b1*Snm1.y[i]+b2*Snm2.y[i]
                                          +b3*Snm3.y[i]+b4*Snm4.y[i]+b5*Snm5.y[i];    
        
        
        return; 
    }
    
    //===================================================================================
    //
    // Functions for the ODE solver (error checking)
    //
    //===================================================================================
    inline double ODE_Solver_error_check(double z, double rtol, double atol, double y, 
                                         double Dy, int mflag, int component, string mess="")
    {
        Dy=fabs(Dy);
        double err_aim=rtol*fabs(y), fac=1.0e+10;
        if(err_aim<=atol)
        { 
            if(mflag>=1) cout << " Time-stepper: absolute error control for ( " << mess 
                              << " / component " << component << " ) at z= " 
                              << z << " and y= " << y << " |Dy|= " << Dy 
                              << " Dy_rel= " << err_aim << " Dy_abs= " << atol << endl;
            err_aim=atol;
        }
        
        if(Dy!=0.0) fac=err_aim/Dy;
        
        return fac;
    }
    
    int ODE_Solver_estimate_error_and_next_possible_step_size(int order, const double zin, 
                                                              double z, double zend, 
                                                              double &Dz_z, const double Dz_z_max, 
                                                              ODE_Solver_data &ODE_Solver_info, 
                                                              const ODE_solver_Solution &S, 
                                                              const ODE_solver_Solution &Sn, 
                                                              const ODE_solver_Solution &Snm1, 
                                                              const ODE_solver_Solution &Snm2, 
                                                              const ODE_solver_Solution &Snm3, 
                                                              const ODE_solver_Solution &Snm4, 
                                                              const ODE_solver_Solution &Snm5, 
                                                              int messflag=1)
    {
        //--------------------------------------------------------
        // even if higher order is used just go with fourth order
        // for error check (added June 15)
        //--------------------------------------------------------        
        int loc_order=(int) min(max(order, 1), 4);
        
        //--------------------------------------------------------
        // S contains best approximation to order 'order'
        // *Stemp will contain best approximation to order 'order-1'
        //--------------------------------------------------------
        int neq=S.y.size();
        ODE_Solver_info.Stemp.z=z;
        ODE_Solver_extrapolate_ynp1(order, ODE_Solver_info.Stemp, 
                                    Sn, Snm1, Snm2, Snm3, Snm4, Snm5); 
        
        //--------------------------------------------------------
        // compute difference
        //--------------------------------------------------------
        for(int i=0; i<neq; i++) ODE_Solver_info.Stemp.y[i]-=S.y[i];
        
        //--------------------------------------------------------
        // estimate errors
        //--------------------------------------------------------
        double power, fac=1.0e+10, dum;
        int max_err_time=-1;
        //--------------------------------------------------------
        if(loc_order<=1) power=1.0;
        else power=1.0/loc_order;
        
        double tol=1.0e-5;
/*
        //--------------------------------------------------------
        // rho
        //--------------------------------------------------------
        dum=ODE_Solver_error_check(z, tol, ODE_Solver_info.abs_vector[0], S.y[0], 
                                   ODE_Solver_info.Stemp.y[0], messflag, 0, "rho");
        if(fac>dum){ fac=dum; max_err_time=0; }
        
        //--------------------------------------------------------
        // Hydrogen X1s
        //--------------------------------------------------------
        dum=ODE_Solver_error_check(z, tol, ODE_Solver_info.abs_vector[1], S.y[1], 
                                   ODE_Solver_info.Stemp.y[1], messflag, 1, "HI 1s");
        if(fac>dum){ fac=dum; max_err_time=1; }
        // X2s
        dum=ODE_Solver_error_check(z, tol, ODE_Solver_info.abs_vector[2], S.y[2], 
                                   ODE_Solver_info.Stemp.y[2], messflag, 2, "HI 2s");
        if(fac>dum){ fac=dum; max_err_time=2; }
        // X2p
        dum=ODE_Solver_error_check(z, tol, ODE_Solver_info.abs_vector[3], S.y[3], 
                                   ODE_Solver_info.Stemp.y[3], messflag, 3, "HI 2p");
        if(fac>dum){ fac=dum; max_err_time=3; }
*/      
        for(int l=0; l<neq; l++)
        {
            if(ODE_Solver_info.rel_vector[l]!=0.0) 
            {
//              dum=ODE_Solver_error_check(z, ODE_Solver_info.rel_vector[l], 
//                                        ODE_Solver_info.abs_vector[l], S.y[l], 
//                                        ODE_Solver_info.Stemp.y[l], messflag, l, "?");
                
                dum=ODE_Solver_error_check(z, tol, ODE_Solver_info.abs_vector[l], S.y[l], 
                                           ODE_Solver_info.Stemp.y[l], messflag, l, "?");
                if(fac>dum){ fac=dum; max_err_time=l; }
            }
        }
        
        if(messflag>=2) cout << " Max error in component " << max_err_time << endl;
        
        //--------------------------------------------------------
        // change step-size; 
        // limit change by different values, depending on success
        //--------------------------------------------------------
        double Dz_z_min=1.0e-10;
        double Dz_z_new=Dz_z*min(100.0, pow(fac, power));
        
        Dz_z_new=min(max(Dz_z_new, Dz_z_min), Dz_z_max);
        //
        if(messflag==1)
        {
            cout << " Estimated change: ";
            cout << " order= " << order << " estimate order= " << loc_order 
                 << " current step-size: " << Dz_z 
                 << " -- new suggested step-size: " << Dz_z_new << " -- Dz= " << zin*Dz_z_new 
                 << " z= " << z 
                 << " zend= " << zend << " Dz= " << z-zend << " zin= " << zin << " Dz= " << zin-z;
            //
            if(Dz_z_new==Dz_z) cout << " -- no change " << endl;
            else if(Dz_z_new==Dz_z_max) cout << " -- maximal step-size reached " << endl;
            else if(Dz_z_new>Dz_z) cout << " -- increase by: " << Dz_z_new/Dz_z << endl;
            else cout << " -- decrease by: " << Dz_z_new/Dz_z << endl;
        }
        
        //--------------------------------------------------------
        // to put limits on the changes in Dz (not too small)
        //--------------------------------------------------------
        if(Dz_z*0.9>Dz_z_new)
        {
            Dz_z=Dz_z_new;
            //
            ODE_Solver_info.n_down=(int) min(ODE_Solver_info.n_down+1, 5);
            
            //--------------------------------------------------------
            // if there were several subsequent decreases of the 
            // step-size take some more drastic action
            //--------------------------------------------------------
            if(ODE_Solver_info.n_down>=3)
            {
                if(ODE_Solver_info.n_down>=5)
                {
                    Dz_z=Dz_z_min;
                    cout << "\n Reset to minimal step-size " << endl;
                }
                else 
                {
                    Dz_z/=pow(2.0, ODE_Solver_info.n_down-1); 
                    cout << "\n Checking to decreasing the step-size by additional factor of " 
                         << pow(2.0, ODE_Solver_info.n_down-1) << endl;
                }
            }
            
            if(Dz_z*0.5>Dz_z_new) return 1;
        }
        else if(Dz_z*1.1<Dz_z_new)
        {
            Dz_z=min(Dz_z_new, Dz_z_max);
            //
            ODE_Solver_info.n_down=0;
            
            return 0;
        }
        
        ODE_Solver_info.n_down=0;
        return 0;
    }
    
    //===================================================================================
    //
    // setting up the equation system for matrix solver
    //
    //===================================================================================
    int ODE_Solver_compute_new_Solution(int order, const double zout, const double reltol, 
                                        ODE_Solver_data &ODE_Solver_info, 
                                        ODE_solver_Solution * &Snew, 
                                        const ODE_solver_Solution &Sn, 
                                        const ODE_solver_Solution &Snm1, 
                                        const ODE_solver_Solution &Snm2, 
                                        const ODE_solver_Solution &Snm3, 
                                        const ODE_solver_Solution &Snm4, 
                                        const ODE_solver_Solution &Snm5)
    {
        //==============================================================
        double zin=Sn.z, Dz=zout-zin, delta=0;
        int neq=Sn.y.size(), converged=0;

        //==============================================================
        // tolerances for Jacobian convergence
        //==============================================================
        double tolJac = reltol/2.0;

        //==============================================================
        int Jac_loops=0;
        int LA_error=0;
        int mess=0;
        
        //==============================================================
        // set initial values for iteration: extrapolate from z--> zout
        //==============================================================
        ODE_Solver_info.Sguessptr->z=Snew->z=zout;
        ODE_Solver_extrapolate_ynp1(order, *ODE_Solver_info.Sguessptr, 
                                    Sn, Snm1, Snm2, Snm3, Snm4, Snm5); 
        
        do{
            Jac_loops++;
            
            if(LA_error!=0)
            {
                ODE_Solver_extrapolate_ynp1(order, *ODE_Solver_info.Sguessptr, 
                                            Sn, Snm1, Snm2, Snm3, Snm4, Snm5); 
                LA_error=0;
                if(mess>=2) cout << " resetting Solver " << endl; 
            }
            
            //--------------------------------------------------------
            // redshift z==zn --> compute f_i(yguess, z)
            //--------------------------------------------------------
            ODE_solver_f(neq, zout, 
                         &ODE_Solver_info.Sguessptr->y[0], 
                         &ODE_Solver_info.Sguessptr->dy[0]);
            
            //--------------------------------------------------------
            // after calling this function, ynew will contain y(yguess)
            // delta is the coefficient in front of h f(y); it depends 
            // on the order that is used to compute things
            //--------------------------------------------------------
            delta=ODE_Solver_compute_ynp1(order, *Snew, *ODE_Solver_info.Sguessptr, 
                                          Sn, Snm1, Snm2, Snm3, Snm4, Snm5);
            
            //--------------------------------------------------------
            // this is rhs of J_F Dx = -F
            //--------------------------------------------------------
            for(int i=0; i<neq; i++) 
                ODE_Solver_info.F[i]=Snew->y[i]-ODE_Solver_info.Sguessptr->y[i]; 
            
            //--------------------------------------------------------
            // compute Jacobian
            //--------------------------------------------------------
            ODE_solver_compute_Jacobian_Matrix(zout, Dz*delta, neq, 
                                               &ODE_Solver_info.Sguessptr->y[0], 
                                               *ODE_Solver_info.Jacobian_ptr); 
            
            //--------------------------------------------------------
            // solve equation for dY
            //--------------------------------------------------------
            LA_error=ODE_Solver_Solve_LA_system(*ODE_Solver_info.Jacobian_ptr, 
                                                ODE_Solver_info.F, 
                                                ODE_Solver_info.dY, tolJac, 1);
                
            //--------------------------------------------------------
            // if nan was produced return main
            //--------------------------------------------------------
            if(LA_error==10) return 10;
            
            //--------------------------------------------------------
            // update current solution
            //--------------------------------------------------------
            if(LA_error==0) for(int i=0; i<neq; i++) 
                ODE_Solver_info.Sguessptr->y[i]+=ODE_Solver_info.dY[i];
            
            //--------------------------------------------------------------------
            // if X1s or rho become unphysical one should restart the solver
            // 22.04.2011: changed this condition to be a bit wider.
            //--------------------------------------------------------------------
            if(ODE_Solver_info.Sguessptr->y[0]>10.0 || ODE_Solver_info.Sguessptr->y[0]<0.0 ||
               ODE_Solver_info.Sguessptr->y[1]>10.0 || ODE_Solver_info.Sguessptr->y[1]<0.0)
                LA_error=1;
            
            //--------------------------------------------------------
            // check convergence of iteration
            //--------------------------------------------------------
            converged=1;
            if(LA_error==0)
            {
                for(int k=0; k<neq; k++) 
                    if( ODE_Solver_info.rel_vector[k]!=0.0 && 
                       (fabs(ODE_Solver_info.dY[k])>=
                        max(ODE_Solver_info.rel_vector[k]
                           *fabs(ODE_Solver_info.Sguessptr->y[k]), ODE_Solver_info.abs_vector[k]) 
                       ) )
                    { 
                        converged=0; 
                        if(mess>=2) cout << " " << k << " " << fabs(ODE_Solver_info.dY[k]) << " " 
                                         << ODE_Solver_info.rel_vector[k]
                                            *fabs(ODE_Solver_info.Sguessptr->y[k]) << " " 
                                         << ODE_Solver_info.abs_vector[k] << endl;
                        break; 
                    }
            }

            if(converged==0)
            {
                if(!(Jac_loops%5))
                { 
                    tolJac/=2.0; 
                    cout << " Tightening Jac error setting at it-# " << Jac_loops << endl; 
                }
                
                if(Jac_loops>25) 
                {   
                    converged=1;
                    if(mess>=2){ cout << " Not very happy :S " << endl; wait_f_r_loc(); } 
                    return 1; 
                }
            }
        } 
        while(converged==0);
        
        if(mess>=1) cout << "\n Number of Jacobian iterations " << Jac_loops 
                         << " needed for propagation from z= " << zin 
                         << " by Dz= " << zin-zout << endl;
        
        //--------------------------------------------------------
        // swap variables, so that Snew contains new solution
        //--------------------------------------------------------
        ODE_Solver_info.Sptr=Snew;
        Snew=ODE_Solver_info.Sguessptr;
        ODE_Solver_info.Sguessptr=ODE_Solver_info.Sptr;
        ODE_Solver_info.Sptr=NULL;
        
        return 0;
    }
    
    
    //===================================================================================
    //
    // do time-step
    //
    //===================================================================================
    int ODE_Solver_Solve_history(double zs, double zend, ODE_solver_Solution &Sz, 
                                 ODE_Solver_data &ODE_Solver_info)
    {   
        if(ODE_Solver_info.Snewptr==NULL)
        {
            cout << " ODE_Solver_Solve_history:: please set up the memory first! Exiting" 
                 << endl;
            exit(0);
        }
        
        //---------------------------------------------------------------------------
        // setting this pointer is only important for the ODE-function routine. All 
        // the other parts of the code do not use global variables, but pass the 
        // memory along.
        //---------------------------------------------------------------------------
        global_ODE_Solver_info_ptr_for_function=&ODE_Solver_info;
        
        //===========================================================================
        // define integration direction
        //===========================================================================
        if(zs>zend) ODE_Solver_info.direction=1;
        else ODE_Solver_info.direction=-1;
        
        //===========================================================================
        // setup
        //===========================================================================
        double zout=zs, zin=zs, Dz_z_old;
        double Dz_z=max(1.0e-8, ODE_Solver_info.Dz_z_last), Dz_z_max;
        int redo_run=0, LAerror=0;
        int maxorder=5;
        int messflag=0;
        
        if(zs!=0.0) Dz_z_max=fabs((zend-zs)/zs);
        else Dz_z_max=fabs(zend-zs);
        
        Dz_z_max=10.0/zs;
        
        //===========================================================================
        // for initial step reduce the accuracy
        //===========================================================================
        // make sure that max order is not set larger
        ODE_Solver_info.order=(int)min(ODE_Solver_info.order, maxorder);
        
        //===========================================================================
        // try some small step for very first integration
        //===========================================================================
        if(ODE_Solver_info.Dz_z_last==0.0) Dz_z=ODE_Solver_info.tolSol/2.0;
        
        do{ 
            if(ODE_Solver_info.direction==1) 
            {
                if(zin<0.0) zout=max(zin*(1.0+Dz_z), zend); 
                else if(zin>ODE_Solver_info.tolSol) zout=max(zin*(1.0-Dz_z), zend); 
                else zout=max(-ODE_Solver_info.tolSol, zend); 
            }
            else
            {
                if(zin<-ODE_Solver_info.tolSol) zout=min(zin*(1.0-Dz_z), zend); 
                else if(zin>0.0) zout=min(zin*(1.0+Dz_z), zend); 
                else zout=min(ODE_Solver_info.tolSol, zend); 
            }
            
            LAerror=ODE_Solver_compute_new_Solution(ODE_Solver_info.order, zout, 
                                                    ODE_Solver_info.tolSol, 
                                                    ODE_Solver_info, 
                                                    ODE_Solver_info.Snewptr, 
                                                    *ODE_Solver_info.Snptr, 
                                                    *ODE_Solver_info.Snm1ptr, 
                                                    *ODE_Solver_info.Snm2ptr, 
                                                    *ODE_Solver_info.Snm3ptr, 
                                                    *ODE_Solver_info.Snm4ptr, 
                                                    *ODE_Solver_info.Snm5ptr);
            
            if(LAerror==1)
            {
                LAerror=0;
                Dz_z/=2.0;
                Dz_z=max(Dz_z, 1.0e-10);
                zout=zin; 
                redo_run=1;
            }
            else if(LAerror==10){ global_ODE_Solver_info_ptr_for_function=NULL; return 10; }
            else if(LAerror==0)
            {
                //--------------------------------------------------------
                // estimate error and next possible stepsize
                //--------------------------------------------------------
                redo_run=0;
                Dz_z_old=Dz_z;
                redo_run=ODE_Solver_estimate_error_and_next_possible_step_size(ODE_Solver_info.order, 
                                                                               zin, zout, zend, Dz_z, 
                                                                               Dz_z_max, ODE_Solver_info, 
                                                                               *ODE_Solver_info.Snewptr, 
                                                                               *ODE_Solver_info.Snptr, 
                                                                               *ODE_Solver_info.Snm1ptr, 
                                                                               *ODE_Solver_info.Snm2ptr, 
                                                                               *ODE_Solver_info.Snm3ptr, 
                                                                               *ODE_Solver_info.Snm4ptr, 
                                                                               *ODE_Solver_info.Snm5ptr, 
                                                                               messflag);
                
                //--------------------------------------------------------
                // accepting the current step if 'redo_run==0'
                //--------------------------------------------------------
                if(redo_run==0)
                {
                    ODE_Solver_info.Sptr=ODE_Solver_info.Snm5ptr;
                    ODE_Solver_info.Snm5ptr=ODE_Solver_info.Snm4ptr;
                    ODE_Solver_info.Snm4ptr=ODE_Solver_info.Snm3ptr;
                    ODE_Solver_info.Snm3ptr=ODE_Solver_info.Snm2ptr;
                    ODE_Solver_info.Snm2ptr=ODE_Solver_info.Snm1ptr;
                    ODE_Solver_info.Snm1ptr=ODE_Solver_info.Snptr;
                    ODE_Solver_info.Snptr=ODE_Solver_info.Snewptr;
                    ODE_Solver_info.Snewptr=ODE_Solver_info.Sptr;
                    ODE_Solver_info.Sptr=NULL;
                    ODE_Solver_info.count++;
                    
                    if(ODE_Solver_info.count>=ODE_Solver_info.order+1)
                    { 
                        ODE_Solver_info.order=(int)min(ODE_Solver_info.order+1, maxorder); 
                        ODE_Solver_info.count=1; 
                    }
                    //
                    zin=zout;
                }
                else{ redo_run=1; ODE_Solver_info.Jac_is_set=0; }
            }
            else{ cout << " Unknown error in LA-solver accurred: " << LAerror << endl; exit(1); }
        }       
        while(ODE_Solver_info.direction*zout>ODE_Solver_info.direction*zend);
        
        //===========================================================================
        // if run was accepted then Snptr contains new solution!!!
        //===========================================================================
        ODE_Solver_info.Dz_z_last=Dz_z;
        //
        Sz.z=ODE_Solver_info.Snptr->z;
        Sz.y=ODE_Solver_info.Snptr->y;
        Sz.dy=ODE_Solver_info.Snptr->dy;
        
        global_ODE_Solver_info_ptr_for_function=NULL;
        return 0;
    }
        
}
