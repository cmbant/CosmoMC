#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "plenslike.h"

void free_qe(qest *qe) {
  int i;

  for (i=0; i < qe->ntrm; i++) {
    free( qe->w12L[i][2] );
    free( qe->w12L[i][1] );
    free( qe->w12L[i][0] );
    free( qe->w12L[i] );
    free( qe->s12L[i] );
  }
  
  free( qe->w12L );
  free( qe->s12L );
}

void fill_qe_resp(int lmax, double *resp, qest *qee, qest *qes, 
		  double *f1, int f1lmax, double *f2, int f2lmax) {
  int l;
  fill_qe_covn( lmax, resp, qee, qes, f1, f1lmax, f2, f2lmax, false, false);
  for (l=0; l<=lmax; l++) {
    resp[l] *= 2.0;
  }
}

void fill_qe_covn(int lmax, double *covn, qest *q12, qest *q34, 
		  double *c1, int c1lmax, double *c2, int c2lmax, 
		  bool switch_34, bool mult_n1_l1l2L) {

  int i, j, tl1, tl2, tl, ip;
  int ngl, tl1min, tl2min, tl1max, tl2max;
  double *glz, *glw, *cl1, *cl2, *clL, *gp1, *gp2, *gp12;

  int i34_0 = switch_34 ? 1 : 0;
  int i34_1 = switch_34 ? 0 : 1;

  tl1max = min3(q12->lmax, q34->lmax, c1lmax);
  tl2max = min3(q12->lmax, q34->lmax, c2lmax);
  ngl  = (tl1max + tl2max + lmax)/2 + 1;

  cl1  = malloc( (tl1max+1)*sizeof(double) );
  cl2  = malloc( (tl2max+1)*sizeof(double) );
  clL  = malloc( (lmax+1)*sizeof(double) );
  glz  = malloc( ngl*sizeof(double) ); 
  glw  = malloc( ngl*sizeof(double) );
  gp1  = malloc( ngl*sizeof(double) ); 
  gp2  = malloc( ngl*sizeof(double) );
  gp12 = malloc( ngl*sizeof(double) );

  memset(covn, 0, (lmax+1)*sizeof(double));
  for (i=0; i < q12->ntrm; i++) {
    for (j=0; j < q34->ntrm; j++) {
      tl1min = max(abs(q12->s12L[i][0]), abs(q34->s12L[j][i34_0]));
      memset(cl1, 0, (tl1max+1)*sizeof(double));
      for (tl1=tl1min; tl1 <= tl1max; tl1++) {
	cl1[tl1] = 
	  ((mult_n1_l1l2L && tl1%2) ? -1. : 1.) * 
	  q12->w12L[i][0][tl1] *
	  q34->w12L[j][i34_0][tl1] *
	  (2.*tl1+1.) * c1[tl1];
      }
      
      tl2min = max(abs(q12->s12L[i][1]), abs(q34->s12L[j][i34_1]));
      memset(cl2, 0, (tl2max+1)*sizeof(double));
      for (tl2=tl2min; tl2 <= tl2max; tl2++) {
	cl2[tl2] = 
	  ((mult_n1_l1l2L && tl2%2) ? -1. : 1.) * 
	  q12->w12L[i][1][tl2] * 
	  q34->w12L[j][i34_1][tl2] *
	  (2.*tl2+1.) * c2[tl2];
      }
      
      init_gauss_legendre_quadrature(ngl, glz, glw);
      wignerd_cf_from_cl(q12->s12L[i][0], q34->s12L[j][i34_0], 1, ngl, tl1max, glz, gp1, cl1);
      wignerd_cf_from_cl(q12->s12L[i][1], q34->s12L[j][i34_1], 1, ngl, tl2max, glz, gp2, cl2);
      
      for (ip=0; ip < ngl; ip++) {
	gp12[ip] = gp1[ip] * gp2[ip];
      }
      
      memset(clL, 0, (lmax+1)*sizeof(double));
      wignerd_cl_from_cf(q12->s12L[i][2], q34->s12L[j][2], 1, ngl, lmax, glz, glw, clL, gp12);
      
      for (tl=0; tl <= lmax; tl++) {
	covn[tl] += 
	  clL[tl] *
	  q12->w12L[i][2][tl] * 
	  q34->w12L[j][2][tl] * 
	  ( (mult_n1_l1l2L && (tl % 2)) ? -1. : 1.) / (32.*M_PI);
      }
    }
  }

  free (cl1); free(cl2); free(clL); free(glz); free(glw); free(gp1); free(gp2); free(gp12);
}

void init_qe_plm( qest *qe, int lmax, double *cltt) {
  int i, l;
  double sq, ln;

  qe->ntrm = 4;
  qe->lmax = lmax;

  qe->s12L = malloc( qe->ntrm*sizeof(int *) );
  qe->w12L = malloc( 4*sizeof(double *) );
  for (i=0; i < qe->ntrm; i++) {
    qe->s12L[i] = malloc( 3*sizeof(int) );
    qe->w12L[i] = malloc( 3*sizeof(double *) );
    qe->w12L[i][0] = malloc( (lmax+1)*sizeof(double) );
    qe->w12L[i][1] = malloc( (lmax+1)*sizeof(double) );
    qe->w12L[i][2] = malloc( (lmax+1)*sizeof(double) );
  }

  qe->s12L[0][0] = 1; qe->s12L[0][1] = 0; qe->s12L[0][2] = 1;
  qe->s12L[1][0] = 0; qe->s12L[1][1] = 1; qe->s12L[1][2] = 1;
  qe->s12L[2][0] = 1; qe->s12L[2][1] = 0; qe->s12L[2][2] = 1;
  qe->s12L[3][0] = 0; qe->s12L[3][1] = 1; qe->s12L[3][2] = 1;

  for (l=0; l<=lmax; l++) {
    sq = sqrt(1.0*l*(l+1.));
    ln = (l%2 ? -1. : 1.);
    qe->w12L[0][0][l] = sq * cltt[l];      qe->w12L[0][1][l] = -0.5;              qe->w12L[0][2][l] = sq;
    qe->w12L[1][0][l] = -0.5;              qe->w12L[1][1][l] = sq * cltt[l];      qe->w12L[1][2][l] = sq;
    qe->w12L[2][0][l] = sq * cltt[l] * ln; qe->w12L[2][1][l] = -0.5 * ln;         qe->w12L[2][2][l] = sq* ln;
    qe->w12L[3][0][l] = -0.5 * ln;         qe->w12L[3][1][l] = sq * cltt[l] * ln; qe->w12L[3][2][l] = sq* ln;
  }
}

void init_qe_slm( qest *qe, int lmax) {
  int i, l;

  qe->ntrm = 1;
  qe->lmax = lmax;

  qe->s12L = malloc( qe->ntrm*sizeof(int *) );
  qe->w12L = malloc( qe->ntrm*sizeof(double *) );
  for (i=0; i < qe->ntrm; i++) {
    qe->s12L[i] = malloc( 3*sizeof(int) );
    qe->w12L[i] = malloc( 3*sizeof(double *) );
    qe->w12L[i][0] = malloc( (lmax+1)*sizeof(double) );
    qe->w12L[i][1] = malloc( (lmax+1)*sizeof(double) );
    qe->w12L[i][2] = malloc( (lmax+1)*sizeof(double) );
  }

  qe->s12L[0][0] = 0; qe->s12L[0][1] = 0; qe->s12L[0][2] = 0;

  for (l=0; l<=lmax; l++) {
    qe->w12L[0][0][l] = 1.0; qe->w12L[0][1][l] = 1.0; qe->w12L[0][2][l] = 1.0;
  }
}

void init_qe_mlm( qest *qe, int lmax, double *cltt, double *fbl1, double *fbl2) {
  int i, l;

  qe->ntrm = 2;
  qe->lmax = lmax;

  qe->s12L = malloc( qe->ntrm*sizeof(int *) );
  qe->w12L = malloc( qe->ntrm*sizeof(double *) );
  for (i=0; i < qe->ntrm; i++) {
    qe->s12L[i] = malloc( 3*sizeof(int) );
    qe->w12L[i] = malloc( 3*sizeof(double *) );
    qe->w12L[i][0] = malloc( (lmax+1)*sizeof(double) );
    qe->w12L[i][1] = malloc( (lmax+1)*sizeof(double) );
    qe->w12L[i][2] = malloc( (lmax+1)*sizeof(double) );
  }

  qe->s12L[0][0] = 0; qe->s12L[0][1] = 0; qe->s12L[0][2] = 0;
  qe->s12L[1][0] = 0; qe->s12L[1][1] = 0; qe->s12L[1][2] = 0;

  for (l=0; l<=lmax; l++) {
    qe->w12L[0][0][l] = cltt[l] / fbl2[l];
    qe->w12L[0][1][l] = fbl1[l];
    qe->w12L[0][2][l] = 1.0;

    qe->w12L[0][0][l] = fbl1[l];
    qe->w12L[0][1][l] = cltt[l] / fbl2[l];
    qe->w12L[0][2][l] = 1.0;
  }
}

void init_qe_nlm( qest *qe, int lmax, double *fbl1, double *fbl2) {
  int i, l;

  qe->ntrm = 1;
  qe->lmax = lmax;

  qe->s12L = malloc( qe->ntrm*sizeof(int *) );
  qe->w12L = malloc( qe->ntrm*sizeof(double *) );
  for (i=0; i < qe->ntrm; i++) {
    qe->s12L[i] = malloc( 3*sizeof(int) );
    qe->w12L[i] = malloc( 3*sizeof(double *) );
    qe->w12L[i][0] = malloc( (lmax+1)*sizeof(double) );
    qe->w12L[i][1] = malloc( (lmax+1)*sizeof(double) );
    qe->w12L[i][2] = malloc( (lmax+1)*sizeof(double) );
  }

  qe->s12L[0][0] = 0; qe->s12L[0][1] = 0; qe->s12L[0][2] = 0;

  for (l=0; l<=lmax; l++) {
    qe->w12L[0][0][l] = fbl1[l];
    qe->w12L[0][1][l] = fbl2[l];
    qe->w12L[0][2][l] = 1.0;
  }
}
