#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "plenslike.h"

void load_plenslike_dat_quad( plenslike_dat_quad *dat, char *tfname ) {
  char line[32767];
  FILE *tf;
  int i, j, l, err, tmp;
  double dmp;
  
  tf = fopen(tfname, "r");
  assert( tf != NULL );    

  // bypass comments
  while (fgets(line, sizeof line, tf)) {
    if (*line == '#') {
      continue;
    } else {
      break;
    }
  }

  // read size header
  err = sscanf(line, "%d", &dat->nbins);
  err = fscanf(tf,   "%d", &dat->lmax);
  err = fscanf(tf,   "%d %d %d %d %d", &dat->lmaxt, &dat->lmax1, &dat->lmax2, &dat->lmax3, &dat->lmax4);
  err = fscanf(tf,   "%lf %lf", &dat->s4hat, &dat->s4std);

  // allocate memory
  dat->bin_lmins     = malloc( dat->nbins * sizeof(int) );
  dat->bin_lmaxs     = malloc( dat->nbins * sizeof(int) );
  dat->bin_vals      = malloc( dat->nbins * sizeof(double) );
  dat->mat_sigma     = malloc( dat->nbins * dat->nbins * sizeof(double) );
  dat->mat_sigma_inv = malloc( dat->nbins * dat->nbins * sizeof(double) );
  dat->clpp_fid      = malloc( (dat->lmax+1)  * sizeof(double) );
  dat->vl_inv        = malloc( (dat->lmax+1)  * sizeof(double) );
  dat->rl_inv        = malloc( (dat->lmax+1)  * sizeof(double) );
  dat->sl_fid        = malloc( (dat->lmax+1)  * sizeof(double) );
  dat->cltt_fid      = malloc( (dat->lmaxt+1) * sizeof(double) );
  dat->bl1n1_fid     = malloc( (dat->lmaxt+1) * sizeof(double) );
  dat->bl2n1_fid     = malloc( (dat->lmaxt+1) * sizeof(double) );
  dat->bl3n1_fid     = malloc( (dat->lmaxt+1) * sizeof(double) );
  dat->bl4n1_fid     = malloc( (dat->lmaxt+1) * sizeof(double) );

  // read bin info
  for (i=0; i<dat->nbins; i++) {
    err = fscanf(tf, "%d %d %d %lf", &tmp, &dat->bin_lmins[i], &dat->bin_lmaxs[i], &dat->bin_vals[i]);
    assert( tmp == i );
    assert( err == 4 );
  }

  // read sigma matrix
  for (i=0; i<dat->nbins; i++) {
    for (j=0; j<dat->nbins; j++) {
      err = fscanf(tf, "%lf", &dat->mat_sigma[i*dat->nbins + j]);
      assert( err == 1 );
    }
  }

  // read sigma inv matrix
  for (i=0; i<dat->nbins; i++) {
    for (j=0; j<dat->nbins; j++) {
      err = fscanf(tf, "%lf", &dat->mat_sigma_inv[i*dat->nbins + j]);
      assert( err == 1 );
    }
  }

  // read spectra
  for (l=0; l<=dat->lmax; l++) {
    err = fscanf(tf, "%lf %lf %lf %lf %lf", &dmp, &dat->clpp_fid[l], &dat->vl_inv[l], &dat->rl_inv[l], &dat->sl_fid[l]);
    assert( err == 5 );
    assert( dmp == l );
  }

  // read spectra
  for (l=0; l<=dat->lmaxt; l++) {
    err = fscanf(tf, "%lf %lf %lf %lf %lf %lf", &dmp, &dat->cltt_fid[l], &dat->bl1n1_fid[l], &dat->bl2n1_fid[l], &dat->bl3n1_fid[l], &dat->bl4n1_fid[l] );
    assert( err == 6 );
    assert( dmp == l );
  }

  // read qe12
  dat->qe12 = malloc( sizeof(qest) );
  err = fscanf(tf, "%d", &dat->qe12->ntrm);
  err = fscanf(tf, "%d", &dat->qe12->lmax);

  dat->qe12->s12L = malloc( dat->qe12->ntrm*sizeof(int *) );
  for (i=0; i < dat->qe12->ntrm; i++) {
    dat->qe12->s12L[i] = malloc( 3*sizeof(int) );
    err = fscanf(tf, "%d %d %d", &dat->qe12->s12L[i][0], &dat->qe12->s12L[i][1], &dat->qe12->s12L[i][2]);
    assert( err = 3 );
  }
  
  dat->qe12->w12L = malloc( dat->qe12->ntrm*sizeof(double **) );
  for (i=0; i < dat->qe12->ntrm; i++) {
    dat->qe12->w12L[i] = malloc( 3*sizeof(double *) );
    dat->qe12->w12L[i][0] = malloc( (dat->qe12->lmax+1)*sizeof(double) );
    dat->qe12->w12L[i][1] = malloc( (dat->qe12->lmax+1)*sizeof(double) );
    dat->qe12->w12L[i][2] = malloc( (dat->qe12->lmax+1)*sizeof(double) );

    for (l=0; l <= dat->qe12->lmax; l++) {
      err = fscanf(tf, "%lf %lf %lf %lf", &dmp, &dat->qe12->w12L[i][0][l], &dat->qe12->w12L[i][1][l], &dat->qe12->w12L[i][2][l] );
      assert( err == 4 );
      assert( dmp == l );
    }
  }
  // --

  // read qe34
  err = fscanf(tf, "%d", &tmp);

  if (tmp == -1) {
    dat->qe34 = dat->qe12;
  } else {
    dat->qe34 = malloc( sizeof(qest) );
    dat->qe34->ntrm = tmp;
    err = fscanf(tf, "%d", &dat->qe34->lmax);
    
    dat->qe34->s12L = malloc( dat->qe34->ntrm*sizeof(int *) );
    for (i=0; i < dat->qe34->ntrm; i++) {
      dat->qe34->s12L[i] = malloc( 3*sizeof(int) );
      err = fscanf(tf, "%d %d %d", &dat->qe34->s12L[i][0], &dat->qe34->s12L[i][1], &dat->qe34->s12L[i][2]);
      assert( err = 2 );
    }
    
    dat->qe34->w12L = malloc( dat->qe34->ntrm*sizeof(double **) );
    for (i=0; i < dat->qe34->ntrm; i++) {
      dat->qe34->w12L[i] = malloc( 3*sizeof(double *) );
      dat->qe34->w12L[i][0] = malloc( (dat->qe34->lmax+1)*sizeof(double) );
      dat->qe34->w12L[i][1] = malloc( (dat->qe34->lmax+1)*sizeof(double) );
      dat->qe34->w12L[i][2] = malloc( (dat->qe34->lmax+1)*sizeof(double) );
      
      for (l=0; l <= dat->qe34->lmax; l++) {
	err = fscanf(tf, "%lf %lf %lf %lf", &dmp, &dat->qe34->w12L[i][0][l], &dat->qe34->w12L[i][1][l], &dat->qe34->w12L[i][2][l] );
	assert( err == 4 );
	assert( dmp == l );
      }
    }
  }
  // --
  
  fclose(tf);
  return;
};

void free_plenslike_dat_quad( plenslike_dat_quad *dat ) {
  free(dat->bin_lmins);
  free(dat->bin_lmaxs);
  free(dat->bin_vals);
  free(dat->mat_sigma);
  free(dat->mat_sigma_inv);
  free(dat->clpp_fid);
  free(dat->vl_inv);
  free(dat->rl_inv);
  free(dat->sl_fid);
  free(dat->cltt_fid);
  free(dat->bl1n1_fid);
  free(dat->bl2n1_fid);
  free(dat->bl3n1_fid);
  free(dat->bl4n1_fid);

  free_qe(dat->qe12);

  if (dat->qe34 != dat->qe12) {
    free_qe(dat->qe34);
    free(dat->qe34);
  }

  free(dat->qe12);
};

void fill_plenslike_quad_bins( plenslike_dat_quad *dat, double *bins, double *clpp ) {
  int i, l;
  double num, den;

  for (i=0; i<dat->nbins; i++) {
    num = 0; den=0;
    for (l=dat->bin_lmins[i]; l<=dat->bin_lmaxs[i]; l++) {
      num += clpp[l] * dat->clpp_fid[l] * dat->vl_inv[l];
      den += dat->clpp_fid[l] * dat->clpp_fid[l] * dat->vl_inv[l];
    }
    bins[i] = num/den;
  }
}

void   fill_quad_resp_pp_cltt( int lmax, double *rl, plenslike_dat_quad *dat, double *cltt ) {
  int lmax_fl = max( max(dat->lmax1, dat->lmax2), max(dat->lmax3, dat->lmax4) );
  qest tqe;
  init_qe_plm( &tqe, max(lmax, lmax_fl), cltt );
  fill_quad_resp_qq_bfid( lmax, rl, dat, &tqe );

  free_qe(&tqe);
}

void   fill_quad_resp_ss_bfid( int lmax, double *rl, plenslike_dat_quad *dat ) {
  int lmax_fl = max( max(dat->lmax1, dat->lmax2), max(dat->lmax3, dat->lmax4) );
  qest tqe;
  init_qe_slm( &tqe, max(lmax, lmax_fl) );
  fill_quad_resp_qq_bfid( lmax, rl, dat, &tqe );

  free_qe(&tqe);
}

void   fill_quad_resp_qq_bfid( int lmax, double *rl, plenslike_dat_quad *dat, qest *qq ) {
  int l, lmax_fl;
  double *fl, *resp12, *resp34;

  lmax_fl = max( max(dat->lmax1, dat->lmax2), max(dat->lmax3, dat->lmax4) );
  fl = malloc( (lmax_fl+1) * sizeof(double) ); for (l=0; l<= lmax_fl; l++) { fl[l] = 1.0; }

  resp12 = malloc( (lmax+1) * sizeof(double) ); memset( resp12, 0, (lmax+1) * sizeof(double) );
  fill_qe_resp(lmax, resp12, dat->qe12, qq, fl, dat->lmax1, fl, dat->lmax2);

  if ((dat->qe12 == dat->qe34) && (dat->lmax1 == dat->lmax3) && (dat->lmax2 == dat->lmax4)) {
    resp34 = resp12;
  } else {
    resp34 = malloc( (lmax+1) * sizeof(double) ); memset( resp34, 0, (lmax+1) * sizeof(double) );
    fill_qe_resp(lmax, resp34, dat->qe34, qq, fl, dat->lmax3, fl, dat->lmax4);
  }

  for (l=0; l<=lmax; l++) {
    rl[l] = resp12[l] * resp34[l];
  }

  if (resp12 != resp34) { free(resp34); }
  free(resp12);
  free(fl);
}

double calc_plenslike_quad( plenslike_dat_quad *dat, double *clpp ) {
  int i, j;
  double ret;
  double *bins = malloc( dat->nbins * sizeof(double) );

  fill_plenslike_quad_bins(dat, bins, clpp);

  ret = 0.0;
  for (i=0; i<dat->nbins; i++) {
    for (j=0; j<dat->nbins; j++) {
      ret += (bins[i] - dat->bin_vals[i]) * (bins[j] - dat->bin_vals[j]) * dat->mat_sigma_inv[i*dat->nbins + j];
    }
  }

  free(bins);

  return -0.5*ret;
};

double calc_plenslike_quad_renorm_cltt( plenslike_dat_quad *dat, double *clpp, double *cltt ) {
  int i, j, l;
  double ret;
  double *bins = malloc( dat->nbins * sizeof(double) );
  double *clpp_renorm = malloc( (dat->lmax+1) * sizeof(double) );
  double *qc_resp = malloc( (dat->lmax+1) * sizeof(double) );

  fill_quad_resp_pp_cltt( dat->lmax, qc_resp, dat, cltt );

  for (l=0; l<=dat->lmax; l++) {
    clpp_renorm[l] = qc_resp[l] * dat->rl_inv[l] * clpp[l];
  }

  fill_plenslike_quad_bins(dat, bins, clpp_renorm);

  ret = 0.0;
  for (i=0; i<dat->nbins; i++) {
    for (j=0; j<dat->nbins; j++) {
      ret += (bins[i] - dat->bin_vals[i]) * (bins[j] - dat->bin_vals[j]) * dat->mat_sigma_inv[i*dat->nbins + j];
    }
  }

  free(bins);
  free(clpp_renorm);
  free(qc_resp);

  return -0.5*ret;
}
