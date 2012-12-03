#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "plenslike.h"

double calc_plenslike_mono( plenslike_dat_mono *dat, double *clpp ) {
  int i, j;
  double ret;
  double *bins = malloc( dat->nbins * sizeof(double) );

  fill_plenslike_mono_bins(dat, bins, clpp);

  ret = 0.0;
  for (i=0; i<dat->nbins; i++) {
    for (j=0; j<dat->nbins; j++) {
      ret += (bins[i] - dat->bin_vals[i]) * (bins[j] - dat->bin_vals[j]) * dat->mat_sigma_inv[i*dat->nbins + j];
    }
  }

  free(bins);

  return -0.5*ret;
};

double calc_plenslike_mono_renorm( plenslike_dat_mono *dat, double *clpp, double *cltt, double *bl ) {
  int i, j, l;
  double ret;
  double *bins = malloc( dat->nbins * sizeof(double) );
  double *clpp_renorm = malloc( (dat->lmax+1) * sizeof(double) );
  double *resp = malloc( (dat->lmax+1) * sizeof(double) );

  fill_qe_plm_resp_plm_mono( dat->lmax, resp, dat->cltt_fid, dat->bl_fid, dat->fl, cltt, bl );

  for (l=0; l<=dat->lmax; l++) {
    clpp_renorm[l] = resp[l] * resp[l] / dat->al_inv[l] / dat->al_inv[l] * clpp[l];
  }

  fill_plenslike_mono_bins(dat, bins, clpp_renorm);

  ret = 0.0;
  for (i=0; i<dat->nbins; i++) {
    for (j=0; j<dat->nbins; j++) {
      ret += (bins[i] - dat->bin_vals[i]) * (bins[j] - dat->bin_vals[j]) * dat->mat_sigma_inv[i*dat->nbins + j];
    }
  }

  free(bins);
  free(clpp_renorm);
  free(resp);

  return -0.5*ret;
}

void fill_plenslike_mono_bins( plenslike_dat_mono *dat, double *bins, double *clpp ) {
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

void fill_qe_plm_resp_plm_mono( int lmax, double *resp, double *cltt_fid, double *bl_fid, double *fl, double *cltt, double *bl) {
  int i, l, ngl;
  double fbl, *glz, *glw;
  double *c11, *c10, *c01, *c00;
  double *z11p, *z11m, *z10p, *z10m, *z01p, *z00p;
  double *ztp, *ztm, *ctp, *ctm;

  ngl = (3*lmax)/2 + 1;
  glz = malloc( ngl*sizeof(double) );
  glw = malloc( ngl*sizeof(double) );

  init_gauss_legendre_quadrature(ngl, glz, glw);

  c11 = malloc( (lmax+1)*sizeof(double) );
  c10 = malloc( (lmax+1)*sizeof(double) );
  c01 = malloc( (lmax+1)*sizeof(double) );
  c00 = malloc( (lmax+1)*sizeof(double) );

  for (l=0; l<=lmax; l++) {
    fbl     = fl[l] * bl[l] / bl_fid[l];
    c11[l] = (2.*l+1.) * l*(l+1.) * cltt_fid[l] * cltt[l] * fbl;
    c10[l] = (2.*l+1.) * sqrt(l*(l+1.)) * cltt_fid[l] * fbl;
    c01[l] = (2.*l+1.) * sqrt(l*(l+1.)) * cltt[l] * fbl;
    c00[l] = (2.*l+1.) * fbl;
  }

  z11p = malloc( ngl*sizeof(double) ); z11m = malloc( ngl*sizeof(double) );
  z10p = malloc( ngl*sizeof(double) ); z10m = malloc( ngl*sizeof(double) );
  z01p = malloc( ngl*sizeof(double) );
  z00p = malloc( ngl*sizeof(double) );

  wignerd_cf_from_cl(1, 1, 1, ngl, lmax, glz, z11p, c11); wignerd_cf_from_cl(1, -1, 1, ngl, lmax, glz, z11m, c11);
  wignerd_cf_from_cl(0, 1, 1, ngl, lmax, glz, z10p, c10); wignerd_cf_from_cl(0, -1, 1, ngl, lmax, glz, z10m, c10);
  wignerd_cf_from_cl(0, 1, 1, ngl, lmax, glz, z01p, c01);
  wignerd_cf_from_cl(0, 0, 1, ngl, lmax, glz, z00p, c00);

  free(c11); free(c00); free(c10); free(c01);

  ztp = malloc( ngl*sizeof(double) ); ztm = malloc( ngl*sizeof(double) );
  for (i=0; i < ngl; i++) {
    ztp[i] = z11p[i] * z00p[i] - z01p[i] * z10p[i];
    ztm[i] = z11m[i] * z00p[i] - z01p[i] * z10m[i];
  }

  free(z11p); free(z00p); free(z10p); free(z01p); free(z11m); free(z10m);
  ctp = malloc( (lmax+1)*sizeof(double) ); ctm = malloc( (lmax+1)*sizeof(double) );

  wignerd_cl_from_cf(1, +1, 1, ngl, lmax, glz, glw, ctp, ztp);
  wignerd_cl_from_cf(1, -1, 1, ngl, lmax, glz, glw, ctm, ztm);

  free(ztp); free(ztm);
  free(glz); free(glw);

  for (l=0; l <= lmax; l++) {
    resp[l] = (ctp[l] + ctm[l]) * l*(l+1.) / (16. * M_PI);
  }

  free(ctp); free(ctm);

  return;
};

void load_plenslike_dat_mono( plenslike_dat_mono *dat, char *tfname ) {
  char line[32767];
  FILE *tf;
  int i, j, err, tmp;
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
  err = fscanf(tf, "%d", &dat->lmax);

  // allocate memory
  dat->bin_lmins     = malloc( dat->nbins * sizeof(int) );
  dat->bin_lmaxs     = malloc( dat->nbins * sizeof(int) );
  dat->bin_vals      = malloc( dat->nbins * sizeof(double) );
  dat->mat_sigma     = malloc( dat->nbins * dat->nbins * sizeof(double) );
  dat->mat_sigma_inv = malloc( dat->nbins * dat->nbins * sizeof(double) );
  dat->clpp_fid      = malloc( (dat->lmax+1) * sizeof(double) );
  dat->cltt_fid      = malloc( (dat->lmax+1) * sizeof(double) );
  dat->bl_fid        = malloc( (dat->lmax+1) * sizeof(double) );
  dat->fl            = malloc( (dat->lmax+1) * sizeof(double) );
  dat->vl_inv        = malloc( (dat->lmax+1) * sizeof(double) );
  dat->al_inv        = malloc( (dat->lmax+1) * sizeof(double) );

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
  for (i=0; i<dat->lmax+1; i++) {
    err = fscanf(tf, "%lf %lf %lf %lf %lf %lf %lf", &dmp, &dat->clpp_fid[i], &dat->cltt_fid[i], &dat->bl_fid[i], &dat->fl[i], &dat->vl_inv[i], &dat->al_inv[i]);
    assert( err == 7 );
    assert( dmp == i );
  }
  
  fclose(tf);
  return;
};

void free_plenslike_dat_mono( plenslike_dat_mono *dat ) {
  free(dat->bin_lmins);
  free(dat->bin_lmaxs);
  free(dat->bin_vals);
  free(dat->mat_sigma);
  free(dat->mat_sigma_inv);
  free(dat->clpp_fid);
  free(dat->cltt_fid);
  free(dat->bl_fid);
  free(dat->fl);
  free(dat->vl_inv);
  free(dat->al_inv);
};
