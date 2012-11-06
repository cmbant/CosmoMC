typedef struct {
  int     nbins;
  int     lmax;
  int    *bin_lmins;
  int    *bin_lmaxs;
  double *bin_vals;
  double *mat_sigma;
  double *mat_sigma_inv;
  double *clpp_fid;
  double *cltt_fid;
  double *bl_fid;
  double *fl;
  double *vl_inv;
  double *al_inv;
} plenslike_dat_mono;

// wignerd.c
void init_gauss_legendre_quadrature(int n, double *x, double *w);
void wignerd_cf_from_cl(int s1, int s2, int nfunc, int ntheta, int lmax, const double *cos_theta, double *out_cf, const double *in_cl);
void wignerd_cl_from_cf(int s1, int s2, int nfunc, int ntheta, int lmax, const double *cos_theta, const double *integration_weights, double *out_cl, const double *in_cf);

// plenslike.c
double calc_plenslike_mono( plenslike_dat_mono *dat, double *clpp );
double calc_plenslike_mono_renorm( plenslike_dat_mono *dat, double *clpp, double *cltt, double *bl );

void fill_qe_plm_resp_plm_mono( int lmax, double *resp, double *cltt_fid, double *bl_fid, double *fl, double *cltt, double *bl);
void fill_plenslike_mono_bins( plenslike_dat_mono *dat, double *clpp, double *bins );

void load_plenslike_dat_mono( plenslike_dat_mono *dat, char *tfname );
void free_plenslike_dat_mono( plenslike_dat_mono *dat );

