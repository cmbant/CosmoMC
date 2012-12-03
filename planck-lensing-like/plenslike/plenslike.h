#include <stdbool.h>

#define min(X,Y) ((X) < (Y) ? (X) : (Y))
#define max(X,Y) ((X) > (Y) ? (X) : (Y))
#define min3(X,Y,Z) ( min( X, min(Y, Z) ) )

typedef struct {
  int     ntrm;
  int     lmax;
  int     **s12L;
  double ***w12L;
} qest;

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

typedef struct {
  int     nbins;
  int     lmax;
  int     lmaxt;
  int     lmax1;
  int     lmax2;
  int     lmax3;
  int     lmax4;
  double  s4hat;
  double  s4std;
  int    *bin_lmins;
  int    *bin_lmaxs;
  double *bin_vals;
  double *mat_sigma;
  double *mat_sigma_inv;
  double *clpp_fid;
  double *vl_inv;
  double *rl_inv;
  double *sl_fid;
  double *cltt_fid;
  double *bl1n1_fid;
  double *bl2n1_fid;
  double *bl3n1_fid;
  double *bl4n1_fid;
  qest   *qe12;
  qest   *qe34;
} plenslike_dat_quad;

// qest.c
void free_qe( qest *qe );

void fill_qe_resp( int lmax, double *resp, qest *qee, qest *qes, 
		   double *f1, int f1lmax, double *f2, int f2lmax );
void fill_qe_covn( int lmax, double *covn, qest *q12, qest *q34, 
		   double *c1, int c1lmax, double *c2, int c2lmax, 
		   bool switch_34, bool mult_n1_l1l2L );

void init_qe_plm( qest *qe, int lmax, double *cltt);
void init_qe_slm( qest *qe, int lmax );
void init_qe_mlm( qest *qe, int lmax, double *cltt, double *fbl1, double *fbl2 );
void init_qe_nlm( qest *qe, int lmax, double *fbl1, double *fbl2 );

// wignerd.c
void init_gauss_legendre_quadrature( int n, double *x, double *w );
void wignerd_cf_from_cl( int s1, int s2, int nfunc, int ntheta, int lmax, const double *cos_theta, 
			 double *out_cf, const double *in_cl );
void wignerd_cl_from_cf( int s1, int s2, int nfunc, int ntheta, int lmax, const double *cos_theta, 
			 const double *integration_weights, double *out_cl, const double *in_cf );

// plenslike_dat_mono.c
void   load_plenslike_dat_mono( plenslike_dat_mono *dat, char *tfname );
void   free_plenslike_dat_mono( plenslike_dat_mono *dat );

void   fill_plenslike_mono_bins( plenslike_dat_mono *dat, double *bins, double *clpp );
void   fill_qe_plm_resp_plm_mono( int lmax, double *resp, double *cltt_fid, double *bl_fid, double *fl, double *cltt, double *bl);

double calc_plenslike_mono( plenslike_dat_mono *dat, double *clpp );
double calc_plenslike_mono_renorm( plenslike_dat_mono *dat, double *clpp, double *cltt, double *bl );

// plenslike_dat_quad.c
void   load_plenslike_dat_quad( plenslike_dat_quad *dat, char *tfname );
void   free_plenslike_dat_quad( plenslike_dat_quad *dat );

void   fill_plenslike_quad_bins( plenslike_dat_quad *dat, double *bins, double *clpp );

void   fill_quad_resp_pp_cltt( int lmax, double *rl, plenslike_dat_quad *dat, double *cltt );
void   fill_quad_resp_ss_bfid( int lmax, double *rl, plenslike_dat_quad *dat );
void   fill_quad_resp_qq_bfid( int lmax, double *rl, plenslike_dat_quad *dat, qest *qq );

double calc_plenslike_quad( plenslike_dat_quad *dat, double *clpp );
double calc_plenslike_quad_renorm_cltt( plenslike_dat_quad *dat, double *clpp, double *cltt );


