%set environment variable getdist_plot_data to your path to plot_data directory
%e.g. export getdist_plot_data=/scratch/aml1005/testchain/plot_data/

p=getdist_defaults();

pl=base_planck_CAMspec;
pl_bao=base_planck_CAMspec_post_BAO;

%getdist_3D(p,pl,pl.omegam,pl.H0,pl.ns, pl_bao);
%getdist_print(p,'omegam-H0-ns_planck_bao',1,1);

getdist_3D(p,pl,pl.omegam,pl.H0,pl.ns, pl_bao);
getdist_print(p,'omegam-H0-ns_planck_bao',1,1);
