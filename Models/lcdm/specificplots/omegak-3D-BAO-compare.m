%set environment variable getdist_plot_data to your path to plot_data directory
%e.g. export getdist_plot_data=/scratch/aml1005/testchain/plot_data/

p=getdist_defaults();

pl=base_omegak_planck_CAMspec;
pl_bao=base_omegak_planck_CAMspec_BAO;

getdist_3D(p,pl,pl.H0,pl.omegam,pl.omegak, pl_bao);

getdist_print(p,'omegak-3D-BAO-compare',1,1);

