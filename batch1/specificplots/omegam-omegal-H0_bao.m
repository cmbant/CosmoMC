%set environment variable getdist_plot_data to your path to plot_data directory
%e.g. export getdist_plot_data=/scratch/aml1005/testchain/plot_data/

p=getdist_defaults();
clf;
pl=base_omegak_planck_CAMspec;
pl_bao=base_omegak_planck_CAMspec_BAO;

%getdist_solid_2D(p,pl_bao,pl.omegam,pl.omegal,'-k','gc');

getdist_3D(p,pl,pl.omegam,pl.omegal,pl.H0,pl_bao);


xlim([0,1]);
ylim([0,1]);
line([1 0],[0 1],'color','k');

getdist_print(p,'omegam-omegal-H0_planck_bao',1,1);

