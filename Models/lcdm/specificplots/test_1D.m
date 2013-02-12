p=getdist_defaults();

pl=base_omegak_planck_CAMspec;
pl_bao=base_omegak_planck_CAMspec_BAO;

getdist_plots_1D(p,pl,{pl.H0,pl.omegam,pl.omegabh2},pl_bao);

getdist_print(p,'test_1D');

