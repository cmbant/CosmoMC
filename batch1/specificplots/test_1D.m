p=getdist_defaults();

pl=base_omegak_planck_CAMspec_lowl_lowLike;
pl_bao=base_omegak_planck_CAMspec_lowl_lowLike_BAO;

getdist_plots_1D(p,pl,{pl.H0,pl.omegam,pl.omegabh2},pl_bao);

getdist_print(p,'test_1D');

