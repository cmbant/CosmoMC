p=getdist_defaults();

pl=base_planck_CAMspec_lowl;
p2=base_planck_CAMspec_lowl_lowLike;
p3=base_planck_CAMspec_lowl_lowLike_highL;
p4=base_planck_CAMspec_lowl_lowLike_post_lensing;

getdist_plots_1D(p,pl,{pl.omegabh2,pl.omegach2,pl.tau,pl.ns,pl.H0,pl.omegam,pl.sigma8},p2,p3,p4);

getdist_print(p,'planck_datasets_1D');

