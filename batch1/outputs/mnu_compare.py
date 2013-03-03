import planckStyle as s
g=s.plotter

labels=None
roots=['base_mnu_planck_lowl_lowLike','base_mnu_planck_lowl_lowLike_lensing','base_mnu_Alens_planck_lowl_lowLike','base_mnu_planck_tauprior','base_planck_lowl_lowLike']
g.plots_1d(roots, legend_labels=labels)

g.exportExtra('mnu_compare.eps')

