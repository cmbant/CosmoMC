#import planckStyle as s
#g=s.plotter

import GetDistPlots
g=GetDistPlots.GetDistPlotter('main/plot_data')

labels=None
roots=['base_mnu_planck_lowl_lowLike','base_mnu_planck_lowl_lowLike_lensing','base_mnu_Alens_planck_lowl_lowLike','base_mnu_Alens_planck_lowl_lowLike_post_lensing','base_mnu_planck_tauprior','base_mnu_planck_tauprior_post_lensing']
g.plots_1d(roots, legend_labels=labels)

g.export('plots/mnu_compare.pdf')

