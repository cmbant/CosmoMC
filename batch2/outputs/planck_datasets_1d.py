import planckStyle as s
from pylab import *


g=s.getSubplotPlotter()

labels=[s.planck,s.lensing,s.WP,s.WPhighL]
roots=['base_planck_lowl','base_planck_lowl_post_lensing','base_planck_lowl_lowLike','base_planck_lowl_lowLike_highL']
g.plots_1d(roots, legend_labels=labels)
g.exportExtra('planck_datasets_1d')

#g.newPlot()
#g.plots_1d(roots,['omegabh2','omegach2','thetastar','ns','omegam','H0','tau','zrei','A','sigma8'],nx=2, legend_labels=labels)
#g.export('planck_datasets_1d_params')

g.newPlot()
g.settings.legend_frac_subplot_margin = 0.15
g.plots_1d(roots,['omegabh2','thetastar','A','tau','omegam','omegach2','ns','sigma8','zrei','H0'],nx=5, legend_ncol=4,legend_labels=labels, share_y=True)
g.export('planck_datasets_1d_params')

g.newPlot()
labels=[s.planck,s.lensing,s.WP]
roots=['base_planck_lowl','base_planck_lowl_post_lensing','base_planck_lowl_lowLike']

#g.plots_1d(roots,['tau','zrei','A','sigma8'],nx=2, legend_labels=labels)
#g.exportExtra('reion_tau_1D')
