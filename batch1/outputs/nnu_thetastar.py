import planckStyle as s
from pylab import *

g=s.getSinglePlotter(ratio=1)

roots = ['base_nnu_yhe_planck_lowl_lowLike_highL', 'base_nnu_planck_lowl_lowLike_highL']


g.plot_2d(roots, param_pair=['nnu','thetastar'], filled=True,lims=[1.0, 6.0, 1.036, 1.047])

g.add_legend([ r'${\Lambda}CDM+N_{\rm eff}+Y_p$',r'${\Lambda}CDM+N_{\rm eff}$'],legend_loc='upper right',colored_text=True);

text(1.2, 1.0365, s.WPhighL, color='#000000')

g.export('neff_thetas')
