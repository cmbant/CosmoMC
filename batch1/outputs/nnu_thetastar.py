import planckStyle as s
from pylab import *

g=s.getSinglePlotter(ratio=1)

roots = ['base_nnu_yhe_planck_lowl_lowLike_highL', 'base_nnu_planck_lowl_lowLike_highL']


g.plot_2d(roots, param_pair=['nnu','thetastar'], filled=True,lims=[1.0, 6.0, 1.036, 1.047])

nnu = g.param_latex_label(roots[0], 'nnu')
yhe = g.param_latex_label(roots[0], 'yheused')

g.add_legend([ s.LCDM + '+'+nnu+'+'+ yhe,s.LCDM+ '+'+nnu],legend_loc='upper right',colored_text=True);

text(1.2, 1.0365, s.WPhighL, color='#000000', fontsize=g.settings.legend_fontsize)

g.export('neff_thetas')
