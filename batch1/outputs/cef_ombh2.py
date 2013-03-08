import planckStyle as s
from pylab import *

g=s.plotter

roots = ['test_H0_alpha_cosmoplanckmass_WMAP9', 'test_H0_alpha_cosmoplanckmass_lowllikeCAMspec_new', 'test_H0_alpha_cosmoplanckmass_lowllikeCAMspec_new_BAO']
params = g.get_param_array(roots[0], ['cef','omegabh2'])
g.settings.setWithSubplotSize(4)
g.settings.axes_fontsize = 24
g.settings.lab_fontsize = 32
g.settings.lw_contour = 0.2
g.alpha_filled_add = 1.0
g.setAxes(params,lims=[0.95, 1.05, 0.020, 0.026], pos=[0.15,0.15,0.60,0.80])
g.add_2d_contours(roots[0], params[0], params[1], filled=True, cols=('#F7BAA6','#E03424'), alpha=0.9)
g.add_2d_contours(roots[1], params[0], params[1], filled=True, cols=('#8CD3F5','#006FED'), alpha=0.7)
g.add_2d_contours(roots[2], params[0], params[1], filled=True, cols=('#3CBA75','#006B30'), alpha=0.7)

text(0.96, 0.0253, r'\textit{WMAP}9', color='#E03424',fontsize=24)
text(0.96, 0.0248, s.WP, color='#006FED',fontsize=24)
text(0.96, 0.0243, s.WP+'+BAO', color='#006B30',fontsize=24)

g.export('cef_ombh2')
