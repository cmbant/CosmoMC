import planckStyle as s
from pylab import *

g=s.plotter

roots = ['base_nnu_yhe_planck_lowl_lowLike_highL', 'base_nnu_planck_lowl_lowLike_highL']
params = g.get_param_array(roots[0], ['nnu','thetastar'])
g.settings.setWithSubplotSize(4)
g.settings.axes_fontsize = 24
g.settings.lab_fontsize = 32
g.settings.lw_contour = 0.2
g.alpha_filled_add = 1.0
g.setAxes(params,lims=[1.0, 6.0, 1.036, 1.047], pos=[0.15,0.15,0.60,0.80])
g.add_2d_contours(roots[0], params[0], params[1], filled=True, cols=('#F7BAA6','#E03424'), alpha=0.9)
g.add_2d_contours(roots[1], params[0], params[1], filled=True, cols=('#8CD3F5','#006FED'), alpha=0.7)

text(3.0, 1.046, r'${\Lambda}CDM+N_{\rm eff}+Y_p$' ,color='#E03424',fontsize=24)
text(3.0, 1.045, r'${\Lambda}CDM+N_{\rm eff}$' ,color='#006FED',fontsize=24)
text(1.2, 1.0365, s.WPhighL, color='#000000',fontsize=24)

g.export('neff_thetas')
