import planckStyle as s
from pylab import *

g=s.plotter

roots = ['base_mnu_Alens_planck_lowl_lowLike_highL', 'base_mnu_Alens_planck_lowl_lowLike_highL_post_lensing']
params = g.get_param_array(roots[0], ['mnu','Alens'])
g.settings.setWithSubplotSize(4)
g.settings.axes_fontsize = 24
g.settings.lab_fontsize = 32
g.settings.lw_contour = 0.2
g.settings.alpha_filled_add = 1.0
g.setAxes(params,lims=[0.0, 2.0, 0.80, 1.80], pos=[0.15,0.15,0.60,0.80])
g.add_2d_contours(roots[0], params[0], params[1], filled=True, cols=('#F7BAA6','#E03424'), alpha=0.9)
g.add_2d_contours(roots[1], params[0], params[1], filled=True, cols=('#8CD3F5','#006FED'), alpha=0.7)

text(0.4, 0.93, s.WPhighL ,color='#E03424',fontsize=24)
text(0.4 ,0.85, s.WPhighLlensing ,color='#006FED',fontsize=24)

g.export('mnu_alens')
