import planckStyle as s
from pylab import *

g=s.getSinglePlotter(ratio=1)

roots = ['base_mnu_Alens_planck_lowl_lowLike_highL', 'base_mnu_Alens_planck_lowl_lowLike_highL_post_lensing']

g.plot_2d(roots, param_pair=['mnu','Alens'], filled=True,lims=[0.0, 2.0, 0.85, 1.85])
g.add_legend([s.WPhighL,s.WPhighLlensing]);

g.export('mnu_alens')
