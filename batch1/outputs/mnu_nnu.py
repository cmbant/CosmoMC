import planckStyle as s
from pylab import *

g=s.getSinglePlotter()

roots = ['base_nnu_mnu_planck_lowl_lowLike_highL','base_nnu_mnu_planck_lowl_lowLike_highL_post_BAO']
params = g.get_param_array(roots[0], ['mnu','nnu'])

g.plot_2d(roots, param_pair=['mnu','nnu'], filled=True,lims=[0.0, 1.0, 2.0, 5.0])

g.add_legend([s.WPhighL,s.WPhighL+'+BAO'],legend_loc='upper right',colored_text=True);


g.export('mnu_nnu')
