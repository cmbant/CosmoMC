import planckStyle as s
from pylab import *

g=s.getSinglePlotter()

roots = ['base_nnu_planck_lowl_lowLike','base_nnu_planck_lowl_lowLike_highL','base_nnu_planck_lowl_lowLike_highL_BAO']

g.plot_2d(roots, param_pair=['ns','nnu'], filled=True,lims=[0.92, 1.021 , 2.2, 4.8])

g.add_y_marker(3.046)
g.add_x_marker(1)

g.add_legend([s.WP,s.WPhighL,s.WPhighL+'+BAO'],legend_loc='upper left',colored_text=True);

xticks([0.92,0.94,0.96,0.98,1.00,1.02])

g.export('ns_nnu_0')
