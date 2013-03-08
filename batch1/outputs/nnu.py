import planckStyle as s
from pylab import *

g=s.plotter
g.settings.line_labels = False
g.settings.setWithSubplotSize(4)

g.settings.tight_layout=True
g.settings.lab_fontsize=32
g.settings.axes_fontsize = 24

g.setAxes(lims=[2.0, 5.0, 0.0, 1.1])

roots=['base_nnu_planck_lowl_lowLike_highL','base_nnu_planck_lowl_lowLike_highL_BAO', 'base_nnu_planck_lowl_lowLike_highL_post_HST', 'base_nnu_planck_lowl_lowLike_highL_BAO_post_HST']

g.add_1d(roots[0],'nnu');
g.add_1d(roots[1],'nnu',color='b');
g.add_1d(roots[2],'nnu',color='r');
g.add_1d(roots[3],'nnu',color='g');

text(2.1,1.02, s.WPhighL, fontsize=24)
text(2.1,0.92, '+BAO', color='b', fontsize=24)
text(2.1,0.82, r'+$H_0$', color='r', fontsize=24)
text(2.1,0.72, '+BAO+$H_0$', color='g', fontsize=24)

xlabel(r'$N_{\rm eff}$',fontsize=32)
ylabel(r'$P/P_{\rm max}$',fontsize=32)

g.export('nnu')
