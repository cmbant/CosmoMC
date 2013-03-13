import planckStyle as s
from pylab import *

g=s.getSinglePlotter()

g.settings.lineM = ['-k', '-b', '-r', '-g', '-m', '-y']

roots=['base_nnu_planck_lowl_lowLike_highL','base_nnu_planck_lowl_lowLike_highL_BAO', 'base_nnu_planck_lowl_lowLike_highL_post_HST', 'base_nnu_planck_lowl_lowLike_highL_BAO_post_HST']

par = g.check_param(roots[0],'nnu')

#g.setAxes([par],lims=[2.0, 5.0, 0.0, 1.1])

g.plot_1d(roots,par);

xlim([2.0,4.8])
sz=9
text(2.1,1.01, s.WPhighL, fontsize=9)
text(2.1,0.92, '+BAO', color='b', fontsize=9)
text(2.1,0.82, r'+$H_0$', color='r', fontsize=9)
text(2.1,0.72, '+BAO+$H_0$', color='g', fontsize=9)

g.export('nnu')
