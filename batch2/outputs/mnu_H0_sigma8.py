import planckStyle as s
from pylab import *

g=s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)

base = 'base_mnu_'
roots=[base+s.defdata_TT]


mnu=[0, 2.3]
H=70.6
sigma=3.3
c='gray'
one=array([1,1])
fill_between(mnu, one*(H-sigma*2), one*(H+sigma*2), facecolor=c, alpha=0.1, edgecolor=c, lw=0)
fill_between(mnu, one*(H-sigma), one*(H+sigma), facecolor=c, alpha=0.15, edgecolor=c, lw=0)

  
g.plot_3d(roots, ['mnu', 'H0', 'sigma8'])

g.add_2d_contours(base+s.defdata_all+'_lensing', 'mnu', 'H0', plotno=0)
g.add_2d_contours(base+s.defdata_all+'_lensing_post_BAO_HST70p6_JLA', 'mnu', 'H0', filled=True, zorder=2, alpha=0.6)

#gca().set_xticks([2, 2.5,3.0,3.5,4])
ylim([51,79])

g.export()
