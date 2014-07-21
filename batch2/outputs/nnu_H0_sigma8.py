import planckStyle as s
from pylab import *

g=s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)

base = 'base_nnu_'
roots=[base+s.defdata_TT]


neff=[1, 5]
H=70.6
sigma=3.3
c='gray'
one=array([1,1])
fill_between(neff, one*(H-sigma*2), one*(H+sigma*2), facecolor=c, alpha=0.1, edgecolor=c, lw=0)
fill_between(neff, one*(H-sigma), one*(H+sigma), facecolor=c, alpha=0.15, edgecolor=c, lw=0)

if False:
#too messy..
    H=72.5
    sigma=2.5
    c='burlywood'
    g.add_y_marker(H+sigma,color=c, ls='-')
    g.add_y_marker(H-sigma,color=c, ls='-')
    g.add_y_marker(H+2*sigma,color=c, ls='-', lw=0.5)
    g.add_y_marker(H-2*sigma,color=c, ls='-', lw=0.5)

    
g.plot_3d(roots, ['nnu', 'H0', 'sigma8'])
g.add_1d_marker(3.046,ls='-')
g.add_1d_marker(3.395)
g.add_1d_marker(4.046)

g.add_2d_contours(base+s.defdata_all+'_post_BAO', 'nnu', 'H0')

gca().set_xticks([2, 2.5,3.0,3.5,4])

g.export()