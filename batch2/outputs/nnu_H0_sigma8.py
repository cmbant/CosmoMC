import planckStyle as s
from pylab import *

g=s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)
roots=['base_nnu_'+s.defdata_TT]


neff=[1, 4.5]
H=70.6
sigma=3.3
c='gray'
one=array([1,1])
fill_between(neff, one*(H-sigma*2), one*(H+sigma*2), facecolor=c, alpha=0.1, edgecolor=c, lw=0)
fill_between(neff, one*(H-sigma), one*(H+sigma), facecolor=c, alpha=0.15, edgecolor=c, lw=0)
    
g.plot_3d(roots, ['nnu', 'H0', 'sigma8'])
g.add_1d_marker(3.046,ls='-')
g.add_1d_marker(3.395)
g.add_1d_marker(4.046)
gca().set_xticks([2, 2.5,3.0,3.5,4])

g.export()
