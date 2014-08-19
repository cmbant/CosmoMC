import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)
base = 'base_nnu_meffsterile_'
roots = [base + s.defdata_TT]

s8 = np.arange(0.5, 1, 0.01)


# plotBounds(galaxygalaxy,'green')
s.plotBounds(s8, s.CFTHlens, 'burlywood')
s.plotBounds(s8, s.PLSZ, 'gray')


g.plot_3d(roots, ['sigma8', 'omegam', 'H0'])

g.add_2d_contours(base + s.defdata_all + '_post_BAO', 'sigma8', 'omegam')
# g.add_2d_contours(base+s.defdata_all+'_post_HST70p6', 'sigma8', 'omegam', plotno=1)

g.export()
