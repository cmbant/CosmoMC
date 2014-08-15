import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)


g.plot_3d('base_lensonly', ['sigma8', 'omegam', 'H0'])

g.add_2d_contours('base_lensonly_HST70p6', 'sigma8', 'omegam', ls='-', color='blue')
g.add_2d_contours('base_lensonly_BAO', 'sigma8', 'omegam', ls='-', color='green')
g.add_2d_contours('base_lensonly_theta', 'sigma8', 'omegam', ls='-', color='red')
g.add_2d_contours('base_' + s.defdata_all, 'sigma8', 'omegam', ls='-', color='black')

g.add_legend(['lensonly+HST', 'lensonly+BAO', 'lensonly+theta', 'Planck'])

xlim([0.6, 1.02])
ylim([0.13, 0.6])

g.export()
