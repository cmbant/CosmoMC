import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)


g.plot_3d('base_mnu_lensonly', ['mnu', 'H0', 'sigma8'])

g.add_2d_contours('base_mnu_lensonly_HST70p6', 'mnu', 'H0', ls='-', color='blue')
g.add_2d_contours('base_mnu_lensonly_BAO', 'mnu', 'H0', ls='-', color='green')
g.add_2d_contours('base_mnu_lensonly_theta', 'mnu', 'H0', ls='-', color='red')
g.add_2d_contours('base_mnu_' + s.defdata_all, 'mnu', 'H0', ls='-', color='black')

g.add_legend(['lensonly+HST', 'lensonly+BAO', 'lensonly+theta', 'Planck'])

# xlim([0.6, 1.02])
# ylim([0.13, 0.6])

g.export()
