import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)

roots = ['base_mnu_lensonly_theta_BAO', 'base_mnu_' + s.defdata_all, 'base_mnu_' + s.defdata_all + '_lensing']

g.plot_3d(roots[0], ['mnu', 'H0', 'sigma8'])
g.add_2d_contours(roots[1], 'mnu', 'H0', ls='-', color='red')
g.add_2d_contours(roots[2], 'mnu', 'H0', ls='-', color='black')


if False:

    g.plot_3d('base_mnu_lensonly', ['mnu', 'H0', 'sigma8'])

    g.add_2d_contours('base_mnu_lensonly_HST70p6', 'mnu', 'H0', ls='-', color='blue')
    g.add_2d_contours('base_mnu_lensonly_BAO', 'mnu', 'H0', ls='-', color='green')
    g.add_2d_contours('base_mnu_lensonly_theta', 'mnu', 'H0', ls='-', color='red')

    g.add_2d_contours('base_mnu_' + s.defdata_all, 'mnu', 'H0', ls='-', color='black')

    g.add_legend(['lensonly+HST', 'lensonly+BAO', 'lensonly+theta', 'Planck'])

g.add_legend(['lensonly+theta+BAO', s.planckall, s.lensing])

# xlim([0.6, 1.02])
# ylim([0.13, 0.6])

g.export()
