import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)


s8 = np.arange(0.5, 1, 0.01)

s.plotBounds(s8, s.planck_lensing)

g.plot_3d('base_lensonly', ['sigma8', 'omegam', 'H0'])

g.add_2d_contours('base_lensonly_HST70p6', 'sigma8', 'omegam', ls='-', color='blue', lw=0.6)
g.add_2d_contours('base_lensonly_BAO', 'sigma8', 'omegam', ls='-', color='green', lw=0.6)
g.add_2d_contours('base_lensonly_theta', 'sigma8', 'omegam', ls='-', color='red', lw=0.6)
g.add_2d_contours('base_' + s.defdata_all, 'sigma8', 'omegam', ls='-', color='black', lw=1.2)

# g.add_2d_contours('base_lensonly_HST70p6_widerns', 'sigma8', 'omegam', ls='--', color='orange')

if False:
    g.add_2d_contours('base_mnu_lensonly', 'sigma8', 'omegam', ls='--', color='black')
    g.add_2d_contours('base_mnu_lensonly_HST70p6', 'sigma8', 'omegam', ls='-', color='blue')
    g.add_2d_contours('base_mnu_lensonly_BAO', 'sigma8', 'omegam', ls='--', color='green')
    g.add_2d_contours('base_mnu_lensonly_theta', 'sigma8', 'omegam', ls='--', color='red')


g.add_legend([s.lensonly + '+' + s.HST, s.lensonly + '+' + s.BAO, s.lensonly + r'+ \theta_{MC}', s.planckall])

xlim([0.6, 1.02])
ylim([0.14, 0.6])

g.export()
