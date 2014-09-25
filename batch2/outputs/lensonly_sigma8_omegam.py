import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)


omm = np.arange(0.1, 0.7, 0.01)

s.plotBounds(omm, s.planck_lensing)

g.plot_3d('base_lensonly', ['omegam', 'sigma8', 'H0'])

g.add_2d_contours('base_lensonly_HST70p6', 'omegam', 'sigma8', ls='-', color='blue', lw=0.6)
g.add_2d_contours('base_lensonly_BAO', 'omegam', 'sigma8', ls='-', color='green', lw=0.6)
g.add_2d_contours('base_lensonly_theta', 'omegam', 'sigma8', ls='-', color='red', lw=0.6)
g.add_2d_contours('base_' + s.defdata_all, 'omegam', 'sigma8', ls='-', color='black', lw=1.2)

# g.add_2d_contours('base_lensonly_HST70p6_widerns', 'omegam','sigma8', ls='--', color='orange')

if False:
    g.add_2d_contours('base_mnu_lensonly', 'omegam', 'sigma8', ls='--', color='black')
    g.add_2d_contours('base_mnu_lensonly_HST70p6', 'omegam', 'sigma8', ls='-', color='blue')
    g.add_2d_contours('base_mnu_lensonly_BAO', 'omegam', 'sigma8', ls='--', color='green')
    g.add_2d_contours('base_mnu_lensonly_theta', 'omegam', 'sigma8', ls='--', color='red')


g.add_legend([s.lensonly + '+' + s.HST, s.lensonly + '+' + s.BAO, s.lensonly + r'+ $\theta_{MC}$', s.planckall])

xlim([0.14, 0.6])
ylim([0.6, 1.02])

g.export()
