import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

omm = np.arange(0.1, 0.7, 0.01)

np.random.seed(343535)

s.plotBounds(omm, s.planck_lensing)

g.plot_3d('base_lensonly', ['omegam', 'sigma8', 'H0'])

g.add_2d_contours('base_lensonly_BAO', 'omegam', 'sigma8', ls='-', color='blue', lw=0.6, alpha=0.4)
g.add_2d_contours('base_lensonly_BAO_theta', 'omegam', 'sigma8', ls='-', color='red', lw=0.6, alpha=0.4)
g.add_2d_contours('base_' + s.defdata, 'omegam', 'sigma8', ls='-', color='black', lw=1.2, alpha=0.5)
# g.add_2d_contours('base_' + s.defdata + '_lensing', 'omegam', 'sigma8', ls='-', color='red')


if False:
    g.add_2d_contours('base_mnu_lensonly', 'omegam', 'sigma8', ls='--', color='black')
    # g.add_2d_contours('base_mnu_lensonly_HST70p6', 'omegam', 'sigma8', ls='-', color='blue')
    g.add_2d_contours('base_mnu_lensonly_BAO', 'omegam', 'sigma8', ls='--', color='green')
    g.add_2d_contours('base_mnu_lensonly_theta', 'omegam', 'sigma8', ls='--', color='red')

# g.add_legend([s.defplanck, s.defplanck + '+lensing'])

g.add_legend([s.lensonly + '+' + s.BAO, s.lensonly + '+' + s.BAO + r'+$\theta_{\rm MC}$', s.defplanck])

xlim([0.15, 0.6])
ylim([0.6, 1.02])

g.export()
