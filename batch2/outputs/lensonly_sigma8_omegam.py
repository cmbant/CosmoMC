import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)


s8 = np.arange(0.5, 1, 0.01)

def lensing(sigma):  # mandelbaum
    return  ((0.572 + 0.019 * sigma) / s8) ** (1 / 0.25)

def plotBounds(data, c):
    fill_between(s8, data(-2), data(2), facecolor=c, alpha=0.15, edgecolor=c, lw=0)
    fill_between(s8, data(-1), data(1), facecolor=c, alpha=0.25, edgecolor=c, lw=0)

plotBounds(lensing, 'grey')

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



g.add_legend(['lensonly+HST', 'lensonly+BAO', 'lensonly+theta', 'Planck'])

xlim([0.6, 1.02])
ylim([0.14, 0.6])

g.export()
