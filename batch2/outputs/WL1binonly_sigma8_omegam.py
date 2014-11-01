import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)


omm = np.arange(0.1, 0.7, 0.01)

s.plotBounds(omm, s.planck_lensing)

g.add_2d_contours('base_WLonly1bin', 'omegam', 'sigma8', ls='-', color='red', lw=0.6)
# g.add_2d_contours('base_lensonly_HST70p6', 'omegam', 'sigma8', ls='-', color='blue', lw=0.6)
g.add_2d_contours('base_WLonly1bin_BAO', 'omegam', 'sigma8', ls='-', color='magenta', lw=0.6)
g.add_2d_contours('base_WLonlyCons1bin', 'omegam', 'sigma8', ls='--', color='blue', lw=0.6)
g.add_2d_contours('base_WLonlyCons1bin_BAO', 'omegam', 'sigma8', ls='-', color='cyan', lw=0.6)

g.add_2d_contours('base_' + s.defdata, 'omegam', 'sigma8', ls='-', color='black', lw=1.2)

g.add_legend([ 'WL 1bin', 'WL 1bin+' + s.BAO, 'WL 1bin (Cons)' , 'WL 1 bin (Cons) + BAO', s.defplanck], colored_text=True)

xlim([0.13, 0.6])
ylim([0.31, 1.1])

g.export()
