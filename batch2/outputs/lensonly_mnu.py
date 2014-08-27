import planckStyle as s
from pylab import *


g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)


s8 = np.arange(0.5, 1, 0.01)

s.plotBounds(s8, s.planck_lensing)

base = 'nnu_meffsterile'

g.plot_3d(g.getRoot(base, 'lensonly'), ['sigma8', 'omegam', 'H0'])

# g.add_2d_contours('base_mnu_lensonly_HST70p6', 'sigma8', 'omegam', ls='-', color='blue', lw=0.6)
# g.add_2d_contours(base + 'lensonly_BAO', 'sigma8', 'omegam', ls='-', color='green', lw=0.6)
g.add_2d_contours(g.getRoot(base, s.defdata_all), 'sigma8', 'omegam', ls='--', color='black', lw=1.2)
# g.add_2d_contours(g.getRoot(base, s.defdata_all + '_BAO'), 'sigma8', 'omegam', ls='--', color='red', lw=1.2)

# g.add_legend([s.lensonly + '+' + s.BAO, s.planckall, s.planckall + '+' + s.BAO], colored_text=True)

# g.add_legend([s.lensonly + '+' + s.BAO, s.planckall, s.planckall + '+' + s.BAO], colored_text=True)

xlim([0.6, 1.02])
ylim([0.14, 0.6])

g.export()
