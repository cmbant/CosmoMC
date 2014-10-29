import planckStyle as s
from pylab import *


g = s.getSinglePlotter()

for tag, pars in zip(['sigma8_omegam', 'omegam_H0'], [['omegam', 'sigma8', 'H0'], ['omegam', 'H0', 'sigma8']]):
    for base in ['', 'mnu', 'nnu', 'nnu_meffsterile']:
        g.newPlot()
        g.make_figure(1, xstretch=1.3)

        if pars[1] == 'sigma8':
            omm = np.arange(0.1, 0.7, 0.01)
            s.plotBounds(omm, s.planck_lensing)

        g.plot_3d(g.getRoot(base, 'lensonly'), pars)

        # g.add_2d_contours('base_mnu_lensonly_HST70p6', 'sigma8', 'omegam', ls='-', color='blue', lw=0.6)
        # g.add_2d_contours(base + 'lensonly_BAO', 'sigma8', 'omegam', ls='-', color='green', lw=0.6)
        g.add_2d_contours(g.getRoot(base, s.defdata), pars[0], pars[1], ls='-', color='black', lw=1.2)
        # g.add_2d_contours(g.getRoot(base, s.defdata_all + '_BAO'), 'sigma8', 'omegam', ls='--', color='red', lw=1.2)

        # g.add_legend([s.lensonly + '+' + s.BAO, s.planckall, s.planckall + '+' + s.BAO], colored_text=True)

        # g.add_legend([s.lensonly + '+' + s.BAO, s.planckall, s.planckall + '+' + s.BAO], colored_text=True)
        if pars[1] == 'sigma8':
            pass
            ylim([0.6, 1.02])
            xlim([0.14, 0.6])
        if base == '': base = 'base'
        g.export('lensonly_' + tag + '_' + base)
