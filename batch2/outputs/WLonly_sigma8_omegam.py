import planckStyle as s
from pylab import *

g = s.getSinglePlotter()



omm = np.arange(0.1, 0.7, 0.01)

s.plotBounds(omm, s.planck_lensing)

for base in ['base', 'base_mnu', 'base_nnu_meffsterile', 'base_nnu_mnu']:

    g.newPlot()
    g.make_figure(1, xstretch=1.3)

    g.plot_3d(base + '_WLonly', ['omegam', 'sigma8', 'H0'])

    # g.add_2d_contours('base_lensonly_HST70p6', 'omegam', 'sigma8', ls='-', color='blue', lw=0.6)
    g.add_2d_contours(base + '_WLonly_BAO', 'omegam', 'sigma8', ls='-', color='darkred', lw=0.6)
    g.add_2d_contours(base + '_WLonly_BAO_theta', 'omegam', 'sigma8', ls='-', color='green', lw=0.6)

    g.add_2d_contours(base + '_WLonlyHeymans', 'omegam', 'sigma8', ls='--', color='blue', lw=0.6)
    g.add_2d_contours(base + '_WLonlyHeymans_BAO', 'omegam', 'sigma8', ls='--', color='darkred', lw=0.6)
    g.add_2d_contours(base + '_WLonlyHeymans_BAO_theta', 'omegam', 'sigma8', ls='--', color='green', lw=0.6)

    g.add_2d_contours(base + '_' + s.defdata, 'omegam', 'sigma8', ls='-', color='black', lw=1.2)

    # g.add_2d_contours('base_lensonly_HST70p6_widerns', 'omegam','sigma8', ls='--', color='orange')


    g.add_legend([  'WL +' + s.BAO, r'WL + \theta+' + s.BAO, 'WL (Heymans)' , 'WL (Heymans) + BAO', r'WL (Heymans) + $\theta$ +BAO', s.defplanck], colored_text=True)

    xlim([0.13, 0.6])
    ylim([0.31, 1.1])

    g.doExport(None, tag=base, adir='outputs/WLonly')
