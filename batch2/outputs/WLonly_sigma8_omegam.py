import planckStyle as s
from pylab import *


omm = np.arange(0.1, 0.7, 0.01)

# for base in ['base', 'base_mnu', 'base_nnu_meffsterile', 'base_nnu_mnu']:
for base in ['base']:

    g = s.getSinglePlotter()

    if base == 'base':
        s.plotBounds(omm, s.planck_lensing)

    g.plot_3d(base + '_WLonlyHeymans', ['omegam', 'sigma8', 'H0'])

    g.add_2d_contours(base + '_WLonlyHeymans_BAO', 'omegam', 'sigma8', ls='-', color='blue', lw=0.6)
    g.add_2d_contours(base + '_WLonlyHeymans_BAO_theta', 'omegam', 'sigma8', ls='-', color='green', lw=0.6)

    g.add_2d_contours(base + '_' + s.defdata, 'omegam', 'sigma8', ls='-', color='black', lw=1.2)

    g.add_legend([  'WL+' + s.BAO, r'WL+$\theta_{\rm MC}$+' + s.BAO, s.defplanck], colored_text=True)

    xlim([0.13, 0.6])
    ylim([0.45, 1.1])

    g.doExport(tag=base, adir='outputs/WLonly')
