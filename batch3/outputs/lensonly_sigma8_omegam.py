import planckStyle as s
import numpy as np
import pylab as plt

omm = np.arange(0.1, 0.7, 0.01)

for p in [0, 1, 2, 3, 4]:
    print(p)

    g = s.getSinglePlotter()

    np.random.seed(343535)

    s.plotBounds(omm, s.planck_lensing)

    g.plot_3d('base_lensing_lenspriors', ['omegam', 'sigma8', 'H0'])

    if p == 0:
        g.add_2d_contours('base_' + s.defdata_all, 'omegam', 'sigma8', ls='-', color='black', alpha=1)

        g.add_2d_contours('base_lensing_lenspriors_BAO', 'omegam', 'sigma8', ls='--', color='black', alpha=0.6)

        g.add_legend([s.planckall, r'$\textit{Planck}$ lensing+BAO'])
    else:
        if p > 1:
            g.add_2d_contours('base_lensing_lenspriors_BAO', 'omegam', 'sigma8', ls='--', color='black', alpha=0.6)
        if p > 2:
            g.add_2d_contours('base_' + s.defdata_all, 'omegam', 'sigma8', ls='-', color='black', alpha=1)
        if p > 3:
            g.add_2d_contours('base_' + s.defdata_all_lensing, 'omegam', 'sigma8', ls='-', color='red', alpha=1)
        if p > 1: g.add_legend([r'$\textit{Planck}$ lensing+BAO', s.planckall, s.planckall + '+lensing'][:p - 1],
                               fontsize=6, align_right=True)

    plt.xlim([0.15, 0.6])
    plt.ylim([0.6, 1.02])

    if p == 0:
        g.export()
    else:
        g.export(tag=str(p))
    g.newPlot()
