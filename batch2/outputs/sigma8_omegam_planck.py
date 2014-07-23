import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)

ranges = [0.73, 0.92, 0.23, 0.37]

for TT in [False, True]:
    g.newPlot()
    if not TT:
        basedat = s.defdata_allNoLowE
        basedatname = s.NoLowLE
        fname = 'Planckall'
    else:
        basedat = s.defdata_TTonly
        basedatname = s.planckTT
        fname = 'PlanckTT'

    roots = ['base_' + basedat,
             'base_' + basedat + '_lensing',
             'base_' + basedat + '_lensing_post_BAO'
             ]
    legends = [basedatname, '+lensing', '+BAO']
    g.plot_2d(roots, param_pair=['sigma8', 'omegam'], filled=True, lims=ranges)

    if TT:
        g.add_2d_contours('base_' + s.defdata_TT, 'sigma8', 'omegam', ls='--', color='magenta')
        g.add_2d_contours('base_' + s.defdata_TTonly + '_post_WMAPtau', 'sigma8', 'omegam', ls='-', color='brown', alpha=0.2)
    else:
        g.add_2d_contours('base_' + s.defdata_all, 'sigma8', 'omegam', ls='--', color='magenta')
        g.add_2d_contours('base_' + s.defdata_allNoLowE + '_post_WMAPtau', 'sigma8', 'omegam', ls='-', color='brown', alpha=0.2)

    g.add_text(basedatname + '+$\\tau(0.07\pm 0.02)$', 0.96, 0.12, color='magenta')
    g.add_text(basedatname + '+$\\tau(0.09\pm 0.013)$', 0.96, 0.06, color='brown', alpha=0.4)

    g.add_legend(legends, legend_loc='upper right', colored_text=True);

    g.export(fname + '_sigma_8-omega_m')

g.newPlot()
roots = ['base_' + s.defdata_TTonly,
                 'base_' + s.defdata_allNoLowE,
                 'base_' + s.defdata_allNoLowE + '_lensing',
                 'base_' + s.defdata_allNoLowE + '_lensing_post_BAO',
                 ]

g.plot_2d(roots, param_pair=['sigma8', 'omegam'], filled=True, lims=ranges)

g.add_2d_contours('base_' + s.defdata_all, 'sigma8', 'omegam', ls='--', color='magenta')
g.add_2d_contours('base_' + s.defdata_TTonly + '_post_WMAPtau', 'sigma8', 'omegam', ls='-', color='brown', alpha=0.2)

legends = [s.planckTT, s.NoLowLE, '+lensing', '+BAO']

g.add_text(s.NoLowLE + '+$\\tau(0.07\pm 0.02)$', 0.96, 0.12, color='magenta')
g.add_text(s.planckTT + '+$\\tau(0.09\pm 0.013)$', 0.96, 0.06, color='brown', alpha=0.4)
g.add_legend(legends, legend_loc='upper right', colored_text=True);
g.export()
