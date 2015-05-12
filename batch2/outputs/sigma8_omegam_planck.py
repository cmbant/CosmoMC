import planckStyle as s
from matplotlib.pyplot import *

g = s.getSinglePlotter()

ranges = [0.246, 0.37, 0.73, 0.965]
pair = ['omegam', 'sigma8']

if False:
    for TT in [False, True]:
        g.newPlot()
        if not TT:
            basedat = s.defdata_allNoLowE
            basedatname = s.NoLowLE
            allname = s.planckall
            fname = 'Planckall'
        else:
            basedat = s.defdata_TTonly
            basedatname = s.planckTT
            fname = 'PlanckTT'
            allname = s.planckTTlowTEB

        roots = [g.getRoot('', basedat),
                 g.getRoot('', basedat + '_lensing'),
                 g.getRoot('', basedat + '_lensing_BAO')]


        legends = [basedatname, '+lensing', '+BAO']
        g.plot_2d(roots, param_pair=pair, filled=True, lims=ranges)

        if TT:
            g.add_2d_contours('base_' + s.defdata_TT, param_pair=pair, ls='--', color='magenta')
    #        g.add_2d_contours('base_' + s.defdata_TTonly + '_post_WMAPtau', param_pair=pair, ls='-', color='brown', alpha=0.2)
        else:
            g.add_2d_contours('base_' + s.defdata_all, param_pair=pair, ls='--', color='magenta')
    #        g.add_2d_contours('base_' + s.defdata_allNoLowE + '_post_WMAPtau', param_pair=pair, ls='-', color='brown', alpha=0.2)

        g.add_text(allname, 0.96, 0.12, color='magenta')
    #    g.add_text(basedatname + '+$\\tau(0.09\pm 0.013)$', 0.96, 0.06, color='brown', alpha=0.4)

        g.add_legend(legends, legend_loc='upper right', colored_text=True)

        g.export(fname + '_sigma_8-omega_m')

g = s.getSinglePlotter(ratio=1)

g.settings.solid_colors[2] = ['slategray']


dataroots = [s.defdata_TTonly,
             s.defdata_allNoLowE,
             s.defdata_allNoLowE + '_lensing',
             s.defdata_allNoLowE + '_lensing_BAO',
             s.defdata_allNoLowE + '_lensing_BAO_zre6p5',
#             s.defdata_allNoLowE + '_reion'
             ]
roots = [g.getRoot('', x) for x in dataroots]

omm = np.arange(0.1, 0.7, 0.01)
s.plotBounds(omm, s.planck_lensing)


g.plot_2d(roots, param_pair=pair, filled=True, lims=ranges)

g.add_2d_contours('base_' + s.defdata_all, param_pair=pair, ls='-', color='olive')
#
g.add_2d_contours(g.getRoot('', s.defdata_all + '_lensing'), param_pair=pair, ls='-', color='darkred')
g.add_2d_contours(g.getRoot('', s.defdata_allNoLowE + '_reion'), param_pair=pair, ls='--', color='midnightblue')

legends = [s.planckTT, s.NoLowLE, r'+lensing', r'+BAO', r'+$z_{\rm re} > 6.5$']

g.add_text(s.planckall, 0.96, 0.18, color='olive')
g.add_text(s.planckall + '+lensing', 0.96, 0.12, color='darkred')
# g.add_text(s.planckall + '+lensing ($z_{\\rm re}>6.5$)', 0.96, 0.06, color='midnightblue', alpha=1)
g.add_text(s.NoLowLE + '+reion prior', 0.96, 0.06, color='midnightblue', alpha=1)

g.add_legend(legends, legend_loc='upper right', colored_text=True, align_right=True)
gca().set_yticks(np.arange(0.75, 1, 0.05))
g.export()
