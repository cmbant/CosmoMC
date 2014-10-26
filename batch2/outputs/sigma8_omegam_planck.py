import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)

ranges = [0.245, 0.37, 0.73, 0.96]
pair = ['omegam', 'sigma8']

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

    g.add_legend(legends, legend_loc='upper right', colored_text=True);

    g.export(fname + '_sigma_8-omega_m')

g.newPlot()
dataroots = [s.defdata_TTonly, s.defdata_allNoLowE, s.defdata_allNoLowE + '_lensing', s.defdata_allNoLowE + '_lensing_BAO']
roots = [g.getRoot('', x) for x in dataroots]

g.plot_2d(roots, param_pair=pair, filled=True, lims=ranges)

g.add_2d_contours('base_' + s.defdata_all, param_pair=pair, ls='-', color='olive')
g.add_2d_contours(g.getRoot('', s.defdata_all + '_lensing_zre6p5'), param_pair=pair, ls='-', color='midnightblue')
# g.add_2d_contours(g.getRoot('', s.defdata_all + '_lensing_BAO'), 'sigma8', 'omegam', ls='-', color='pink')

# g.add_2d_contours('base_' + s.defdata_TTonly + '_post_WMAPtau', param_pair=pair, ls='--', color='brown', alpha=0.2)

# g.add_2d_contours('base_' + s.defdata_allNoLowE + '_lowtau', 'sigma8', 'omegam', color='red', filled=True)


legends = [s.planckTT, s.NoLowLE, r'+lensing', r'+BAO']

g.add_text(s.planckall, 0.96, 0.18, color='olive')
g.add_text(s.planckall + '+lensing', 0.96, 0.12, color='midnightblue')

# g.add_text(s.planckTT + '+$\\tau(0.09\pm 0.013)$', 0.96, 0.06, color='brown', alpha=0.4)
g.add_legend(legends, legend_loc='upper right', colored_text=True);
g.export()
