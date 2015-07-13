import planckStyle as s
from pylab import *

g = s.getSinglePlotter()


ranges = [0, 22, 0.76, 0.93]
pair = ['zrei', 'sigma8']

g.newPlot()
g.make_figure(1, xstretch=1.3)


dataroots = [s.defdata_TTonly, s.defdata_allNoLowE, s.defdata_allNoLowE + '_lensing', s.defdata_allNoLowE + '_lensing_BAO']
roots = [g.getRoot('', x) for x in dataroots]

g.plot_2d(roots, param_pair=pair, filled=True, lims=ranges)


legends = [s.planckTT, s.NoLowLE, r'+lensing', r'+BAO']
g.add_x_marker(6.5, ls='-')

c = 'gray'
one = array([1, 1])
fill_between([0, 6.5], one * 0.7, one * 1, facecolor=c, alpha=0.1, edgecolor=c, lw=0)

# g.add_2d_contours(g.getRoot('', s.defdata_TTonly + '_reion_BAO'), param_pair=pair, ls='-', color='red')

# g.add_text(s.planckall, 0.96, 0.18, color='olive')
# g.add_text(s.planckall + '+lensing', 0.96, 0.12, color='midnightblue')

g.add_legend(legends, legend_loc='lower right', colored_text=True)
g.export()
