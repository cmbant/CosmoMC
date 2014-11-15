import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)

base = 'base_nnu_'
roots = [base + s.defdata]


neff = [1, 5]
H = 70.6
sigma = 3.3
c = 'gray'
one = array([1, 1])
fill_between(neff, one * (H - sigma * 2), one * (H + sigma * 2), facecolor=c, alpha=0.1, edgecolor=c, lw=0)
fill_between(neff, one * (H - sigma), one * (H + sigma), facecolor=c, alpha=0.15, edgecolor=c, lw=0)

if False:
# too messy..
    H = 72.5
    sigma = 2.5
    c = 'burlywood'
    g.add_y_marker(H + sigma, color=c, ls='-')
    g.add_y_marker(H - sigma, color=c, ls='-')
    g.add_y_marker(H + 2 * sigma, color=c, ls='-', lw=0.5)
    g.add_y_marker(H - 2 * sigma, color=c, ls='-', lw=0.5)


g.plot_3d(roots, ['nnu', 'H0', 'sigma8'])
norm = 3.046
g.add_1d_marker(norm, ls='-')
g.add_1d_marker(norm + 0.39)
g.add_1d_marker(norm + 0.57)
g.add_1d_marker(norm + 1)

g.add_2d_contours(g.getRoot('nnu', s.defdata_all + '_BAO'), 'nnu', 'H0', color='black')
# g.add_2d_contours(g.getRoot('nnu', s.defdata_all + '_abundances'), 'nnu', 'H0', color='red', ls='--')
# g.add_2d_contours(g.getRoot('nnu', s.defdata + '_H073p9'), 'nnu', 'H0', color='magenta', ls='-')

g.add_legend([s.planckall + '+BAO', s.planckall + '+$Y_P$'], legend_loc='upper left', colored_text=True)

gca().set_xticks([2, 2.5, 3.0, 3.5, 4])

g.export()
