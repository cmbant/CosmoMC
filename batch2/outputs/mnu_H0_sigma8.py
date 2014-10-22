import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)

base = 'base_mnu_'
roots = [base + s.defdata]


mnu = [0, 2.3]
H, sigma = s.H0_gpe
c = 'gray'
one = array([1, 1])
fill_between(mnu, one * (H - sigma * 2), one * (H + sigma * 2), facecolor=c, alpha=0.1, edgecolor=c, lw=0)
fill_between(mnu, one * (H - sigma), one * (H + sigma), facecolor=c, alpha=0.15, edgecolor=c, lw=0)

if False:
    H, sigma = s.H0_high
    c = 'burlywood'
    plot(mnu, one * (H - sigma), color=c, ls='--')
    plot(mnu, one * (H + sigma), color=c, ls='--')
    plot(mnu, one * (H - 2 * sigma), color=c, ls=':')
    plot(mnu, one * (H + 2 * sigma), color=c, ls=':')


g.plot_3d(roots, ['mnu', 'H0', 'sigma8'])

root = g.getRoot('mnu', s.defdata + '_lensing')
g.add_2d_contours(root, 'mnu', 'H0', plotno=0)

root = g.getRoot('mnu', s.defdata + '_BAO_HST70p6_JLA_lensing')
g.add_2d_contours(root, 'mnu', 'H0', filled=True, zorder=2, alpha=0.6)

# gca().set_xticks([2, 2.5,3.0,3.5,4])
ylim([51, 79])

g.export()
