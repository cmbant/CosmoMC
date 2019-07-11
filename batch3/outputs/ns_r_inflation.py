# Uses first-order slow-roll results, with naive N definition, as in inflation paper

import planckStyle as s
from pylab import *
import numpy as np

BK = 'BK15'


def N_r_ns(r, ns):
    return (r - 16) / (8. * ns - 8 + r) / 2.


def r_ns(ns, p):
    return 8. * (1 - ns) * p / (2. + p)


def r_ns_2(ns, p):
    return -4. / 3 * (
            32 * ns - 12 * ns * p - 28 - 6 * p - 3 * p ** 2 + 3 * p ** 2 * ns ** 2 + 18 * p * ns ** 2 - 4 * ns ** 2) * p / (
                   (2. + p) ** 3)


def ns_N(N, p):  # first order
    return (4 * N - p - 4) / (4. * N + p)


def ns_N_2(N, p):  # solved used above expresion for r_ns, but naive expression for N
    return 0.5 * (-12 * p ** 2 + 96 * N + 80 * N * p + 12 * N * p ** 3 + 32 * p + 96 * N * p ** 2 - 2 * sqrt(
        384. * p - 1008 * p ** 2 - 2304 * p ** 3 + 384 * N * p + 4608 * N * p ** 2 + 6912 * N * p ** 3 + 3936 * p ** 4 * N + 936 * p ** 5 * N + 2304 * N ** 2 + 4608 * N ** 2 * p + 3456 * N ** 2 * p ** 2 + 1152 * N ** 2 * p ** 3 + 144 * N ** 2 * p ** 4 - 1464 * p ** 4 - 360 * p ** 5 + 72 * N * p ** 6 - 27 * p ** 6)) / (
                   6. * N * p ** 3 + 36 * N * p ** 2 - 8 * N * p - 3 * p ** 3 - 18 * p ** 2 + 4 * p)


# return (-12.*p ** 2 + 96 * N + 80 * N * p + 12 * N * p ** 3 + 32 * p + 96 * N * p ** 2 + 2 * sqrt(384. * p - 1008 * p ** 2 - 2304 * p ** 3 + 384 * N * p + 4608 * N * p ** 2 + 6912 * N * p ** 3 + 3936 * p ** 4 * N + 936 * p ** 5 * N + 2304 * N ** 2 + 4608 * N ** 2 * p + 3456 * N ** 2 * p ** 2 + 1152 * N ** 2 * p ** 3 + 144 * N ** 2 * p ** 4 - 1464 * p ** 4 - 360 * p ** 5 + 72 * N * p ** 6 - 27 * p ** 6)) / (6. * N * p ** 3 + 36 * N * p ** 2 - 8 * N * p - 3 * p ** 3 - 18 * p ** 2 + 4 * p)


g = s.getSinglePlotter(ratio=1)

roots = ['base_r_' + s.defdata_TTTEEE,
         'base_r_' + s.defdata_TTTEEE + '_lensing',
         'base_r_plikHM_TTTEEE_lowl_lowE_%s_lensing_post_BAO' % BK,
         ]
colors = [g.settings.solid_colors[1], g.settings.solid_colors[3], g.settings.solid_colors[0]]

g.plot_2d(roots, ['ns', 'r02'], filled=True, colors=colors)

g.add_2d_contours('base_r_CamSpecHM_TTTEEE_lowl_lowE_%s_post_BAO_lensing' % BK, 'ns', 'r02', ls='--', color='gray')
# g.add_2d_contours('base_r_plikHM_TT_lowl_lowE_BK15_post_BAO_lensing', 'ns', 'r02', ls='-.', color='darkcyan')

ns = np.arange(0.91, 1.02, 0.0005)
r = np.arange(0, 0.34, 0.002)
ns, r = np.meshgrid(ns, r)

# this is the first order result
N = N_r_ns(r, ns)
N[r > 8 * (1 - ns)] = 100

CS = contour(ns, r, N, origin='lower', levels=[50, 60], colors='k', linestyles=':', linewidths=0.3,
             extent=[0.95, 1.01, 0, 0.25])

for x, y, lab in zip([0.954, 0.9585], [0.2, 0.211], ['N=50', 'N=60']):
    plt.text(x, y, lab, size=7, rotation=-62, color='k',
             ha="center", va="center", bbox=dict(ec='1', fc='1', alpha=0))

ns = arange(0.9, 1.1, 0.001)  #
plot(ns, r_ns(ns, 1), ls='-', color='k', lw=1, alpha=0.8)

plt.text(0.954, 0.13, 'Convex', size=7, rotation=-30, color='k',
         ha="center", va="center")
plt.text(0.954, 0.116, 'Concave', size=7, rotation=-30, color='k',
         ha="center", va="center")

modcol = 'red'

text(0.962, 0.155, r'$\phi^2$', fontsize=9, color=modcol, bbox=dict(ec='1', fc='1', alpha=0.8), zorder=-2)

text(0.971, 0.081, r'$\phi$', fontsize=9, color=modcol)

for p in [1, 2]:
    ns = arange(0.93, 1.1, 0.002)  #

    if p != 1: plot(ns, r_ns(ns, p), ls='-', color='black', lw=0.6, alpha=0.4)

    #    print p, ns_N(50, p), ns_N(60, p)
    ns = arange(ns_N(50, p), ns_N(60, p), 0.0005)  #
    plot(ns, r_ns(ns, p), ls='-', color=modcol, lw=1.2, alpha=1)

labels = [s.planckall, s.planckall + '+lensing', '+%s+BAO' % BK]

leg = g.add_legend(labels, colored_text=True, align_right=True)

xlim([0.945, 1])
ylim([0, 0.26])

gca().set_xticks([0.95, 0.96, 0.97, 0.98, 0.99, 1.0])
gca().set_yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25])

g.export('ns_r_inflation', tag=BK)
