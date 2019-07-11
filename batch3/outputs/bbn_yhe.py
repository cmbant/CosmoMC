import planckStyle as s
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, r'c:\work\dist\git\camb')
from camb.bbn import BBN_table_interpolator

BBNstandard = BBN_table_interpolator('PArthENoPE_880.2_standard.dat')
# BBN theoretical error
sigma_yp_theo = 0.0003
# resolution of the theoretical BBN curve (number of omega_b values)
num_ob = 50
# omegab range in the plot
ob_min = 0.019
ob_max = 0.025
# yhe range in the plot
yp_min = 0.175
yp_max = 0.28
# helium data: Aver et al. 2015
aver_mean = 0.2449
aver_sigma = 0.004
# helium data: Serenelli and Basu 2010
sere_minus = 0.294
sere_plus = yp_max

sere_b = np.zeros(2, dtype='float')
sere_y1 = np.zeros(2, dtype='float')
sere_y2 = np.zeros(2, dtype='float')
sere_b[0] = ob_min
sere_b[1] = ob_max
sere_y1[0] = sere_minus
sere_y1[1] = sere_minus
sere_y2[0] = sere_plus
sere_y2[1] = sere_plus

labels = [s.planckall, s.planckall + "+lensing+BAO"]
datatag = [s.defdata_all, s.defdata_all + '_lensing_BAO']

########### ombh2 -Yhe #############

g = s.getSinglePlotter()

colors = g.settings.solid_colors[3:0:-1]
del colors[1]

bbn_b = np.arange(ob_min, ob_max + 0.1, (ob_max - ob_min) / num_ob)
bbn_y = np.array([BBNstandard.Y_p(x, 0) for x in bbn_b])
bbn_y1 = bbn_y - sigma_yp_theo
bbn_y2 = bbn_y + sigma_yp_theo

g.add_y_bands(aver_mean, aver_sigma)

# plt.fill_between(sere_b, sere_y1, yp_max, alpha=0.07, color='gray')
# plt.plot(sere_b, sere_y1, alpha=0.2, color='gray', linestyle='-')
plt.text(0.0193, 0.249, "Aver et al. (2015)", fontsize=7.)
# plt.text(0.0183, 0.325, "Excluded by Serenelli \& Basu (2010)", fontsize=6.5)

bbn_y1 = bbn_y - 2 * sigma_yp_theo
bbn_y2 = bbn_y + 2 * sigma_yp_theo
plt.fill_between(bbn_b, bbn_y1, bbn_y2, alpha=0.4, color='green', lw=0, zorder=10)
bbn_y1 = bbn_y - sigma_yp_theo
bbn_y2 = bbn_y + sigma_yp_theo
plt.fill_between(bbn_b, bbn_y1, bbn_y2, alpha=0.9, color='green', lw=0, zorder=11)

# plt.plot(bbn_b, bbn_y1, color='green', linestyle='solid')
# plt.plot(bbn_b, bbn_y2, color='green', linestyle='solid')


roots = [g.getRoot('yhe', d) for d in datatag]

g.settings.legend_fontsize = 8
g.plot_2d(roots, 'omegabh2', 'YpBBN', filled=True, lims=[ob_min + 0.0001, ob_max, yp_min, yp_max])
g.add_legend(labels, legend_loc='lower left', colored_text=False)
# plt.gca().set_yticks([0.2, 0.25, 0.3])

plt.gca().annotate('Standard BBN',
                   xy=(0.0242, 0.249),
                   xycoords='data',
                   xytext=(-35, -30),
                   textcoords='offset points',
                   arrowprops=dict(arrowstyle="->",
                                   connectionstyle="arc3,rad=.2"),
                   fontsize=8.
                   )
g.export()

########### Neff -Yhe #############
g = s.getSinglePlotter()

N_min = 0.01
N_max = 5
Neff = np.arange(N_min, N_max + 0.1, 0.1)

Nrange = [N_min, N_max]

g.add_y_bands(aver_mean, aver_sigma)

plt.fill_between(Nrange, Neff[-1], sere_y1, alpha=0.07, color='gray')
plt.plot(Nrange, sere_y1, alpha=0.2, color='gray', linestyle='-')
plt.text(0.17, 0.242, "Aver et al. (2015)", fontsize=6)
plt.text(0.17, 0.337, "Excluded by Serenelli \& Basu (2010)", fontsize=6)

roots = [g.getRoot('nnu_yhe', d) for d in datatag]
# roots += ['base_nnu_yhe_' + s.defdata_all + '_Aver15']

g.plot_2d(roots, 'nnu', 'YpBBN', filled=True, lims=[0, N_max, yp_min, yp_max])

g.add_2d_contours('base_nnu_yhe_' + s.defdata_all + '_Aver15_post_BAO_lensing', 'nnu', 'YpBBN', filled=False)

ombh2mean = 0.0224
bbn_y = np.array([BBNstandard.Y_p(ombh2mean, x - 3.046) for x in Neff])
bbn_y1 = bbn_y - 2 * sigma_yp_theo
bbn_y2 = bbn_y + 2 * sigma_yp_theo
plt.fill_between(Neff, bbn_y1, bbn_y2, alpha=0.4, color='green', lw=0)
bbn_y1 = bbn_y - sigma_yp_theo
bbn_y2 = bbn_y + sigma_yp_theo
plt.fill_between(Neff, bbn_y1, bbn_y2, alpha=0.9, color='green', lw=0)

# plt.plot(Neff, bbn_y1, color='green', linestyle='solid')
# plt.plot(Neff, bbn_y2, color='green', linestyle='solid')

labels = labels[:1] + ['+lensing+BAO']

g.add_legend(labels, legend_loc='lower left', colored_text=True, fontsize=8)
g.add_x_marker(3.046)
plt.gca().set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
# g.rotate_yticklabels()

plt.gca().annotate('Standard BBN\n' + r'($\Omega_{\rm b} h^2=0.0224$)',
                   xy=(4.5, 0.262),
                   xycoords='data',
                   xytext=(-44, 30),
                   textcoords='offset points',
                   arrowprops=dict(arrowstyle="->",
                                   connectionstyle="arc3,rad=.2"),
                   fontsize=6.5
                   )

g.export(tag='neff')
