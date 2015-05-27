import planckStyle as s
import numpy as np
import matplotlib.pyplot as plt

def yhe_to_ypBBN(Yp):
    m_proton = 1.672621637e-27
    m_H = 1.673575e-27
    not4 = 3.9715
    m_He = m_H * not4
    return -4 * m_H * Yp / (Yp * m_He - 4 * Yp * m_H - m_He)


def addYpFit(samples):
    p = samples.getParams()
    yBBN = yhe_to_ypBBN(p.yhe)
    YpBBN = samples.addDerived(yBBN, name='YpFit', label=r'Y_P^{\rm BBN}')
    samples.updateChainBaseStatistics()
    return YpBBN


def addYpFitRoots(roots):
    for root in roots:
        samples = g.sampleAnalyser.samplesForRoot(root)
        Yp = addYpFit(samples)
        samples.getMargeStats()
        print root, Yp.mean, Yp.limits[1]


# Parthenelope fits, from Julien 3 Dec 2014
def ypBBN(omegab, dneff, taun=880.3):
    return (
               0.2311 + 0.9502 * omegab - 11.27 * omegab * omegab
               + dneff * (0.01356 + 0.008581 * omegab - 0.1810 * omegab * omegab)
               + dneff * dneff * (-0.0009795 - 0.001370 * omegab + 0.01746 * omegab * omegab)
           ) * pow(taun / 880.3, 0.728)


def dh_fit(omegab, dneff, taun=880.3):
    return (
               18.754 - 1534.4 * omegab + 48656. * omegab * omegab - 552670. * omegab * omegab * omegab
               + dneff * (2.4914 - 208.11 * omegab + 6760.9 * omegab * omegab - 78007. * omegab * omegab * omegab)
               + dneff * dneff * (
                   0.012907 - 1.3653 * omegab + 37.388 * omegab * omegab - 267.78 * omegab * omegab * omegab)
           ) * pow(taun / 880.3, 0.418) * 1e-5


# 68% range for omega_b from cmb (TT+lowP+BAO)
# omegab_minus = 0.02228 - 0.00020 #0.02178
# omegab_plus  = 0.02228 + 0.00020 #0.02232

# BBN theoretical error
sigma_yp_theo = 0.0003
# resolution of the theoretical BBN curve (number of omega_b values)
num_ob = 50
# omegab range in the plot
ob_min = 0.018
ob_max = 0.026
# yhe range in the plot
yp_min = 0.15
yp_max = 0.35
# helium data: Aver et al. 2013, arXiv:1309.0047
aver_mean = 0.2465
aver_sigma = 0.0097
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

labels = [s.defplanck, s.defplanck + "+BAO", s.planckall]
datatag = [s.defdata, s.defdata + '_BAO', s.defdata_all]

########### ombh2 -Yhe #############

g = s.getSinglePlotter()

colors = g.settings.solid_colors[3:0:-1]
del colors[1]

bbn_b = np.arange(ob_min, ob_max + 0.1, (ob_max - ob_min) / num_ob)
bbn_y = ypBBN(bbn_b, 0)
bbn_y1 = bbn_y - sigma_yp_theo
bbn_y2 = bbn_y + sigma_yp_theo

g.add_y_bands(aver_mean, aver_sigma, xlim=sere_b)

plt.fill_between(sere_b, sere_y1, yp_max, alpha=0.07, color='gray')
plt.plot(sere_b, sere_y1, alpha=0.2, color='gray', linestyle='-')
plt.text(0.0183, 0.249, "Aver et al. (2013)", fontsize=7.)
plt.text(0.0183, 0.325, "Excluded by Serenelli \& Basu (2010)", fontsize=7.)

plt.fill_between(bbn_b, bbn_y1, bbn_y2, alpha=0.9, color='green')
plt.plot(bbn_b, bbn_y1, color='green', linestyle='solid')
plt.plot(bbn_b, bbn_y2, color='green', linestyle='solid')

roots = [g.getRoot('yhe', d) for d in datatag]

g.settings.legend_fontsize = 8
g.plot_2d(roots, 'omegabh2', 'YpBBN', filled=True, lims=[ob_min + 0.0001, ob_max, yp_min, yp_max], colors=colors)
g.add_legend(labels, legend_loc='lower left', colored_text=False)
plt.gca().set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])

plt.gca().annotate('Standard BBN',
                   xy=(0.025, 0.25),
                   xycoords='data',
                   xytext=(-60, -30),
                   textcoords='offset points',
                   arrowprops=dict(arrowstyle="->",
                                   connectionstyle="arc3,rad=.2"),
                   fontsize=9.
                   )
g.export()


########### Neff -Yhe #############
g = s.getSinglePlotter()

N_min = 0.01
N_max = 5
Neff = np.arange(N_min, N_max + 0.1, 0.1)

Nrange = [N_min, N_max]

g.add_y_bands(aver_mean, aver_sigma, xlim=Nrange)

plt.fill_between(Nrange, Neff[-1], sere_y1, alpha=0.07, color='gray')
plt.plot(Nrange, sere_y1, alpha=0.2, color='gray', linestyle='-')
plt.text(0.17, 0.245, "Aver et al. (2013)", fontsize=7.)
plt.text(0.17, 0.337, "Excluded by Serenelli \& Basu (2010)", fontsize=7.)

roots = [g.getRoot('nnu_yhe', d) for d in datatag]

g.plot_2d(roots, 'nnu', 'YpBBN', filled=True, lims=[0, N_max, yp_min, yp_max], colors=colors)

bbn_y = ypBBN(0.02222, Neff - 3.046)
bbn_y1 = bbn_y - sigma_yp_theo
bbn_y2 = bbn_y + sigma_yp_theo

plt.fill_between(Neff, bbn_y1, bbn_y2, alpha=0.9, color='green')
plt.plot(Neff, bbn_y1, color='green', linestyle='solid')
plt.plot(Neff, bbn_y2, color='green', linestyle='solid')
g.add_legend(labels, legend_loc='lower left', colored_text=False, fontsize=8)
g.add_x_marker(3.046)
plt.gca().set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
# g.rotate_yticklabels()

plt.gca().annotate('Standard BBN',
                   xy=(4.5, 0.262),
                   xycoords='data',
                   xytext=(-40, 30),
                   textcoords='offset points',
                   arrowprops=dict(arrowstyle="->",
                                   connectionstyle="arc3,rad=.2"),
                   fontsize=8.
                   )

g.export(tag='neff')
