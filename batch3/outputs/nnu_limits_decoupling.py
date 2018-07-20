import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
import pylab as plt

# Use MeV units
MeV = 1
GeV = 1000
m_e = 0.511 * MeV
m_mu = 105.658 * MeV
m_t = 172.44 * GeV
m_b = 4.18 * GeV
m_H = 126 * GeV
m_W = 80.5 * GeV
m_Z = 91.18 * GeV
Kelvin_to_MeV = 8.61733e-11 * MeV
TCMB = 2.7255 * Kelvin_to_MeV

# Standard d.o.f.
g_quarks = 12 * 6 * 7. / 8
g_gluons = 16
g_leptons = (3 * 4 + 3 * 2) * 7. / 8
g_vector_bosons = 3 * 3
g_photon = 2
g_Higgs = 1
g_SM = g_quarks + g_gluons + g_leptons + g_vector_bosons + g_photon + g_Higgs

# Data from Borsanyi et al 2016 (T in MeV) arXiv:1606.07494
# Note I added fake data point at end to try to keep intepolation spline fairly flat there
# Asymptotic behavioue for high T not known (eventually converge to expected g_SM, but only at much higher temperature?)

log10T_nodes = np.array([0.00, 0.50, 1.00, 1.25, 1.60, 2.00, 2.15, 2.20,
                         2.40, 2.50, 3.00, 4.00, 4.30, 4.60, 5.00, 5.45, 5.5])  # MeV up to 281 GeV
geff_nodes = np.array([10.71, 10.74, 10.76, 11.09, 13.68, 17.61, 24.07, 29.84,
                       47.83, 53.04, 73.48, 83.10, 85.56, 91.97, 102.17, 104.98, 105])
geff_gs_nodes = np.array([1.00228, 1.00029, 1.00048, 1.00505, 1.02159, 1.02324, 1.05423, 1.07578,
                          1.06118, 1.04690, 1.01778, 1.00123, 1.00389, 1.00887, 1.00750, 1.00023, 1])

gs_nodes = geff_nodes / geff_gs_nodes

_mx = 100


def fermi_dirac_rho(p, m, T):
    E = np.sqrt(m ** 2 + p ** 2)
    if E / T > _mx:
        return 0
    else:
        return p ** 2 * E / (np.exp(E / T) + 1) / (2 * np.pi ** 2)


def bose_einstein_rho(p, m, T):
    E = np.sqrt(m ** 2 + p ** 2)
    if E / T > _mx:
        return 0
    else:
        return p ** 2 * E / (np.exp(E / T) - 1) / (2 * np.pi ** 2)


def fermi_dirac_P(p, m, T):
    E = np.sqrt(m ** 2 + p ** 2)
    if E / T > _mx:
        return 0
    else:
        return p ** 4 / E / (np.exp(E / T) + 1) / (2 * np.pi ** 2) / 3


def bose_einstein_P(p, m, T):
    E = np.sqrt(m ** 2 + p ** 2)
    if E / T > _mx:
        return 0
    else:
        return p ** 4 / E / (np.exp(E / T) - 1) / (2 * np.pi ** 2) / 3


def spline(x, y):
    return UnivariateSpline(x, y, s=0, ext=3)


# Load numerical results with errors from 1803.01038 appendix C
# URL is http://member.ipmu.jp/satoshi.shirai/standardmodel2018.dat
import requests
import tempfile
import os

filename = os.path.join(tempfile.gettempdir(), 'standardmodel2018.dat')
if not os.path.exists(filename):
    r = requests.get('http://member.ipmu.jp/satoshi.shirai/standardmodel2018.dat', allow_redirects=True)
    open(filename, 'w').write(r.content)
table = np.loadtxt(filename)
table_gstar = spline(table[:, 0], table[:, 1])
table_gstar_err = spline(table[:, 0], table[:, 2])

# Do some standard thermal history through electron-positron annhilation
Tgamma = np.logspace(-4, 2, 1000)
T_init = Tgamma[-1]
T_nu0 = (4. / 11) ** (1. / 3) * TCMB

a_init = T_nu0 / T_init  # start at same temperature
rho = np.zeros(Tgamma.size)
P = np.zeros(Tgamma.size)

for i, T in enumerate(Tgamma):
    rho[i], err = quad(fermi_dirac_rho, 0, np.inf, args=(m_e, T))
    P[i], err = quad(fermi_dirac_P, 0, np.inf, args=(m_e, T))

massless_boson_rho = Tgamma ** 4 * np.pi ** 2 / 30

# electrons, positrons and photon
s = 4 * (rho + P) / Tgamma + 2 * 4. / 3 * massless_boson_rho / Tgamma
rhotot = 4 * rho + 2 * massless_boson_rho

a = a_init * (s[-1] / s) ** (1. / 3)
z = 1 / a - 1
Tnu = T_nu0 / a

a_Tgamma = spline(Tgamma, a)
Tgamma_a = spline(a[::-1], Tgamma[::-1])
Tnu_a = spline(a[::-1], Tnu[::-1])

# maybe should use 3.045, using this for consistency
Neff_standard = 3.046

# splines after neutrino decoupling, slightly fudging 0.046 factor to give right result at late times
gstar_s = spline(Tgamma, s / Tgamma ** 3 * 45. / 2 / np.pi ** 2 + 7. / 8 * Neff_standard * 2 * (Tnu / Tgamma) ** 3)
gstar = spline(Tgamma, rhotot / Tgamma ** 4 * 30 / np.pi ** 2 + 7. / 8 * Neff_standard * 2 * (Tnu / Tgamma) ** 4)
gstar_s_nonu = spline(Tgamma, s / Tgamma ** 3 * 45. / 2 / np.pi ** 2)
gstar_P = spline(Tgamma, (s * Tgamma - rhotot) / Tgamma ** 4 * 30 / np.pi ** 2 + 7. / 8 * Neff_standard * 2 * (
        Tnu / Tgamma) ** 4 / 3)

# Where to match spline to the standard thermal at <1MeV
logTmatch = -0.1
Tcut = np.unique(Tgamma[np.log10(Tgamma) < logTmatch])
logTcut = np.log10(Tcut)
print 'Match at T=%s, g=%s for match to g=%s' % (10. ** logTmatch, gstar_s(Tcut)[-1], gs_nodes[0])
print 'max GeV', 10. ** log10T_nodes[-2] / 1000

# join results
full_T = np.concatenate((np.log10(Tcut), log10T_nodes))
maxT = 10. ** (log10T_nodes[-2])
gs_full = spline(full_T, np.concatenate((gstar_s(Tcut), gs_nodes)))
gs_full_therm = spline(full_T, np.concatenate((gstar_s_nonu(Tcut), gs_nodes)))

geff_full = spline(full_T, np.concatenate((gstar(Tcut), geff_nodes)))
w_full = spline(full_T, np.concatenate((gstar_P(Tcut) / gstar(Tcut), 4. / 3 * gs_nodes / geff_nodes - 1)))

log10a_full = - full_T + 1. / 3 * np.log10(gs_full(np.log10(TCMB)) / gs_full(full_T)) + np.log10(TCMB)
# a*T/TCMB as function of log10(a)
aT_a_full = spline(log10a_full[::-1], (gs_full(full_T[::-1]) / gs_full(np.log10(TCMB))) ** (-1. / 3))


def aTnu_a_full(loga):
    Tgam = aT_a_full(loga)
    Tgam[Tgam > T_nu0 / TCMB] = T_nu0 / TCMB
    return Tgam


# gstar at zero of things that were in equilibrium at given time
def gs_full_therm0(logTd):
    res = np.empty(logTd.shape)
    res[logTd <= logTcut[-1]] = 2
    res[logTd > logTcut[-1]] = gstar_s(TCMB)
    return res


gextra = 1
Textra_0 = (gs_full_therm0(full_T) / gs_full_therm(full_T)) ** (1. / 3) * TCMB
Neff = spline(full_T, gextra * Textra_0 ** 4 / (2 * 7 / 8. * (T_nu0) ** 4))

# Make the plot

import planckStyle as s

g = s.getSinglePlotter()

# Get planck 1-tail limit:
root = dataroot = s.defdata_all + '_lensing_BAO'
samps = g.getSamples('nnu', dataroot)
print('using %s' % dataroot)
print('Two tail:', samps.getLatex(params=['nnu'], limit=2))
p = samps.getParams()
samps.filter(p.nnu >= Neff_standard)
samps.setRanges({'nnu': [Neff_standard, None]})
samps.updateChainBaseStatistics()
limits = [samps.confidence('nnu', 0.32, upper=True) - Neff_standard,
          samps.confidence('nnu', 0.05, upper=True) - Neff_standard]
print('One tail limits:', limits)

log10T = np.linspace(-2, log10T_nodes[-1] + 1, 500)
T = 10. ** log10T


def mark_masses(ax):
    for m in [m_e, m_mu, m_b, m_W, m_t]:
        ax.axvline(m, color='gray', ls='--', lw=1)
    # These choices of bounds are fairly random:
    ax.axvspan(0.65, 2.8, color='gray', alpha=0.1, linewidth=0)
    ax.axvspan(130, 400, color='gray', alpha=0.1, linewidth=0)


# What to plot
wants = ['g', 'Neff']
plot_gs = False
do_data = 'Planck'
data_limit = {}
do_expts = False

data_limit['Planck'] = limits
data_limit['SO'] = [0.06, 0.12]
data_limit['S4'] = [0.03, 0.06]
data_limit['0.1 muK-arcmin'] = [0.014, 0.028]

fig, axs = plt.subplots(len(wants), 1, sharex='col', figsize=(5.5, 5.5/8*3.5 * len(wants)))
axs = np.atleast_1d(axs)

for i, (ax, want) in enumerate(zip(axs.reshape(-1), wants)):
    if want == 'g':
        if plot_gs: ax.semilogx(T, gs_full(log10T))
        ax.fill_between(T, table_gstar(T / 1000.) - table_gstar_err(T / 1000.),
                        table_gstar(T / 1000.) + table_gstar_err(T / 1000.),
                        color='C0', alpha=0.1)
        ax.semilogx(T, geff_full(log10T))
        if plot_gs:
            ax.legend(['$g_s$', '$g_*$'])
            ax.set_ylabel('$g$')
        else:
            ax.set_ylabel('$g_*$')
        ax.set_xlim([T[0], T[-1]]);
        ax.set_ylim([0, 110])
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.text(1.25, 78, 'Neutrino decoupling', rotation=-90, verticalalignment='top', horizontalalignment='center')
        ax.text(230, 78, 'QCD', rotation=0, fontsize=7, verticalalignment='top', horizontalalignment='center')

        ax2 = ax.twinx()
        # without Higgs, W,Z, t
        g_most = 5. / 6 * g_quarks + g_gluons + g_leptons + g_photon
        labs = [3.38, 10.75, g_most, g_SM]  # 3.38 allows for 3.046 neutrinos; 86.25 is all minus higgs, W,Z,top quark
        for lab in labs:
            ax.axhline(lab, color='gray', alpha=0.3)
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks(labs)
        ax2.get_yaxis().set_tick_params(which='both', direction='in', length=0)
        ax2.set_yticklabels(['%.5g' % lab for lab in labs])

    elif want == 'w':
        ax.semilogx(T, w_full(log10T))
        ax.set_ylabel(r'$w = P/\rho$')
        ax.axhline(1 / 3., color='gray', alpha=0.3)
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.set_ylim([0.2, 0.37])

    elif want == 'Neff':

        ax.loglog(T, Neff(log10T))
        ax.loglog(T, Neff(log10T) * 2)
        ax.loglog(T, Neff(log10T) * 7. / 8 * 2)
        ax.loglog(T, Neff(log10T) * 7. / 8 * 4)
        ax.axhline(1, color='gray', alpha=0.4)
        ax.axhline(0.57, color='gray', alpha=0.4)
        Neff_min = ((gs_full_therm0(full_T[-1]) / g_SM) ** (1. / 3) * TCMB) ** 4 / (2 * 7 / 8. * (T_nu0) ** 4)
        ax.legend(['Boson (g=1)', 'Boson (g=2)', 'Fermion (g=2)', 'Fermion (g=4)'], fontsize=7, loc='upper right')
        ax.set_ylim([1.1e-2, 20])
        ax.set_ylabel(r'$\Delta N_{\rm eff}$')
        ax.set_yticks([0.1, 1, 10])
        ax.set_yticklabels(['0.1', '1', '10'])
        if do_expts:
            for ix, (ex, limits) in enumerate(data_limit.items()):
                ax.axhline(limits[1], alpha=0.6, color='C%s' % (ix + 4), ls=':')
                ax.text(0.02, limits[1], ex + r' $2\sigma$', color='C%s' % (ix + 4),
                        verticalalignment='center',
                        bbox=dict(facecolor='white', ec='white', alpha=1, boxstyle="square,pad=0.25"))

        elif do_data:
            ax.axhspan(data_limit[do_data][1],20, color='gold', alpha=0.2)
            ax.axhspan(data_limit[do_data][0],20, color='gold', alpha=0.1)
#            ax.axhspan(0, data_limit[do_data][1], color='gold', alpha=0.1)
#            ax.axhspan(0, data_limit[do_data][0], color='gold', alpha=0.2)

        ax2 = ax.twinx()
        labs = [Neff_min, 7. / 8 * 2 * Neff_min, 7. / 8 * 4 * Neff_min, 0.57, 1]
        ax2.set_yscale('log')
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks(labs)
        ax2.get_yaxis().set_tick_params(which='both', direction='in', length=0)
        ax2.set_yticklabels(['%.2g' % lab for lab in labs])
        [t.set_color(x) for (x, t) in zip(['C0', 'C2', 'C3', 'C0', 'C2'], ax2.yaxis.get_ticklabels())]

    if i == 0:
        labs = [m_e, m_mu, m_b, m_W, m_t]
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xscale('log')
        ax2.set_xticks(labs)
        ax2.get_xaxis().set_tick_params(which='both', direction='in', length=0)
        ax2.set_xticklabels(['$m_e$', r'$m_\mu$', '$m_b$', '$m_W$', '$m_t$'])

    if i == len(wants) - 1:
        ax.set_xlabel(r'$T_\gamma$ [MeV]')

    mark_masses(ax)

plt.subplots_adjust(hspace=0)
plt.savefig('../../outputs/nnu_limits_decoupling.pdf', bbox_inches='tight')
