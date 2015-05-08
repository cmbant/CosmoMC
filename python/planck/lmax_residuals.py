import os, sys
from numpy import *
import matplotlib.pyplot as plt
from getdist import covmat, chains, types
import matplotlib.animation as animation

if False:
    d = loadtxt(r'z://base_plikHM_TT_lowTEB_lensing_50k_lenspotentialCls.dat')
    d2 = loadtxt(r'z://base_plikHM_TT_lowTEB_lensing_lenspotentialCls.dat')
#   d3 = loadtxt(r'z://base_plikHM_TT_lowTEB_lensing_10k_lenspotentialCls.dat')

    plt.plot(d2[:2000, 0], 1 - d2[:2000, 5] / d[:2000, 5], '-r')
#    plot(d3[:2000, 0], 1 - d3[:2000, 5] / d[:2000, 5], '-b')

    if False:
        mxs = [30, 18, 14, 12]
        for mx in mxs:
            d2 = loadtxt(r'z:\base_plikHM_TT_lowTEB_lensing_%uk_lenspotentialCls.dat' % mx)

            # plot(d[:, 0], d[:, 5])
            plt.plot(d2[:2000, 0], 1 - d2[:2000, 5] / d[:2000, 5], '-r')
            # savefig(r'z:\bestfit_cl_compare.pdf', bbox_inches='tight')
    plt.xlim([2, 400])
#    legend(mxs)
    plt.ylim([-0.02, 0.02])
    plt.show()
    sys.exit()


def plotFit(d, df=None, dref=None, plikLite=False, capt='', leg=True):
    plt.plot(d[:, 0], (d[:, 1] - LCDM[:, 1]), '-b')
    if dref is not None:
        plt.plot(dref[:, 0], (dref[:, 1] - LCDM[:, 1]), ':b')
    foreground = df is not None
    if foreground:
        plt.plot(df[:, 0], (df[:, 1] - LCDMf[:, 1]), ':c')
        plt.plot(df[:, 0], (df[:, 2] - LCDMf[:, 2]), '-r')
        plt.plot(df[:, 0], (df[:, 3] - LCDMf[:, 3]), '-g')
        plt.plot(df[:, 0], (df[:, 4] - LCDMf[:, 4]), '-m')

    plt.ylim([-80, 80])
    plt.axhline(0, color='k')
    # plot(d[:, 0], (d[:, 1] - d2[:, 1]) / (np.sqrt(2 * d[:, 1] ** 2 / (2 * d[:, 0] + 1))))
    if not plikLite:
        datapoints = loadtxt(r'C:\tmp\Planck\final_Nov14\residual_v9.10CMH_15_Nov_2014.dat')
        pts = datapoints[:, 3]
        sig = datapoints[:, 5]
    else:
        datapoints = loadtxt(r'C:\Users\antony\Dropbox\Planck\final_Nov14\residuals\PlikLite_v17spectra_all.txt')
        pts = datapoints[:, 1]
        sig = datapoints[:, 2]

    Lbin = datapoints[:, 0]
    for i, L in enumerate(Lbin):
        ix = np.where(d[:, 0] == L)
        pts[i] -= LCDM[ ix, 1]

    plt.errorbar(Lbin, pts, sig, fmt='.', color='k')
    plt.xlim([2, 2500])
    plt.xlabel('$l$')
    if capt: plt.text(1000, -70, capt)
    if foreground and leg:
        legs = [r'$\Delta D_l$']
        if dref is not None: legs += [r'$\Delta D_l^{\rm base}$']
        legs += [r'$\Delta D_l^{100}$', r'$\Delta D_l^{143}$', r'$\Delta D_l^{143\times 217}$', r'$\Delta D_l^{217}$']
        plt.legend(legs, ncol=2)
    # savefig(r'C:\Users\antony\Dropbox\Planck\final_Nov14\residuals\plikHM_TT_lowTEB_vs_Alens.pdf', bbox_inches='tight')
    # savefig(r'C:\Users\antony\Dropbox\Planck\final_Nov14\residuals\plikHM_TT_lowTEB_vs_Alens_plikLite.pdf', bbox_inches='tight')

# d = loadtxt(r'C:\tmp\Planck\final_Nov14\base_Alens\plikHM_TT_lowTEB\base_Alens_plikHM_TT_lowTEB.minimum.theory_cl')[:2499, :]
# df = loadtxt(r'C:\tmp\Planck\final_Nov14\base_Alens\plikHM_TT_lowTEB\base_Alens_plikHM_TT_lowTEB.minimum.plik_foregrounds')[:2499, :]
# d = loadtxt(r'C:\tmp\Planck\final_Nov14\base_CamSpecHM_TT_lowTEB_lmax800.minimum.theory_cl')[:2499, :]
# df = loadtxt(r'C:\tmp\Planck\final_Nov14\base_CamSpecHM_TT_lowTEB_lmax800.minimum.camspec_foregrounds')[:2499, :]

lmx = 2500
dat = 'CamSpecHM_TT_lowl_tau07_lmax' + str(lmx)
basepath = r'C:\tmp\Planck\final_Nov14\cuts\base'
root = basepath + os.sep + dat + os.sep + 'base_' + dat + '.minimum'
LCDM = loadtxt(root + '.theory_cl')[:2499, :]
LCDMf = loadtxt(root + '.camspec_foregrounds')[:2499, :]

# LCDM = loadtxt(r'C:\tmp\Planck\final_Nov14\base\plikHM_TT_lowTEB\base_plikHM_TT_lowTEB.minimum.theory_cl')[:2499, :]
# 'LCDMf = loadtxt(r'C:\tmp\Planck\final_Nov14\base\plikHM_TT_lowTEB\base_plikHM_TT_lowTEB.minimum.plik_foregrounds')[:2499, :]

if False:
    root = 'z://beam217.minimum'
    d = loadtxt(root + '.theory_cl')[:2499, :]
    df = loadtxt(root + '.camspec_foregrounds')[:2499, :]
    plotFit(d, df, False, leg=True)

if True:
    i = 0
    separate = False
    extparam = ['Alens']
    base = "_".join(['base'] + extparam)
    if not separate:
        lmaxs = range(800, 2501, 100)
        f, axs = plt.subplots(4, 4, figsize=(34, 30))
        axs = axs.reshape(-1)
    else:
        lmaxs = range(700, 2501, 50)
        f = plt.figure()
    plt.rc('text', usetex=True)
    def drawFor(i):
        lmx = lmaxs[i]
        if separate:
            plt.clf()
        else:
            if i >= len(axs): return
            plt.sca(axs[i])
        dat = 'CamSpecHM_TT_lowl_tau07_lmax' + str(lmx)
        path = r'C:\tmp\Planck\final_Nov14\cuts' + os.sep + base
        root = path + os.sep + dat + os.sep + base + '_' + dat + '.minimum'
        d = loadtxt(root + '.theory_cl')[:2499, :]
        df = loadtxt(root + '.camspec_foregrounds')[:2499, :]
        bestfit = types.BestFit(root)
        H0 = bestfit.parWithName('H0').best_fit
        f217 = bestfit.parWithName('f2000_217').best_fit
        label = r'$l_{\rm max}$: %u, $H_0$: %.2f, $f^{217}_{2000}$: %.1f ' % (lmx, H0, f217)
        if len(extparam):
            baseroot = basepath + os.sep + dat + os.sep + 'base_' + dat + '.minimum'
            basefit = types.BestFit(baseroot)
            delta = (bestfit.logLike - basefit.logLike) * 2
            delta_lowl = bestfit.parWithName('chi2_lowl').best_fit - basefit.parWithName('chi2_lowl').best_fit
            delta_CamSpec = bestfit.parWithName('chi2_CamSpec').best_fit - basefit.parWithName('chi2_CamSpec').best_fit
            tit = r' $\Delta\chi^2 = %.1f$, $\Delta\chi^2_{\rm lowl}=%.1f$, $\Delta\chi^2_{\rm CamSpec}=%.1f $' % (delta, delta_lowl, delta_CamSpec)
            based = loadtxt(baseroot + '.theory_cl')[:2499, :]
            plotFit(d, df, plikLite=False, dref=based, capt=label, leg=separate)
            for p in extparam:
                par = bestfit.parWithName(p)
                tit = '$%s = %.2f$; ' % (par.label, par.best_fit) + tit
            plt.title(tit, fontsize=10)
        else:
            plotFit(d, df, plikLite=False, capt=label, leg=separate)

        plt.axvline(lmx, color='gray', ls=':')
        if separate: plt.savefig(r'C:\tmp\Planck\final_Nov14\cuts\%s_Lmax%u.pdf' % (base, lmx), bbox_inches='tight')

    if separate:
        ani = animation.FuncAnimation(f, drawFor, range(len(lmaxs)), interval=800, blit=False, repeat_delay=2000)
        ani.save(r'C:\tmp\Planck\final_Nov14\cuts' + os.sep + base + '_lmax_fit.mp4')
    else:
        for i in range(len(lmaxs)):
            drawFor(i)
        plt.savefig(r'C:\tmp\Planck\final_Nov14\cuts' + os.sep + base + '_Lmax_panels.pdf', bbox_inches='tight')
plt.show()

sys.exit()

d = np.loadtxt(r'C:\Users\antony\Dropbox\Planck\final_Oct14\lensing\PostBorn\spectra_planck_uncorrected.data')
d2 = np.loadtxt(r'C:\Users\antony\Dropbox\Planck\final_Oct14\lensing\PostBorn\spectra_planck.data')
my = np.loadtxt(r'c:\work\dist\git\cosmomcplanck\camb\base_plikHM_TT_lowTEB_lensing_lenspotentialCls.dat')

# plot(d[:, 0], d[:, 1])

var = 0
cum = my[:, 5] * 0
for i in range(1, 6000):
    L = my[i, 0]
    var += (L + 0.5) * my[i, 5] / 4  # * 3.1415 * 2 / 4
    cum[i] = var

print 'rms = ', var
gca().set_xscale('log')

L = my[:, 0]
plot(L, cum)
savefig(r'C:\Users\antony\Dropbox\Planck\final_Oct14\lensing\PostBorn\rms_kappa_with_L.pdf', bbox_inches='tight')
show()
sys.exit()

plt.plot(L, L ** 2 * (d2[:, 1]) / (2 * 3.1415))
plot(L, L ** 2 * (d2[:, 3]) / (2 * 3.1415))

# plot(my[:, 0], my[:, 5] * (2 * 3.14159) / 4, '-r')
xlim([20, 2000])
ylim([1e-8, 0.01])
# ylim([0, 0.2])
gca().set_xscale('log')
gca().set_yscale('log')
# savefig(r'C:\Users\antony\Dropbox\Planck\final_Oct14\lensing\PostBorn\Curl.pdf', bbox_inches='tight')
show()
sys.exit()


c = chains.loadGridChain("C://tmp/Planck/clik9.0/", 'base', 'lensonly', ignore_frac=0.3)

p = c.getParams()

derived = p.A * (p.omegam ** 0.6 * p.H0 / 100) ** 2.3
print c.mean(derived), c.std(derived), c.std(derived) / c.mean(derived)

sys.exit()

cov = covmat.CovMat(os.path.dirname(os.path.abspath(__file__)) + r'/../planck_covmats/base_TTTEEE_lowTEB_plik.covmat')

cov.plot()
show()

sys.exit()

cov = os.path.dirname(os.path.abspath(__file__)) + r'/../planck_covmats/base_TT_lensing_lowTEB_plik.covmat'

with open(cov) as f:
    header = f.readline().split()[1:]
d = loadtxt(cov)

d = linalg.inv(d)

for i in range(len(d[:, 0])):
    print sqrt(1 / d[i, i]), header[i]

