import planckStyle as s
from getdist import inifile
from paramgrid import batchjob
from pylab import *
from scipy.interpolate import UnivariateSpline

# Uses CAMB 1.0+
sys.path.insert(0, r'C:\Work\Dist\git\camb')
from cosmomc_to_camb import get_camb_params
import camb

redshifts = [0.01, 0.05, 0.106, 0.15, 0.2, 0.25, 0.3, 0.34, 0.38, 0.42, 0.46, 0.51, 0.56, 0.61, 0.65, 0.72, 0.85, 1.,
             1.2, 1.5, 2., 2.33, 2.5,
             3.]

updated2019 = True

plot_DV = False  # if false plot D_V or D_M depending on measurement

reds = np.asarray(redshifts)

g = s.getSinglePlotter()

datavar = '_lensing'
samples = g.sampleAnalyser.samplesForRoot('base_plikHM_TTTEEE_lowl_lowE' + datavar)
variant = None  # e.g. to add another set of bands in H(z) 'base_nnu_plikHM_TTTEEE_lowl_lowE'

datasets = ['sdss_6DF_bao.dataset', 'sdss_MGS_bao.dataset', 'sdss_DR14_quasar_bao.dataset']
names = ['6DFGS', 'SDSS\nMGS', 'SDSS quasars']

dataredshifts = np.zeros(len(datasets))
datapoints = np.zeros(len(datasets))
dataerrs = np.zeros(len(datasets))
colors = ['g', 'm', 'r', 'darkred']

for i, dat in enumerate(datasets):
    if '.dataset in dat':
        ini = inifile.IniFile(batchjob.getCodeRootPath() + 'data/' + dat)
        datapoints[i], dataerrs[i] = [float(f) for f in ini.split('bao_measurement')]
        dataredshifts[i] = ini.float('zeff')
        rescale = ini.float('rs_rescale', 1.)
        tp = ini.string('measurement_type')
        if tp == 'rs_over_DV':
            dataerrs[i] = -0.5 / (datapoints[i] + dataerrs[i]) + 0.5 / (datapoints[i] - dataerrs[i])
            datapoints[i] = 1 / datapoints[i]
        elif tp != 'DV_over_rs':
            raise Exception('error')
        datapoints[i] *= rescale
        dataerrs[i] *= rescale
    print(dataredshifts[i], datapoints[i], dataerrs[i])


def GetBackgroundFuncs(samples):
    ixs = samples.randomSingleSamples_indices()[::40]
    DMs = np.zeros((len(ixs), len(redshifts)))
    Hs = np.zeros(DMs.shape)
    rsDV = np.zeros(DMs.shape)
    for i, ix in enumerate(ixs):
        print(i, ix)
        dic = samples.getParamSampleDict(ix)
        pars = get_camb_params(dic)
        pars.z_outputs = redshifts
        results = camb.get_background(pars)
        bao = results.get_background_outputs()
        rsDV[i, :] = 1 / bao[:, 0]
        DMs[i, :] = bao[:, 2] * (1 + reds)
        Hs[i, :] = bao[:, 1]

    Hmeans = np.zeros(len(redshifts))
    Herrs = np.zeros(len(redshifts))
    DMmeans = np.zeros(len(redshifts))
    DMerrs = np.zeros(len(redshifts))
    for i, z in enumerate(redshifts):
        Hmeans[i] = np.mean(Hs[:, i]) / (1 + z)
        Herrs[i] = np.std(Hs[:, i]) / (1 + z)
        DMmeans[i] = np.mean(DMs[:, i])
        DMerrs[i] = np.std(DMs[:, i])

    Hinterp = UnivariateSpline([0] + redshifts, [samples.mean('H0')] + list(Hmeans), s=0)
    DMinterp = UnivariateSpline([0] + redshifts, [0] + list(DMmeans), s=0)
    Herrinterp = UnivariateSpline([0] + redshifts, [samples.std('H0')] + list(Herrs), s=0)
    DMerrinterp = UnivariateSpline([0] + redshifts, [0] + list(DMerrs), s=0)
    return Hinterp, Herrinterp, DMinterp, DMerrinterp, rsDV


Hinterp, Herrinterp, DMinterp, DMerrinterp, rsDV = GetBackgroundFuncs(samples)

rdrag = samples.mean('rdrag')
c = 299792
if True:  # H(z) plot

    redplot = np.arange(0, 3, 0.05)
    if variant:
        samples2 = g.sampleAnalyser.samplesForRoot(variant)
        Hinterp2, Herrinterp2, _, _, _ = GetBackgroundFuncs(samples2)
        s.plotBands(redplot, Hinterp2(redplot), Herrinterp2(redplot), color='blue')

    Hmeans2 = Hinterp(redplot)
    Herrs2 = Herrinterp(redplot)
    s.plotBands(redplot, Hmeans2, Herrs2)

    if updated2019:
        # Aganthe https://arxiv.org/abs/1904.03400
        z = 2.34
        sdssdata = [c / (8.86 * rdrag) / (1 + z)]
        sdsserr = [c / (8.86 * rdrag) * 0.29 / 8.86 / (1 + z)]
        print('Ly-alpha H:', sdssdata, sdsserr)
        sdssredshifts = [z-0.007]
        err = errorbar(sdssredshifts, sdssdata, sdsserr, fmt='.', color='orange', marker='x',
                       markersize=3)
        err[-1][0].set_linestyle('--')
        text(2.33, 63, r'DR14 Ly-$\alpha$', color='orange', horizontalalignment='center')
        # joint with cross-correlation https://arxiv.org/abs/1904.03430
        z = 2.34
        sdssdata = [c / (9.0 * rdrag) / (1 + z)]
        sdsserr = [c / (9.0 * rdrag) * 0.22 / 9.0 / (1 + z)]
        sdssredshifts = [z + 0.007]
    else:
        # Bautista https://arxiv.org/abs/1702.00176
        z = 2.33
        sdssredshifts = [z]
        sdssdata = [c / (9.07 * rdrag) / (1 + z)]
        sdsserr = [c / (9.07 * rdrag) * 0.31 / 9.07 / (1 + z)]
        print('Ly-alpha H:', sdssdata, sdsserr)
        err = errorbar(sdssredshifts, sdssdata, sdsserr, fmt='.', color='orange', marker='x',
                       markersize=3)
        err[-1][0].set_linestyle('--')
        text(2.33, 63, r'BOSS Ly-$\alpha$', color='orange', horizontalalignment='center')
        # arXiv: du Mas des Bourboux 1708.02225
        z = 2.4
        sdssredshifts = [z]
        sdssdata = [c / (8.94 * rdrag) / (1 + z)]
        sdsserr = [c / (8.94 * rdrag) * 0.22 / 8.94 / (1 + z)]
    print('Ly-alpha H joint:', sdssdata, sdsserr)
    errorbar(sdssredshifts, sdssdata, sdsserr, fmt='.', color='gold', marker='x',
             markersize=3)
    #    text(2.33, 63, r'BOSS Ly-$\alpha$', color='orange', horizontalalignment='center')

    pts = np.loadtxt(batchjob.getCodeRootPath() + 'data/DR12/sdss_DR12Consensus_bao.dat', usecols=[0, 1])
    cov = np.loadtxt(batchjob.getCodeRootPath() + 'data/DR12/BAO_consensus_covtot_dM_Hz.txt')
    rdfid = 147.78
    sdssredshifts = pts[1::2, 0]
    sdssdata = pts[1::2, 1] * rdfid / rdrag / (1 + sdssredshifts)
    sdsserr = np.sqrt(np.diag(cov)[1::2]) * rdfid / rdrag / (1 + sdssredshifts)
    print('sdss DR12 H:', sdssdata, sdsserr)
    errorbar(sdssredshifts, sdssdata, sdsserr, fmt='.', color='orangered', marker='^',
             markersize=3)
    text(0.5, 63.5, 'BOSS DR12', color='orangered', horizontalalignment='center')
    if variant:
        rdrag2 = samples2.mean('rdrag')
        sdssdata = pts[1::2, 1] * rdfid / rdrag2 / (1 + sdssredshifts)
        sdsserr = np.sqrt(np.diag(cov)[1::2]) * rdfid / rdrag2 / (1 + sdssredshifts)
        errorbar(sdssredshifts + 0.05, sdssdata, sdsserr, fmt='.', color='darkred', marker='.',
                 markersize=3)

    # convert Dv/rs to H using Planck values of DM and rdrag
    if False:
        z = 0.15
        sdssredshifts = [z]
        sdssdata = [DMinterp(z) ** 2 * z / (4.465666824 * rdrag) ** 3 * c / (1 + z)]
        sdsserr = [DMinterp(z) ** 2 * z / ((4.465666824 - .168135) * rdrag) ** 3 * c / (1 + z) - sdssdata[0]]
        print('sdss MGS H:', sdssdata, sdsserr)
        errorbar(sdssredshifts, sdssdata, sdsserr, fmt='.', color='m', marker='s',
                 markersize=3)
        text(0.1, 56, 'SDSS MGS', color='m', horizontalalignment='left')

    if updated2019:
        H0red = [0.0075]
        H0red = [0]
        H0data = [74.03]
        H0err = [1.42]
        errorbar(H0red, H0data, H0err, fmt='.', color='blue', marker='.', markersize=3)
        text(0.06, 74, 'Riess et al. (2019)', color='blue', horizontalalignment='left')
    else:
        H0red = [0.0075]
        H0red = [0]
        H0data = [73.45]
        H0err = [1.66]
        errorbar(H0red, H0data, H0err, fmt='.', color='blue', marker='.', markersize=3)
        text(0.06, 73, 'Riess et al. (2018)', color='blue', horizontalalignment='left')

    if False:
        # DR14 quasar http://arxiv.org/pdf/1801.03043
        rdfid = 147.78
        H0red = np.array([0.978, 1.23, 1.526, 1.944])
        H0data = np.array([113.72, 131.44, 148.11, 172.63]) * rdfid / rdrag / (1 + H0red)
        H0err = np.array([14.63, 12.42, 12.75, 14.79]) * rdfid / rdrag / (1 + H0red)
        errorbar(H0red, H0data, H0err, fmt='.', color='darkgreen', marker='o', markersize=3)
        text(1.05, 52, 'DR14 Quasars', color='darkgreen', horizontalalignment='left')
    else:
        # DR14 quasar http://arxiv.org/abs/1801.02689, single bin
        # Note three bins are all systematically lower than one bin. Hector explained he thinks this is due
        # to non-Gaussian tails. Hence plotting 1 bin probably fairly visual comparison with Planck.
        rdfid = 147.78
        H0red = np.array([1.52])
        if False:  # Gil-Marin
            H0data = np.array([162]) * rdfid / rdrag / (1 + H0red)
            H0err = np.array([12]) * rdfid / rdrag / (1 + H0red)
            text(1.2, 58, 'DR14 quasars', color='darkgreen', horizontalalignment='left')
        else:  # Zarouk 1801.03062
            H0data = np.array([159]) * rdfid / rdrag / (1 + H0red)
            H0err = np.atleast_2d([13, 12]).T * rdfid / rdrag / (1 + H0red)
            text(1.2, 57, 'DR14 quasars', color='darkgreen', horizontalalignment='left')
        errorbar(H0red, H0data, H0err, fmt='.', color='darkgreen', marker='o', markersize=3)

    xlim([-0.015, 2.7])
    # axvline(0,color='k',lw=0.5, zorder=-10)
    ylim([55, 77])
    yticks(np.arange(54, 76.1, 2))
    xlabel('$z$')
    ylabel(r'$H(z)/(1+z)\, [{\rm km}\,{\rm s}^{-1}\,{\rm Mpc}^{-1}]$')
    g.export(tag=datavar.replace('_', '') + '_Hz')

g = s.getSinglePlotter()
rsDVmeans = np.zeros(len(redshifts))
rsDVerrs = np.zeros(len(redshifts))

for i in range(len(redshifts)):
    rsDVmeans[i] = np.mean(rsDV[:, i])
    rsDVerrs[i] = np.std(rsDV[:, i])

planckvals = np.interp(dataredshifts, redshifts, rsDVmeans)
datapoints /= planckvals
dataerrs /= planckvals

print(planckvals)
print(datapoints)
print(dataerrs)

s.plotBands(redshifts, 1, rsDVerrs / rsDVmeans)
# plt.plot(redshifts, 1+DMerrinterp(redshifts)/DMinterp(redshifts), color='r', lw=0.5)

offsets = [0.93, 1.08, 0.94]
for i, (form, col, name, offset) in enumerate(
        zip(['*', 's', 'o'], ['g', 'm', 'red'], names, offsets)):
    errorbar(dataredshifts[i:i + 1], datapoints[i:i + 1], dataerrs[i:i + 1], fmt='.', marker=form, color=col,
             markersize=3,
             markeredgecolor=col)
    text(max(0.2, dataredshifts[i]), offset, name, color=col, horizontalalignment='center')

rs_fid_wig = 148.6
wigredshifts = [0.44, 0.6, 0.73]
wigdata = np.array([1716., 2221., 2516.]) / rs_fid_wig
wigerr = np.array([83., 101., 86.]) / rs_fid_wig
print('wig:', wigdata)
print('wig err:', wigerr)
planckvals = np.interp(wigredshifts, redshifts, rsDVmeans)
errorbar(wigredshifts, wigdata / planckvals, wigerr / planckvals, fmt='.', color='#006FED', marker='o', markersize=1.5)
text(0.72, 1.08, 'WiggleZ', color='#006FED', horizontalalignment='center')

if False:  # Two less correlated points
    sdssredshifts = [0.32, 0.57]
    sdssdata = [1270 / 147.78, 2033 / 147.78]
    sdsserr = [14 / 147.78, 21 / 147.78]
    planckvals = np.interp(sdssredshifts, redshifts, rsDVmeans)
    print('sdss DR12:', sdssdata, sdsserr, planckvals)
    errorbar(sdssredshifts, sdssdata / planckvals, sdsserr / planckvals, fmt='.', color='orangered', marker='^',
             markersize=3)
    text(0.5, 0.95, 'BOSS\nDR12', color='orangered', horizontalalignment='center')
else:  # three correlated points
    pts = np.loadtxt(batchjob.getCodeRootPath() + 'data/DR12/BAO_consensus_results_dV_FAP.dat', usecols=[0, 1])
    cov = np.loadtxt(batchjob.getCodeRootPath() + 'data/DR12/BAO_consensus_covtot_dV_FAP.txt')
    sdssredshifts = pts[::2, 0]
    sdssdata = pts[::2, 1]
    sdsserr = np.sqrt(np.diag(cov)[::2])
    planckvals = np.interp(sdssredshifts, redshifts, rsDVmeans)
    print('sdss DR12:', sdssdata, sdsserr, planckvals)
    errorbar(sdssredshifts, sdssdata / planckvals, sdsserr / planckvals, fmt='.', color='orangered', marker='^',
             markersize=3)
    text(0.5, 0.94, 'BOSS\nDR12', color='orangered', horizontalalignment='center')

# Ly-alpha
# https://arxiv.org/abs/1702.00176
# https://arxiv.org/abs/1708.02225
# Take numbers of abstract of 1708.02225 and naively combine
redshift = 2.4
if plot_DV:
    data = [(2.4 * 8.94 * 36.6 ** 2) ** (1. / 3)]
    err = [
        ((2.4 * (8.94 + 0.22) * (36.6 + 1.2) ** 2) ** (1. / 3) - (2.4 * (8.94 - 0.22) * (36.6 - 1.2) ** 2) ** (
                1. / 3)) / 2]
    planckvals = np.interp([redshift], redshifts, rsDVmeans)
    print('Ly-alpha:', data, err, planckvals)
else:
    if updated2019:
        # Ly-alpha and cross-correlation: https://arxiv.org/abs/1904.03430 (quoted to different precision in the two papers)
        redshift = 2.34
        data = np.array([36.981])
        err = np.atleast_2d([1.18,1.26])
    else:
        # DM https: // arxiv.org / abs / 1708.02225, Bourboux
        data = np.array([36.6])
        err = np.array([1.2])
    planckvals = DMinterp(redshift) / rdrag
    print('Ly-alpha:', data, err, planckvals)
    errorbar([redshift], data / planckvals, err / planckvals, fmt='.', color='orange', marker='H', markersize=3)
if updated2019:
    text(2.00, 0.973, r'DR14 Ly-$\alpha$/'+'\n\\quad quasars ($D_M$)', color='orange', horizontalalignment='left', fontsize=7.5)
else:
    text(2.23, 0.97, r'BOSS Ly-$\alpha$ ($D_M$)', color='orange', horizontalalignment='center')

# DES, convert from DM using PLanck H
# https://arxiv.org/abs/1712.06209
redshift = 0.81
if not plot_DV:
    planckvals = [DMinterp(redshift) / rdrag]
    data = np.array([10.75 * (1 + redshift)])
    err = np.array([0.43 * (1 + redshift)])
    print('DES:', data, err, planckvals)
else:
    planckvals = np.interp([redshift], redshifts, rsDVmeans)
    pl_err = DMerrinterp(0.81)

    data = planckvals * (10.75 * (1 + 0.81) * rdrag / DMinterp(0.81)) ** (2. / 3)
    err = (planckvals * ((10.75 + 0.43) * (1 + 0.81) * rdrag / (DMinterp(0.81) - pl_err)) ** (2. / 3)
           - planckvals * ((10.75 - 0.43) * (1 + 0.81) * rdrag / (DMinterp(0.81) + pl_err)) ** (2. / 3)) / 2
    print('DES:', sdssdata, sdsserr, planckvals)
errorbar([redshift], data / planckvals, err / planckvals, fmt='.', color='darkblue', marker='v',
         markersize=4)
text(0.82, 1.03, r'DES ($D_M$)', color='darkblue', horizontalalignment='left')

# eBoss 0.72 DR14 LRG
# https://arxiv.org/abs/1712.08064
# Quote h in units???? assume typo
sdssredshifts = [0.72]
sdssdata = [2353 / 147.78]
sdsserr = [62 / 147.78]
planckvals = np.interp(sdssredshifts, redshifts, rsDVmeans)
print('DR14 LRG:', sdssdata, sdsserr, planckvals)
errorbar(sdssredshifts, sdssdata / planckvals, sdsserr / planckvals, fmt='.', color='darkcyan', marker='x',
         markersize=4)
text(0.81, 0.91, r'DR14 LRG', color='darkcyan', horizontalalignment='center')

# 6Df arXiv:1803.01746
sdssredshifts = [0.122]
sdssdata = [539 / 147.5]
sdsserr = [17 / 147.5]
planckvals = np.interp(sdssredshifts, redshifts, rsDVmeans)
print('6DF:', sdssdata, sdsserr, planckvals)
err = errorbar(sdssredshifts, sdssdata / planckvals, sdsserr / planckvals, fmt='.', color='m', marker='o', mfc='C2',
               mec='C2',
               markersize=2, lw=1)
err[-1][0].set_linestyle('--')

xlim([0.01, 2.7])
ylim([0.88, 1.11])
xlabel('$z$')
if plot_DV:
    ylabel(r'$(D_{\rm V}/r_{\rm drag})/(D_{\rm V}/r_{\rm drag})_{\rm \it Planck}$')
else:
    ylabel(r'$(D/r_{\rm drag})/(D/r_{\rm drag})_{\rm \it Planck}$')
gca().set_yticks(np.arange(0.9, 1.11, 0.05))
g.export(tag=datavar.replace('_', ''))
