import planckStyle as s
from getdist import inifile
from paramgrid import batchjob
from pylab import *
from scipy.interpolate import UnivariateSpline
import os
import tempfile
import pickle

sys.path.insert(0, 'C:\Work\Dist\git\camb')
from cosmomc_to_camb import get_camb_params
import camb

max_z=1.7
redshifts = [0., 0.01, 0.02, 0.05, 0.106, 0.15, 0.2, 0.25, 0.3, 0.32, 0.34, 0.38, 0.4, 0.44, 0.48, 0.51, 0.54, 0.57,
             0.60, 0.61, 0.65, 0.73, 0.85, 1., 1.2, 1.5, 2., 2.33, 2.5, 3.]

g = s.getSinglePlotter()

datavar = '_lensing'
root = 'base_plikHM_TTTEEE_lowl_lowE' + datavar
samples = g.sampleAnalyser.samplesForRoot(root)

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

DR12 = np.loadtxt(batchjob.getCodeRootPath() + 'data/DR12/final_consensus_results_dM_Hz_fsig.dat',
                  usecols=[0, 1])
DR12cov = np.loadtxt(batchjob.getCodeRootPath() + 'data/DR12/final_consensus_covtot_dM_Hz_fsig.txt')

p = samples.getParams()

reds = np.asarray(redshifts)

ommh2 = samples.mean('omegamh2')
# Need to get f sigma(8) of z using CAMB. This takes at least some minutes, so cache.

cachename = os.path.join(tempfile.gettempdir(), 'planck2018' + str(hash((samples.jobItem.chainRoot, os.path.getmtime(
    samples.jobItem.chainRoot + '_1.txt')) + tuple(redshifts)))
                         + '.fsig_evolve')
if os.path.isfile(cachename):
    with open(cachename, 'rb') as inp:
        ixs, f8s, Hs, DMs, FAPs = pickle.load(inp)
else:
    camb.set_z_outputs(redshifts)
    ixs = samples.randomSingleSamples_indices()[::4]
    DMs = np.zeros((len(ixs), len(redshifts)))
    Hs = np.zeros(DMs.shape)
    FAPs = np.zeros(DMs.shape)
    f8s = np.zeros(DMs.shape)
    for i, ix in enumerate(ixs):
        print(i, ix)
        dic = samples.getParamSampleDict(ix)
        pars = get_camb_params(dic)
        pars.set_matter_power(redshifts, kmax=2)
        results = camb.get_results(pars)
        bao = results.get_background_outputs()
        DMs[i, :] = bao[:, 2] * (1 + reds)
        Hs[i, :] = bao[:, 1]
        FAPs[i, :] = bao[:, 3]
        f8s[i, :] = results.get_fsigma8()[::-1]
        assert (abs(dic['fsigma8z038'] / f8s[i, redshifts.index(0.38)] - 1) < 0.001)

    with open(cachename, 'wb') as output:
        pickle.dump([ixs, f8s, Hs, DMs, FAPs], output, pickle.HIGHEST_PROTOCOL)

fsigmeans = np.zeros(len(redshifts))
fsigerrs = np.zeros(len(redshifts))
FAPmeans = np.zeros(len(redshifts))
FAPerrs = np.zeros(len(redshifts))
for i, z in enumerate(redshifts):
    fsigmeans[i] = np.mean(f8s[:, i])
    fsigerrs[i] = np.std(f8s[:, i])
    FAPmeans[i] = np.mean(FAPs[:, i])
    FAPerrs[i] = np.std(FAPs[:, i])
fsiginterp = UnivariateSpline(redshifts, fsigmeans, s=0)
fsigerrsinterp = UnivariateSpline(redshifts, fsigerrs, s=0)

for conditional in [False, True]:
    g = s.getSinglePlotter()

    redplot = np.arange(0, max_z, 0.02)
    means = fsiginterp(redplot)
    errs2 = fsigerrsinterp(redplot)
    fill_between(redplot, means - 2 * errs2, means + 2 * errs2, color='gray', alpha=0.15, lw=0.2)
    fill_between(redplot, means - errs2, means + errs2, color='gray', alpha=0.5, lw=0.2)
    plot(redplot, means, color='k', lw=0.5)

    name = '6dFGS'
    dataredshifts = [0.067]
    datapoints = [0.423]
    dataerr = [0.055]
    errorbar(dataredshifts, datapoints, dataerr, fmt='.', marker='*', color='g', markeredgecolor='g',
             markersize=3)
    text(0.005, 0.33, name, color='g', horizontalalignment='left', fontsize=5)

    # Howlett  arXiv:1409.3238
    name = 'SDSS MGS'
    dataredshifts = [0.15]
    if conditional:
        datapoints = [0.44]
        dataerr =[[0.12], [0.16]]
    else:
        datapoints = [0.53]
        dataerr = [0.19]
    errorbar(dataredshifts, datapoints, dataerr, fmt='.', marker='s', color='m', markeredgecolor='m',
             markersize=3)
    text(0.07,0.75, name, color='m', horizontalalignment='left')

    #https://arxiv.org/abs/1310.2820
    #Errors in https: // arxiv.org / abs / 1102.1014 are tighter, but conditional on WMAP
    name = 'SDSS LRG'
    dataredshifts = [0.3]
    datapoints = [0.49]
    dataerr = [sqrt((0.08) ** 2 + 0.04 ** 2)]
    errorbar(dataredshifts, datapoints, dataerr, fmt='.', marker='x', color='c', markeredgecolor='c',
             markersize=3)
    text(0.22,0.6, name, color='c', horizontalalignment='left')

    # SN Ia and galaxy velocities
    # https://arxiv.org/abs/1611.09862
    name = '6dFGS\n+SnIa'
    dataredshifts += [0.02]
    datapoints += [0.428]
    dataerr += [0.0465]
    errorbar(dataredshifts, datapoints, dataerr, fmt='.', marker='^', color='darkcyan', markeredgecolor='darkcyan',
             markersize=3)
    text(0.003, 0.5, name, color='darkcyan', horizontalalignment='left', fontsize=5)

     #GAMA https://arxiv.org/abs/1309.5556
    dataredshifts = [0.18,0.38-0.008]
    datapoints = [0.36, 0.44]
    dataerr  = [0.09, 0.06]
    col= 'darkred'
    errorbar(dataredshifts, datapoints, dataerr, fmt='.', marker='*', color=col, markeredgecolor=col,
             markersize=3)
    text(0.3, 0.3, 'GAMA', color=col, horizontalalignment='center')


    # Vipers Also two-point updated version in https://arxiv.org/abs/1612.05645
    # Note 1708.00026 shifts upper point
    dataredshifts = [0.6-0.008,0.86]
    datapoints = [0.55, 0.4]
    dataerr  = [0.12, 0.11]
    col= 'olive'
    errorbar(dataredshifts, datapoints, dataerr, fmt='.', marker='*', color=col, markeredgecolor=col,
             markersize=3)
    text(0.63, 0.54, 'VIPERS', color=col, horizontalalignment='left')


    # DR12
    dataredshifts = DR12[::3, 0]  # [0.38, 0.51, 0.61]
    if not conditional:
        datapoints = DR12[2::3, 1]
        dataerr = np.sqrt(np.diag(DR12cov))[2::3]
    else:
        data = DR12[:, 1]
        invcov = np.linalg.inv(DR12cov)
        fixedix = [0, 1, 3, 4, 6, 7]
        wantix = [2, 5, 8]
        modelix = [redshifts.index(x) for x in dataredshifts if x in redshifts]
        assert (len(modelix) == 3)
        model = np.zeros(9)
        rdfid = 147.78
        condcov = np.linalg.inv(invcov[np.ix_(wantix, wantix)])
        covav = np.zeros((3, 3))
        meanav = np.zeros(3)
        for i, ix in enumerate(ixs):
            rdrag = samples.getParamSampleDict(ix)['rdrag']
            model[::3] = DMs[i, modelix] * rdfid / rdrag
            model[1::3] = Hs[i, modelix] * rdrag / rdfid
            deltamean = condcov.dot(invcov[np.ix_(wantix, fixedix)]).dot(
                model[fixedix] - data[fixedix])
            meanav += deltamean
            covav += np.outer(deltamean, deltamean)
        meanav /= len(ixs)
        covav = condcov + covav / len(ixs) - np.outer(meanav, meanav)
        print(condcov, covav, meanav)
        datapoints = data[wantix] - meanav
        dataerr = np.sqrt(np.diag(covav))
    errorbar(dataredshifts, datapoints, dataerr, fmt='.', color='orangered')
    text(0.4, 0.52, 'BOSS\nDR12', color='orangered', horizontalalignment='left', fontsize=6)

    # #Wigglez 1204.3674
    wigredshifts = [0.44, 0.6, 0.73]
    wigdata = np.array([0.413, 0.390, 0.437])
    wigerr = np.array([0.08, 0.063, 0.072])
    if conditional:
        WigFAPdata = np.array([0.482, 0.650, 0.865])
        # covariance for FAP FAP FAP fis8 fsig8 fsig8
        if True:
            modelix = [redshifts.index(x) for x in wigredshifts if x in redshifts]
            assert (len(modelix) == 3)
            wigcov = np.array([[2.401, 1.350, 0.000, 2.862, 1.080, 0.000],
                               [1.350, 2.809, 1.934, 1.611, 2.471, 1.641],
                               [0, 1.934, 5.329, 0, 1.978, 4.468],
                               [2.862, 1.611, 0, 6.400, 2.570, 0.],
                               [1.080, 2.471, 1.978, 2.570, 3.969, 2.540],
                               [0, 1.641, 4.468, 0, 2.540, 5.184]
                               ]) * 1e-3
            cinv = inv(wigcov)
            cov_fsig_inv = cinv[3:, 3:]
            cov_fsig = inv(cov_fsig_inv)

            if True:
                meanav = np.zeros(3)
                covav = np.zeros((3, 3))
                for i, ix in enumerate(ixs):
                    deltamean = cov_fsig.dot(cinv[3:, 0:3]).dot(FAPs[i, modelix] - WigFAPdata)
                    meanav += deltamean
                    covav += np.outer(deltamean, deltamean)
                meanav /= len(ixs)
                covav = cov_fsig + covav / len(ixs) - np.outer(meanav, meanav)
                cond_wigdata = wigdata - meanav
                dataerr = np.sqrt(np.diag(covav))
                print ('wig1', cond_wigdata, dataerr)
            else:  # agrees well (correction from Planck scatter is small anyway)
                WigFAP = np.interp(wigredshifts, redshifts, FAPmeans)
                cond_wigdata = wigdata - cov_fsig.dot(cinv[3:, 0:3]).dot(WigFAP - WigFAPdata)
                WigFAPerr = np.interp(wigredshifts, redshifts, FAPerrs)
                FAPtheoryErr = cov_fsig.dot(cinv[3:, 0:3]).dot(WigFAPerr)
                diagerr = sqrt(np.diagonal(cov_fsig))
                dataerr = sqrt(diagerr ** 2 + FAPtheoryErr ** 2)
                print ('wig2', cond_wigdata, sqrt(diagerr ** 2 + FAPtheoryErr ** 2))

        errorbar(wigredshifts, cond_wigdata, dataerr, fmt='.', color='#006FED')
        text(0.53, 0.37, 'WiggleZ', color='#006FED', horizontalalignment='left')
    else:
        errorbar(wigredshifts, wigdata, wigerr, fmt='.', color='#006FED')
        text(0.51, 0.27, 'WiggleZ', color='#006FED', horizontalalignment='left')

#SDSS quasars https://arxiv.org/abs/1801.02689
    if not conditional:
#        errorbar([1.52], [0.42], [0.076], fmt='o', color='orange')
        errorbar([1.52], [0.426], [0.077], fmt='o', color='orange')
        text(1.58, 0.31, 'DR14 quasars', color='orange', horizontalalignment='right')
        xlim([0, 1.6])
    else:
        xlim([0, 0.9])

#FastSound
    errorbar([1.4], [0.482], [0.116], fmt='+', color='darkblue')
    text(1.4, 0.62, 'FastSound', color='darkblue', horizontalalignment='center')
    xlim([0, 1.6])

    if conditional:
        ylim([0.2, 0.8])
    else:
        ylim([0.2, 0.8])
    xlabel('$z$')
    ylabel(r'$f\sigma_8$')

    if conditional:
        g.export(tag='conditional' + datavar)
    else:
        g.export(tag=datavar)
