import planckStyle as s
from paramgrid import batchjob
import GetDistPlots
import pylab as plt
import numpy as np

g = s.getSubplotPlotter(subplot_size=4)

rd_fid = 147.78

FAP = True
if FAP:
    cov = np.loadtxt(batchjob.getCodeRootPath() + 'data/DR12/final_consensus_covtot_dV_FAP_fsig.txt')
    pts = np.loadtxt(batchjob.getCodeRootPath() + 'data/DR12/final_consensus_results_dV_FAP_fsig.dat', usecols=[0, 1])
else:
    cov = np.loadtxt(batchjob.getCodeRootPath() + 'data/DR12/final_consensus_covtot_dM_Hz_fsig.txt')
    pts = np.loadtxt(batchjob.getCodeRootPath() + 'data/DR12/sdss_DR12Consensus_final.dat', usecols=[0, 1])

pnames = ['DM038', 'Hubble038', 'fsigma8z038', 'DM051', 'Hubble051', 'fsigma8z051', 'DM061', 'Hubble061', 'fsigma8z061']
redshifts = pts[:, 0]
data = pts[:, 1]
planckmeans = []


def BAOdensity(p1, p2, marge=True):
    err = np.sqrt(cov[p1, p1])
    DAv = np.arange(data[p1] - 4 * err, data[p1] + 4 * err, err / 100)
    err = np.sqrt(cov[p2, p2])
    Hv = np.arange(data[p2] - 4 * err, data[p2] + 4 * err, err / 100)
    DA, H = np.meshgrid(DAv, Hv)
    if marge:
        v1 = data[p1]
        v2 = data[p2]
        mcov = cov[np.ix_([p1, p2], [p1, p2])]
        invcov = np.linalg.inv(mcov)
        like = (DA - v1) ** 2 * invcov[0, 0] + 2 * (DA - v1) * (H - v2) * invcov[0, 1] + (H - v2) ** 2 * invcov[
            1, 1]
    else:
        # Not done yet
        invcov = np.linalg.inv(condcov[np.ix_([p1, p2], [p1, p2])])
        v1 = conddata[p1]
        v2 = conddata[p2]
        like = (DA - v1) ** 2 * invcov[0, 0] + 2 * (DA - v1) * (H - v2) * invcov[0, 1] + (H - v2) ** 2 * invcov[
            1, 1]

    density = GetDistPlots.Density2D(DAv, Hv, np.exp(-like / 2))
    density.contours = np.exp(-np.array([1.509, 2.4477]) ** 2 / 2)
    return density


rootTT = 'base_plikHM_TT_lowl_lowE'
root = 'base_plikHM_TTTEEE_lowl_lowE'
rootL = 'base_plikHM_TTTEEE_lowl_lowE_lensing'

c = 299792.458


def makeNewSamples(root):
    samples = g.sampleAnalyser.samplesForRoot(root)
    p = samples.getParams()
    pars = []
    for p1, p2, p3 in zip(range(0, 9, 3), range(1, 9, 3), range(2, 9, 3)):
        if not FAP:
            DM = getattr(p, pnames[p1]) * rd_fid / p.rdrag
            rsH = getattr(p, pnames[p2]) * p.rdrag / rd_fid

            par1 = samples.addDerived(DM, name=pnames[p1] + 'n',
                                      label=r'D_{\rm M}(%.2f) (r_{\mathrm{drag}}^{\rm fid}/r_{\mathrm{drag}})\,[\rm{Mpc}]' %
                                            redshifts[p1])
            par2 = samples.addDerived(rsH, name=pnames[p2] + 'n',
                                      label=r'H(%.2f) (r_{\mathrm{drag}}/r_{\mathrm{drag}}^{\rm fid})\, [{\rm km} \,{\rm s}^{-1}{\rm Mpc}^{-1}]' % (
                                          redshifts[p2]))
        else:
            AP = getattr(p, pnames[p1]) * getattr(p, pnames[p2]) / c
            par2 = samples.addDerived(AP, name=pnames[p1] + 'FAP',
                                      label=r'F_{\rm AP(%.2f)}' % redshifts[p1])
            par1 = par2

        pars += [par1, par2, samples.paramNames.parWithName(pnames[p3])]
        #
    samples.updateBaseStatistics()
    return samples, pars


samples, pars = makeNewSamples(root)
samplesT, parsTT = makeNewSamples(rootTT)
samplesL, parsL = makeNewSamples(rootL)


def makeConditionalData(samples):
    invcov = np.linalg.inv(cov)
    fixedix = [0, 3, 6]  # DV/rs
    wantix = [p for p in range(9) if not p in fixedix]
    ixs = samples.randomSingleSamples_indices()[::4]
    meanav = np.zeros(6)
    covav = np.zeros((6, 6))
    theory = np.zeros(3)
    cov_FAP_fsig = np.linalg.inv(invcov[np.ix_(wantix, wantix)])
    for i, ix in enumerate(ixs):
        pars = samples.getParamSampleDict(ix)
        for j, (z, tag) in enumerate(zip([0.38, 0.51, 0.61], ['038', '051', '061'])):
            DM = pars['DM' + tag]
            H = pars['Hubble' + tag] / c
            theory[j] = pow(DM ** 2 / H * z, 1. / 3) / pars['rdrag']  # DV/rd
        deltamean = cov_FAP_fsig.dot(invcov[np.ix_(wantix, fixedix)]).dot(theory - data[fixedix])
        meanav += deltamean
        covav += np.outer(deltamean, deltamean)
    meanav /= len(ixs)
    covav = cov_FAP_fsig + covav / len(ixs) - np.outer(meanav, meanav)
    condcov = cov.copy()
    condcov[np.ix_(wantix, wantix)] = covav
    conddata = data.copy()
    conddata[wantix] = data[wantix] - meanav
    return conddata, condcov


conddata, condcov = makeConditionalData(samples)

for p1, p2, p3, z in zip(pars[::3], pars[1::3], pars[2::3], redshifts[::3]):
    print('p1 %s' % z, samples.mean(p1), samples.std(p1))
    print('p2 %s' % z, samples.mean(p2), samples.std(p2))
    print('fs8 %s' % z, samples.mean(p3), samples.std(p3))
    planckmeans += [samples.mean(p1), samples.mean(p2), samples.mean(p3)]

    # print samples.PCA([p1.name, p2.name], 'LLL', p1.name)


def doplot(ix1, ix2, ax=None):
    density = BAOdensity(ix1, ix2, marge=True)
    g.add_2d_contours(root, pars[ix1], pars[ix2], filled=True, density=density, ax=ax)

    g.add_2d_contours(rootTT, parsTT[ix1], parsTT[ix2], filled=True, ax=ax, plotno=1)
 #   g.add_2d_contours(root, pars[ix1], pars[ix2], filled=True, ax=ax, plotno=2)
    g.add_2d_contours(rootL, parsL[ix1], parsL[ix2], filled=True, ax=ax, plotno=3)

    density = BAOdensity(ix1, ix2, marge=False)
    g.add_2d_contours(root, pars[ix1], pars[ix2], filled=False, density=density, ls='--', ax=ax)

    #   g.add_legend([r'DR12 ($z_{\rm eff}=%s$)' % redshifts[ix1]], legend_loc='upper left', ax=ax)
    g.add_text(r'$z_{\rm eff}=%s$' % redshifts[ix1], fontsize=16)
    # g.setAxes([pars[ix1], pars[ix2]], lims=[1350, 1500, 85, 110])
    g.setAxes([pars[ix1], pars[ix2]], lims=[density.x[0], density.x[-1], density.y[0], density.y[-1]], ax=ax)


g.make_figure(nx=3)
for i in range(0, 9, 3):
    if not FAP:
        doplot(i, i + 2, g._subplot_number(i // 3))
    else:
        doplot(i + 1, i + 2, g._subplot_number(i // 3))

        # g.export(tag='z%s' % (i // 2 + 1))
# g.newPlot()
g.settings.legend_frac_subplot_margin = 0.15
g.finish_plot(legend_labels=['SDSS DR12', s.planckTT, s.lensingall], legend_ncol=3)

g.export()
