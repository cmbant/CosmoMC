import planckStyle as s
from paramgrid import batchjob
import GetDistPlots
import pylab as plt
import numpy as np

g = s.getSubplotPlotter(subplot_size=4)

rd_fid = 147.78

cov = np.loadtxt(batchjob.getCodeRootPath() + 'data/DR12/BAO_consensus_covtot_dM_Hz.txt')

pts = np.loadtxt(batchjob.getCodeRootPath() + 'data/DR12/sdss_DR12Consensus_bao.dat', usecols=[0, 1])

pnames = ['DM038', 'Hubble038', 'DM051', 'Hubble051', 'DM061', 'Hubble061']
redshifts = pts[:, 0]
data = pts[:, 1]
planckmeans = []


def BAOdensity(p1, p2, marge=True):
    err = np.sqrt(cov[p1, p1])
    DAv = np.arange(data[p1] - 4 * err, data[p1] + 4 * err, 4)
    err = np.sqrt(cov[p2, p2])
    Hv = np.arange(data[p2] - 4 * err, data[p2] + 4 * err, 0.3)
    DA, H = np.meshgrid(DAv, Hv)
    v1 = data[p1]
    v2 = data[p2]
    if marge:
        mcov = cov[np.ix_([p1, p2], [p1, p2])]
        invcov = np.linalg.inv(mcov)
        like = (DA - v1) ** 2 * invcov[0, 0] + 2 * (DA - v1) * (H - v2) * invcov[0, 1] + (H - v2) ** 2 * invcov[
            1, 1]
    else:
        # Not done yet
        means = np.array(planckmeans)
        invcov = np.linalg.inv(cov)

        fixedix = [p for p in range(0, 6) if not p in [p1, p2]]
        wantix = [p1, p2]
        deltamean = cov[np.ix_(wantix, wantix)].dot(invcov[np.ix_(wantix, fixedix)]).dot(means[fixedix] - data[fixedix])
        print('Data, conditional correction', [v1, v2], deltamean)
        v1 -= deltamean[0]
        v2 -= deltamean[1]
        invcov = invcov[np.ix_(wantix, wantix)]
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
    for p1, p2 in zip(range(0, 6, 2), range(1, 6, 2)):
        DM = getattr(p, pnames[p1])
        H = getattr(p, pnames[p2])
        Da = DM * rd_fid / p.rdrag
        DMpar = samples.addDerived(Da, name=pnames[p1] + 'n',
                                   label=r'D_{\rm M}(%.2f) (r_{\mathrm{drag}}^{\rm fid}/r_{\mathrm{drag}})\,[\rm{Mpc}]' %
                                         redshifts[p1])
        rsH = H * p.rdrag / rd_fid
        rsHpar = samples.addDerived(rsH, name=pnames[p2] + 'n',
                                    label=r'H(%.2f) (r_{\mathrm{drag}}/r_{\mathrm{drag}}^{\rm fid})\, [{\rm km} \,{\rm s}^{-1}{\rm Mpc}^{-1}]' % (
                                        redshifts[p2]))
        dV = pow(DM ** 2 / (H/c) * redshifts[p2], 1. / 3) / p.rdrag
        samples.addDerived(dV, name='DV'+ pnames[p2][-3:],
                                    label=r'dV/rd' )
        pars += [DMpar, rsHpar]
        #
    samples.updateBaseStatistics()
    return samples, pars


samples, pars = makeNewSamples(root)
samplesT, parsTT = makeNewSamples(rootTT)
samplesL, parsL = makeNewSamples(rootL)

for p1, p2, z in zip(pars[::2], pars[1::2], redshifts[::2]):
    print('DM %s' % z, samples.mean(p1), samples.std(p1))
    print('rsH %s' % z, samples.mean(p2), samples.std(p2))
    planckmeans += [samples.mean(p1), samples.mean(p2)]

    # print samples.PCA([p1.name, p2.name], 'LLL', p1.name)


def doplot(ix1, ix2, ax=None):
    density = BAOdensity(ix1, ix2, marge=True)
    g.add_2d_contours(root, pars[ix1], pars[ix2], filled=True, density=density, ax=ax)
    #    density = BAOdensity(ix1, ix2, marge=False)
    #    g.add_2d_contours(root, pars[ix1], pars[ix2], filled=False, density=density, ls='--', ax=ax)

    g.settings.scatter_size = 2
    #    g.add_3d_scatter(root, [pars[ix1], pars[ix2], 'omegach2'], alpha=0.07, extra_thin=1, ax=ax, color_bar = ix1==4)

    #    g.add_2d_scatter(rootTT, parsTT[ix1], parsTT[ix2], alpha=0.07, extra_thin=1, ax=ax, color ='g')

    #    g.add_2d_scatter(rootL, parsL[ix1], parsL[ix2], alpha=0.07, extra_thin=1, ax=ax, color ='r')
    g.add_2d_scatter(rootTT, parsTT[ix1], parsTT[ix2], alpha=0.07, extra_thin=1, ax=ax, color='g')

    #    g.add_2d_scatter(root, pars[ix1], pars[ix2], alpha=0.07, extra_thin=1, ax=ax, color =g.settings.solid_colors[1][1])
    g.add_2d_scatter(rootL, parsL[ix1], parsL[ix2], alpha=0.07, extra_thin=1, ax=ax, color='r')

    g.add_legend([r'DR12 ($z_{\rm eff}=%s$)' % redshifts[ix1]], legend_loc='upper left', ax=ax)
    # g.setAxes([pars[ix1], pars[ix2]], lims=[1350, 1500, 85, 110])
    g.setAxes([pars[ix1], pars[ix2]], lims=[density.x[0], density.x[-1], density.y[0], density.y[-1]], ax=ax)


g.make_figure(nx=3)
for i in range(0, 6, 2):
    doplot(i, i + 1, g._subplot_number(i // 2))
    # g.export(tag='z%s' % (i // 2 + 1))
#    g.newPlot()
g.finish_plot()
g.export(tag='TT-vs-TTTEEElens')
