import planckStyle as s
from paramgrid import batchjob
import GetDistPlots
from pylab import *

g = s.getSinglePlotter()

alpha_npoints = 280
rd_fid = 149.28
H_fid = 93.558
DA_fid = 1359.72

DAbar = 1421
DAerr = 20
Hbar = 96.8
Herr = 3.4
corr = 0.539


def BAOdensityG():
    cov = np.zeros((2, 2))
    cov[0, 0] = DAerr ** 2
    cov[1, 1] = Herr ** 2
    cov[0, 1] = corr * DAerr * Herr
    cov[1, 0] = cov[0, 1]
    invcov = inv(cov)

    Hv = np.arange(85, 110, 0.3)
    DAv = np.arange(1000, 1500, 4)

    H, DA = np.meshgrid(Hv, DAv)
    like = (DA - DAbar) ** 2 * invcov[0, 0] + 2 * (DA - DAbar) * (H - Hbar) * invcov[0, 1] + (H - Hbar) ** 2 * invcov[
        1, 1]

    density = GetDistPlots.Density2D(Hv, DAv, exp(-like / 2))
    density.contours = exp(-np.array([1.509, 2.4477]) ** 2 / 2)
    return density


def BAOdensity(prob_file):
    d = loadtxt(prob_file)
    ix = 0
    prob = np.zeros((alpha_npoints, alpha_npoints))
    alpha_perp = np.zeros(alpha_npoints)
    alpha_pl = np.zeros(alpha_npoints)
    for i in range(alpha_npoints):
        for j in range(alpha_npoints):
            alpha_perp[i] = d[ix, 0]
            alpha_pl[j] = d[ix, 1]
            prob[j, i] = d[ix, 2]
            ix += 1
    prob = prob / np.max(prob)
    return alpha_perp, alpha_pl, prob


alpha_perp, alpha_pl, prob = BAOdensity(batchjob.getCodeRootPath() + 'data/sdss_DR11CMASS_consensus.dat')

densityG = BAOdensityG()

perp = alpha_perp * DA_fid
para = H_fid / alpha_pl

density = GetDistPlots.Density2D(perp, para, prob)
density.contours = exp(-np.array([1.509, 2.4477]) ** 2 / 2)

root = 'base_plikHM_TT_lowTEB_lensing'

c = 29979.2458


def makeNew(samples):
    p = samples.getParams()
    rsH = p.Hubble057 * p.rdrag / rd_fid
    rsHpar = samples.addDerived(rsH, name='rsH',
                                label=r'H(0.57) (r_{\mathrm{drag}}/r_{\mathrm{drag}}^{\rm fid})\, [{\rm km} \,{\rm s}^{-1}{\rm Mpc}^{-1}]')
    Da = p.DA057 * rd_fid / p.rdrag
    Dapar = samples.addDerived(Da, name='Da',
                               label=r'D_A(0.57) (r_{\mathrm{drag}}^{\rm fid}/r_{\mathrm{drag}})\,[\rm{Mpc}]')
    # comb = (Da / 1410) * (rsH / 91.7) ** 1.72
    Hphys = rsH * rd_fid / c
    Daphys = Da / rd_fid
    comb = (Daphys / 9.385) * (Hphys / 0.4582) ** 1.7

    samples.updateBaseStatistics()
    print 'Da ', samples.mean(Da), samples.std(Da)
    print 'rsH ', samples.mean(rsH), samples.std(rsH)
    print 'Daphys ', samples.mean(Daphys), samples.std(Daphys)
    print 'Hphys ', samples.mean(Hphys), samples.std(Hphys)

    print samples.mean(comb), samples.std(comb)
    print samples.PCA(['Da', 'rsH'], 'LLL', 'Da')

    return rsHpar, Dapar


samples = g.sampleAnalyser.samplesForRoot(root)
rsH, Da = makeNew(samples)

# print samples.get2DContourLevels(prob)
# print density.contours

g.add_2d_contours(root, 'Da', 'rsH', filled=True, density=density)
# g.add_2d_contours(root, 'Da', 'rsH', filled=False, density=densityG)
g.settings.scatter_size = 2
g.add_3d_scatter(root, ['Da', 'rsH', 'omegach2'], alpha=0.07, extra_thin=1)
g.add_legend(['BOSS CMASS'], legend_loc='upper left')
g.setAxes([Da, rsH], lims=[1350, 1500, 85, 110])

g.export()

