import planckStyle as s
import os
from pylab import *

g = s.getSinglePlotter()

g.settings.param_names_for_labels = 'clik_Hunits.paramnames'

outdir = ''
roots = ['base_nnu_meffsterile_' + s.defdata_all + '_post_lensing']

if False:
    samples = g.sampleAnalyser.samplesForRoot('base_nnu_meffsterile_' + s.defdata_all + '_lensing_BAO')
    p = samples.getParams()
    mphys = p.meffsterile / (p.nnu - 3.046) ** 0.75
    samples.filter(mphys < 2)
    samples.updateChainBaseStatistics()
    p = samples.getParams()
    print 'meffsterile < %s'%samples.confidence(p.meffsterile, 0.05, upper=True)
    print 'nnu < %s'%samples.confidence(p.nnu, 0.05, upper=True)

mmax = 1.7
nmax = 3.8
g.plot_3d(roots, ['meffsterile', 'nnu', 'H0'])

m = np.arange(0, mmax, 0.01)
N = np.arange(3.047, nmax, 0.01)
N, m = np.meshgrid(N, m)
z = m / (N - 3.046) ** 0.75
CS = contour(m, N, z, origin='lower', levels=[0.5, 1, 2, 5], colors='k', linestyles='--', linewidths=0.3,
             extent=[0, 3, 3.046, nmax])
clabel(CS, CS.levels, inline=True, fmt='%1.1f', fontsize=7)

z = m / (N - 3.046)
CS = contour(m, N, z, origin='lower', levels=[0.5, 1, 2, 5], colors='gray', linestyles=':', linewidths=0.3,
             extent=[0, 3, 3.046, nmax])

m = np.arange(0, mmax * 1.1, 0.01)
fill_between(m, 0, (m / 10) ** (4. / 3) + 3.046, color='gray')

# g.add_2d_contours('base_nnu_meffsterile_' + s.defdata_all + '_BAO','meffsterile', 'nnu')

xlim([0, mmax])
ylim([3.046, nmax])

g.export()
