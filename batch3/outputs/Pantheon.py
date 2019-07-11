from __future__ import absolute_import
from __future__ import print_function
import pylab as plt
import sys

sys.path.insert(0, r'c:\work\dist\git\camb')
import camb
from cosmomc_to_camb import get_camb_params
import planckStyle as s
import numpy as np
from planck import SN
import os

g = s.getSinglePlotter()

like = SN.SN_likelihood(os.path.join(os.path.dirname(__file__), r'../../data/Pantheon/full_long18.dataset'))
JLA = SN.SN_likelihood(os.path.join(os.path.dirname(__file__), r'../../data/jla.dataset'), marginalize=False)

common = []
for name in like.names:
    common.append(name in JLA.names or 'SDSS' + name in JLA.names or 'sn' + name in JLA.names)
common = np.array(common, dtype=np.bool)
print(like.nsn, np.sum(common), like.nsn - np.sum(common))

redshifts = np.logspace(-2, 1, 1000)
samples = g.sampleAnalyser.samplesForRoot('base_plikHM_TTTEEE_lowl_lowE_lensing')
ixs = samples.randomSingleSamples_indices()
dists = np.zeros((len(ixs), len(redshifts)))
sndists = np.zeros((len(ixs), like.nsn))
for i, ix in enumerate(ixs):
    dic = samples.getParamSampleDict(ix)
    camb_pars = get_camb_params(dic)
    results = camb.get_background(camb_pars)
    dists[i, :] = 5 * np.log10((1 + redshifts) ** 2 * results.angular_diameter_distance(redshifts))
    sndists[i, :] = 5 * np.log10((1 + like.zcmb) ** 2 * results.angular_diameter_distance(like.zcmb))

paramdic = g.bestfit('base_plikHM_TTTEEE_lowl_lowE_lensing').getParamDict()
camb_pars = get_camb_params(paramdic)
results = camb.get_background(camb_pars)

invvars = 1.0 / like.pre_vars
wtval = np.sum(invvars)

offset = 5 * np.log10(1e-5)
lumdists = 5 * np.log10((1 + like.zcmb) * (1+like.zhel) * results.angular_diameter_distance(like.zcmb))

redshifts = np.logspace(-2, 1, 1000)
d = results.angular_diameter_distance(redshifts)
theory = 5 * np.log10((1 + redshifts) ** 2 * d)
planck_means = np.zeros(redshifts.shape)
planck_err = planck_means.copy()
plotdists = dists.copy()
for i in range(dists.shape[0]):
    # normalize optimally as though SN points
    estimated_scriptm = np.sum((sndists[i, :] - lumdists) * invvars) / wtval
    plotdists[i, :] -= estimated_scriptm

for i in range(len(planck_means)):
    planck_means[i] = np.mean(plotdists[:, i]) - theory[i]
    planck_err[i] = np.std(plotdists[:, i])

estimated_scriptm = np.sum((like.mag - lumdists) * invvars) / wtval
m = like.mag - estimated_scriptm - offset
pred = lumdists - offset

cov = np.linalg.inv(like.invcov)
diagerrs = np.sqrt(np.diag(cov))

##Binned
ix = np.argsort(like.zcmb)

mins = []
maxs = []
Cinvd = like.invcov.dot(m - pred)
zinv = like.invcov.dot(np.log(like.zcmb))

nbin = 20
di = like.nsn // nbin + 1
bins = [ix[i:i + di] for i in range(0, like.nsn, di)]
bins = bins[:-1] + [bins[-1][:-10], bins[-1][-10:]]
nbin = len(bins)
Csmall = np.zeros((nbin, nbin))
x = np.zeros(nbin)
zx = np.zeros(nbin)

for i, ixs in enumerate(bins):
    for j, ixs2 in enumerate(bins):
        Csmall[i, j] = np.sum(like.invcov[np.ix_(ixs, ixs2)])
    x[i] = np.sum(Cinvd[ixs])
    zx[i] = np.sum(zinv[ixs])

smallcov = np.linalg.inv(Csmall)
bandpowers = smallcov.dot(x)
zbands = np.exp(smallcov.dot(zx))
errs = np.sqrt(np.diag(smallcov))

for ixs in bins:
    mins += [np.min(like.zcmb[ixs])]
    maxs += [np.max(like.zcmb[ixs])]
mins = np.array(mins)
maxs = np.array(maxs)

fig, axs = plt.subplots(2, 1, figsize=(5, 4), sharex='col', gridspec_kw={'height_ratios': [3, 2]})
ax = axs[0]
s.plotBands(redshifts, planck_means, planck_err)
ax.errorbar(like.zcmb[common], m[common] - pred[common], diagerrs[common], fmt='.', lw=0.3, markersize=2,
            label='JLA and Pantheon')
ax.set_xscale('log')
common = np.array([not x for x in common])
ax.errorbar(like.zcmb[common], m[common] - pred[common], diagerrs[common], fmt='.', lw=0.3, markersize=2,
            label='Pantheon only')
ax.axhline(0, color='k', lw=1)
ax.set_xlim(0.01, 3)
ax.legend(loc='lower right');
ax.set_yticks([-0.8, -0.4, 0, 0.4, 0.8])
ax.tick_params(axis='both', which='major', labelsize=10)

ax = axs[1]
ax.errorbar(zbands, bandpowers, errs, [zbands - mins, maxs - zbands], fmt='.', marker='o', color='C2', markersize=4)
ax.set_xscale('log')
ax.axhline(0, color='k', lw=1)
ax.set_xlabel('$z$', fontsize=13)
#fig.text(0.01, 0.5, r'$\mu -\mathcal{M} - \mu_{\rm{Planck}}$', va='center', rotation='vertical', fontsize=12)
fig.text(0.01, 0.5, r'$\mu - \mu_{\rm{Planck}}$', va='center', rotation='vertical', fontsize=13)
ax.set_xticks([0.01, 0.1, 1])
ax.set_xticklabels(['0.01', '0.1', '1'])
ax.set_yticks([-0.1, -0.05, 0, 0.05])
ax.tick_params(axis='both', which='major', labelsize=10)
plt.subplots_adjust(hspace=0)

plt.savefig('../../outputs/Pantheon.pdf', bbox_inches='tight')
