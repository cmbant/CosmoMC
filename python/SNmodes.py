import chains
import numpy as np
from pylab import *

rootdir = r'C:\tmp\Planck\SNmodes'

lmax = 800
lmin = 800

params = ['theta', 'omegabh2', 'omegamh2', 'ns', 'tau', 'clamp']
# params = ['tau', 'logA', 'ns', 'omegabh2', 'omegach2', 'thetastar']

highL = chains.loadGridChain(rootdir, 'base', 'v97_TT_tau07_alllmin' + str(lmax + 1), ignore_frac=0.3)
lowl = chains.loadGridChain(rootdir, 'base', 'v97_TT_tau07_lowl_alllmax' + str(lmin), ignore_frac=0.3)

allL = chains.loadGridChain(rootdir, 'base', 'v97_TT_tau07_lowl', ignore_frac=0.3)
pol = chains.loadGridChain(rootdir, 'base', 'v97_all_tau07_lowl', ignore_frac=0.3)

allcov = allL.cov(params)

noise = pol.cov(params)
if False:
    i = params.index('tau')
    N = inv(noise)
    N[i, i] -= 1 / 0.02 ** 2
    noise = inv(N)

w, U = allL.getSignalToNoise(params, noise)
print params

for i in range(len(params)):
    U[i, :] *= np.sqrt(np.diagonal(allcov))
    print w[i], U[i, :]


print 'SVD'
for i in range(len(params)):
    U[i, :] *= sqrt(w[i])

X = np.dot(U.T, U)
w2, U2 = np.linalg.eigh(X)
print params

for i in range(len(params)):
    print w2[i], U2[i, :]


# u, s, v = np.linalg.svd(U[:, :], full_matrices=False)
# print s
# print params
# print v

print 'combined..'
U1 = U[-1, :]
U2 = U[-2, :]
U2p = U2 - U1 * U2[0] / U1[0]
print U2p




