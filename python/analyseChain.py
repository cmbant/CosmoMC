import chains
import sys, numpy as np

rootdir = r'C://tmp/Planck/final_nominal/PLA/'

chains = chains.loadGridChain(rootdir, 'base', 'planck_lowl_lowLike_highL', 'BAO')

p = chains.getParams()

# chains.reweight_plus_logLike((p.sigma8 - 0.786) ** 2 / 0.031 ** 2 / 2)
p = chains.getParams()
derived = p.H0 * p.age * 0.00102269032

# derived = p.sigma8 * p.omegam ** 0.6
print 'mean, err = ', chains.mean(derived), chains.std(derived)

print '95% limits: ', chains.twoTailLimits(derived, 0.95)
print '95% upper: ', chains.confidence(derived, 0.05, upper=True)

# chain.setParams(sys.modules[__name__])

# H0, s8, Omm, omb, tau, omegach2, omegabh2, omegam = chain.valuesForParam(['H0', 'sigma8', 'omegam', 'omegabh2', 'tau', 'omegach2', 'omegabh2', 'omegam'])
# mnu, nnu = chain.valuesForParam(['meffsterile', 'nnu'])

# theta = chain.valuesForParam(['thetastar'])


# h = H0 / 100

# derived = Omm * h ** 3.2 / omb ** 0.54

# derived = sigma8 * omegam ** 0.6  # * np.exp(-tau)


# print chain.mean(nnu), r'\pm', chain.std(nnu)
# print chain.confidence(nnu, 0.05, upper=True), chain.confidence(mnu, 0.05, upper=True)
# print chain.confidence(mnu, 0.05, upper=False)

