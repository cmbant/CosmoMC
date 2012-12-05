# settings for grid of planck+BAO runs where importance sampling does not work well

extparams = [['omegak'], ['nnu']]

# dataset names
planck = 'planck_CAMspec'
highL = 'highL'
BAO = 'BAO'
lowl = 'lowl'
lowLike = 'lowLike'
lensing = 'lensing'


datasets = []
# lists of dataset names to combine, with corresponding sets of inis to include
datasets.append([[planck, lowl, lowLike, BAO], ['CAMspec_defaults.ini', 'lowl.ini', 'lowLike.ini', 'BAO.ini']])
datasets.append([[planck, lowl, lowLike, highL, BAO], ['CAMspec_ACTSPT_defaults.ini', 'lowl.ini', 'lowLike.ini', 'BAO.ini']])

importanceRuns = []

# priors and widths for parameters which are varied
params = dict()
params['mnu'] = '0 0 5 0.1 0.03'
params['omegak'] = '0 -0.3 0.3 0.001 0.001'
params['w'] = '-1 -3 -0.3 0.02 0.02'
params['nnu'] = '3.046 0 10 0.05 0.05'
params['nrun'] = '0 -1 1 0.001 0.001'
params['r'] = '0 0 2 0.03 0.03'
params['Alens'] = '1 0 10 0.05 0.05'
params['yhe'] = '0.245 0.1 0.5 0.006 0.006'
params['alpha1'] = '0 -1 1 0.0003 0.0003'
params['deltazrei'] = '0.5 0.1 3 0.3 0.3'
params['wa'] = '0 -2 2 0.3 0.3'

skip = []
skip.append('base_omegak_mnu_planck_CAMspec_highL')

# if covmats are unreliable, so start learning ASAP
newCovmat = True

start_at_bestfit = True

# ini files you want to base each set of runs on
defaults = ['common_batch1.ini']

importanceDefaults = ['importance_sampling.ini']

covrenames = dict()
covrenames['planck_CAMspec_lowl_lowLike'] = 'planck_CAMspec'
