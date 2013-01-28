# settings for grid of planck+BAO runs where importance sampling does not work well

extparams = [['omegak'], ['nnu'], ['mnu', 'omegak'], ['w'], ['w', 'wa']]

# dataset names
planck = 'planck_CAMspec'
highL = 'highL'
BAO = 'BAO'
lowl = 'lowl'
lowLike = 'lowLike'
lensing = 'lensing'
SNLS = 'SNLS'
HST = 'HST'

datasets = []
# lists of dataset names to combine, with corresponding sets of inis to include
datasets.append([[planck, lowl, lowLike, BAO], ['CAMspec_defaults.ini', 'lowl.ini', 'lowLike.ini', 'BAO.ini']])
datasets.append([[planck, lowl, lowLike, highL, BAO], ['CAMspec_ACTSPT_defaults.ini', 'lowl.ini', 'lowLike.ini', 'BAO.ini']])
datasets.append([[planck, lowl, lowLike, SNLS], ['CAMspec_defaults.ini', 'lowl.ini', 'lowLike.ini', 'SNLS.ini']])

importanceRuns = []
importanceRuns.append([lensing, ['lensing.ini']])
importanceRuns.append([BAO, ['BAO.ini']])
importanceRuns.append([HST, ['HST.ini']])
# importanceRuns.append([SNLS, ['SNLS.ini']])

skip = ['HST_HST', 'BAO_BAO', 'SNLS_SNLS']

# if covmats are unreliable, so start learning ASAP
newCovmat = True

start_at_bestfit = False

# ini files you want to base each set of runs on
defaults = ['common_batch1.ini']

importanceDefaults = ['importance_sampling.ini']

covrenames = dict()
covrenames['lowl_lowLike_highL'] = 'lowl_lowLike'
covrenames['lowl_BAO'] = 'lowl_lowLike_BAO'
covrenames['SNLS'] = 'BAO'

covrenames['w_wa_'] = 'w_'
