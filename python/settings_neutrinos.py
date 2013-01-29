# settings for grid of 2-extra-parameter neutrino runs

extparams = [['nnu', 'yhe'], ['nnu', 'mnu'], ['nnu', 'meffsterile'], ['mnu', 'omegak'], ['nrun', 'r', 'omegak']]

# dataset names
planck = 'planck_CAMspec'
highL = 'highL'
BAO = 'BAO'
lowl = 'lowl'
lowLike = 'lowLike'
lensing = 'lensing'
HST = 'HST'


datasets = []
# lists of dataset names to combine, with corresponding sets of inis to include
datasets.append([[planck, lowl, lowLike], ['CAMspec_defaults.ini', 'lowl.ini', 'lowLike.ini']])
datasets.append([[planck, lowl, lowLike, highL], ['CAMspec_ACTSPT_defaults.ini', 'lowl.ini', 'lowLike.ini']])

importanceRuns = []
importanceRuns.append([lensing, ['lensing.ini']])
importanceRuns.append([BAO, ['BAO.ini']])
importanceRuns.append([HST, ['HST.ini']])


skip = []

start_at_bestfit = False

# ini files you want to base each set of runs on
defaults = ['common_batch1.ini']

importanceDefaults = ['importance_sampling.ini']


covrenames = dict()
covrenames['lowl_lowLike_highL'] = 'lowl_lowLike'
covrenames['lowl_lowLike_highL_BAO'] = 'lowl_lowLike'
covrenames['_nrun_r'] = ''
