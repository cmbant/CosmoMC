# sample settings for a particular grid run

# sets of parameters to vary in addition to baseline
extparams = [[], ['omegak'], ['mnu'], ['nrun', 'r'], ['r'], ['nnu'], ['nrun'], ['Alens'], ['alpha1']]

# dataset names
planck = 'planck_CAMspec'
lowl = 'lowl'
lowLike = 'lowLike'
lensing = 'lensing'
highL = 'highL'
WMAP = 'WMAP'
BAO = 'BAO'
HST = 'HST'

datasets = []
# lists of dataset names to combine, with corresponding sets of inis to include
datasets.append([[planck, lowl, lowLike], ['CAMspec_defaults.ini', 'lowl.ini', 'lowLike.ini']])
datasets.append([[planck, lowl], ['CAMspec_defaults.ini', 'lowl.ini']])
datasets.append([[planck, lowl, lowLike, highL], ['CAMspec_ACTSPT_defaults.ini', 'lowl.ini', 'lowLike.ini']])
datasets.append([[WMAP], ['WMAP.ini']])

class importanceFilterPlanck:
    def wantImportance(self, jobItem):
        return planck in jobItem.data_set[0]

class importanceFilterAcc:
    def wantImportance(self, jobItem):
        return lowLike in jobItem.data_set[0]

# add importance name tags, and list of specific .ini files to include (in batch1/)
importanceRuns = []
importanceRuns.append([lensing, ['lensing.ini'], importanceFilterPlanck()])
importanceRuns.append([BAO, ['BAO.ini']])
importanceRuns.append([HST, ['HST.ini']])
# importanceRuns.append(['acc', ['accuracy.ini'], importanceFilterAcc()])
# importanceRuns.append(['v61N', ['v61N.ini'], importanceFilterAcc()])
# importanceRuns.append(['lensing_acc', ['lensing_acc.ini'], importanceFilterAcc()])

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
skip.append('WMAP_lensing')

# if covmats are unreliable, so start learning ASAP
newCovmat = True

start_at_bestfit = True

# try to match run to exisitng covmat
covrenames = dict()
covrenames['_BAO'] = '_post_BAO'
covrenames['_HST'] = '_post_HST'
covrenames['_lowl'] = '_lowl_lowLike'
covrenames['_lowl_lowLike'] = ''
covrenames['planck_CAMspec_'] = 'planck_CAMspec_lowl_'
covrenames['_lensing'] = ''
covrenames['_alpha1'] = ''
covrenames['_r'] = ''

# ini files you want to base each set of runs on
defaults = ['common_batch1.ini']

importanceDefaults = ['importance_sampling.ini']
