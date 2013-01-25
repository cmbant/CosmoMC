# sample settings for a particular grid run

# sets of parameters to vary in addition to baseline
extparams = [[], ['omegak'], ['mnu'], ['nrun', 'r'], ['r'], ['nnu'], ['nrun'], ['Alens'], ['w'], ['yhe'], ['alpha1']]

# dataset names
planck = 'planck_CAMspec'
lowl = 'lowl'
lowLike = 'lowLike'
lensing = 'lensing'
highL = 'highL'
WMAP = 'WMAP'
BAO = 'BAO'
HST = 'HST'
SNLS = 'SNLS'

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
importanceRuns.append([SNLS, ['SNLS.ini']])

# importanceRuns.append(['acc', ['accuracy.ini'], importanceFilterAcc()])
# importanceRuns.append(['v61N', ['v61N.ini'], importanceFilterAcc()])
# importanceRuns.append(['lensing_acc', ['lensing_acc.ini'], importanceFilterAcc()])


skip = []
skip.append('WMAP_lensing')

# if covmats are unreliable, so start learning ASAP
newCovmat = True

start_at_bestfit = False

# try to match run to exisitng covmat
covrenames = dict()
covrenames['_BAO'] = '_post_BAO'
covrenames['_HST'] = '_post_HST'
covrenames['_lowl'] = '_lowl_lowLike'
covrenames['_lowl_lowLike_highL'] = '_lowl_lowLike'
covrenames['_lensing'] = ''
covrenames['_alpha1'] = ''
covrenames['_r'] = ''

# ini files you want to base each set of runs on
defaults = ['common_batch1.ini']

importanceDefaults = ['importance_sampling.ini']
