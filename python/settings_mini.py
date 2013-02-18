# sample settings for a particular grid run

# sets of parameters to vary in addition to baseline
extparams = [[], ['r'], ['nnu'], ['nrun'], ['Alens'], ['yhe']]


newCovmats = True

# dataset names
planck = 'planck_CAMspec_lmax1000'
lowl = 'lowl49'
lowLike = 'lowLike'
lensing = 'lensing'
highL = 'highL'
WMAP = 'WMAP'
BAO = 'BAO'
HST = 'HST'
SNLS = 'SNLS'
Union = 'Union2'

datasets = []
# lists of dataset names to combine, with corresponding sets of inis to include
datasets.append([[planck, lowl, lowLike], ['CAMspec_lmax1000_defaults.ini', 'lowl49.ini', 'lowLike.ini']])
datasets.append([[planck, lowl], ['CAMspec_lmax_1000_defaults.ini', 'lowl49.ini']])
datasets.append([[planck, lowl, lowLike, highL], ['CAMspec_lmax1000_ACTSPT_defaults.ini', 'lowl49.ini', 'lowLike.ini']])

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

skip = []


start_at_bestfit = False

# try to match run to exisitng covmat
covrenames = dict()
covrenames['_lowl49'] = '_lowl'
covrenames['_lmax1000'] = ''
covrenames['_lmax1000_lowl49'] = '_lowl'

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
