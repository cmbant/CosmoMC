# sample settings for a particular grid run

newCovmats = True
start_at_bestfit = False

# dataset names
planck = 'planck'
plik = 'plik'
lowl = 'lowl'
lowLike = 'lowLike'
lensing = 'lensing'
highL = 'highL'
WMAP = 'WMAP'
BAO = 'BAO'
HST = 'HST'
SNLS = 'SNLS'
Union = 'Union2'

# set up groups of parameters and data sets
class group:pass


g = group()
g.params = [[], ['mnu'], ['nnu'], ['nrun'], ['Alens'], ['yhe']]

g.datasets = []
# lists of dataset names to combine, with corresponding sets of inis to include
# g.datasets.append([[planck, lowl, lowLike], ['CAMspec_nonclik.ini', 'lowl.ini', 'lowLike.ini']])
# g.datasets.append([[planck, lowl, lowLike, highL], ['CAMspec_ACTSPT_nonclik.ini', 'lowl.ini', 'lowLike.ini']])
camSpec = ['CAMspec_nonclik.ini', 'lowl.ini', 'lowLike.ini']
camSpechighL = ['CAMspec_ACTSPT_nonclik.ini', 'lowl.ini', 'lowLike.ini']

g.datasets.append([[plik, lowl, lowLike], ['PLik_CAMspec_defaults.ini', 'lowl.ini', 'lowLike.ini']])
g.datasets.append([[plik, lowl, lowLike, highL], ['PLik_CAMspec_ACTSPT_defaults.ini', 'lowl.ini', 'lowLike.ini']])
g.datasets.append([[planck, 'lmax2000', lowl, lowLike], ['CAMspec_lmax2000.ini'] + camSpec])
g.datasets.append([[planck, 'no217auto', lowl, lowLike], ['CAMspec_no217auto.ini'] + camSpec])
g.datasets.append([[planck, 'lmin1200', lowl, lowLike], ['CAMspec_lmin1200.ini'] + camSpec])
g.datasets.append([[planck, 'v61N', lowl, lowLike], ['CAMspec_v61N.ini'] + camSpec])
g.datasets.append([[planck, 'lmax2000', lowl, lowLike, highL], ['CAMspec_lmax2000.ini'] + camSpechighL])
g.datasets.append([[planck, 'no217auto', lowl, lowLike, highL], ['CAMspec_no217auto.ini'] + camSpechighL])
g.datasets.append([[planck, 'lmin1200', lowl, lowLike, highL], ['CAMspec_lmin1200.ini'] + camSpechighL])

# g.datasets.append([[planck, 'v61N', lowl, lowLike, highL], ['CAMspec_v61N.ini'] + camSpechighL])


class importanceFilterPlanck:
    def wantImportance(self, jobItem):
        return planck in jobItem.data_set[0]

class importanceFilterAcc:
    def wantImportance(self, jobItem):
        return lowLike in jobItem.data_set[0]

# add importance name tags, and list of specific .ini files to include (in batch1/)
g.importanceRuns = []
# g.importanceRuns.append([[lensing], ['lensing.ini'], importanceFilterPlanck()])
# g.importanceRuns.append([[BAO], ['BAO.ini']])

groups = [g]

# try to match run to exisitng covmat
covrenames = [['plik', 'planck'], ['_lmax2000', ''], ['_no217auto', ''], ['_lmin1200', ''], ['_v61N', '']]

ini_dir = 'batch1/'

# ini files you want to base each set of runs on
defaults = ['common_batch1.ini']

importanceDefaults = ['importance_sampling.ini']
