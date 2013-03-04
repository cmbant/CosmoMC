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
g.datasets.append([[plik, lowl, lowLike], ['PLik_CAMspec_defaults.ini', 'lowl.ini', 'lowLike.ini']])
g.datasets.append([[plik, lowl, lowLike, highL], ['PLik_CAMspec_ACTSPT_defaults.ini', 'lowl.ini', 'lowLike.ini']])
g.datasets.append([['planck_lmax2000', lowl, lowLike], ['CAMspec_lmax2000.ini', 'lowl.ini', 'lowLike.ini']])
g.datasets.append([['planck_no217auto', lowl, lowLike], ['CAMspec_no217auto.ini', 'lowl.ini', 'lowLike.ini']])


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
covrenames = [['plik', 'planck'], ['planck_lmax2000', 'planck'], ['planck_no217auto', 'planck']]

# ini files you want to base each set of runs on
defaults = ['common_batch1.ini']

importanceDefaults = ['importance_sampling.ini']
