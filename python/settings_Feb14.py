# provisional updated (reduced) grid settings for full mission
# New BBN fitting built in
# New BAO -> DR11
import re, batchJob

ini_dir = 'batch2/'

defaults = ['common.ini']

importanceDefaults = ['importance_sampling.ini']

# dataset names
planck_vars = ['nonclik_v65CS', 'nonclik_v85F']
planck_pol = ['nonclik_pol']

cuts = ['CAMspec_lmax2000', 'CAMspec_217lmax2000', 'CAMspec_no217auto']

lowl = 'lowl'
lowLike = 'lowLike'
lensing = 'lensing'
highL = 'highL'
WMAP = 'WMAP'
BAO = 'BAO'
HST = 'HST'
JLA = 'JLA'

BAOdata = 'BAODR11.ini'

Camspec = 'CAMspec_defaults.ini'

basedata = [Camspec, lowl, lowLike]

planck_sets = [ [[planck, lowl, lowLike], [planck] + basedata] for planck in planck_vars]
planck_sets += [batchJob.dataSet([planck, lowl, lowLike]) for planck in planck_pol]


basedata = ['nonclik_v85F'] + basedata
planck_cuts = [ [[cut, lowl, lowLike], [cut] + basedata] for cut in cuts]


start_at_bestfit = False
newCovmats = True

# Importance sampling settings

class importanceFilterLensing:
    def wantImportance(self, jobItem):
        return planck in jobItem.data_set.names and (not'omegak' in jobItem.param_set or (len(jobItem.param_set) == 1))

class importanceFilterNotOmegakLowl:
    def wantImportance(self, jobItem):
        return not ('omegak' in jobItem.param_set and jobItem.datatag == planck + '_' + lowl)


post_lensing = [[lensing], ['lensing.ini'], importanceFilterLensing()]
post_BAO = [[BAO], [BAOdata], importanceFilterNotOmegakLowl()]
post_HST = [[HST], ['HST.ini'], importanceFilterNotOmegakLowl()]
post_SN = [[JLA], ['JLA_marge.ini'], importanceFilterNotOmegakLowl()]
# set up groups of parameters and data sets
class group:pass

groups = []

g = group()
# sets of parameters to vary in addition to baseline
g.params = [[], ['nnu'], ['Alens']]

# lists of dataset names to combine, with corresponding sets of inis to include
g.datasets = planck_sets

# add importance name tags, and list of specific .ini files to include (in batch1/)
g.importanceRuns = [post_BAO, post_SN]
g.groupName = 'main'
groups.append(g)

g = group()
# sets of parameters to vary in addition to baseline
g.params = [[], ['Alens']]

# lists of dataset names to combine, with corresponding sets of inis to include
g.datasets = planck_cuts

# add importance name tags, and list of specific .ini files to include (in batch1/)
g.importanceRuns = []
g.groupName = 'cuts'
groups.append(g)



skip = ['base_nnu_meffsterile_planck_lowl_lowLike']



def covRenamer(name):
    renamed = re.sub(r'nonclik_([^_]*)', 'planck', name, re.I)
    renamed = re.sub(r'CAMspec_([^_]*)', 'planck', renamed, re.I)
    if renamed == name: return[]
    else: return [renamed]


# try to match run to exisitng covmat
covrenames = []
covrenames.append([planck, 'planck'])
covrenames.append(['planck', 'planck_CAMspec'])
covrenames.append(['tauprior', 'lowl_lowLike'])
covrenames.append(['_lensing', '_post_lensing'])
covrenames.append(['_BAO', '_post_BAO'])
covrenames.append(['_HST', '_post_HST'])
covrenames.append(['_SNLS', '_post_SNLS'])

