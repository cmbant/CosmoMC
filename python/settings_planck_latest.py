# provisional updated (reduced) grid settings for full mission
# New BBN fitting built in
# New BAO -> DR11

ini_dir = 'batch2/'

defaults = ['common.ini']

importanceDefaults = ['importance_sampling.ini']

# dataset names
planck = 'v65F'
lowl = 'lowl'
lowLike = 'lowLike'
lensing = 'lensing'
highL = 'highL'
WMAP = 'WMAP'
BAO = 'BAO'
HST = 'HST'
SNLS = 'SNLS'

BAOdata = 'BAODR11.ini'

Camspec = 'CAMspec_v5F_nonclik.ini'
CamspecHighL = 'CAMspec_ACTSPT_v5F_nonclik.ini'

planck_lowl_lowLike = [[planck, lowl, lowLike], [Camspec, 'lowl.ini', 'lowLike.ini']]
planck_lowl_lowLike_highL = [[planck, lowl, lowLike, highL], [CamspecHighL, 'lowl.ini', 'lowLike.ini']]
planck_lowl = [[planck, lowl], [Camspec, 'lowl.ini']]
WMAP9 = [[WMAP], ['WMAP.ini']]

planck_lowl_lowLike_BAO = [[planck, lowl, lowLike, BAO], [ Camspec, 'lowl.ini', 'lowLike.ini', BAOdata]]
planck_lowl_lowLike_highL_BAO = [[planck, lowl, lowLike, highL, BAO], [CamspecHighL, 'lowl.ini', 'lowLike.ini', BAOdata]]
planck_lowl_lowLike_SNLS = [[planck, lowl, lowLike, SNLS], [Camspec, 'lowl.ini', 'lowLike.ini', 'SNLS.ini']]
planck_lowl_lowLike_HST = [[planck, lowl, lowLike, HST], [Camspec, 'lowl.ini', 'lowLike.ini', 'HST.ini']]
planck_lowl_lowLike_lensing = [[planck, lowl, lowLike, lensing], [Camspec, 'lowl.ini', 'lowLike.ini', 'lensing.ini']]
planck_lowl_lowLike_highL_lensing = [[planck, lowl, lowLike, highL, lensing], [CamspecHighL, 'lowl.ini', 'lowLike.ini', 'lensing.ini']]
planck_tauprior = [[planck, 'tauprior'], [Camspec, 'tauprior.ini']]
planck_tauprior_highL = [[planck, 'tauprior', highL], [CamspecHighL, 'tauprior.ini']]


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
post_SNLS = [[SNLS], ['SNLS_marge.ini'], importanceFilterNotOmegakLowl()]
# set up groups of parameters and data sets
class group:pass

groups = []

g1 = group()
# sets of parameters to vary in addition to baseline
g1.params = [[], ['omegak'], ['mnu'], ['r'], ['nnu'], ['nrun'], ['Alens'], ['yhe']]

# lists of dataset names to combine, with corresponding sets of inis to include
g1.datasets = [planck_lowl_lowLike, planck_lowl_lowLike_highL]

# add importance name tags, and list of specific .ini files to include (in batch1/)
g1.importanceRuns = [post_BAO, post_HST]
g1.groupName = 'main'
groups.append(g1)

g2 = group()
# lists of dataset names to combine, with corresponding sets of inis to include
g2.params = [['nnu', 'yhe'], ['nnu', 'mnu'], ['mnu', 'Alens']]
# todo: prior on m_phys for ['nnu', 'meffsterile'],
g2.datasets = [planck_lowl_lowLike, planck_lowl_lowLike_highL]
g2.importanceRuns = [post_BAO, post_HST]
g2.groupName = 'ext'
groups.append(g2)

g3 = group()
g3.params = [['omegak']]
g3.datasets = [planck_lowl_lowLike_BAO, planck_lowl_lowLike_highL_BAO]
g3.importanceRuns = [post_lensing , post_HST]
g3.groupName = 'geom'
groups.append(g3)

g4 = group()
g4.params = [[]]
g4.datasets = [planck_lowl]
g4.importanceRuns = [post_BAO, post_HST]
g4.groupName = 'planckonly'
groups.append(g4)


skip = ['base_nnu_meffsterile_planck_lowl_lowLike']

# try to match run to exisitng covmat
covrenames = []
covrenames.append([planck, 'planck'])
covrenames.append(['planck', 'planck_CAMspec'])
covrenames.append(['tauprior', 'lowl_lowLike'])
covrenames.append(['_lensing', '_post_lensing'])
covrenames.append(['_BAO', '_post_BAO'])
covrenames.append(['_HST', '_post_HST'])
covrenames.append(['_SNLS', '_post_SNLS'])

