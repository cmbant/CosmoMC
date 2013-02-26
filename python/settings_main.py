# sample settings for a particular grid run


# ini files you want to base each set of runs on
defaults = ['common_batch1.ini']

importanceDefaults = ['importance_sampling.ini']

# dataset names
planck = 'planck'
lowl = 'lowl'
lowLike = 'lowLike'
lensing = 'lensing'
highL = 'highL'
WMAP = 'WMAP'
BAO = 'BAO'
HST = 'HST'
SNLS = 'SNLS'
Union = 'Union2'

Camspec = 'CAMspec_nonclik.ini'
CamspecHighL = 'CAMspec_ACTSPT_nonclik.ini'

planck_lowl_lowLike = [[planck, lowl, lowLike], [Camspec, 'lowl.ini', 'lowLike.ini']]
planck_lowl_lowLike_highL = [[planck, lowl, lowLike, highL], [CamspecHighL, 'lowl.ini', 'lowLike.ini']]
planck_lowl = [[planck, lowl], [Camspec, 'lowl.ini']]
WMAP9 = [[WMAP], ['WMAP.ini']]

planck_lowl_lowLike_BAO = [[planck, lowl, lowLike, BAO], [Camspec, 'lowl.ini', 'lowLike.ini', 'BAO.ini']]
planck_lowl_lowLike_highL_BAO = [[planck, lowl, lowLike, highL, BAO], [CamspecHighL, 'lowl.ini', 'lowLike.ini', 'BAO.ini']]
planck_lowl_lowLike_SNLS = [[planck, lowl, lowLike, SNLS], [Camspec, 'lowl.ini', 'lowLike.ini', 'SNLS.ini']]
planck_lowl_lowLike_Union = [[planck, lowl, lowLike, Union], [Camspec, 'lowl.ini', 'lowLike.ini', 'Union.ini']]
planck_lowl_lowLike_HST = [[planck, lowl, lowLike, HST], [Camspec, 'lowl.ini', 'lowLike.ini', 'HST.ini']]


start_at_bestfit = False
newCovmats = True

# Importance sampling settings

class importanceFilterPlanck:
    def wantImportance(self, jobItem):
        return planck in jobItem.dataname_set

class importanceFilterNotOmegakLowl:
    def wantImportance(self, jobItem):
        return not ('omegak' in jobItem.param_set and jobItem.datatag == planck + '_' + lowl)


post_lensing = [[lensing], ['lensing.ini'], importanceFilterPlanck()]
post_BAO = [[BAO], ['BAO.ini'], importanceFilterNotOmegakLowl()]
post_HST = [[HST], ['HST.ini'], importanceFilterNotOmegakLowl()]
post_SNLS = [[SNLS], ['SNLS_marge.ini'], importanceFilterNotOmegakLowl()]
post_Union = [[Union], ['Union.ini'], importanceFilterNotOmegakLowl()]

# set up groups of parameters and data sets
class group:pass


g1 = group()
# sets of parameters to vary in addition to baseline
g1.params = [[], ['omegak'], ['mnu'], ['nrun', 'r'], ['r'], ['nnu'], ['nrun'], ['Alens'], ['w'], ['yhe'], ['alpha1']]

# lists of dataset names to combine, with corresponding sets of inis to include
g1.datasets = [planck_lowl_lowLike, planck_lowl_lowLike_highL, planck_lowl, WMAP9]

# add importance name tags, and list of specific .ini files to include (in batch1/)
g1.importanceRuns = [post_lensing , post_BAO, post_HST, post_SNLS, post_Union]
g1.groupName = 'main'


g2 = group()
# lists of dataset names to combine, with corresponding sets of inis to include
g2.params = [['nnu', 'yhe'], ['nnu', 'mnu'], ['nnu', 'meffsterile'], ['mnu', 'omegak'], ['mnu', 'Alens'], ['nrun', 'r', 'omegak']]
g2.datasets = [planck_lowl_lowLike, planck_lowl_lowLike_highL]
g2.importanceRuns = [post_lensing , post_BAO, post_HST]
g2.groupName = 'ext'


g3 = group()
g3.params = [['omegak'], ['nnu'], ['w'], ['w', 'wa'], ['nrun', 'r', 'omegak'], ['mnu', 'omegak']]
g3.datasets = [planck_lowl_lowLike_BAO, planck_lowl_lowLike_highL_BAO, planck_lowl_lowLike_SNLS, planck_lowl_lowLike_Union, planck_lowl_lowLike_HST]
g3.importanceRuns = [post_lensing , post_BAO, post_HST]
g3.groupName = 'geom'

groups = [g1, g2, g3]

# try to match run to exisitng covmat
covrenames = dict()
covrenames['planck'] = 'planck_CAMspec'
covrenames['_BAO'] = '_post_BAO'
covrenames['_HST'] = '_post_HST'
covrenames['_lowl'] = '_lowl_lowLike'
covrenames['_lowl_lowLike_highL'] = '_lowl_lowLike'
covrenames['_lensing'] = ''
covrenames['_alpha1'] = ''
covrenames['_r'] = ''
covrenames['lowl_lowLike_highL'] = 'lowl_lowLike'
covrenames['lowl_BAO'] = 'lowl_lowLike_BAO'
covrenames['SNLS'] = 'BAO'
covrenames['Union2'] = 'SNLS'
covrenames['Union2'] = 'BAO'
covrenames['HST'] = 'BAO'
covrenames['w_wa_'] = 'w_'
covrenames['lowl_lowLike_highL_BAO'] = 'lowl_lowLike'
covrenames['mnu_Alens'] = 'mnu'
covrenames['_nrun_r'] = ''
