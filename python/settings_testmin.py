# sample settings for a particular grid run


params = dict()
params['mnu'] = '0.7 0 5 0.1 0.03'
params['omegak'] = '-0.01 -0.3 0.3 0.001 0.001'  # starting exactly on flat seems to confuse minimizer
params['w'] = '-1.2 -3 -0.3 0.02 0.02'
params['nnu'] = '3.8 0.05 10 0.05 0.05'
params['nrun'] = '-0.01 -1 1 0.001 0.001'
params['r'] = '0.1 0 2 0.03 0.03'
params['Alens'] = '1.2 0 10 0.05 0.05'
params['yhe'] = '0.26 0.1 0.5 0.006 0.006'
params['wa'] = '-0.2 -2 2 0.3 0.3'
params['meffsterile'] = '0.2 0 3 0.1 0.03'

extra_opts = dict()
extra_opts['param[omegabh2]'] = '0.022 0.005 0.1 0.001 0.001'
extra_opts['param[tau]'] = '0.095 0.01 0.8 0.003 0.003'
extra_opts['param[logA]'] = '3.04 2.7 4 0.002 0.002'
extra_opts['param[logA]'] = '3.04 2.7 4 0.002 0.002'
extra_opts['minimization_points_factor'] = 8
# ini files you want to base each set of runs on
defaults = ['common_batch1.ini']

importanceDefaults = ['importance_sampling.ini']

# dataset names
planck = 'testmin'
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
planck_lowl_lowLike_lensing = [[planck, lowl, lowLike, lensing], [Camspec, 'lowl.ini', 'lowLike.ini', 'lensing.ini']]
planck_lowl_lowLike_highL_lensing = [[planck, lowl, lowLike, highL, lensing], [CamspecHighL, 'lowl.ini', 'lowLike.ini', 'lensing.ini']]
planck_tauprior = [[planck, 'tauprior'], [Camspec, 'tauprior.ini']]
planck_tauprior_highL = [[planck, 'tauprior', highL], [CamspecHighL, 'tauprior.ini']]


start_at_bestfit = False
newCovmats = True

# Importance sampling settings

class importanceFilterLensing:
    def wantImportance(self, jobItem):
        return planck in jobItem.dataname_set and (not'omegak' in jobItem.param_set or (len(jobItem.param_set) == 1))

class importanceFilterNotOmegakLowl:
    def wantImportance(self, jobItem):
        return not ('omegak' in jobItem.param_set and jobItem.datatag == planck + '_' + lowl)


post_lensing = [[lensing], ['lensing.ini'], importanceFilterLensing()]
post_BAO = [[BAO], ['BAO.ini'], importanceFilterNotOmegakLowl()]
post_HST = [[HST], ['HST.ini'], importanceFilterNotOmegakLowl()]
post_SNLS = [[SNLS], ['SNLS_marge.ini'], importanceFilterNotOmegakLowl()]
post_Union = [[Union], ['Union.ini'], importanceFilterNotOmegakLowl()]

# set up groups of parameters and data sets
class group:pass

groups = []

g1 = group()
# sets of parameters to vary in addition to baseline
g1.params = [[], ['omegak'], ['mnu'], ['nrun', 'r'], ['r'], ['nnu'], ['nrun'], ['Alens'], ['w'], ['yhe'], ['alpha1']]

# lists of dataset names to combine, with corresponding sets of inis to include
g1.datasets = [planck_lowl_lowLike, planck_lowl_lowLike_highL, planck_lowl]

# add importance name tags, and list of specific .ini files to include (in batch1/)
g1.importanceRuns = [post_lensing , post_BAO, post_HST, post_SNLS, post_Union]
g1.groupName = 'main'
groups.append(g1)

g2 = group()
# lists of dataset names to combine, with corresponding sets of inis to include
g2.params = [['nnu', 'yhe'], ['nnu', 'mnu'], ['nnu', 'meffsterile'], ['mnu', 'omegak'], ['mnu', 'Alens'], ['nrun', 'r', 'omegak']]
g2.datasets = [planck_lowl_lowLike, planck_lowl_lowLike_highL]
g2.importanceRuns = [post_lensing , post_BAO, post_HST]
g2.groupName = 'ext'
groups.append(g2)

g3 = group()
g3.params = [['omegak'], ['nnu'], ['w'], ['w', 'wa'], ['nrun', 'r', 'omegak'], ['mnu', 'omegak']]
g3.datasets = [planck_lowl_lowLike_BAO, planck_lowl_lowLike_highL_BAO, planck_lowl_lowLike_SNLS, planck_lowl_lowLike_Union, planck_lowl_lowLike_HST]
g3.importanceRuns = [post_lensing , post_BAO, post_HST]
g3.groupName = 'geom'
groups.append(g3)


g4 = group()
g4.params = [['mnu'], ['omegak']]
g4.datasets = [planck_lowl_lowLike_lensing, planck_lowl_lowLike_highL_lensing]
g4.importanceRuns = [post_BAO, post_HST, post_SNLS]
g4.groupName = 'lensing'
groups.append(g4)

g = group()
g.params = [[], ['mnu'], ['Alens']]
g.datasets = [planck_tauprior, planck_tauprior_highL]
g.importanceRuns = [post_lensing]
g.groupName = 'tauprior'
groups.append(g)


skip = ['base_nnu_meffsterile_planck_lowl_lowLike']

# try to match run to exisitng covmat
covrenames = []
covrenames.append(['testmin', 'planck'])
covrenames.append(['tauprior', 'lowl_lowLike'])
covrenames.append(['_lensing', '_post_lensing'])
covrenames.append(['_BAO', '_post_BAO'])
covrenames.append(['_HST', '_post_HST'])
covrenames.append(['_SNLS', '_post_SNLS'])
covrenames.append(['_Union2', '_post_Union2'])
