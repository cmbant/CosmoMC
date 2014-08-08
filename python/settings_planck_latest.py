# provisional updated (reduced) grid settings for full mission
# New BBN fitting built in
# New BAO -> DR11
import batchJob, copy, re

ini_dir = 'batch2/'

defaults = ['common.ini']

importanceDefaults = ['importance_sampling.ini']

# dataset names
lowl = 'lowl'
lowLike = 'lowLike'
lensing = 'lensing'
highL = 'highL'
WMAP = 'WMAP'
BAO = 'BAO'
HST = 'HST70p6'
JLA = 'JLA'

BAOdata = 'BAODR11.ini'
HSTdata = 'HST_GPE70p6.ini'

Camspec = 'CAMspec_defaults.ini'
highL = 'highL'
lowl = 'lowl'
lowTEB = 'lowTEB'
# dataset names
tauprior = {'prior[tau]':'0.07 0.02'}
tauname = 'tau07'
WMAPtau = {'prior[tau]':'0.09 0.013'}


camspec_detsets = ['nonclik_v97F.ini']
camspec_CS = ['nonclik_v97CS.ini']


variant_tag = ['TTTEEE', 'TT']  # ['TTTEEE', 'TT', 'TE', 'EE']
variants = variant_tag

planck_highL_sets = []
planck_vars = ['v97CS']
planck_ini = ['CAMspec_%s.ini']
planck_base = [camspec_CS]
planck_covmats = [None]

if True:
    planck_vars += ['plik']
    planck_ini += ['plik_dx11c_%s_v13.ini']
    planck_base += [[]]
    planck_covmats += ['planck_covmats/plik_dx11c_%s_v13.covmat']
for planck, ini, base, cov in zip(planck_vars, planck_ini, planck_base, planck_covmats):
    for name, var in zip(variant_tag, variants):
        if cov is not None:
            covM = cov % (var)
        else: covM = None
        planck_highL_sets.append(batchJob.dataSet([planck , name], base + [ ini % (var)], covmat=covM))


WMAP9 = [[WMAP], ['WMAP.ini']]

start_at_bestfit = False
newCovmats = False

# Importance sampling settings

class importanceFilterLensing:
    def wantImportance(self, jobItem):
        return planck in jobItem.data_set.names and (not'omegak' in jobItem.param_set or (len(jobItem.param_set) == 1))

class importanceFilterNotOmegakLowl:
    def wantImportance(self, jobItem):
        return not ('omegak' in jobItem.param_set and jobItem.datatag == planck + '_' + lowl)


post_lensing = [[lensing], ['lensing.ini'], importanceFilterLensing()]
post_BAO = [[BAO], [BAOdata], importanceFilterNotOmegakLowl()]
post_HST = [[HST], [HSTdata], importanceFilterNotOmegakLowl()]
post_JLA = [[JLA], ['JLA_marge.ini'], importanceFilterNotOmegakLowl()]
post_nonCMB = [[BAO, HST, JLA], [BAOdata, HSTdata, 'JLA_marge.ini'], importanceFilterNotOmegakLowl()]
post_all = [[lensing, BAO, HST, JLA], [lensing, BAOdata, HSTdata, 'JLA_marge.ini'], importanceFilterNotOmegakLowl()]
post_WP = [[ 'WMAPtau'], [WMAPtau, {'redo_no_new_data':'T'}]]

# set up groups of parameters and data sets

groups = []

g = batchJob.jobGroup('main')
# Main group with just tau prior

g.datasets = copy.deepcopy(planck_highL_sets)
for d in g.datasets:
    d.add(lowTEB)

g.params = [[], ['omegak'], ['mnu'], ['r'], ['nnu'], ['nrun'], ['Alens'], ['yhe']]
g.importanceRuns = [post_BAO, post_JLA, post_lensing, post_HST, post_all]
groups.append(g)


g2 = batchJob.jobGroup('ext')
g2.datasets = copy.deepcopy(g.datasets)
g2.params = [ ['nnu', 'meffsterile'], ['nnu', 'mnu'], ['nnu', 'yhe']]
g2.importanceRuns = [post_BAO, post_JLA, post_HST, post_nonCMB]
groups.append(g2)

g3 = batchJob.jobGroup('geom')
g3.params = [['omegak'], [ 'mnu']]
g3.datasets = []
for d in copy.deepcopy(g.datasets):
    d.add(BAO, BAOdata)
    g3.datasets.append(d)
for d in copy.deepcopy(g.datasets):
    d.add(BAO, BAOdata)
    d.add(HST, HSTdata)
    d.add(JLA)
    g3.datasets.append(d)

g3.importanceRuns = [post_lensing]
groups.append(g3)

# Not needed, just included first time by mistkae..
# g4 = batchJob.jobGroup('planckonlynolow')
# g4.params = [[]]
# g4.datasets = copy.deepcopy(CS)
# g4.importanceRuns = [post_BAO, post_JLA, post_lensing]
# groups.append(g4)

g5 = batchJob.jobGroup('nopoltau')
g5.params = [[]]
g5.datasets = copy.deepcopy(planck_highL_sets)
for d in g5.datasets:
    d.add(lowl)
for d in copy.deepcopy(g5.datasets):
    d.add(lensing)
    g5.datasets.append(d)

g5.importanceRuns = [post_BAO, post_nonCMB, post_WP]
groups.append(g5)

g6 = batchJob.jobGroup('lensing')
g6.datasets = copy.deepcopy(g.datasets)
for d in g6.datasets:
    d.add(lensing)
    d.add(None, {'redo_theory':'F'})

g6.params = [[], ['omegak'], ['mnu'], ['nnu', 'meffsterile'], ['nnu', 'mnu'], ['Alens']]
g6.importanceRuns = []
groups.append(g6)


skip = []

# try to match run to exisitng covmat
covrenames = []
for planck in planck_vars:
    covrenames.append([planck, 'planck'])
covrenames.append(['planck', 'planck_CAMspec'])
covrenames.append(['tauprior', 'lowl_lowLike'])
covrenames.append(['_lensing', '_post_lensing'])
covrenames.append(['_BAO', '_post_BAO'])
covrenames.append(['_HST', '_post_HST'])
covrenames.append(['_JLA', '_post_SNLS'])

def covRenamer(name):
    renamed = re.sub(r'_v.*_highL', '_planck_lowl_lowLike_highL', name, re.I)
    renamed = re.sub(r'_v.*', '_planck_lowl_lowLike', renamed, re.I)
    if renamed == name: return[]
    else: return [renamed]

