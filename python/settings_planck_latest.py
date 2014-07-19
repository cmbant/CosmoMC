# provisional updated (reduced) grid settings for full mission
# New BBN fitting built in
# New BAO -> DR11
import batchJob, copy, re

ini_dir = 'batch2/'

defaults = ['common.ini']

importanceDefaults = ['importance_sampling.ini']

# dataset names
planck = 'v97CS'
lowl = 'lowl'
lowLike = 'lowLike'
lensing = 'lensing'
highL = 'highL'
WMAP = 'WMAP'
BAO = 'BAO'
HST = 'HST70p6'
JLA = 'JLA'

BAOdata = 'BAODR11.ini'


Camspec = 'CAMspec_defaults.ini'
highL = 'highL'
lowl = 'lowl'
# dataset names
tauprior = {'prior[tau]':'0.07 0.02'}
tauname = 'tau07'
WMAPtau = {'prior[tau]':'0.09 0.013'}

varTE = {'param[calTE]': '1 0.1 2 0.005 0.005'}
varEE = {'param[calEE]': '1 0.1 2 0.01 0.01'}


TT = {'want_spec':'T T T T F F'}
EE = {'want_spec':'F F F F F T'}
TE = {'want_spec':'F F F F T F'}
TEEE = {'want_spec':'F F F F T T'}
TTTE = {'want_spec':'T T T T T F'}
full = {'want_spec':'T T T T T T'}

freecal = 'freecal.ini'
freecalEE = {'param[calEE]':'1 0.1 2 0.01 0.01', 'prior[calEE]':'1 1'}
freecalTE = {'param[calTE]':'1 0.1 2 0.005 0.005', 'prior[calTE]': '1 1'}

planck_detsets = [freecal, 'nonclik_v97F.ini', Camspec]
planck_CS = [freecal, 'nonclik_v97CS.ini', Camspec]


detsets = []
CS = []
for name, datasets, planck_vars in zip(['v97CS'], [CS], [planck_CS]):
    datasets.append(batchJob.dataSet([name , 'TT'], [TT] + planck_vars))
#    datasets.append(batchJob.dataSet([name , 'TE'], [TE, varTE, freecalTE] + planck_vars))
#    datasets.append(batchJob.dataSet([name , 'EE'], [EE, varEE, freecalEE] + planck_vars))
    datasets.append(batchJob.dataSet([name, 'all'], [full, varTE, varEE, freecalTE, freecalEE] + planck_vars))


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
post_HST = [[HST], ['HST_GPE70p6.ini'], importanceFilterNotOmegakLowl()]
post_JLA = [[JLA], ['JLA_marge.ini'], importanceFilterNotOmegakLowl()]
# set up groups of parameters and data sets

groups = []

g = batchJob.jobGroup('main')
# Main group with just tau prior

g.datasets = copy.deepcopy(CS)
for d in g.datasets:
    d.add(tauname, tauprior)
    d.add(lowl)

g.params = [[], ['omegak'], ['mnu'], ['r'], ['nnu'], ['nrun'], ['Alens'], ['yhe']]
g.importanceRuns = [post_BAO, post_JLA, post_lensing, post_HST]
groups.append(g)


g2 = batchJob.jobGroup('ext')
g2.datasets = copy.deepcopy(g.datasets)
g2.params = [ ['nnu', 'meffsterile'], ['nnu', 'mnu'], ['nnu', 'yhe']]
g2.importanceRuns = [post_BAO, post_JLA, post_lensing, post_HST]
groups.append(g2)

g3 = batchJob.jobGroup('geom')
g3.params = [['omegak']]
g3.datasets = []
for d in copy.deepcopy(g.datasets):
    d.add(BAO, BAOdata)
    g3.datasets.append(d)
for d in copy.deepcopy(g.datasets):
    d.add(JLA)
    g3.datasets.append(d)
for d in copy.deepcopy(g.datasets):
    d.add(lensing)
    g3.datasets.append(d)

g3.importanceRuns = [post_BAO , post_JLA]
groups.append(g3)

g4 = batchJob.jobGroup('planckonlynolow')
g4.params = [[]]
g4.datasets = copy.deepcopy(CS)
g4.importanceRuns = [post_BAO, post_JLA, post_lensing]
groups.append(g4)

g5 = batchJob.jobGroup('planckonly')
g5.params = [[]]
g5.datasets = copy.deepcopy(CS)
for d in g5.datasets:
    d.add(lowl)
g5.importanceRuns = [post_BAO, post_JLA, post_lensing]
groups.append(g5)



skip = []

# try to match run to exisitng covmat
covrenames = []
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

