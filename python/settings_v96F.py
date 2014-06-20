import re, batchJob, copy

# Look at Alens=1 for fixed lensing template amplitude

ini_dir = 'batch2/'

defaults = ['common.ini']

importanceDefaults = ['importance_sampling.ini']

Camspec = 'CAMspec_defaults.ini'
# dataset names
tauprior = {'prior[tau]':'0.06 0.02'}
WMAPtau = {'prior[tau]':'0.09 0.013'}

TT = {'want_spec':'T T T T F F', 'pre_marged':'F'}
EE = {'want_spec':'F F F F F T', 'pre_marged':'F'}
TE = {'want_spec':'F F F F T F', 'pre_marged':'F'}
TEEE = {'want_spec':'F F F F T T', 'pre_marged':'F'}
TTTE = {'want_spec':'T T T T T F', 'pre_marged':'F'}

full = {'want_spec':'T T T T T T', 'pre_marged':'T'}

freecal = {'prior[cal0]':' 1 0.05', 'prior[cal2]':'1 0.05'}

planck_vars = ['pico.ini', 'nonclik_v96F.ini', Camspec, {'indep_sample':0} ]

wig1800_217 = {'param[wig2_217]':'0 -50 50 3 3'}
wig1800_143 = {'param[wig2_143]':'0 -50 50 3 3'}
wig1460_217 = {'param[wig1_217]':'0 -50 50 3 3'}

datasets = []

datasets.append(batchJob.dataSet('v96', [full] + planck_vars))
datasets.append(batchJob.dataSet('v96TE', [TE] + planck_vars))
datasets.append(batchJob.dataSet('v96EE', [EE] + planck_vars))
datasets.append(batchJob.dataSet('v96TT', [TT] + planck_vars))
datasets.append(batchJob.dataSet('v96TEEE', [TEEE] + planck_vars))
datasets.append(batchJob.dataSet('v96TTTE', [TTTE] + planck_vars))


covmat = 'planck_covmats/base_planck_lowl_lowLike.covmat'

start_at_bestfit = False
newCovmats = False

groups = []

g = batchJob.jobGroup('main')
g.datasets = []
for d in copy.deepcopy(datasets):
    d.add(None, tauprior)
    g.datasets.append(d)

for d in copy.deepcopy(datasets):
    d.add(None, tauprior)
    d.add('lowl', ['lowl.ini'])
    g.datasets.append(d)

g.params = [[], ['Alens']]

groups.append(g)


g = batchJob.jobGroup('commlowl')
g.datasets = []

for d in copy.deepcopy(datasets):
    d.add('commlowl', ['commlowl.ini'])
    g.datasets.append(d)

g.params = [[]]
groups.append(g)



g = batchJob.jobGroup('freecal')
g.datasets = copy.deepcopy(datasets)
for d in g.datasets:
    d.add(None, tauprior)
    d.add('freecal', freecal)

groups.append(g)




g = batchJob.jobGroup('WMAPtau')
g.datasets = copy.deepcopy(datasets)
for d in g.datasets:
    d.add('WMAPtau', WMAPtau)
g.params = [[], ['Alens']]

#    d.addFirst('freecal', freecal)
# sets of parameters to vary in addition to baseline

groups.append(g)



