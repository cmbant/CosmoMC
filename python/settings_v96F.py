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

freecal = {'prior[cal0]':' 1 1', 'prior[cal2]':'1 1'}
freecalpol = {'prior[cal0]':' 1 1', 'prior[cal2]':'1 1', 'param[calTE]':'1 0.1 2 0.005 0.005', 'param[calEE]':'1 0.1 2 0.01 0.01'}


planck_detsets = ['pico.ini', 'nonclik_v96F.ini', Camspec, {'indep_sample':0} ]
planck_CS = ['pico.ini', 'nonclik_v96CS.ini', Camspec, {'indep_sample':0} ]

wig1800_217 = {'param[wig2_217]':'0 -50 50 3 3'}
wig1800_143 = {'param[wig2_143]':'0 -50 50 3 3'}
wig1460_217 = {'param[wig1_217]':'0 -50 50 3 3'}

v96detsets = []
v96CS = []
for name, datasets, planck_vars in zip(['v96', 'v96CS'], [v96detsets, v96CS], [planck_detsets, planck_CS]):
    datasets.append(batchJob.dataSet(name, [full] + planck_vars))
    datasets.append(batchJob.dataSet(name + 'TE', [TE] + planck_vars))
    datasets.append(batchJob.dataSet(name + 'EE', [EE] + planck_vars))
    datasets.append(batchJob.dataSet(name + 'TT', [TT] + planck_vars))
    datasets.append(batchJob.dataSet(name + 'TEEE', [TEEE] + planck_vars))
    datasets.append(batchJob.dataSet(name + 'TTTE', [TTTE] + planck_vars))


covmat = 'planck_covmats/base_planck_lowl_lowLike.covmat'

start_at_bestfit = False
newCovmats = False

groups = []

g = batchJob.jobGroup('main')
g.datasets = []
for d in copy.deepcopy(v96detsets):
    d.add(None, tauprior)
    g.datasets.append(d)

for d in copy.deepcopy(v96detsets):
    d.add(None, tauprior)
    d.add('lowl', ['lowl.ini'])
    g.datasets.append(d)
g.params = [[], ['Alens'], ['nnu'], ['mnu'], ['r'], ['nrun'], ['yhe']]

groups.append(g)

g = batchJob.jobGroup('plik')
g.datasets = []
plik_vars = ['pico.ini', {'indep_sample':0} ]
g.datasets.append(batchJob.dataSet('plikTT', plik_vars + ['plik_dx11c_TT_v11.ini'], covmat='planck_covmats/plik_dx11c_TT_v11.covmat'))
g.datasets.append(batchJob.dataSet('plikTE', plik_vars + ['plik_dx11c_TE_v11.ini'], covmat='planck_covmats/plik_dx11c_TE_v11.covmat'))
g.datasets.append(batchJob.dataSet('plikEE', plik_vars + ['plik_dx11c_EE_v11.ini'], covmat='planck_covmats/plik_dx11c_EE_v11.covmat'))
g.datasets.append(batchJob.dataSet('plikTTTEEE', plik_vars + ['plik_dx11c_TTTEEE_v11.ini'], covmat='planck_covmats/plik_dx11c_TTTEEE_v11.covmat'))
for d in g.datasets:
    d.add(None, tauprior)

groups.append(g)


g = batchJob.jobGroup('commlowl')
g.datasets = []

for d in copy.deepcopy(v96detsets):
    d.add('commlowl', ['commlowl.ini'])
    g.datasets.append(d)

g.params = [[]]
groups.append(g)

g = batchJob.jobGroup('freecal')
g.datasets = copy.deepcopy(v96detsets)
for d in g.datasets:
    d.add(None, tauprior)
    d.add('freecal', freecal)
groups.append(g)

g = batchJob.jobGroup('freecalCS')
g.datasets = copy.deepcopy(v96CS)
for d in g.datasets:
    d.add(None, tauprior)
    d.add('freecal', freecal)
groups.append(g)


g = batchJob.jobGroup('freecalpol')
g.datasets = [batchJob.dataSet('v96', [full] + planck_vars)]
for d in g.datasets:
    d.add(None, tauprior)
    d.add('freecalpol', freecalpol)
groups.append(g)


g = batchJob.jobGroup('WMAPtau')
g.datasets = copy.deepcopy(v96detsets)
for d in g.datasets:
    d.add('WMAPtau', WMAPtau)
g.params = [[], ['Alens']]

groups.append(g)


g = batchJob.jobGroup('old')
g.datasets = []
for name, x in zip(['v62TN', 'v65F', 'v85F'], ['nonclik_v62TN.ini', 'nonclik.ini', 'nonclik_v85F.ini']):
    g.datasets.append(batchJob.dataSet(name, ['pico.ini', x, Camspec, {'indep_sample':0} ]))
g.datasets.append(batchJob.dataSet('v96F', [TT] + planck_vars))

for d in g.datasets:
    d.add('freecal', freecal)


groups.append(g)



