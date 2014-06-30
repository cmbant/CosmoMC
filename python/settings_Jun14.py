import re, batchJob, copy

# Look at Alens=1 for fixed lensing template amplitude

ini_dir = 'batch2/'

defaults = ['common.ini']

importanceDefaults = ['importance_sampling.ini']

override_defaults = ['pico.ini']
extra_opts = {'indep_sample':0, 'checkpoint':'F'}

Camspec = 'CAMspec_defaults.ini'
highL = 'highL'
lowl = 'lowl'
# dataset names
tauprior = {'prior[tau]':'0.07 0.02'}
tauname = 'tau07'
WMAPtau = {'prior[tau]':'0.09 0.013'}

TT = {'want_spec':'T T T T F F', 'pre_marged':'F'}
EE = {'want_spec':'F F F F F T', 'pre_marged':'F'}
TE = {'want_spec':'F F F F T F', 'pre_marged':'F'}
TEEE = {'want_spec':'F F F F T T', 'pre_marged':'F'}
TTTE = {'want_spec':'T T T T T F', 'pre_marged':'F'}
full = {'want_spec':'T T T T T T', 'pre_marged':'F'}

freecal = {'prior[cal0]':'1 1', 'prior[cal2]':'1 1'}
freecalEE = {'param[calEE]':'1 0.1 2 0.01 0.01', 'prior[calEE]':'1 1'}
freecalTE = {'param[calTE]':'1 0.1 2 0.005 0.005', 'prior[calTE]': '1 1'}


planck_detsets = [freecal, 'nonclik_v96F.ini', Camspec]
planck_CS = [freecal, 'nonclik_v96CS.ini', Camspec]

# not using these checks yet
wig1800_217 = {'param[wig2_217]':'0 -50 50 3 3'}
wig1800_143 = {'param[wig2_143]':'0 -50 50 3 3'}
wig1460_217 = {'param[wig1_217]':'0 -50 50 3 3'}

v96detsets = []
v96CS = []
for name, datasets, planck_vars in zip(['v96', 'v96CS'], [v96detsets, v96CS], [planck_detsets, planck_CS]):
    datasets.append(batchJob.dataSet([name , 'TT'], [TT] + planck_vars))
    datasets.append(batchJob.dataSet([name , 'TE'], [TE, freecalTE] + planck_vars))
    datasets.append(batchJob.dataSet([name , 'EE'], [EE, freecalEE] + planck_vars))
    datasets.append(batchJob.dataSet([name, 'all'], [full, freecalTE, freecalEE] + planck_vars))
#    datasets.append(batchJob.dataSet(name + 'TEEE', [TEEE] + planck_vars))
#    datasets.append(batchJob.dataSet(name + 'TTTE', [TTTE] + planck_vars))

plik = []
plik.append(batchJob.dataSet(['plik', 'TT'], ['plik_dx11c_TT_v11.ini'], covmat='planck_covmats/plik_dx11c_TT_v11.covmat'))
plik.append(batchJob.dataSet(['plik', 'TE'], ['plik_dx11c_TE_v11.ini'], covmat='planck_covmats/plik_dx11c_TE_v11.covmat'))
plik.append(batchJob.dataSet(['plik', 'EE'], ['plik_dx11c_EE_v11.ini'], covmat='planck_covmats/plik_dx11c_EE_v11.covmat'))


start_at_bestfit = False
newCovmats = True

groups = []

g = batchJob.jobGroup('main')
# Main group with just tau prior

g.datasets = copy.deepcopy(v96detsets) + copy.deepcopy(v96CS) + copy.deepcopy(plik)

for d in g.datasets:
    d.add(tauname, tauprior)
g.params = [[], ['Alens'], ['nnu'], ['mnu'], ['nrun'], ['yhe']]
groups.append(g)


# lowl, same as main but will all + lowl
g = copy.deepcopy(g)
g.groupName = 'lowl'
for d in g.datasets:
    d.add('lowl', ['lowl.ini'])
groups.append(g)


# group to see the effect of tau prior
g = batchJob.jobGroup('WMAPtau')
g.datasets = copy.deepcopy(v96detsets)
for d in g.datasets:
    d.add('WMAPtau', WMAPtau)
g.params = [[], ['Alens']]
groups.append(g)

g = batchJob.jobGroup('highL')
g.datasets = [copy.deepcopy(v96detsets[0])] + [copy.deepcopy(plik[0])]
for d in g.datasets:
    d.add(tauname, tauprior)
    d.add(lowl)
    d.add(highL)
g.params = [[], ['Alens']]
groups.append(g)


g = batchJob.jobGroup('old')
g.datasets = []
for name, x in zip(['v62TN', 'v65F', 'v85F'], ['nonclik_v62TN.ini', 'nonclik.ini', 'nonclik_v85F.ini']):
    g.datasets.append(batchJob.dataSet(name, [x, Camspec]))

for d in g.datasets:
    d.addEnd(tauname, tauprior)
groups.append(g)



def covRenamer(name):
    renamed = re.sub(r'_v.*_highL', '_planck_lowl_lowLike_highL', name, re.I)
    renamed = re.sub(r'_v.*', '_planck_lowl_lowLike', renamed, re.I)
    if renamed == name: return[]
    else: return [renamed]
