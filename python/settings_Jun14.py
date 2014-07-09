import re, batchJob, copy

# Look at Alens=1 for fixed lensing template amplitude

ini_dir = 'batch2/'

defaults = ['common.ini']

importanceDefaults = ['importance_sampling.ini']

override_defaults = ['pico.ini']
extra_opts = {'indep_sample':0, 'checkpoint':'F', 'pre_marged':'F'}

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

TT100_217 = {'want_spec':'T F T F F F', 'param[cal2]':'0.992', 'param[aps143]':'0', 'param[acib143]':'0', 'param[psr]':'1', 'param[cibr]':'1'}
TT100_143 = {'want_spec':'T T F F F F', 'param[cal2]':'0.992', 'param[aps217]':'0', 'param[acib217]':'0', 'param[psr]':'1', 'param[cibr]':'1'}
no217auto = {'want_spec':'T T F T F F'}


freecal = 'freecal.ini'
freecalEE = {'param[calEE]':'1 0.1 2 0.01 0.01', 'prior[calEE]':'1 1'}
freecalTE = {'param[calTE]':'1 0.1 2 0.005 0.005', 'prior[calTE]': '1 1'}


planck_detsets = [freecal, 'nonclik_v97F.ini', Camspec]
planck_CS = [freecal, 'nonclik_v97CS.ini', Camspec]

# not using these checks yet
wig1800_217 = {'param[wig2_217]':'0 -50 50 3 3'}
wig1800_143 = {'param[wig2_143]':'0 -50 50 3 3'}
wig1460_217 = {'param[wig1_217]':'0 -50 50 3 3'}

detsets = []
CS = []
for name, datasets, planck_vars in zip(['v97', 'v97CS'], [detsets, CS], [planck_detsets, planck_CS]):
    datasets.append(batchJob.dataSet([name , 'TT'], [TT] + planck_vars))
    datasets.append(batchJob.dataSet([name , 'TE'], [TE, varTE, freecalTE] + planck_vars))
    datasets.append(batchJob.dataSet([name , 'EE'], [EE, varEE, freecalEE] + planck_vars))
    datasets.append(batchJob.dataSet([name, 'all'], [full, varTE, varEE, freecalTE, freecalEE] + planck_vars))
#    datasets.append(batchJob.dataSet(name + 'TEEE', [TEEE] + planck_vars))
#    datasets.append(batchJob.dataSet(name + 'TTTE', [TTTE] + planck_vars))

plik = []
plik.append(batchJob.dataSet(['plik', 'TT'], ['plik_dx11c_TT_v12.ini'], covmat='planck_covmats/plik_dx11c_TT_v12.covmat'))
plik.append(batchJob.dataSet(['plik', 'TE'], ['plik_dx11c_TE_v12.ini'], covmat='planck_covmats/plik_dx11c_TE_v12.covmat'))
plik.append(batchJob.dataSet(['plik', 'EE'], ['plik_dx11c_EE_v12.ini'], covmat='planck_covmats/plik_dx11c_EE_v12.covmat'))
plik.append(batchJob.dataSet(['plik', 'all'], ['plik_dx11c_TTTEEE_v12.ini'], covmat='planck_covmats/plik_dx11c_TTTEEE_v12.covmat'))


start_at_bestfit = False
newCovmats = True

groups = []

g = batchJob.jobGroup('main')
# Main group with just tau prior

g.datasets = copy.deepcopy(detsets) + copy.deepcopy(CS) + copy.deepcopy(plik)

for d in g.datasets:
    d.add(tauname, tauprior)
g.params = [[], ['Alens'], ['nnu'], ['mnu'], ['nrun'], ['yhe']]
groups.append(g)


# lowl, same as main but will all + lowl
g = copy.deepcopy(g)
g.groupName = 'lowl'
g.datasets = [d for d in g.datasets if ('TT' in d.names or 'all' in d.names)]
for d in g.datasets:
    d.add('lowl', ['lowl.ini'])
groups.append(g)


# group to see the effect of tau prior
g = batchJob.jobGroup('WMAPtau')
g.datasets = copy.deepcopy(detsets)
for d in g.datasets:
    d.add('WMAPtau', WMAPtau)
g.params = [[], ['Alens']]
groups.append(g)

g = batchJob.jobGroup('highL')
g.datasets = [copy.deepcopy(detsets[0])] + [copy.deepcopy(plik[0])]
g.datasets = [d for d in g.datasets if ('TT' in d.names or 'all' in d.names)]
for d in g.datasets:
    d.add(tauname, tauprior)
    d.add(lowl)
    d.add(highL)
g.params = [[], ['Alens']]
groups.append(g)


g = batchJob.jobGroup('old')
g.datasets = []
for name, x in zip(['v62TN', 'v65F', 'v85F', 'v96F'],
                    ['nonclik_v62TN.ini', 'nonclik.ini', 'nonclik_v85F.ini', 'nonclik_v96F.ini']):
    g.datasets.append(batchJob.dataSet(name, [x, Camspec]))

for d in g.datasets:
    d.addEnd(tauname, tauprior)
groups.append(g)


g = batchJob.jobGroup('lmax')
g.datasets = []
for lmax in range(550, 2600, 150):
    sets = copy.deepcopy(detsets) + copy.deepcopy(CS)
    for d in sets:
        d.add(tauname, tauprior)
        d.add('lmax' + str(lmax), {'camspec_lmax': (str(lmax) + ' ') * 6})
    g.datasets += sets

g.params = [[], ['yhe']]
groups.append(g)


g = batchJob.jobGroup('channels')
datasets = []
for name, planck_vars in zip(['v97', 'v97CS'], [planck_detsets, planck_CS]):
    for namecut, cutvars in zip(['no143', 'no217', 'no217auto'], [TT100_217, TT100_143, no217auto]):
        datasets.append(batchJob.dataSet([name , 'TT', namecut], [TT, cutvars] + planck_vars))
g.datasets = []
for lmax in range(550, 2550, 75):
    sets = copy.deepcopy(datasets)
    for d in sets:
        d.add(tauname, tauprior)
        d.add('lmax' + str(lmax), {'camspec_lmax': (str(lmax) + ' ') * 6})
    g.datasets += sets
g.params = [[], ['Alens']]
groups.append(g)

g = batchJob.jobGroup('lmin')
datasets = []
lmins = [800, 1200]
for name, planck_vars in zip(['v97', 'v97CS'], [planck_detsets, planck_CS]):
    for namecut, cutvars in zip(['no143', 'no217', 'no217auto'], [TT100_217, TT100_143, no217auto]):
        datasets.append(batchJob.dataSet([name , 'TT', namecut], [TT, cutvars] + planck_vars))
g.datasets = []
for lmin in lmins:
    sets = copy.deepcopy(datasets)
    for d in sets:
        d.add(tauname, tauprior)
        d.add('lmin' + str(lmin), {'param[cal0]':'0.9997', 'camspec_lmin': '2500 ' + (str(lmin) + ' ') * 5})
    g.datasets += sets
g.params = [[]]
groups.append(g)



def covRenamer(name):
    renamed = re.sub(r'_v.*_highL', '_planck_lowl_lowLike_highL', name, re.I)
    renamed = re.sub(r'_v.*', '_planck_lowl_lowLike', renamed, re.I)
    if renamed == name: return[]
    else: return [renamed]
