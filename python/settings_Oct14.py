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
tauprior = {'prior[tau]':'0.07 0.015'}
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

TT100_217 = {'want_spec':'T F T F F F', 'param[cal2]':'0.995', 'param[aps143]':'0', 'param[psr]':'1', 'param[cibr]':'1'}
TT100_143 = {'want_spec':'T T F F F F', 'param[cal2]':'0.995', 'param[aps217]':'0', 'param[psr]':'1', 'param[cibr]':'1'}
no217auto = {'want_spec':'T T F T F F'}


freecal = 'freecal.ini'


CamSpecVars = ['v910F', 'v910CMH']
planck_detsets = [freecal, 'nonclik_detsets.ini']
planck_CS = [freecal, 'nonclik.ini']

# not using these checks yet
wig1800_217 = {'param[wig2_217]':'0 -50 50 3 3'}
wig1800_143 = {'param[wig2_143]':'0 -50 50 3 3'}
wig1460_217 = {'param[wig1_217]':'0 -50 50 3 3'}

detsets = []
CS = []
for name, datasets, planck_vars in zip(CamSpecVars, [detsets, CS], [planck_detsets, planck_CS]):
    datasets.append(batchJob.dataSet([name , 'TT'], planck_vars + ['CAMspec_TT.ini']))

detsetsTT = copy.copy(detsets)
CSTT = copy.copy(CS)

for name, datasets, planck_vars in zip(CamSpecVars, [detsets, CS], [planck_detsets, planck_CS]):
    datasets.append(batchJob.dataSet([name , 'TE'], planck_vars + ['CAMspec_TE.ini']))
    datasets.append(batchJob.dataSet([name , 'EE'], planck_vars + ['CAMspec_EE.ini']))
    datasets.append(batchJob.dataSet([name , 'TTTEEE'], planck_vars + ['CAMspec_TTTEEE.ini']))


plikHM = []
plikHM.append(batchJob.dataSet(['plikHM', 'TT'], ['plik_dx11dr2_HM_v14_TT.ini'], covmat='planck_covmats/plik_dx11dr2_DS_v14_TT.covmat'))
plikHM.append(batchJob.dataSet(['plikHM', 'TE'], ['plik_dx11dr2_HM_v14_TE.ini'], covmat='planck_covmats/plik_dx11dr2_DS_v14_TE.covmat'))
plikHM.append(batchJob.dataSet(['plikHM', 'EE'], ['plik_dx11dr2_HM_v14_EE.ini'], covmat='planck_covmats/plik_dx11dr2_DS_v14_EE.covmat'))
plikHM.append(batchJob.dataSet(['plikHM', 'TTTEEE'], ['plik_dx11dr2_HM_v14_TTTEEE.ini'], covmat='planck_covmats/plik_dx11dr2_DS_v14_TTTEEE.covmat'))

plikDS = []
plikDS.append(batchJob.dataSet(['plikDS', 'TT'], ['plik_dx11dr2_DS_v14_TT.ini'], covmat='planck_covmats/plik_dx11dr2_HM_v14_TT.covmat'))
plikDS.append(batchJob.dataSet(['plikDS', 'TE'], ['plik_dx11dr2_DS_v14_TE.ini'], covmat='planck_covmats/plik_dx11dr2_HM_v14_TE.covmat'))
plikDS.append(batchJob.dataSet(['plikDS', 'EE'], ['plik_dx11dr2_DS_v14_EE.ini'], covmat='planck_covmats/plik_dx11dr2_HM_v14_EE.covmat'))
plikDS.append(batchJob.dataSet(['plikDS', 'TTTEEE'], ['plik_dx11dr2_DS_v14_TTTEEE.ini'], covmat='planck_covmats/plik_dx11dr2_HM_v14_TTTEEE.covmat'))

plik = plikHM + plikDS

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
g.datasets = [d for d in g.datasets if ('TT' in d.names or 'TTTEEE' in d.names)]
for d in g.datasets:
    d.add('lowl', ['lowl.ini'])
groups.append(g)


g = batchJob.jobGroup('lmax')
g.datasets = []
for lmax in range(550, 2600, 150):
    sets = copy.deepcopy(detsets) + copy.deepcopy(CS)
    for d in sets:
        d.add(tauname, tauprior)
        d.add('lmax' + str(lmax), {'camspec_lmax': (str(lmax) + ' ') * 6})
    g.datasets += sets

g.params = [[], ['yhe'], ['nnu']]
groups.append(g)

g = batchJob.jobGroup('PCA')
g.datasets = []
for lmax in range(600, 1601, 200):
    datasets = detsetsTT + CSTT
    sets = copy.deepcopy(datasets)
    for d in sets:
        d.add(tauname, tauprior)
        d.add('lowl', ['lowl.ini'])
        d.add('alllmax' + str(lmax), {'camspec_lmax': (str(lmax) + ' ') * 6})
    g.datasets += sets
    sets = copy.deepcopy(datasets)  # + copy.deepcopy(CS)
    for d in sets:
        d.add(tauname, tauprior)
        d.add('alllmin' + str(lmax + 1), {'camspec_lmin': (str(lmax + 1) + ' ') * 6})
    g.datasets += sets
g.params = [[]]
groups.append(g)


g = batchJob.jobGroup('channels')
chopdatasets = []
for aset in detsetsTT + CSTT:
# for name, planck_vars in zip(CamSpecVars, [planck_detsets, planck_CS]):
    for namecut, cutvars in zip(['no143', 'no217', 'no217auto'], [TT100_217, TT100_143, no217auto]):
        d = copy.deepcopy(aset)
        d.add(namecut, cutvars)
        datasets.append(d)
g.datasets = []
for lmax in range(550, 2550, 150):
    sets = copy.deepcopy(chopdatasets)
    for d in sets:
        d.add(tauname, tauprior)
        d.add('lmax' + str(lmax), {'camspec_lmax': (str(lmax) + ' ') * 6})
    g.datasets += sets
g.params = [[], ['Alens'], ['nnu']]
groups.append(g)

g = batchJob.jobGroup('lmin')
lmins = [800, 1200]
g.datasets = []
for lmin in lmins:
    sets = copy.deepcopy(chopdatasets)
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
