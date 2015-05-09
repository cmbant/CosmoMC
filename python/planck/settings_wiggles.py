from __future__ import absolute_import
import copy

from paramgrid import batchjob


ini_dir = 'batch2/'

defaults = ['common.ini']

importanceDefaults = ['importance_sampling.ini']

Camspec = 'CAMspec_defaults.ini'
# dataset names
tauprior = {'prior[tau]': '0.06 0.02'}
WMAPtau = {'prior[tau]': '0.09 0.013'}

no217auto = {'want_spec': 'T T F T'}
freecal = {'prior[cal0]': ' 1 0.05', 'prior[cal2]': '1 0.05'}

planck_vars = ['pico.ini', 'nonclik_v85F.ini', Camspec, 'lowl.ini', {'indep_sample': 0}]

wig1800_217 = {'param[wig2_217]': '0 -50 50 3 3'}
wig1800_143 = {'param[wig2_143]': '0 -50 50 3 3'}
wig1460_217 = {'param[wig1_217]': '0 -50 50 3 3'}

datasets = []

datasets.append(batchjob.dataSet('v85F', planck_vars))
datasets.append(batchjob.dataSet('no217auto', [no217auto] + planck_vars))
datasets.append(batchjob.dataSet('wig1800_217', planck_vars + [wig1800_217]))
datasets.append(batchjob.dataSet('wig1800_143', planck_vars + [wig1800_143]))
datasets.append(batchjob.dataSet('wig1460_217', planck_vars + [wig1460_217]))
datasets.append(batchjob.dataSet('wig1800_143_217', planck_vars + [wig1800_217, wig1800_143]))
datasets.append(batchjob.dataSet('wig_both', planck_vars + [wig1800_217, wig1800_143, wig1460_217]))

covmat = 'planck_covmats/base_planck_lowl_lowLike.covmat'

start_at_bestfit = False
newCovmats = True

groups = []

g = batchjob.jobGroup('main')
g.datasets = copy.deepcopy(datasets)
for d in g.datasets:
    d.add(None, tauprior)

groups.append(g)

g = batchjob.jobGroup('freecal')
g.datasets = copy.deepcopy(datasets)
for d in g.datasets:
    d.add(None, tauprior)
    d.add('freecal', freecal)

groups.append(g)

g = batchjob.jobGroup('WMAPtau')
g.datasets = copy.deepcopy(datasets)
for d in g.datasets:
    d.add('WMAPtau', WMAPtau)
    d.addFirst('freecal', freecal)
# sets of parameters to vary in addition to baseline

groups.append(g)



