import re, batchJob, copy

ini_dir = 'batch2/'

defaults = ['common.ini']

importanceDefaults = ['importance_sampling.ini']

Camspec = 'CAMspec_defaults.ini'
# dataset names
tauprior = {'prior[tau]':'0.06 0.02'}
no217auto = {'want_spec':'T T F T'}
freecal = {'prior[cal0]':' 1 0.05', 'prior[cal2]':'1 0.05'}


planck_vars = ['pico.ini', 'nonclik_v85F.ini', Camspec, 'lowl.ini', tauprior, {'indep_sample':0} ]

wig1800_217 = {'param[wig2_217]':'0 -50 50 3 3'}
wig1800_143 = {'param[wig2_143]':'0 -50 50 3 3'}
wig1460_217 = {'param[wig1_217]':'0 -50 50 3 3'}

datasets = []

datasets.append(batchJob.dataSet('v85F', planck_vars))
datasets.append(batchJob.dataSet('no217auto', [no217auto] + planck_vars))
datasets.append(batchJob.dataSet('wig1800_217', planck_vars + [wig1800_217]))
datasets.append(batchJob.dataSet('wig1800_143', planck_vars + [wig1800_143]))
datasets.append(batchJob.dataSet('wig1460_217', planck_vars + [wig1460_217]))
datasets.append(batchJob.dataSet('wig1800_143_217', planck_vars + [wig1800_217, wig1800_143]))
datasets.append(batchJob.dataSet('wig_both', planck_vars + [wig1800_217, wig1800_143, wig1460_217]))

covmat = 'planck_covmats/base_planck_lowl_lowLike.covmat'

start_at_bestfit = False
newCovmats = True

# set up groups of parameters and data sets
class group:pass

groups = []

g = group()
# sets of parameters to vary in addition to baseline
g.params = [[]]

# lists of dataset names to combine, with corresponding sets of inis to include
g.datasets = datasets

g.importanceRuns = []
g.groupName = 'main'
groups.append(g)


g = group()
datasets = copy.deepcopy(datasets)
for d in datasets:
    d.names = ['freecal'] + set.names
    d.params = [freecal] + set.params
# sets of parameters to vary in addition to baseline
g.params = [[]]

# lists of dataset names to combine, with corresponding sets of inis to include
g.datasets = datasets

g.importanceRuns = []
g.groupName = 'freecal'
groups.append(g)
