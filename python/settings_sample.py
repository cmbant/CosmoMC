# sample settings for a particular grid run

from paramgrid import batchjob

# Directory to find .ini files
ini_dir = 'batch2/'

# directory to look for existing covariance matrices
cov_dir = 'planck_covmats/'

# ini files you want to base each set of runs on
defaults = ['common.ini']
importanceDefaults = ['importance_sampling.ini']

# set up list of groups of parameters and data sets
groups = []

# make first group of runs (all parameter variations with all data combinations)
g = batchjob.jobGroup('main')

g.params = [[], ['mnu'], ['nnu']]

g.datasets = []

# lists of dataset names to combine, with corresponding sets of inis to include
g.datasets.append(batchjob.dataSet(['plikHM', 'TT', 'lowTEB'], ['plik_dx11dr2_HM_v18_TT.ini', 'lowTEB.ini']))
g.datasets.append(batchjob.dataSet(['plikHM', 'TT', 'lowTEB', 'lensing'],
                                   ['plik_dx11dr2_HM_v18_TT.ini', 'lowTEB.ini', 'lensing.ini']))


# add importance name tags, and list of specific .ini files to include (in batch1/)
g.importanceRuns = []
g.importanceRuns.append([['BAO'], ['BAO.ini']])

groups.append(g)

# ranges for parameters when they are varied (can delete params if you just want to use defaults)
params = dict()
params['w'] = '-0.99 -3. 1 0.02 0.02'
params['wa'] = '0 -3 2 0.05 0.05'
params['mnu'] = '0.02 0 5 0.1 0.03'
params['omegak'] = '-0.0008 -0.3 0.3 0.001 0.001'  # starting exactly on flat seems to confuse minimizer
params['nnu'] = '3.046 0.05 10 0.05 0.05'
params['nrun'] = '0 -1 1 0.005 0.001'
params['r'] = '0 0 3 0.03 0.03'
params['Alens'] = '1 0 10 0.05 0.05'
params['yhe'] = '0.245 0.1 0.5 0.006 0.006'
params['alpha1'] = '0 -1 1 0.0003 0.0003'
params['meffsterile'] = '0.1 0 3 0.1 0.03'


# extra parameters that are set only when specific parameters are varied. Can deleted to get defaults.
param_extra_opts = {
    'mnu': {'num_massive_neutrinos': 3},
    'meffsterile': {'param[mnu]': '0.06', 'param[nnu]': '3.1 3.046 10 0.05 0.05', 'num_massive_neutrinos': 1,
                    'accuracy_level': 1.2},
    'yhe': {'bbn_consistency': False},
    'r': {'compute_tensors': True},
    'nt': {'inflation_consistency': False, 'lmax_tensor': 1000}
}


# try to match each new run name to exisitng covmat (covariance matrix for efficient exploration)
# e.g. try to map to get name without particular data combinations
covWithoutNameOrder = ['lensing', 'BAO']
# or try replacing various names (these are standard for the provided planck_covmats)
covNameMappings = {'plikHM': 'plik', 'plikLite': 'plik'}
# for mapping to nominal mission names try
# covNameMappings = {'plikHM':'planck','TT':'', 'lowTEB':'lowLike'}

