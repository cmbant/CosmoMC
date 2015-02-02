# sample settings for a particular grid run

import batchJob

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
g = batchJob.jobGroup('main')

g.params = [[], ['mnu'], ['nnu']]

g.datasets = []

# lists of dataset names to combine, with corresponding sets of inis to include
g.datasets.append(batchJob.dataSet(['plikHM', 'TT', 'lowTEB'], ['plik_dx11dr2_HM_v18_TT.ini', 'lowTEB.ini']))
g.datasets.append(batchJob.dataSet(['plikHM', 'TT', 'lowTEB', 'lensing'], ['plik_dx11dr2_HM_v18_TT.ini', 'lowTEB.ini', 'lensing.ini']))


# add importance name tags, and list of specific .ini files to include (in batch1/)
g.importanceRuns = []
g.importanceRuns.append([['BAO'], ['BAO.ini']])

groups.append(g)

# try to match each new run name to exisitng covmat
# e.g. get name without particular data combinations
covWithoutNameOrder = ['lensing', 'BAO']
# or try replacing various names (these are standard for provided planck_covmats)
covNameMappings = {'plikHM':'plik', 'plikLite':'plik'}
# for mapping to nominal mission names try
# covNameMappings = {'plikHM':'planck','TT':'', 'lowTEB':'lowLike'}

