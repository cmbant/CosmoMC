# sample settings for a particular grid run

# sets of parameters to vary in addition to baseline
extparams = [[], ['omegak', 'mnu'], ['nrun', 'r'], ['omegak'], ['nnu'], ['nrun']]

# dataset names
planck = 'planck_CAMspec'
highL = 'highL'
BAO = 'BAO'
HST = 'HST'

datasets = []
# lists of dataset names to combine, with corresponding sets of inis to include
datasets.append([[planck], ['CAMspec_defaults.ini']])
datasets.append([[planck, highL], ['CAMspec_ACTSPT_defaults.ini']])

# add importance name tags, and list of specific .ini files to include (in batch1/)
importanceRuns = []
importanceRuns.append([BAO, ['BAO.ini']])
importanceRuns.append([HST, ['HST.ini']])

# priors and widths for parameters which are varied
params = dict()
params['mnu'] = '0 0 5 0.1 0.03'
params['omegak'] = '0 -0.3 0.3 0.001 0.001'
params['w'] = '-1 -3 -0.3 0.02 0.02'
params['nnu'] = '3.046 0 10 0.05 0.05'
params['nrun'] = '0 -1 1 0.001 0.001'
params['r'] = '0 0 2 0.03 0.03'
params['Alens'] = '0 0 10 0.05 0.05'
params['yhe'] = '0.245 0.1 0.5 0.006 0.006'
params['alpha1'] = '0 -1 1 0.0003 0.0003'
params['deltazrei'] = '0.5 0.1 3 0.3 0.3'
params['wa'] = '0 -2 2 0.3 0.3'

skip = []
skip.append('base_omegak_mnu_planck_CAMspec_highL')

# if covmats are unreliable, so start learning ASAP
newCovmat = True

start_at_bestfit = True

# ini files you want to base each set of runs on
defaults = ['common_batch1.ini']

importanceDefaults = ['importance_sampling.ini']
