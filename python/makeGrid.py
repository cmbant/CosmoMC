import os, sys, batchJob, iniFile


if len(sys.argv) < 2:
    print 'Usage: python/makeGrid.py new_directory_for_outputs [action]'
    sys.exit()


batchPath = os.path.abspath(sys.argv[1]) + os.sep

# 0: chains, 1: importance sampling, 2: best-fit, 3: best-fit and Hessian
cosmomcAction = 0
if len(sys.argv) > 2: cosmomcAction = int(sys.argv[2])


batch = batchJob.batchJob(batchPath)

# sets of parameters to vary in addition to baseline
batch.extparams = [[], ['omegak', 'mnu'], ['nrun', 'r'], ['omegak'], ['nnu'], ['nrun']]

# dataset names
planck = 'planck_CAMspec'
highL = 'highL'


batch.datasets = []
# lists of dataset names to combine, with corresponding sets of inis to include
batch.datasets.append([[planck], ['CAMspec_defaults.ini']])
batch.datasets.append([[planck, highL], ['CAMspec_ACTSPT_defaults.ini']])

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

# if covmats are unreliable, so start learning ASAP
newCovmat = True

# ini files you want to base each set of runs on
defaults = ['common_batch1.ini']

batch.makeDirectories()
batch.save()

for jobItem in batch.items():
        jobItem.makeChainPath()
        ini = iniFile.iniFile()

        for param in jobItem.param_set:
            ini.params['param[' + param + ']'] = params[param]
        for deffile in defaults:
            ini.defaults.append(batch.commonPath + deffile)
        for iniitem in jobItem.data_set[1]:
            ini.includes.append(batch.commonPath + iniitem)

        if 'mnu' in jobItem.param_set:
            ini.params['num_massive_neutrinos'] = 3
        if 'nnu' in jobItem.param_set:
            if ('mnu' in jobItem.param_set): raise Exception('no support for nnu and mnu')
            ini.params['param[mnu]'] = '0 0 0 0 0'
        if 'r' in jobItem.param_set:
            ini.params['compute_tensors'] = True

        ini.params['file_root'] = jobItem.chainRoot
#        f.append('propose_matrix =' + jobItem.datatag + '.covmat')

        if newCovmat:
            ini.params['MPI_Max_R_ProposeUpdate'] = 20
            ini.params['start_at_bestfit'] = True

        ini.params['action'] = cosmomcAction
        ini.saveFile(jobItem.iniFile())

comment = 'Done... to run do: python python/runbatch.py ' + batchPath + 'iniFiles'
if cosmomcAction == 3 or cosmomcAction == 2: comment += ' 0'
print comment
