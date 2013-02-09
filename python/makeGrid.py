import os, sys, batchJob, iniFile


if len(sys.argv) < 2:
    print 'Usage: python/makeGrid.py new_directory_for_outputs grid_settings_python_file'
    print 'e.g. python/makeGrid.py /scratch/aml1005/planckgrids/testchain settings_testchain'
    sys.exit()


batchPath = os.path.abspath(sys.argv[1]) + os.sep

# 0: chains, 1: importance sampling, 2: best-fit, 3: best-fit and Hessian
cosmomcAction = 0

batch = batchJob.batchJob(batchPath)

settings = __import__(sys.argv[2])

# priors and widths for parameters which are varied
if not hasattr(settings, 'params'):
    params = dict()
    params['mnu'] = '0 0 5 0.1 0.03'
    params['omegak'] = '-0.0008 -0.3 0.3 0.001 0.001'  # starting exactly on flat seems to confuse minimizer
    params['w'] = '-1 -3 -0.3 0.02 0.02'
    params['nnu'] = '3.046 0.05 10 0.05 0.05'
    params['nrun'] = '0 -1 1 0.001 0.001'
    params['r'] = '0 0 2 0.03 0.03'
    params['Alens'] = '1 0 10 0.05 0.05'
    params['yhe'] = '0.245 0.1 0.5 0.006 0.006'
    params['alpha1'] = '0 -1 1 0.0003 0.0003'
    params['deltazrei'] = '0.5 0.1 3 0.3 0.3'
    params['wa'] = '0 -2 2 0.3 0.3'
    params['meffsterile'] = '0 0 5 0.1 0.03'
    settings.params = params

batch.datasets = settings.datasets
batch.extparams = settings.extparams
batch.skip = settings.skip
batch.makeItems(settings.importanceRuns)
batch.makeDirectories()
batch.save()

for jobItem in batch.items(wantSubItems=False):

        jobItem.makeChainPath()
        ini = iniFile.iniFile()

        for param in jobItem.param_set:
            ini.params['param[' + param + ']'] = settings.params[param]
        for iniitem in jobItem.data_set[1]:
            ini.defaults.append(batch.commonPath + iniitem)
        for deffile in settings.defaults:
            ini.defaults.append(batch.commonPath + deffile)

        if 'mnu' in jobItem.param_set:
            ini.params['num_massive_neutrinos'] = 3
        if 'meffsterile' in jobItem.param_set:
            ini.params['param[mnu]'] = '0.06 0.06 0.06 0 0'
            ini.params['param[nnu]'] = '3.1 3.046 10 0.05 0.05'
            ini.params['num_massive_neutrinos'] = 1
        if 'yhe' in jobItem.param_set:
            ini.params['bbn_consistency'] = False
        if 'r' in jobItem.param_set:
            ini.params['compute_tensors'] = True

        ini.params['file_root'] = jobItem.chainRoot

        covmat = batch.basePath + 'planck_covmats/' + jobItem.name + '.covmat'
        if os.path.exists(covmat):
            ini.params['propose_matrix'] = covmat
        else:
            hasCov = False
            ini.params['MPI_Max_R_ProposeUpdate'] = 20
            for new, old in settings.covrenames.items():
                covmat = batch.basePath + 'planck_covmats/' + jobItem.name.replace(new, old) + '.covmat'
                if os.path.exists(covmat):
                    ini.params['propose_matrix'] = covmat
                    hasCov = True
                    break

        ini.params['start_at_bestfit'] = settings.start_at_bestfit
        ini.params['action'] = cosmomcAction
        ini.saveFile(jobItem.iniFile())
        if not settings.start_at_bestfit:
            ini.params['action'] = 2
            variant = '_minimize'
            ini.saveFile(jobItem.iniFile(variant))


# add ini files for importance sampling runs
        for imp in jobItem.importanceJobs():
            if batch.hasName(imp.name.replace('_post', '')): raise Exception('importance sampling something you already have?')
            for minimize in (False, True):
                ini = iniFile.iniFile()
                for inc in imp.importanceSettings:
                    ini.includes.append(batch.commonPath + inc)
                if cosmomcAction == 0 and not minimize:
                    for deffile in settings.importanceDefaults:
                        ini.defaults.append(batch.commonPath + deffile)
                    ini.params['redo_outroot'] = imp.chainRoot
                    ini.params['action'] = 1
                else:
                    ini.params['file_root'] = imp.chainRoot
                if minimize:
                    ini.params['action'] = 2
                    variant = '_minimize'
                else: variant = ''
                ini.defaults.append(jobItem.iniFile())
                ini.saveFile(imp.iniFile(variant))
                if cosmomcAction != 0: break


comment = 'Done... to run do: python python/runbatch.py ' + batchPath + ' --noimportance'
print comment
if not settings.start_at_bestfit:
    comment = '....... for best fits: python python/runbatch.py ' + batchPath + ' --minimize'
    print comment
print ''
comment = 'for importance sampled: python python/runbatch.py ' + batchPath + ' --importance'
print comment
comment = 'for best-fit for importance sampled: python python/runbatch.py ' + batchPath + ' --importance_minimize'
print comment
