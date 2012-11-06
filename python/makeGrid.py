import os, sys, batchJob, iniFile


if len(sys.argv) < 3:
    print 'Usage: python/makeGrid.py new_directory_for_outputs grid_settings_python_file [action]'
    print 'e.g. python/makeGrid.py /scratch/aml1005/planckgrids/testchain settings_testchain'
    sys.exit()


batchPath = os.path.abspath(sys.argv[1]) + os.sep

# 0: chains, 1: importance sampling, 2: best-fit, 3: best-fit and Hessian
cosmomcAction = 0
if len(sys.argv) > 3: cosmomcAction = int(sys.argv[3])

if cosmomcAction == 1: print 'actually you don''t need to do this, use action=0 and configure importanceRuns'

batch = batchJob.batchJob(batchPath)

settings = __import__(sys.argv[2])

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
        if 'nnu' in jobItem.param_set:
            if ('mnu' in jobItem.param_set): raise Exception('no support for nnu and mnu')
            ini.params['param[mnu]'] = '0 0 0 0 0'
        if 'r' in jobItem.param_set:
            ini.params['compute_tensors'] = True

        ini.params['file_root'] = jobItem.chainRoot

        covmat = batch.basePath + 'planck_covmats/' + jobItem.name + '.covmat'
        if os.path.exists(covmat):
            ini.params['propose_matrix'] = covmat
        elif settings.newCovmat:
            ini.params['MPI_Max_R_ProposeUpdate'] = 20

        ini.params['start_at_bestfit'] = settings.start_at_bestfit

        ini.params['action'] = cosmomcAction
        ini.saveFile(jobItem.iniFile())

# add ini files for importance sampling runs
        for imp in jobItem.importanceJobs():
            if batch.hasName(imp.name.replace('_post', '')): raise Exception('importance sampling something you already have?')
            ini = iniFile.iniFile()
            for inc in imp.importaceSettings:
                ini.includes.append(batch.commonPath + inc)
            if cosmomcAction == 0:
                for deffile in settings.importanceDefaults:
                    ini.defaults.append(batch.commonPath + deffile)
                ini.params['redo_outroot'] = imp.chainRoot
                ini.params['action'] = 1
            else:
                ini.params['file_root'] = imp.chainRoot
            ini.defaults.append(jobItem.iniFile())
            ini.saveFile(imp.iniFile())


comment = 'Done... to run do: python python/runbatch.py ' + batchPath
if cosmomcAction == 3 or cosmomcAction == 2: comment += ' --nodes 0'
print comment
comment = 'for importance sampled: python python/runbatch.py ' + batchPath + ' --importance'
if cosmomcAction == 3 or cosmomcAction == 2: comment += ' --nodes 0'
print comment

