import os, iniFile, batchJobArgs

def checkDir(fname):
    if not os.path.exists(fname): os.makedirs(fname)


Opts = batchJobArgs.batchArgs('Run getdist over the grid of models', notExist=True)
Opts.parser.add_argument('--plots', action='store_true')
Opts.parser.add_argument('--norun', action='store_true')
Opts.parser.add_argument('--noplots', action='store_true')
Opts.parser.add_argument('--specific', action='store_true')


(batch, args) = Opts.parseForBatch()

base_ini = 'getdist_common_batch1.ini'

matlab = 'matlab'

plot_types = ['.m', '_2D.m', '_3D.m']
# '_tri.m' is very slow for so many


data_dir = batch.batchPath + 'plot_data' + os.sep
ini_dir = batch.batchPath + 'getdist' + os.sep

checkDir(data_dir)
checkDir(ini_dir)

if not args.plots and not args.specific:
    for jobItem in Opts.filteredBatchItems():
            ini = iniFile.iniFile()
            ini.params['file_root'] = jobItem.chainRoot
            checkDir(jobItem.distPath)
            ini.params['out_dir'] = jobItem.distPath
            ini.params['plot_data_dir'] = data_dir
            custom_plot = batch.commonPath + 'plots' + os.sep + jobItem.paramtag + '.ini'
            custom_plot2 = batch.commonPath + 'plots' + os.sep + jobItem.name + '.ini'
            if os.path.exists(custom_plot2):
                ini.includes.append(custom_plot2)
            elif os.path.exists(custom_plot):
                ini.includes.append(custom_plot)
            elif len(jobItem.param_set) > 0:
                ini.params['plot_2D_param'] = jobItem.param_set[0]
            ini.defaults.append(batch.commonPath + base_ini)
            if jobItem.isImportanceJob:
                ini.params['compare_num'] = 1
                ini.params['compare1'] = jobItem.parent.chainRoot
            fname = ini_dir + jobItem.name + '.ini'
            ini.saveFile(fname)
            if not args.norun and (not args.notexist or not jobItem.getDistExists()):
                if jobItem.chainExists():
                    print "running: " + fname
                    os.system('./getdist ' + fname)
                else: print "Chains do not exist yet: " + jobItem.chainRoot


if not args.norun and not args.noplots:
    if not args.specific:
        cat_cmd = 'cat '
        for jobItem in Opts.filteredBatchItems():
                os.chdir(jobItem.distPath)
                for tp in plot_types:
                    fname = jobItem.distRoot + tp
                    print fname
                    if os.path.exists(fname):
                        cat_cmd = cat_cmd + ' ' + fname
        if len(cat_cmd) > 5: os.system(cat_cmd + '|' + matlab)

    if args.name is None:
        specficPath = batch.batchPath + 'specific_plots' + os.sep
        checkDir(specficPath)
        os.chdir(specficPath)
        os.system('export getdist_plot_data=' + data_dir + ';export MATLABPATH=' + batch.basePath + 'mscripts' + os.sep
                    + ';cat ' + batch.commonPath + 'specificplots' + os.sep + '*.m | ' + matlab)

