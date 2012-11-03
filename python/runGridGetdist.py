import os, iniFile, batchJobArgs

def checkDir(fname):
    if not os.path.exists(fname): os.makedirs(fname)


Opts = batchJobArgs.batchArgs('Run getdist over the grid of models', importance=True)
Opts.parser.add_argument('--plots', action='store_true')
Opts.parser.add_argument('--norun', action='store_true')
Opts.parser.add_argument('--noimportance', action='store_true')

(batch, args) = Opts.parseForBatch()

base_ini = 'getdist_common_batch1.ini'

matlab = 'matlab'

plot_types = ['.m', '_2D.m', '_3D.m']
# '_tri.m' is very slow for so many


data_dir = batch.batchPath + 'plot_data' + os.sep
ini_dir = batch.batchPath + 'getdist' + os.sep

checkDir(data_dir)
checkDir(ini_dir)

if not args.plots:
    for jobItem in batch.items():
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
        if args.importance is None:
            fname = ini_dir + jobItem.name + '.ini'
            ini.saveFile(fname)
            if not args.norun:
                print "running: " + fname
                os.system('./getdist ' + fname)

# add ini files for importance sampling runs
        if not args.noimportance:
            for imp in jobItem.importanceJobs():
                if Opts.wantImportance(imp.importanceTag):
                    ini.params['file_root'] = imp.chainRoot
                    fname = ini_dir + imp.name + '.ini'
                    ini.params['compare_num'] = 1
                    ini.params['compare1'] = jobItem.chainRoot
                    ini.saveFile(fname)
                    if os.path.exists(ini.params['file_root'] + '_1.txt') and not args.norun:
                        print "running: " + fname
                        os.system('./getdist ' + fname)


if not args.norun:
    cat_cmd = 'cat '
    for jobItem in batch.items(wantImportance=not args.noimportance):
        if Opts.jobItemWanted(jobItem):
            os.chdir(jobItem.distPath)
            for tp in plot_types:
                fname = jobItem.distRoot + tp
                print fname
                if os.path.exists(fname):
                    cat_cmd = cat_cmd + ' ' + fname
    if len(cat_cmd) > 5: os.system(cat_cmd + '|' + matlab)
