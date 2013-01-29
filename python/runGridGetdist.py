import os, iniFile, batchJobArgs

def checkDir(fname):
    if not os.path.exists(fname): os.makedirs(fname)


Opts = batchJobArgs.batchArgs('Run getdist over the grid of models', notExist=True)
Opts.parser.add_argument('--plots', action='store_true')
Opts.parser.add_argument('--norun', action='store_true')
Opts.parser.add_argument('--noplots', action='store_true')
Opts.parser.add_argument('--specific', action='store_true')
Opts.parser.add_argument('--compare_only', action='store_true')
Opts.parser.add_argument('--settings', default='getdist_settings')


(batch, args) = Opts.parseForBatch()

base_ini = 'getdist_common_batch1.ini'

matlab = 'matlab'

plot_types = ['.m', '_2D.m', '_3D.m']
# '_tri.m' is very slow for so many


data_dir = batch.batchPath + 'plot_data' + os.sep
ini_dir = batch.batchPath + 'getdist' + os.sep

checkDir(data_dir)
checkDir(ini_dir)

settings = __import__(args.settings)

if not args.plots and not args.specific:
    for compare in (([None], [])[args.compare_only] + settings.compare_datatag):
        for jobItem in Opts.filteredBatchItems():
            if compare is  None or jobItem.datatag == compare.datatag:
                ini = iniFile.iniFile()
                ini.params['file_root'] = jobItem.chainRoot
                checkDir(jobItem.distPath)
                ini.params['out_dir'] = jobItem.distPath
                ini.params['plot_data_dir'] = data_dir
                if 'meffsterile' in jobItem.param_set: ini.params['limits[nnu]'] = '3.046 N'
                custom_plot = batch.commonPath + 'plots' + os.sep + jobItem.paramtag + '.ini'
                custom_plot2 = batch.commonPath + 'plots' + os.sep + jobItem.name + '.ini'
                if os.path.exists(custom_plot2):
                    ini.includes.append(custom_plot2)
                elif os.path.exists(custom_plot):
                    ini.includes.append(custom_plot)
                elif len(jobItem.param_set) > 0:
                    ini.params['plot_2D_param'] = jobItem.param_set[0]
                ini.defaults.append(batch.commonPath + base_ini)
                tag = ''
                if not compare is None:
                    tag = compare.tag
                    ini.params['plots_only'] = True
                    ini.params['out_root'] = jobItem.name + tag
                    ini.params['compare_num'] = len(compare.compares)
                    for i in range(1, len(compare.compares) + 1):
                        ini.params['compare' + str(i)] = jobItem.name.replace(compare.datatag, compare.compares[i - 1])
                    if not hasattr(jobItem, 'compareRoots'): jobItem.compareRoots = []
                    jobItem.compareRoots += [jobItem.distRoot + tag]
                elif jobItem.isImportanceJob:
                    ini.params['compare_num'] = 1
                    ini.params['compare1'] = jobItem.parent.chainRoot
                fname = ini_dir + jobItem.name + tag + '.ini'
                ini.saveFile(fname)
                if not args.norun and (not args.notexist or not jobItem.getDistExists()):
                    if jobItem.chainExists():
                        print "running: " + fname
                        os.system('./getdist ' + fname)
                    else: print "Chains do not exist yet: " + jobItem.chainRoot



if not args.norun and not args.noplots or args.compare_only:
    if not args.specific:
        cat_cmd = 'cat '
        for jobItem in Opts.filteredBatchItems():
                os.chdir(jobItem.distPath)
                for tp in plot_types:
                    if not args.compare_only:
                        fname = jobItem.distRoot + tp
                        print fname
                        if os.path.exists(fname):
                            cat_cmd = cat_cmd + ' ' + fname
                    if hasattr(jobItem, 'compareRoots'):
                        for root in jobItem.compareRoots:
                            fname = root + tp
                            print fname
                            if os.path.exists(fname):
                                cat_cmd = cat_cmd + ' ' + fname

        if len(cat_cmd) > 5: os.system(cat_cmd + '|' + matlab)

    if args.name is None and not args.compare_only:
        specficPath = batch.batchPath + 'specific_plots' + os.sep
        checkDir(specficPath)
        os.chdir(specficPath)
        os.system('export getdist_plot_data=' + data_dir + ';export MATLABPATH=' + batch.basePath + 'mscripts' + os.sep
                    + ';cat ' + batch.commonPath + 'specificplots' + os.sep + '*.m | ' + matlab)

