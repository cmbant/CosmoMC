import os, sys, batchJob, iniFile

def checkDir(fname):
    if not os.path.exists(fname): os.makedirs(fname)


if len(sys.argv) < 2:
    print 'Usage: python/runGridGetdist.py directory_with_outputs [-plots, -imponly]'
    sys.exit()

base_ini = 'getdist_common_batch1.ini'

matlab = 'matlab'

batchPath = os.path.abspath(sys.argv[1]) + os.sep
batch = batchJob.readobject(batchPath)


data_dir = batchPath + 'plot_data' + os.sep
ini_dir = batchPath + 'getdist' + os.sep

checkDir(data_dir)
checkDir(ini_dir)

opt = ''
if len(sys.argv) >= 3: opt = sys.argv[2]

if opt != '-plots':
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
        fname = ini_dir + jobItem.name + '.ini'
        ini.saveFile(fname)
        if opt != '-imponly' and opt != '-norun':
            print "running: " + fname
            os.system('./getdist ' + fname)

# add ini files for importance sampling runs
        for imp in jobItem.importanceJobs():
            ini.params['file_root'] = imp.chainRoot
            fname = ini_dir + imp.name + '.ini'
            ini.params['compare_num'] = 1
            ini.params['compare1'] = jobItem.chainRoot
            ini.saveFile(fname)
            if os.path.exists(ini.params['file_root'] + '_1.txt') and opt != '-norun':
                print "running: " + fname
                os.system('./getdist ' + fname)


if opt != '-norun':

    plot_types = ['.m', '_2D.m', '_3D.m']
    # '_tri.m' is very slow for so many

    if opt != '-imponly':
        cat_cmd = 'cat '
        for jobItem in batch.items():
            os.chdir(jobItem.distPath)
            for tp in plot_types:
                fname = jobItem.distPath + jobItem.name + tp
                print fname
                if os.path.exists(fname):
                    cat_cmd = cat_cmd + ' ' + fname
        if len(cat_cmd) > 5: os.system(cat_cmd + '|' + matlab)

    cat_cmd = 'cat '
    for jobItem in batch.items():
        os.chdir(jobItem.distPath)
        for imp in jobItem.importanceRuns:
            tag = '_post_' + imp
            for tp in plot_types:
                fname = jobItem.distPath + jobItem.name + tag + tp
                print fname
                if os.path.exists(fname):
                    cat_cmd = cat_cmd + ' ' + fname
    if len(cat_cmd) > 5: os.system(cat_cmd + '|' + matlab)


