import os, sys, batchJob, iniFile

def checkDir(fname):
    if not os.path.exists(fname): os.makedirs(fname)


if len(sys.argv) < 2:
    print 'Usage: python/runGridGetdist.py directory_with_outputs [-plots]'

base_ini = 'getdist_common_batch1.ini'

matlab = 'matlab'

batchPath = os.path.abspath(sys.argv[1]) + os.sep
batch = batchJob.readobject(batchPath + 'batch.pyobj')


data_dir = batchPath + 'plot_data' + os.sep
ini_dir = batchPath + 'getdist' + os.sep

checkDir(data_dir)
checkDir(ini_dir)

if len(sys.argv) < 3 or sys.argv[2] != '-plots':
    for jobItem in batch.items():
        ini = iniFile.iniFile()
        ini.params['file_root'] = jobItem.chainRoot
        checkDir(jobItem.distPath)
        ini.params['out_dir'] = jobItem.distPath
        ini.params['plot_data_dir'] = data_dir
        custom_plot = batch.commonPath + 'plots' + os.sep + jobItem.paramtag + '.ini'
        if os.path.exists(custom_plot):
            ini.includes.append(custom_plot)
        elif len(jobItem.param_set) > 0:
            ini.params['plot_2D_param'] = jobItem.param_set[0]
        ini.defaults.append(batch.commonPath + base_ini)
        fname = ini_dir + jobItem.name + '.ini'
        ini.saveFile(fname)
        print "running: " + fname
        os.system('./getdist ' + fname)


plot_types = ['.m', '_2D.m', '_3D.m', '_tri.m']

cat_cmd = 'cat '
for jobItem in batch.items():
    os.chdir(jobItem.distPath)
    for tp in plot_types:
        fname = jobItem.distPath + jobItem.name + tp
        print fname
        if os.path.exists(fname):
            cat_cmd = cat_cmd + ' ' + fname
os.system(cat_cmd + '|' + matlab)



