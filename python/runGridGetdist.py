import os, sys, batchJob, iniFile

def checkDir(fname):
    if not os.path.exists(fname): os.makedirs(fname)


if len(sys.argv) < 2:
    print 'Usage: python/runGridGetdist.py directory_with_outputs'

base_ini = 'getdist_common_batch1.ini'

batchPath = os.path.abspath(sys.argv[1]) + os.sep
batch = batchJob.readobject(batchPath + 'batch.pyobj')

ini = iniFile.iniFile()

data_dir = batchPath + 'plot_data' + os.sep
ini_dir = batchPath + 'getdist' + os.sep

checkDir(data_dir)
checkDir(ini_dir)

for jobItem in batch.items():
    ini.params['file_root'] = jobItem.chainRoot
    out_dir = jobItem.chainPath + 'dist/'
    checkDir(out_dir)
    ini.params['out_dir'] = out_dir
    ini.params['plot_data_dir'] = data_dir
    ini.defaults.append(batch.commonPath + base_ini)
    fname = ini_dir + jobItem.name + '.ini'
    ini.saveFile(fname)
    print "running: " + fname
    os.system('./getdist ' + fname)





