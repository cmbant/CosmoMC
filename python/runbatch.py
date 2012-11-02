#!/usr/bin/env python
import sys, os, batchJob
try: import argparse
except:
    print 'use module load to load python 2.7'
    sys.exit()

parser = argparse.ArgumentParser(description='Submit jobs to run chains or importance sample')

parser.add_argument('batchPath')

parser.add_argument('--nodes', type=int, default=2)

parser.add_argument('--importance', nargs='*', default=None)

args = parser.parse_args()


# if len(sys.argv) < 2:
#    print 'Usage: python/runbatch.py batch_directory [chains/importance] [num_nodes]'
#    sys.exit()

subScript = 'runMPI_HPCS.pl'

batch = batchJob.readobject(args.batchPath)

def submitJob(ini):
        command = 'perl ' + subScript + ' ' + ini + ' ' + str(args.nodes)
        print 'Submitting...' + command
        os.system(command)


for jobItem in batch.items(wantSubItems=False):
    if args.importance == None:
        submitJob(jobItem.iniFile())
    else:
        for imp in jobItem.importanceJobs():
            if (len(args.importance) == 0 or imp.importanceTag in args.importance):
                submitJob(imp.postIniFile())



# iniDir = os.path.abspath(sys.argv[1]) + os.sep
# dirList = os.listdir(iniDir)
# for fname in dirList:
#   command = 'perl ' + subScript + ' ' + iniDir + fname + ' ' + noOfMpiNodes
#    print 'Submitting...' + command
#    os.system(command)
