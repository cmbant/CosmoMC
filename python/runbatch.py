#!/usr/bin/env python
import os, batchJobArgs

Opts = batchJobArgs.batchArgs('Submit jobs to run chains or importance sample', importance=True)
Opts.parser.add_argument('--nodes', type=int, default=2)

(batch, args) = Opts.parseForBatch()


# if len(sys.argv) < 2:
#    print 'Usage: python/runbatch.py batch_directory [chains/importance] [num_nodes]'
#    sys.exit()

subScript = 'runMPI_HPCS.pl'



def submitJob(ini):
        command = 'perl ' + subScript + ' ' + ini + ' ' + str(args.nodes)
        print 'Submitting...' + command
        os.system(command)


for jobItem in batch.items(wantSubItems=False):
    if args.importance == None:
        submitJob(jobItem.iniFile())
    else:
        for imp in jobItem.importanceJobs():
            if Opts.wantImportance(imp.importanceTag): submitJob(imp.postIniFile())



# iniDir = os.path.abspath(sys.argv[1]) + os.sep
# dirList = os.listdir(iniDir)
# for fname in dirList:
#   command = 'perl ' + subScript + ' ' + iniDir + fname + ' ' + noOfMpiNodes
#    print 'Submitting...' + command
#    os.system(command)
