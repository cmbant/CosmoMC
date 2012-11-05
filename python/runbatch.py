#!/usr/bin/env python
import os, batchJobArgs

Opts = batchJobArgs.batchArgs('Submit jobs to run chains or importance sample')
Opts.parser.add_argument('--nodes', type=int, default=2)

(batch, args) = Opts.parseForBatch()

subScript = 'runMPI_HPCS.pl'

def submitJob(ini):
        command = 'perl ' + subScript + ' ' + ini + ' ' + str(args.nodes)
        print 'Submitting...' + command
        os.system(command)


for jobItem in Opts.filteredBatchItems(wantSubItems=False):
        submitJob(jobItem.iniFile())
