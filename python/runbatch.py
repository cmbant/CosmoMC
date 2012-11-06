#!/usr/bin/env python
import os, batchJobArgs

Opts = batchJobArgs.batchArgs('Submit jobs to run chains or importance sample')
Opts.parser.add_argument('--nodes', type=int, default=2)
Opts.parser.add_argument('--script', default='runMPI_HPCS.pl')

(batch, args) = Opts.parseForBatch()


def submitJob(ini):
        command = 'perl ' + args.script + ' ' + ini + ' ' + str(args.nodes)
        print 'Submitting...' + command
        os.system(command)


for jobItem in Opts.filteredBatchItems(wantSubItems=False):
        submitJob(jobItem.iniFile())
