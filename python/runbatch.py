#!/usr/bin/env python
import os, batchJobArgs

Opts = batchJobArgs.batchArgs('Submit jobs to run chains or importance sample')
Opts.parser.add_argument('--nodes', type=int, default=2)
Opts.parser.add_argument('--script', default='runMPI_HPCS.pl')
Opts.parser.add_argument('--notexist', action='store_true')
Opts.parser.add_argument('--dryrun', action='store_true')
Opts.parser.add_argument('--subitems', action='store_true')

(batch, args) = Opts.parseForBatch()


def submitJob(ini):
        command = 'perl ' + args.script + ' ' + ini + ' ' + str(args.nodes)
        print 'Submitting...' + command
        if not args.dryrun: os.system(command)


for jobItem in Opts.filteredBatchItems(wantSubItems=args.subitems):
        if not args.notexist or not os.path.exists(jobItem.chainRoot + '_1.txt'):
            submitJob(jobItem.iniFile())
