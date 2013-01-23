#!/usr/bin/env python
import os, batchJobArgs, sys

Opts = batchJobArgs.batchArgs('Submit jobs to run chains or importance sample', notExist=True, converge=True)
Opts.parser.add_argument('--nodes', type=int, default=2)
Opts.parser.add_argument('--script', default='runMPI_HPCS.pl')

Opts.parser.add_argument('--dryrun', action='store_true')
Opts.parser.add_argument('--subitems', action='store_true')
Opts.parser.add_argument('--importance_minimize', action='store_true')


(batch, args) = Opts.parseForBatch()

variant = ''
if args.importance_minimize:
    variant = '_minimize'
    if args.importance is None: args.importance = []
    args.nodes = 0

def submitJob(ini):
        command = 'perl ' + args.script + ' ' + ini + ' ' + str(args.nodes)
        print 'Submitting...' + command
        if not args.dryrun: os.system(command)


for jobItem in Opts.filteredBatchItems(wantSubItems=args.subitems):
        if not args.notexist or args.importance_minimize and not jobItem.chainMinimumExists() or not args.importance_minimize and not jobItem.chainExists():
            if args.converge == 0 or jobItem.hasConvergeBetterThan(args.converge): submitJob(jobItem.iniFile(variant))
