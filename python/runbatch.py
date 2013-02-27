#!/usr/bin/env python
import os, batchJobArgs

Opts = batchJobArgs.batchArgs('Submit jobs to run chains or importance sample', notExist=True, converge=True)
Opts.parser.add_argument('--nodes', type=int, default=2)
Opts.parser.add_argument('--script', default='runMPI_HPCS.pl')

Opts.parser.add_argument('--dryrun', action='store_true')
Opts.parser.add_argument('--subitems', action='store_true')
Opts.parser.add_argument('--minimize', action='store_true')
Opts.parser.add_argument('--importance_minimize', action='store_true')
Opts.parser.add_argument('--minimize_failed', action='store_true')
Opts.parser.add_argument('--checkpoint_run', action='store_true')


(batch, args) = Opts.parseForBatch()


variant = ''
if args.importance_minimize:
    variant = '_minimize'
    if args.importance is None: args.importance = []
    args.nodes = 0

if args.minimize:
    args.noimportance = True
    variant = '_minimize'
    args.nodes = 0

if args.importance is None: args.noimportance = True


def submitJob(ini):
        command = 'perl ' + args.script + ' ' + ini + ' ' + str(args.nodes)
        if not args.dryrun:
            print 'Submitting...' + command
            os.system(command)
        else: print '...' + command

for jobItem in Opts.filteredBatchItems(wantSubItems=args.subitems):
        if (not args.notexist or (args.importance_minimize or args.minimize) and not jobItem.chainMinimumExists()
       or not args.importance_minimize and not jobItem.chainExists()) and (not args.minimize_failed or not jobItem.chainMinimumConverged()):
            if args.converge == 0 or jobItem.hasConvergeBetterThan(args.converge):
                if not args.checkpoint_run or jobItem.wantCheckpointContinue() and jobItem.notRunning(): submitJob(jobItem.iniFile(variant))
