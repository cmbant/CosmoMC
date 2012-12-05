#!/usr/bin/env python
import os, batchJobArgs, sys

Opts = batchJobArgs.batchArgs('Submit jobs to run chains or importance sample', notExist=True)
Opts.parser.add_argument('--nodes', type=int, default=2)
Opts.parser.add_argument('--script', default='runMPI_HPCS.pl')

Opts.parser.add_argument('--dryrun', action='store_true')
Opts.parser.add_argument('--subitems', action='store_true')
Opts.parser.add_argument('--importance_minimize', action='store_true')


(batch, args) = Opts.parseForBatch()

variant = ''
if args.importance_minimize:
    variant = '_minimize'
    if args.importance is None:
        print 'can only minimize for importance runs'
        sys.exit()

def submitJob(ini):
        command = 'perl ' + args.script + ' ' + ini + ' ' + str(args.nodes)
        print 'Submitting...' + command
        if not args.dryrun: os.system(command)


for jobItem in Opts.filteredBatchItems(wantSubItems=args.subitems):
        if not args.notexist or not jobItem.chainExists():
            submitJob(jobItem.iniFile(variant))
