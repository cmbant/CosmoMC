#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import hashlib
import os
from paramgrid import batchjob_args, jobqueue

Opts = batchjob_args.batchArgs('Submit jobs to run chains or importance sample', notExist=True, notall=True,
                               converge=True)

jobqueue.addArguments(Opts.parser, combinedJobs=True)

Opts.parser.add_argument('--subitems', action='store_true', help='include sub-grid items')
Opts.parser.add_argument('--not_queued', action='store_true')
Opts.parser.add_argument('--filters', action='store_true',
                         help='run any python importance filters on grid (no submission)')
Opts.parser.add_argument('--minimize', action='store_true', help='Run minimization jobs')
Opts.parser.add_argument('--importance_minimize', action='store_true',
                         help='Run minimization jobs for chains that are importance sampled')
Opts.parser.add_argument('--minimize_failed', action='store_true', help='run where minimization previously failed')
Opts.parser.add_argument('--checkpoint_run', nargs='?', default=None, const=0, type=float,
                         help='run if stopped and not finished; if optional value given then only run chains with convergence worse than the given value')
Opts.parser.add_argument('--importance_ready', action='store_true', help='where parent chain has converged and stopped')
Opts.parser.add_argument('--importance_changed', action='store_true',
                         help='run importance jobs where the parent chain has changed since last run')
Opts.parser.add_argument('--parent_converge', type=float, default=0,
                         help='minimum R-1 convergence for importance job parent')
Opts.parser.add_argument('--parent_stopped', action='store_true', help='only run if parent chain is not still running')
Opts.parser.add_argument('--chain_exists', action='store_true', help='Only run if chains already exist')

(batch, args) = Opts.parseForBatch()

if args.not_queued:
    print('Getting queued names...')
    queued = jobqueue.queue_job_names(args.batchPath, queued=True, running=True)


def notQueued(name):
    if args.minimize: name += '_minimize'
    return not name in queued


variant = ''
if args.importance_minimize:
    variant = '_minimize'
    if args.importance is None: args.importance = []

if args.minimize:
    args.noimportance = True
    variant = '_minimize'

if args.importance is None:
    if args.importance_changed or args.importance_ready:
        args.importance = []
    else:
        if args.filters:
            args.importance = []
        else:
            args.noimportance = True

isMinimize = args.importance_minimize or args.minimize

if args.combineOneJobName:
    print('Combining multiple (hopefully fast) into single job script: ' + args.combineOneJobName)

iniFiles = []

if not args.filters: jobqueue.checkArguments(**args.__dict__)


def jobName():
    s = "-".join([os.path.basename(ini) for ini in iniFiles])
    if len(iniFiles) < 2 or len(s) < 70: return s
    base = os.path.basename(iniFiles[0])
    if len(base) > 70: base = base[:70]
    return base + '__' + hashlib.md5(s).hexdigest()[:16]


def submitJob(ini):
    global iniFiles
    ini = ini.replace('.ini', '')
    if not args.dryrun:
        print('Submitting...' + ini)
    else:
        print('... ' + ini)
    iniFiles.append(ini)
    if args.combineOneJobName: return
    if len(iniFiles) >= args.runsPerJob:
        if args.runsPerJob > 1: print('--> jobName: ', jobName())
        jobqueue.submitJob(jobName(), iniFiles, **args.__dict__)
        iniFiles = []


for jobItem in Opts.filteredBatchItems(wantSubItems=args.subitems):
    if ((not args.notexist or isMinimize and not jobItem.chainMinimumExists()
         or not isMinimize and not jobItem.chainExists()) and (
                not args.minimize_failed or not jobItem.chainMinimumConverged())
        and (isMinimize or args.notall is None or not jobItem.allChainExists(args.notall))) \
            and (not isMinimize or getattr(jobItem, 'want_minimize', True)):
        if not args.parent_converge or not jobItem.isImportanceJob or jobItem.parent.hasConvergeBetterThan(
                args.parent_converge):
            if args.converge == 0 or not jobItem.hasConvergeBetterThan(args.converge, returnNotExist=True):
                if args.checkpoint_run is None or jobItem.wantCheckpointContinue(
                        args.checkpoint_run) and jobItem.notRunning():
                    if (not jobItem.isImportanceJob or isMinimize
                        or (args.importance_ready and jobItem.parent.chainFinished()
                            or not args.importance_ready and (args.filters or jobItem.parent.chainExists()))
                        and (not args.importance_changed or jobItem.parentChanged())
                        and (not args.parent_stopped or jobItem.parent.notRunning())) \
                            and (not args.chain_exists or jobItem.chainExists()):
                        if jobItem.isImportanceJob and hasattr(jobItem, 'importanceFilter'):
                            if args.filters:
                                if jobItem.parent.chainExists():
                                    print('Filtering for... %s' % jobItem.name)
                                    jobItem.importanceFilter.filter(batch, jobItem)
                                else:
                                    print("parent chains don't exist: %s" % (jobItem.name))
                        elif not args.filters:
                            if not args.not_queued or notQueued(jobItem.name):
                                submitJob(jobItem.iniFile(variant))

if len(iniFiles) > 0:
    if args.runsPerJob > 1: print('--> jobName: ', jobName())
    jobqueue.submitJob(args.combineOneJobName or jobName(), iniFiles, sequential=args.combineOneJobName is not None,
                       **args.__dict__)
