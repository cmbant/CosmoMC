#!/usr/bin/env python
import hashlib, os, batchJobArgs, jobQueue

Opts = batchJobArgs.batchArgs('Submit jobs to run chains or importance sample', notExist=True, notall=True, converge=True)

jobQueue.addArguments(Opts.parser, combinedJobs=True)

Opts.parser.add_argument('--subitems', action='store_true')
Opts.parser.add_argument('--minimize', action='store_true')
Opts.parser.add_argument('--importance_minimize', action='store_true')
Opts.parser.add_argument('--minimize_failed', action='store_true')
Opts.parser.add_argument('--checkpoint_run', action='store_true')
Opts.parser.add_argument('--importance_ready', action='store_true')
Opts.parser.add_argument('--not_queued', action='store_true')


(batch, args) = Opts.parseForBatch()

if args.not_queued:
    print 'Getting queued names...'
    queued = jobQueue.queue_job_names(args.batchPath)

def notQueued(name):
    for job in queued:
        if name in job:
#            print 'Already running:', name
            return False
    return True


variant = ''
if args.importance_minimize:
    variant = '_minimize'
    if args.importance is None: args.importance = []

if args.minimize:
    args.noimportance = True
    variant = '_minimize'

if args.importance is None: args.noimportance = True

isMinimize = args.importance_minimize or args.minimize

if args.combineOneJobName:
    print 'Combining multiple (hopefully fast) into single job script: ' + args.combineOneJobName

iniFiles = []

jobQueue.checkArguments(**args.__dict__)

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
            print 'Submitting...' + ini
        else: print '... ' + ini
        iniFiles.append(ini)
        if args.combineOneJobName: return
        if len(iniFiles) >= args.runsPerJob:
            if args.runsPerJob > 1: print '--> jobName: ', jobName()
            jobQueue.submitJob(jobName(), iniFiles, **args.__dict__)
            iniFiles = []

for jobItem in Opts.filteredBatchItems(wantSubItems=args.subitems):
    if ((not args.notexist or isMinimize and not jobItem.chainMinimumExists()
       or not isMinimize and not jobItem.chainExists()) and (not args.minimize_failed or not jobItem.chainMinimumConverged())
          and (isMinimize or args.notall is None or not jobItem.allChainExists(args.notall))):
            if args.converge == 0 or jobItem.hasConvergeBetterThan(args.converge):
                if not args.checkpoint_run or jobItem.wantCheckpointContinue() and jobItem.notRunning():
                    if not jobItem.isImportanceJob or  (args.importance_ready and jobItem.parent.chainFinished()
                                                        or not args.importance_ready and jobItem.parent.chainExists()):
                        if not args.not_queued or notQueued(jobItem.name):
                            submitJob(jobItem.iniFile(variant))


if len(iniFiles) > 0:
    if args.runsPerJob > 1: print '--> jobName: ', jobName()
    jobQueue.submitJob(args.combineOneJobName or jobName(), iniFiles, sequential=args.combineOneJobName is not None, **args.__dict__)
