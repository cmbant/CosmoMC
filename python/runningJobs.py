from __future__ import absolute_import
from __future__ import print_function
from paramgrid import batchjob_args, jobqueue

Opts = batchjob_args.batchArgs(
    'List details of running or queued jobs; gives job stats, then current R-1 and job/chain names', importance=True,
    batchPathOptional=True)

group = Opts.parser.add_mutually_exclusive_group()
group.add_argument('--queued', action='store_true')
group.add_argument('--running', action='store_true')

(batch, args) = Opts.parseForBatch()

if batch:
    items = [jobItem for jobItem in Opts.filteredBatchItems()]
    batchNames = set([jobItem.name for jobItem in items] + [jobItem.name + '_minimize' for jobItem in items])
else:
    batchNames= set()

ids, jobNames, nameslist, infos = jobqueue.queue_job_details(args.batchPath, running=not args.queued,
                                                             queued=not args.running)
for jobId, jobName, names, info in zip(ids, jobNames, nameslist, infos):
    if batchNames.intersection(set(names)):
        stats = dict()
        if batch:
            for name in names:
                for jobItem in items:
                    if jobItem.name == name:
                        R = jobItem.convergeStat()[0]
                        if R: stats[name] = "%6.3f" % R
                        break
        R = stats.get(jobName) or ' ' * 6
        print(info + ' |', R, jobName)
        if len(names) > 1:
            for name in names:
                R = stats.get(name) or ' ' * 6
                print('    >> ', R, name)

