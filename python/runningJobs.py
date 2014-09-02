import jobQueue, argparse

parser = argparse.ArgumentParser(description='List details of running or queued jobs')

parser.add_argument('batchPath', help='directory containing a grid, or empty to see non-batch jobs', nargs='?')
group = parser.add_mutually_exclusive_group()
group.add_argument('--queued', action='store_true')
group.add_argument('--running', action='store_true')

args = parser.parse_args()

ids, jobNames, nameslist, infos = jobQueue.queue_job_details(args.batchPath, running=not args.queued, queued=not args.running)
for jobId, jobName, names, info in zip(ids, jobNames, nameslist, infos):
    print  info, jobName
    if len(names) > 1:
        for name in names:
            print ' >> ', name
