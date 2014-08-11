import jobQueue, batchJobArgs

Opts = batchJobArgs.batchArgs('Delete running or queued jobs', importance=True, batchPathOptional=False)

group = Opts.parser.add_mutually_exclusive_group()
group.add_argument('--queued', action='store_true')
group.add_argument('--running', action='store_true')

Opts.parser.add_argument('--delete_id_min', type=int)
Opts.parser.add_argument('--delete_id_range', nargs=2, type=int)
Opts.parser.add_argument('--delete_ids', nargs='+', type=int)

Opts.parser.add_argument('--confirm', action='store_true')


(batch, args) = Opts.parseForBatch()


if args.delete_id_range is not None:
    jobQueue.deleteJobs(args.batchPath, jobId_minmax=args.delete_id_range, confirm=args.confirm)
if args.delete_id_min is not None:
    jobQueue.deleteJobs(args.batchPath, jobId_min=args.delete_id_min, confirm=args.confirm)
elif args.delete_ids is not None:
    jobQueue.deleteJobs(args.batchPath, args.delete_ids, confirm=args.confirm)
else:
    items = [jobItem for jobItem in Opts.filteredBatchItems()]
    batchNames = set([jobItem.name for jobItem in items])
    jobQueue.deleteJobs(args.batchPath, rootNames=batchNames, confirm=args.confirm)

if not args.confirm: print 'jobs not actually deleted: add --confirm to really cancel them'
