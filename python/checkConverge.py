import os, fnmatch, batchJobArgs

Opts = batchJobArgs.batchArgs('Find chains which have failed or not converged.', importance=True, converge=True)

Opts.parser.add_argument('--exist', action='store_true')

(batch, args) = Opts.parseForBatch()

notExist = []
converge = []

for jobItem in Opts.filteredBatchItems():
    if not jobItem.chainExists():
        notExist.append(jobItem)
    elif args.converge != 0 and not jobItem.hasConvergeBetterThan(args.converge, returnNotExist=True):
        converge.append(jobItem)

print 'Checking batch:'
if not args.exist and len(notExist) > 0:
    print 'Not exist...'
    for jobItem in notExist:
        print '...', jobItem.chainRoot

print 'Converge check...'
for jobItem in converge:
    print '...', jobItem.chainRoot, jobItem.R()


