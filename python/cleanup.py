import os, fnmatch, batchJobArgs

Opts = batchJobArgs.batchArgs('delete failed chains, files etc.', importance=True, converge=True)

Opts.parser.add_argument('--dist', action='store_true')
Opts.parser.add_argument('--ext', nargs='+', default=['*'])
Opts.parser.add_argument('--confirm', action='store_true')

(batch, args) = Opts.parseForBatch()


args.ext = ['.' + ext for ext in args.ext] + ['._*' + ext for ext in args.ext]
for jobItem in Opts.filteredBatchItems():
    if args.converge == 0 or not jobItem.hasConvergeBetterThan(args.converge, returnNotExist=True):
        dirs = [jobItem.chainPath]
        if args.dist: dirs = []
        dirs += [jobItem.distPath]
        for adir in dirs:
            for f in os.listdir(adir):
                for ext in args.ext:
                    if fnmatch.fnmatch(f, jobItem.name + ext):
                        fname = adir + f
                        if os.path.exists(fname):
                            print fname
                            if args. confirm: os.remove(fname)

if not args.confirm: print 'Files not actually deleted: add --confirm to delete'
