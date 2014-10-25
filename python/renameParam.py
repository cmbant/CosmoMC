import os, batchJobArgs, paramNames

Opts = batchJobArgs.batchArgs('rename parameter in all .paramnames files in grid', importance=True)
Opts.parser.add_argument('--old_new', nargs='+', help="list of oldname newname oldname2 newname2...")
Opts.parser.add_argument('--labelNames', default='clik_latex.paramnames', help=".paramnames file for new param labels")
Opts.parser.add_argument('--map_file', help="file with rows of oldname newname label")
Opts.parser.add_argument('--confirm', action='store_true', help="true to replace .paramnames files")


(batch, args) = Opts.parseForBatch()

if args.old_new and len(args.old_new) < 2: raise Exception('Must have at least one pair of parameters to rename')

if args.labelNames:
    labels = paramNames.paramNames(args.labelNames)
else:
    labels = None

if args.map_file:
    mapper = dict()
    with open(args.map_file) as f:
        for line in f:
            if line.strip():
                old, new, label = [s.strip() for s in line.split(None, 2)]
                mapper[old] = (new, label)
else:
    mapper = None


for jobItem in Opts.filteredBatchItems():
    name = jobItem.chainRoot + '.paramnames'
    if os.path.exists(name):
        names = paramNames.paramNames(name)
        has = False
        if mapper:
            for p in names.names:
                new = mapper.get(p.name, None)
                if new:
                    p.name, p.label = new
                    has = True
        if args.old_new:
            for old, new in zip(args.old_new[::2], args.old_new[1::2]):
                p = names.parWithName(old)
                if p:
                    has = True
                    p.name = new
        if labels:
            for p in names.names:
                plab = labels.parWithName(p.name)
                if plab:
                    has = True
                    p.label = plab.label
        if has:
            print jobItem.chainRoot
            if args.confirm: names.saveAsText(name)

if args.confirm:
    print 'Done. Re-run getdist to update getdist outputs.'
else:
    print '... run with --confirm to actually replace .paramname files'
