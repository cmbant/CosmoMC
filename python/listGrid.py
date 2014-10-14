import batchJobArgs


Opts = batchJobArgs.batchArgs('List items in a grid', importance=True, converge=True, notExist=True)
Opts.parser.add_argument('--exists', action='store_true', help='chain must exist')


(batch, args) = Opts.parseForBatch()
items = Opts.sortedParamtagDict(chainExist=args.exists)

for paramtag, parambatch in items:
    for jobItem in parambatch:
        print jobItem.name