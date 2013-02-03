import batchJobArgs, ResultObjs


Opts = batchJobArgs.batchArgs('Compare parameter constraints over set of models')
Opts.parser.add_argument('--params', nargs='+', default=None)
Opts.parser.add_argument('--compare', nargs='+', default=None)
Opts.parser.add_argument('--nobestfits', action='store_true')
Opts.parser.add_argument('--single_extparam', action='store_true')
Opts.parser.add_argument('--sigma', type=int, default=2)

(batch, args) = Opts.parseForBatch()
formatter = ResultObjs.numberFormatter()

table = dict()

for jobItem in Opts.filteredBatchItems():
    if (args.compare is None or jobItem.matchesDatatag(args.compare)) and (not args.single_extparam or
                                len(jobItem.param_set) == 1  and len(set(args.params).intersection(jobItem.param_set)) > 0):
        jobItem.loadJobItemResults(paramNameFile=None, bestfit=not args.nobestfits, noconverge=True, silent=True)
        if jobItem.result_marge is not None:
            results = []
            for par in args.params:
                texValues = jobItem.result_marge.texValues(formatter, par, sigma=args.sigma)
                if texValues is not None:
                    if not jobItem.paramtag in table: table[jobItem.paramtag] = dict()
                    dataTable = table[jobItem.paramtag]
                    dataName = jobItem.datatag.replace('_post', '')
                    if not jobItem.datatag in dataTable: dataTable[dataName] = dict()
                    dataTable[dataName][par] = texValues


def sortData(batch, datatags):
    items = []
    for tag in datatags:
        items += [item for item in batch if item == tag]
    return items

for paramtag in table:
    print paramtag
    dataTable = table[paramtag]
    if args.compare is not None: dataTable = sortData(dataTable, args.compare)
    for datatag in dataTable:
        print ' ', table[paramtag][datatag], datatag


