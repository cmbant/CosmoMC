import batchJobArgs, ResultObjs


Opts = batchJobArgs.batchArgs('Compare parameter constraints over set of models')
Opts.parser.add_argument('--params', nargs='+', default=None)
Opts.parser.add_argument('--compare', nargs='+', default=None)
Opts.parser.add_argument('--nobestfits', action='store_true')
Opts.parser.add_argument('--single_extparam', action='store_true')

(batch, args) = Opts.parseForBatch()
formatter = ResultObjs.numberFormatter()

table = dict()

for jobItem in Opts.filteredBatchItems():
    if (args.compare is None or jobItem.datatag in args.compare) and (not args.single_extparam or
                                len(jobItem.param_set) == 1  and len(set(args.params).intersection(jobItem.param_set)) > 0):
        jobItem.loadJobItemResults(paramNameFile=None, bestfit=not args.nobestfits, noconverge=True, silent=True)
        if jobItem.result_marge is not None:
            results = []
            for par in args.params:
                texValues = jobItem.result_marge.texValues(formatter, par)
                if texValues is not None:
                    if not jobItem.paramtag in table: table[jobItem.paramtag] = dict()
                    dataTable = table[jobItem.paramtag]
                    if not jobItem.datatag in dataTable: dataTable[jobItem.datatag] = dict()
                    dataTable[jobItem.datatag][par] = texValues
    #            param = jobItem.result_marge.parWithName(par)
    #            if param is not None:
    #                results += [par] + [param.mean]
    #                results.append(param.limits[1])

#            if len(results) > 0: print results, jobItem.name

for paramtag in table:
    if args.compare is None: print paramtag
    else: print paramtag, " ".join(table[paramtag])
    for datatag in table[paramtag]:
        if args.compare is None: name = datatag
        else: name = ''
        print ' ', table[paramtag][datatag], name


