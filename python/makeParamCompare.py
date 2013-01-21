import os, batchJobArgs


Opts = batchJobArgs.batchArgs('Compare parameter constraints over set of models')
Opts.parser.add_argument('--compare', nargs='+', default=None)

(batch, args) = Opts.parseForBatch()

for jobItem in Opts.filteredBatchItems():
    jobItem.loadJobItemResults(paramNameFile=None, bestfit=False, noconverge=True, silent=True)
    if jobItem.result_marge is not None:
        results = []
        for par in args.compare:
            param = jobItem.result_marge.parWithName(par)
            if param is not None:
                results += [par] + [param.mean]
                results.append(param.limits[1])

        if len(results) > 0: print results, jobItem.name


