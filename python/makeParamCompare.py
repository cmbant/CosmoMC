import os, batchJobArgs


Opts = batchJobArgs.batchArgs('Compare parameter constraints over set of models')
Opts.parser.add_argument('--compare', nargs='+', default=None)

(batch, args) = Opts.parseForBatch()

for jobItem in Opts.filteredBatchItems():
    jobItem.loadJobItemResults(paramNameFile=None, bestfit=True)
    m = jobItem.result_marge
    results = []
    for par in args.compare:
        param = m.parWithName(par)
        if param is not None:
            results += [par] + [param.mean]
            results.append(param.limits[1])

    print jobItem.name, results


