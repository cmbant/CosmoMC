import os , batchJobArgs, ResultObjs, paramNames


Opts = batchJobArgs.batchArgs('Compare parameter constraints over set of models')
Opts.parser.add_argument('--params', nargs='+')
Opts.parser.add_argument('--compare', nargs='+', default=None)
Opts.parser.add_argument('--nobestfits', action='store_true')
Opts.parser.add_argument('--single_extparam', action='store_true')
Opts.parser.add_argument('--sigma', type=int, default=2)
Opts.parser.add_argument('--latex_filename', default=None)
Opts.parser.add_argument('--paramNameFile', default='clik_latex.paramnames')

(batch, args) = Opts.parseForBatch()
formatter = ResultObjs.numberFormatter()

names = paramNames.paramNames(args.paramNameFile)

table = dict()
paramtag_for_param = dict()
for par in args.params:
    paramtag_for_param[par] = []

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
                    if not jobItem.paramtag in paramtag_for_param[par]: paramtag_for_param[par].append(jobItem.paramtag)
                    dataName = jobItem.datatag.replace('_post', '')
                    if not jobItem.datatag in dataTable: dataTable[dataName] = dict()
                    dataTable[dataName][par] = texValues

def textAsColumn(txt, latex=True, bold=False):
        wid = len(txt)
        if latex:
            wid += 2
            if bold: wid += 11
        res = txt + ' ' * max(0, 28 - wid)
        if latex:
            res = '$' + res + '$'
        return res


def sortData(batch, datatags):
    items = []
    for tag in datatags:
        items += [item for item in batch if item == tag]
    return items

lines = []

for i, par in enumerate(args.params):
    for paramtag in paramtag_for_param[par]:
        print paramtag
        dataTable = table[paramtag]
        if args.compare is not None: dataTable = sortData(dataTable, args.compare)
        cols = ['$' + names.parWithName(par, True).label + '$']
        for datatag in dataTable:
            if args.latex_filename is not None:
                for par in args.params:
                    if par in table[paramtag][datatag]:
                        res = table[paramtag][datatag][par]
                        if len(res) > 1: cols.append(textAsColumn(res[1]))
                        cols.append(textAsColumn(res[0], latex=res[0] != '---'))
            else: print ' ', table[paramtag][datatag], datatag
        lines.append(" & ".join(cols) + '\\\\')


if args.latex_filename is not None:
        if args.latex_filename .find('.') < 0: args.latex_filename += '.tex'
        (outdir, outname) = os.path.split(args.latex_filename)
        if not os.path.exists(outdir): os.makedirs(os.path.dirname(outdir))
        ResultObjs.textFile(lines).write(args.latex_filename)
