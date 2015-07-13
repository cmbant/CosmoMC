from __future__ import absolute_import
from __future__ import print_function
import os
from paramgrid import batchjob_args
from getdist import types, paramnames


Opts = batchjob_args.batchArgs('Compare parameter constraints over set of models')
Opts.parser.add_argument('--params', nargs='+')
Opts.parser.add_argument('--chain_name_params', nargs='+')

Opts.parser.add_argument('--compare', nargs='+', default=None)
Opts.parser.add_argument('--nobestfits', action='store_true')
Opts.parser.add_argument('--single_extparam', action='store_true')
Opts.parser.add_argument('--limit', type=int, default=2)
Opts.parser.add_argument('--latex_filename', default=None)
Opts.parser.add_argument('--mathColumns', action='store_true')
Opts.parser.add_argument('--endline', default='\\cr')
Opts.parser.add_argument('--paramNameFile', default='clik_latex.paramnames')

(batch, args) = Opts.parseForBatch()
formatter = types.TableFormatter()

names = paramnames.ParamNames(args.paramNameFile)

if args.chain_name_params is None: args.chain_name_params = args.params

if args.compare: args.compare = [batch.normalizeDataTag(dat) for dat in args.compare]

table = dict()
paramtag_for_param = dict()
for par in args.params:
    paramtag_for_param[par] = []

for jobItem in Opts.filteredBatchItems():
    if (args.compare is None or jobItem.matchesDatatag(args.compare)) and (not args.single_extparam or
                                                                                       len(
                                                                                               jobItem.param_set) == 1 and jobItem.hasParam(
                                                                                       args.chain_name_params)):
        jobItem.loadJobItemResults(paramNameFile=None, bestfit=not args.nobestfits, noconverge=True, silent=True)
        if jobItem.result_marge is not None:
            results = []
            for par in args.params:
                texValues = jobItem.result_marge.texValues(formatter, par, limit=args.limit)
                if texValues is not None:
                    if not jobItem.paramtag in table: table[jobItem.paramtag] = dict()
                    dataTable = table[jobItem.paramtag]
                    if not jobItem.paramtag in paramtag_for_param[par]: paramtag_for_param[par].append(jobItem.paramtag)
                    if not jobItem.normed_data in dataTable: dataTable[jobItem.normed_data] = dict()
                    dataTable[jobItem.normed_data][par] = texValues


def makeMath(txt):
    if args.mathColumns:
        return '$' + txt + '$'
    else:
        return txt


def textAsColumn(txt, latex=True):
    if latex:
        return makeMath(txt)
    else:
        return txt


def sortData(batch, datatags):
    items = []
    for tag in datatags:
        items += [item for item in batch if item == tag]
    return items


lines = []

for i, par in enumerate(args.params):
    for paramtag in paramtag_for_param[par]:
        dataTable = table[paramtag]
        print(paramtag, len(dataTable))
        if args.compare is not None: dataTable = sortData(dataTable, args.compare)
        cols = [makeMath(names.parWithName(par, True).label)]
        for datatag in dataTable:
            if args.latex_filename is not None:
                for par2 in args.params:
                    if par2 in table[paramtag][datatag]:
                        res = table[paramtag][datatag][par2]
                        if len(res) > 1: cols.append(textAsColumn(res[1]))
                        cols.append(textAsColumn(res[0], latex=res[0] != '---'))
            else:
                print(' ', table[paramtag][datatag], datatag)
        lines.append(" & ".join(cols) + args.endline)

if args.latex_filename is not None:
    if args.latex_filename.find('.') < 0: args.latex_filename += '.tex'
    (outdir, outname) = os.path.split(args.latex_filename)
    if not os.path.exists(outdir): os.makedirs(os.path.dirname(outdir))
    types.TextFile(lines).write(args.latex_filename)
