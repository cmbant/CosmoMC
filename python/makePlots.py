from __future__ import absolute_import
from __future__ import print_function
import os
from paramgrid import batchjob_args


Opts = batchjob_args.batchArgs('Make plots from getdist outputs', importance=True, converge=True, plots=True)
Opts.parser.add_argument('out_dir', help='directory to put the produced plots in')

Opts.parser.add_argument('--compare_data', nargs='+', default=None,
                         help='data tags to compare for each parameter combination (data1_data2)')
Opts.parser.add_argument('--compare_importance', nargs='*', default=None)
Opts.parser.add_argument('--compare_paramtag', nargs='+', default=None,
                         help='list of parameter tags to compare for each data combination (param1_param2)')
Opts.parser.add_argument('--compare_alldata', action='store_true',
                         help='compare all data combinations for each parameter combination')
Opts.parser.add_argument('--compare_replacing', nargs='*', default=None,
                         help='compare results replacing data combination x with y, z..(datavar1 datavar2 ...)')

Opts.parser.add_argument('--legend_labels', default=None, nargs='+',
                         help='labels to replace full chain names in legend')
Opts.parser.add_argument('--D2_param', default=None, help='x-parameter for 2D plots')
Opts.parser.add_argument('--D2_y_params', nargs='+', default=None, help='list of y parameter names for 2D plots')
Opts.parser.add_argument('--filled', action='store_true', help='for 2D plots, output filled contours')

Opts.parser.add_argument('--tri_params', nargs='+', default=None, help='list of parameters for triangle plots')

Opts.parser.add_argument('--legend_ncol', type=int, default=None, help='numnber of columns to draw legends')
Opts.parser.add_argument('--allhave', action='store_true', help='only include plots where all combinations exist')
Opts.parser.add_argument('--outtag', default=None, help='tag to add to output filenames to distinguish output')

(batch, args, g) = Opts.parseForBatch()

outdir = args.out_dir
if not os.path.exists(outdir): os.makedirs(outdir)
outdir = os.path.abspath(outdir) + os.sep


def doplot(roots):
    ncol = args.legend_ncol or (1, 2)[len(roots) > 2]
    if args.D2_param is not None:
        g.plots_2d(roots, args.D2_param, params2=args.D2_y_params, nx=args.nx, legend_labels=args.legend_labels,
                   filled=args.filled, legend_ncol=ncol)
    elif args.tri_params is not None:
        g.triangle_plot(roots, args.tri_params, legend_labels=args.legend_labels, filled_compare=args.filled)
    else:
        g.plots_1d(roots, paramList=args.paramList, nx=args.nx, legend_labels=args.legend_labels, legend_ncol=ncol)


def comparePlot(jobItems):
    doplot([jobItem.name for jobItem in jobItems])


tp = ''
if args.data is not None and len(args.data) == 1 and args.compare_alldata is not None: tp = '_' + args.data[0] + tp
if args.outtag is not None: tp = '_' + args.outtag + tp
if args.D2_param is not None: tp = '_' + args.D2_param + '_2D'
if args.tri_params is not None: tp = '_tri'

if args.compare_replacing is not None:
    if args.data is None:
        args.data = args.compare_replacing
    else:
        args.data += args.compare_replacing

items = Opts.sortedParamtagDict()

for paramtag, parambatch in items:
    g.newPlot()
    if args.compare_replacing is not None:
        print('comparing changes in data for: ' + paramtag)
        origCompare = [item for item in parambatch if args.compare_replacing[0] in item.data_set.names]
        if len(origCompare) == 0:
            print('..None')
            continue
        else:
            for jobItem in origCompare:
                compares = [jobItem.datatag]
                for replace in args.compare_replacing[1:]:
                    compares.append(
                        batch.normalizeDataTag(jobItem.data_set.tagReplacing(args.compare_replacing[0], replace)))
                compares = Opts.filterForDataCompare(parambatch, compares)
                if len(compares) == 1 or args.allhave and len(compares) != len(args.compare_replacing): continue
                print('comparing: ', [i.name for i in compares])
                comparePlot(compares)
                for ext in args.outputs: g.export(
                    outdir + jobItem.name + '-vs-' + "-".join(args.compare_replacing[1:]) + tp + '.' + ext)
    elif args.compare_data is not None or args.compare_alldata:
        print('comparing data combinations for: ' + paramtag)
        if args.compare_alldata:
            compares = parambatch
        else:
            compares = Opts.filterForDataCompare(parambatch, args.compare_data)
        if len(compares) == 0:
            print('..None')
            continue
        if not args.compare_alldata and args.allhave and len(compares) != len(args.compare_data):
            print('..not all, skipping')
            continue
        else:
            comparePlot(compares)
        for ext in args.outputs: g.export(outdir + paramtag + tp + '.' + ext)
    elif args.compare_importance is not None:
        for jobItem in parambatch:
            if not jobItem.isImportanceJob:
                print('plotting: ' + jobItem.name)
                roots = [jobItem.name]
                for imp in jobItem.importanceItems:
                    if len(args.compare_importance) == 0 or imp.importanceTag in args.compare_importance: roots.append(
                        imp.name)
                doplot(roots)
                for ext in args.outputs: g.export(outdir + jobItem.name + tp + '.' + ext)
    elif args.compare_paramtag is not None:
        for jobItem in parambatch:
            if not jobItem.paramtag in args.compare_paramtag:
                output = jobItem.name + '-vs-' + "-".join(args.compare_paramtag)
                print('plotting: ' + output)
                roots = [batch.normed_name_item(tag + '_' + jobItem.normed_data, wantImportance=True) for tag in
                         args.compare_paramtag]
                roots = [jobItem.name] + [root.name for root in roots if root is not None]
                if len(roots) > 1:
                    doplot(roots)
                    for ext in args.outputs: g.export(outdir + output + tp + '.' + ext)
    else:
        for jobItem in parambatch:
            print('plotting: ' + jobItem.name)
            doplot([jobItem.name])
            for ext in args.outputs: g.export(outdir + jobItem.name + tp + '.' + ext)



