import os, batchJobArgs, paramNames, GetDistPlots


Opts = batchJobArgs.batchArgs('Make plots from getdist outputs', importance=True, converge=True)
Opts.parser.add_argument('out_dir')

Opts.parser.add_argument('--plot_data', default=None)
Opts.parser.add_argument('--paramNameFile', default='clik_latex.paramnames')
Opts.parser.add_argument('--paramList', default=None)
Opts.parser.add_argument('--compare_data', nargs='+', default=None)
Opts.parser.add_argument('--compare_importance', nargs='*', default=None)
Opts.parser.add_argument('--compare_paramtag', nargs='+', default=None)
Opts.parser.add_argument('--nx', default=None)
Opts.parser.add_argument('--legend_labels', default=None)
Opts.parser.add_argument('--outputs', nargs='+', default=['pdf'])

(batch, args) = Opts.parseForBatch()

if args.paramList is not None: args.paramList = paramNames.paramNames(args.paramList)

outdir = args.out_dir
if not os.path.exists(outdir): os.makedirs(outdir)
outdir = os.path.abspath(outdir) + os.sep

if args.plot_data is None: data = batch.batchPath + '/plot_data'
else: data = args.plot_data

g = GetDistPlots.GetDistPlotter(data)


def legendLabels(jobItems):
    return [jobItem.name for jobItem in jobItems]

def plot1d(roots, legend_labels=None):
    if legend_labels is None: legend_labels = args.legend_labels
    else: legend_labels = roots
    g.plots_1d(roots, paramList=args.paramList, nx=args.nx, legend_labels=legend_labels)

def comparePlot(jobItems, titles=None):
    roots = [jobItem.name for jobItem in jobItems]
    plot1d(roots, legend_labels=legendLabels(jobItems))


items = Opts.sortedParamtagDict()

for paramtag, parambatch in items:
    g.newPlot()
    if not args.compare_data is None:
        print 'comparing: ' + paramtag
        compares = Opts.filterForDataCompare(parambatch, args.compare_data)
        if len(compares) == 0:
            print '..None'
            continue
        else: comparePlot(compares)
        for ext in args.outputs: g.export(outdir + paramtag + '.' + ext)
    elif args.compare_importance is not None:
        for jobItem in parambatch:
            if not jobItem.isImportanceJob:
                print 'plotting: ' + jobItem.name
                roots = [jobItem.name]
                for imp in jobItem.importanceItems:
                    if len(args.compare_importance) == 0 or imp.importanceTag in args.compare_importance: roots.append(imp.name)
                plot1d(roots)
                for ext in args.outputs: g.export(outdir + jobItem.name + '.' + ext)
    elif args.compare_paramtag is not None:
            for jobItem in parambatch:
                if not jobItem.paramtag in args.compare_paramtag:
                    output = jobItem.name + '-vs-' + "-".join(args.compare_paramtag)
                    print 'plotting: ' + output
                    roots = [batch.normed_name_item(tag + '_' + jobItem.normed_data, wantImportance=True) for tag in args.compare_paramtag]
                    roots = [jobItem.name] + [root.name for root in roots if root is not None]
                    if len(roots) > 1:
                        plot1d(roots)
                        for ext in args.outputs: g.export(outdir + output + '.' + ext)
    else:
        for jobItem in parambatch:
            print 'plotting: ' + jobItem.name
            plot1d(jobItem.name)
            for ext in args.outputs: g.export(outdir + jobItem.name + '.' + ext)



