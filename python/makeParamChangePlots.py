import os, batchJobArgs, ResultObjs, paramNames, planckStyle, GetDistPlots

from pylab import *


Opts = batchJobArgs.batchArgs('Compare how parameters means and errors vary', importance=True, converge=True)
Opts.parser.add_argument('fname', help='filename for the produced plots in')

# this is just for the latex labelsm set None to use those in chain .paramnames
Opts.parser.add_argument('--paramNameFile', default='clik_latex.paramnames', help=".paramnames file for custom labels for parameters")

Opts.parser.add_argument('--paramList', default=None, help=".paramnames file listing specific parameters to include (only)")
Opts.parser.add_argument('--compare', nargs='+', default=None)

Opts.parser.add_argument('--size_inch', type=float, default=None, help='output subplot size in inches')
Opts.parser.add_argument('--nx', default=None, help='number of plots per row')
Opts.parser.add_argument('--outputs', nargs='+', default=['pdf'], help='output file type (default: pdf)')
Opts.parser.add_argument('--bands_sigma', type=float, default=None)


(batch, args) = Opts.parseForBatch()

if args.paramList is not None: args.paramList = paramNames.paramNames(args.paramList)

g = GetDistPlots.GetDistPlotter('main/plot_data')
if args.size_inch is not None: g.settings.setWithSubplotSize(args.size_inch)


items = Opts.sortedParamtagDict(chainExist=True)

lmaxs = range(700, 2600, 150)
for paramtag, parambatch in items:
    datanames = []
    parnames = []
    for compi, comp in enumerate(args.compare):
        datanames = [comp.replace('XX', str(lmax)) for lmax in lmaxs]
        print datanames
        compares = Opts.filterForDataCompare(parambatch, datanames)
        if not compares:continue
        for jobItem in compares:
            jobItem.loadJobItemResults(paramNameFile=args.paramNameFile)
            if not parnames:
                parnames = jobItem.result_marge.names
                plot_col, plot_row = g.make_figure(len(parnames), nx=args.nx, xstretch=1.5)
            print jobItem.name

        for i, param in enumerate(parnames):
            if compares[0].result_marge.parWithName(param.name):
                subplot(plot_row, plot_col, i + 1)
                means = [jobItem.result_marge.parWithName(param.name).mean for jobItem in compares]
                if args.bands_sigma is not None:
                    sigmas = [jobItem.result_marge.parWithName(param.name).err for jobItem in compares]
                    delta = args.bands_sigma * array([sqrt(abs(sigma ** 2 - sigmas[-1] ** 2)) for sigma in sigmas])
                    fill_between(lmaxs, means[-1] + delta, means[-1] - delta, facecolor=[1, 1, 0.8], alpha=1, edgecolor=[1, 1, 0.8], lw=0)
                plot(lmaxs, means, **g.get_line_styles(compi))


#            g.set_ylabel(param)
    if parnames:
        for i, param in enumerate(parnames):
            subplot(plot_row, plot_col, i + 1)
            text(0.95, 0.8, '$' + param.label + '$', horizontalalignment='right', verticalalignment='center', transform=gca().transAxes, fontsize=18)
        g.finish_plot(g.default_legend_labels(None, args.compare), legend_ncol=2)
        for ext in args.outputs: g.export(paramtag + '_' + args.fname + '.' + ext)
