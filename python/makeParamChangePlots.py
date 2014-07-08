import os, batchJobArgs, ResultObjs, paramNames, planckStyle, GetDistPlots

from pylab import *


Opts = batchJobArgs.batchArgs('Compare how parameters means and errors vary', importance=True, converge=True)
Opts.parser.add_argument('fname', help='filename root for the produced plots')

# this is just for the latex labelsm set None to use those in chain .paramnames
Opts.parser.add_argument('--paramNameFile', default='clik_latex.paramnames', help=".paramnames file for custom labels for parameters")

Opts.parser.add_argument('--paramList', default=None, help=".paramnames file listing specific parameters to include (only)")
Opts.parser.add_argument('--compare', nargs='+', default=None)

Opts.parser.add_argument('--size_inch', type=float, default=None, help='output subplot size in inches')
Opts.parser.add_argument('--nx', default=None, help='number of plots per row')
Opts.parser.add_argument('--outputs', nargs='+', default=['pdf'], help='output file type (default: pdf)')
Opts.parser.add_argument('--bands_sigma', type=float, nargs='+', default=None)
Opts.parser.add_argument('--sigma_from_left', action='store_true')
Opts.parser.add_argument('--sigma_between', type=int, nargs=2, default=None,
                         help='zero-based indices of the compare items to use for plotting error bands. '
                         + ' e.g. 0 1 uses difference between sigmas^2 between data 1 at each L and data 0 at one end ')


(batch, args) = Opts.parseForBatch()

if args.paramList is not None: args.paramList = paramNames.paramNames(args.paramList)

g = GetDistPlots.GetDistPlotter('main/plot_data')
if args.size_inch is not None: g.settings.setWithSubplotSize(args.size_inch)


items = Opts.sortedParamtagDict(chainExist=True)

fill_colors = [[1, 1, 0.8], [1, 0.9, 0.6], [1, 0.8, 0.2]]

if args.sigma_from_left: compare_index = 0
else: compare_index = -1

lmaxs = range(700, 2501, 150)
# lmaxs = range(700, 2501, 75)

for paramtag, parambatch in items:
    print 'Doing paramtag: ' + paramtag + '...'
    datanames = []
    parnames = []
    par_sigmas = dict()
    clf()
    for parse in [True, False]:
        for compi, comp in enumerate(args.compare):
            datanames = [comp.replace('XX', str(lmax)) for lmax in lmaxs]
            print datanames
            compares = Opts.filterForDataCompare(parambatch, datanames, getDistExists=True)
            if not compares:continue
            if len(compares) != len(datanames):
                print 'not all chain results exist:' + comp
                print 'missing -->', [d for d in datanames if not d in [c.datatag for c in compares]]
                continue
            if parse:
                for jobItem in compares:
                        jobItem.loadJobItemResults(paramNameFile=args.paramNameFile)
                        if not parnames:
                            parnames = jobItem.result_marge.names
                        elif args.sigma_between:
                            parnames = [p for p in parnames if jobItem.result_marge.parWithName(p.name)]
            elif compi == 0: plot_col, plot_row = g.make_figure(len(parnames), nx=args.nx, xstretch=1.5)

            for i, param in enumerate(parnames):
                if compares[0].result_marge.parWithName(param.name):
                    if parse:
                        if args.bands_sigma is not None:
                            par_sigmas[ str(compi) + '_' + param.name] = [jobItem.result_marge.parWithName(param.name).err for jobItem in compares]
                    else:
                        means = [jobItem.result_marge.parWithName(param.name).mean for jobItem in compares]
                        subplot(plot_row, plot_col, i + 1)
                        if args.bands_sigma is not None and compi == 0:
                            if args.sigma_between is not None:
                                sigmas_ref = par_sigmas[str(args.sigma_between[0]) + '_' + param.name ]
                                sigmas = par_sigmas[str(args.sigma_between[1]) + '_' + param.name ]
                            else:
                                sigmas_ref = par_sigmas['0_' + param.name]
                                sigmas = sigmas_ref
                            for i in range(len(args.bands_sigma)):
                                delta = args.bands_sigma[-i - 1] * array([sqrt(abs(sigma ** 2 - sigmas_ref[compare_index] ** 2)) for sigma in sigmas])
                                c = fill_colors[i]
                                fill_between(lmaxs, means[compare_index] + delta, means[compare_index] - delta, facecolor=c, alpha=1, edgecolor=c, lw=0)
                        plot(lmaxs, means, **g.get_line_styles(compi))

    if parnames:
        for i, param in enumerate(parnames):
            subplot(plot_row, plot_col, i + 1)
            text(0.95, 0.8, '$' + param.label + '$', horizontalalignment='right', verticalalignment='center', transform=gca().transAxes, fontsize=18)
            xlim(lmaxs[0], lmaxs[-1])
        g.finish_plot(g.default_legend_labels(None, [paramtag + '_' + p for p in  args.compare]), legend_ncol=2)
        for ext in args.outputs: g.export(paramtag + '_' + args.fname + '.' + ext)
