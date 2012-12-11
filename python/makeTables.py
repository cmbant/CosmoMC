import os, batchJobArgs, ResultObjs


Opts = batchJobArgs.batchArgs('Make pdf tables from latex generated from getdist outputs', importance=True)
Opts.parser.add_argument('latex_filename')
Opts.parser.add_argument('--bestfitonly', action='store_true')
Opts.parser.add_argument('--bestfit', action='store_true', default=True)

# this is just for the latex labelsm set None to use those in chain .paramnames
Opts.parser.add_argument('--paramNameFile', default='clik_latex.paramnames')
Opts.parser.add_argument('--columns', type=int, nargs=1, default=3)
Opts.parser.add_argument('--compare', nargs='+', default=None)
Opts.parser.add_argument('--portrait', action='store_true')


(batch, args) = Opts.parseForBatch()

outfile = args.latex_filename
if outfile.find('.') < 0: outfile += '.tex'

lines = []
lines.append('\\documentclass[10pt]{article}')
lines.append('\\usepackage{fullpage}')
lines.append('\\usepackage[pdftex]{hyperref}')
if args.portrait: lines.append('\\usepackage[paperheight=15in,margin=0.8in]{geometry}')
else: lines.append('\\usepackage[landscape,margin=0.8in]{geometry}')

lines.append('\\renewcommand{\\arraystretch}{1.5}')
lines.append('\\begin{document}')
lines.append('\\tableofcontents')

def texEscapeText(string):
    return string.replace('_', '{\\textunderscore}')

def loadJobItemResults(jobItem):
    bf_file = jobItem.chainRoot + '.minimum'
    marge_root = jobItem.distRoot
    jobItem.result_converge = None
    jobItem.result_marge = None
    jobItem.result_bestfit = None
    if os.path.exists(bf_file):
        jobItem.result_bestfit = ResultObjs.bestFit(bf_file, args.paramNameFile)
    if not args.bestfitonly:
        if os.path.exists(marge_root + '.margestats'):
            jobItem.result_converge = ResultObjs.convergeStats(marge_root + '.converge')
            jobItem.result_marge = ResultObjs.margeStats(marge_root + '.margestats', args.paramNameFile)
            if not jobItem.result_bestfit is None and args.bestfit: jobItem.result_marge.addBestFit(jobItem.result_bestfit)
        else: print 'missing: ' + marge_root


def paramResultTable(jobItem):
    tableLines = []
    caption = ''
    loadJobItemResults(jobItem)
    bf = jobItem.result_bestfit
    if not bf is None:
        caption += ' Best-fit $\\chi^2_{\\rm eff} = ' + ('%.2f' % (jobItem.result_bestfit.logLike * 2)) + '$'
    if args.bestfitonly:
        tableLines += ResultObjs.resultTable(args.columns, [bf]).lines
    else:
        if not jobItem.result_converge is None: caption += '; R-1 =' + jobItem.result_converge.worstR()
        if not jobItem.result_marge is None: tableLines += ResultObjs.resultTable(args.columns, [jobItem.result_marge]).lines
    tableLines.append('')
    tableLines.append(caption)
    if not bf is None:
        tableLines.append('')
        tableLines.append('$\chi^2_{\\rm eff}$:')
        for kind, vals in bf.sortedChiSquareds():
            tableLines.append(kind + ' - ')
            for (name, chisq) in vals:
                tableLines.append('  ' + texEscapeText(name) + ': ' + ('%.2f' % chisq) + ' ')
    return tableLines

def compareTable(jobItems):
    for jobItem in jobItems:
        loadJobItemResults(jobItem)
        print jobItem.name
    return ResultObjs.resultTable(1, [jobItem.result_marge for jobItem in jobItems if jobItem.result_marge is not None],
                                   titles=[jobItem.datatag for jobItem in jobItems if jobItem.result_marge is not None]).lines

def filterBatchData(batch, datatags):
    return [jobItem for jobItem in batch if jobItem.datatag in datatags]

items = dict()
for jobItem in Opts.filteredBatchItems():
    if not jobItem.paramtag in items: items[jobItem.paramtag] = []
    items[jobItem.paramtag].append(jobItem)
items = sorted(items.iteritems())

for paramtag, parambatch in items:
    lines.append('\\section{ ' + texEscapeText("+".join(parambatch[0].param_set)) + '}')
    if not args.compare is None:
        compares = filterBatchData(parambatch, args.compare)
        if len(compares) > 0: lines += compareTable(compares)
        else: print 'no matches for compare'
    else:
        for jobItem in parambatch:
            lines.append('\\subsection{ ' + texEscapeText(jobItem.name) + '}')
            tableLines = paramResultTable(jobItem)
            ResultObjs.textFile(tableLines).write(jobItem.distRoot + '.tex')
            lines += tableLines
    lines.append('\\newpage')

lines.append('\\end{document}')

ResultObjs.textFile(lines).write(outfile)

print 'Now converting to PDF...'
os.system('pdflatex ' + outfile)
# #again to get table of contents
os.system('pdflatex ' + outfile)

