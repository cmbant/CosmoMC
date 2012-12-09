import os, batchJobArgs, ResultObjs


Opts = batchJobArgs.batchArgs('Make pdf tables from latex generated from getdist outputs', importance=True)
Opts.parser.add_argument('latex_filename')
Opts.parser.add_argument('--bestfitonly', action='store_true')
Opts.parser.add_argument('--bestfit', action='store_true', default=True)

# this is just for the latex labelsm set None to use those in chain .paramnames
Opts.parser.add_argument('--paramNameFile', default='clik_latex.paramnames')
Opts.parser.add_argument('--columns', type=int, nargs=1, default=3)


(batch, args) = Opts.parseForBatch()

outfile = args.latex_filename
if outfile.find('.') < 0: outfile += '.tex'

lines = []
lines.append('\\documentclass[10pt]{article}')
lines.append('\\usepackage{fullpage}')
lines.append('\\usepackage[pdftex]{hyperref}')

lines.append('\\usepackage[landscape,margin=0.8in]{geometry}')
lines.append('\\renewcommand{\\arraystretch}{1.5}')
lines.append('\\begin{document}')
lines.append('\\tableofcontents')

def texEscapeText(string):
    return string.replace('_', '{\\textunderscore}')

def paramResultTable(bf_file, marge_file):
    tableLines = []
    caption = ''
    if bf_file != '' and os.path.exists(bf_file + '.minimum'):
        bf = ResultObjs.bestFit(bf_file + '.minimum', args.paramNameFile)
        caption += ' Best-fit $\\chi^2_{\\rm eff} = ' + ('%.2f' % (bf.logLike * 2)) + '$'
    else: bf = None
    if args.bestfitonly:
        tableLines.append(bf.resultTable(args.columns).tableTex())
    else:
        if os.path.exists(marge_file + '.margestats'):
            converge = ResultObjs.convergeStats(marge_file + '.converge')
            caption += '; R-1 =' + converge.worstR()
            marge = ResultObjs.margeStats(marge_file + '.margestats', args.paramNameFile)
            if not bf is None and args.bestfit: marge.addBestFit(bf)
            tableLines.append(marge.resultTable(args.columns).tableTex())
        else: print 'missing: ' + marge_file
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

items = dict()
for jobItem in Opts.filteredBatchItems():
    if not jobItem.paramtag in items: items[jobItem.paramtag] = []
    items[jobItem.paramtag].append(jobItem)
items = sorted(items.iteritems())

lastTag = ''
for paramtag, parambatch in items:
    for jobItem in parambatch:
        if paramtag != lastTag:
            lastTag = paramtag
            lines.append('\\section{ ' + texEscapeText("+".join(jobItem.param_set)) + '}')
        lines.append('\\subsection{ ' + texEscapeText(jobItem.name) + '}')
        tableLines = paramResultTable(jobItem.chainRoot, jobItem.distRoot)
        ResultObjs.textFile(tableLines).write(jobItem.distRoot + '.tex')
        lines = lines + tableLines
        lines.append('\\newpage')

lines.append('\\end{document}')

ResultObjs.textFile(lines).write(outfile)

print 'Now converting to PDF...'
os.system('pdflatex ' + outfile)
# #again to get table of contents
os.system('pdflatex ' + outfile)

