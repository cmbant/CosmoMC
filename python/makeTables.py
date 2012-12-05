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
lines.append('\\usepackage[landscape,margin=0.8in]{geometry}')
lines.append('\\renewcommand{\\arraystretch}{1.5}')
lines.append('\\begin{document}')

def texEscapeText(string):
    return string.replace('_', '{\\textunderscore}')

def addResultTable(caption, bf_file, marge_file):
    if bf_file != '' and os.path.exists(bf_file + '.minimum'):
        bf = ResultObjs.bestFit(bf_file + '.minimum', args.paramNameFile)
        caption += ' (Best-fit $\\chi^2_{\\rm eff} = ' + ('%.2f' % (bf.logLike * 2)) + '$)'
    else: bf = None
    if args.bestfitonly:
        lines.append(bf.resultTable(args.columns, caption).tableTex())
    else:
        if os.path.exists(marge_file + '.margestats'):
            converge = ResultObjs.convergeStats(marge_file + '.converge')
            caption += '; R-1 =' + converge.worstR()
            marge = ResultObjs.margeStats(marge_file + '.margestats', args.paramNameFile)
            if not bf is None and args.bestfit: marge.addBestFit(bf)
            lines.append(marge.resultTable(args.columns, caption).tableTex())
        else: print 'missing: ' + marge_file
    if not bf is None:
        lines.append('$\chi^2_{\\rm eff}$:')
        for kind, vals in bf.sortedChiSquareds():
            lines.append(kind + ' - ')
            for (name, chisq) in vals:
                lines.append('  ' + texEscapeText(name) + ': ' + ('%.2f' % chisq) + ' ')
    lines.append('\\newpage')


for jobItem in Opts.filteredBatchItems():
        lines.append('%' + jobItem.name)
        caption = texEscapeText(jobItem.name)
        addResultTable(caption, jobItem.chainRoot, jobItem.distRoot)

lines.append('\\end{document}')

textFileHandle = open(outfile, 'w')
textFileHandle.write("\n".join(lines))
textFileHandle.close()

print 'Now converting to PDF...'
os.system('pdflatex ' + outfile)
