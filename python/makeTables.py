import os, batchJobArgs, ResultObjs, paramNames, planckStyle


Opts = batchJobArgs.batchArgs('Make pdf tables from latex generated from getdist outputs', importance=True, converge=True)
Opts.parser.add_argument('latex_filename')
Opts.parser.add_argument('--limit', type=int, default=2)
Opts.parser.add_argument('--bestfitonly', action='store_true')
Opts.parser.add_argument('--nobestfit', action='store_true')

# this is just for the latex labelsm set None to use those in chain .paramnames
Opts.parser.add_argument('--paramNameFile', default='clik_latex.paramnames')
Opts.parser.add_argument('--paramList', default=None)
Opts.parser.add_argument('--blockEndParams', default=None)
Opts.parser.add_argument('--columns', type=int, nargs=1, default=3)
Opts.parser.add_argument('--compare', nargs='+', default=None)
Opts.parser.add_argument('--titles', default=None)  # for compare plots
Opts.parser.add_argument('--forpaper', action='store_true')
Opts.parser.add_argument('--separate_tex', action='store_true')
Opts.parser.add_argument('--header_tex', default=None)
Opts.parser.add_argument('--height', default="8in")
Opts.parser.add_argument('--width', default="10in")


(batch, args) = Opts.parseForBatch()

if args.blockEndParams is not None: args.blockEndParams = args.blockEndParams.split(';')
outfile = args.latex_filename

if args.paramList is not None: args.paramList = paramNames.paramNames(args.paramList)

if args.forpaper: formatter = planckStyle.planckStyleTableFormatter()
else: formatter = None

lines = []
if not args.forpaper:
    lines.append('\\documentclass[10pt]{article}')
    lines.append('\\usepackage{fullpage}')
    lines.append('\\usepackage[pdftex]{hyperref}')
    lines.append('\\usepackage[paperheight=' + args.height + ',paperwidth=' + args.width + ',margin=0.8in]{geometry}')
    lines.append('\\renewcommand{\\arraystretch}{1.5}')
    lines.append('\\begin{document}')
    if args.header_tex is not None:
        lines.append(open(args.header_tex, 'r').read())
    lines.append('\\tableofcontents')

def texEscapeText(string):
    return string.replace('_', '{\\textunderscore}')

def getTableLines(content):
    return ResultObjs.resultTable(args.columns, [content], blockEndParams=args.blockEndParams,
                         formatter=formatter, paramList=args.paramList, limit=args.limit).lines

def paramResultTable(jobItem):
    tableLines = []
    caption = ''
    jobItem.loadJobItemResults(paramNameFile=args.paramNameFile, bestfit=not args.nobestfit, bestfitonly=args.bestfitonly)
    bf = jobItem.result_bestfit
    if not bf is None:
        caption += ' Best-fit $\\chi^2_{\\rm eff} = ' + ('%.2f' % (bf.logLike * 2)) + '$'
    if args.bestfitonly:
        if bf is not None: tableLines += getTableLines(bf)
    else:
        if jobItem.result_converge is not None: caption += '; R-1 =' + jobItem.result_converge.worstR()
        if jobItem.result_marge is not None: tableLines += getTableLines(jobItem.result_marge)
    tableLines.append('')
    if not args.forpaper: tableLines.append(caption)
    if not bf is None and not args.forpaper:
        tableLines.append('')
        tableLines.append('$\chi^2_{\\rm eff}$:')
        for kind, vals in bf.sortedChiSquareds():
            tableLines.append(kind + ' - ')
            for (name, chisq) in vals:
                tableLines.append('  ' + texEscapeText(name) + ': ' + ('%.2f' % chisq) + ' ')
    return tableLines

def compareTable(jobItems, titles=None):
    for jobItem in jobItems:
        jobItem.loadJobItemResults(paramNameFile=args.paramNameFile, bestfit=not args.nobestfit, bestfitonly=args.bestfitonly)
        print jobItem.name
    if titles is None: titles = [jobItem.datatag for jobItem in jobItems if jobItem.result_marge is not None]
    else: titles = titles.split(';')
    return ResultObjs.resultTable(1, [jobItem.result_marge for jobItem in jobItems if jobItem.result_marge is not None],
               formatter=formatter, limit=args.limit, titles=titles, blockEndParams=args.blockEndParams, paramList=args.paramList).lines


items = Opts.sortedParamtagDict(chainExist=not args.bestfitonly)

for paramtag, parambatch in items:
    if not args.forpaper:
        if len(parambatch[0].param_set) == 0: paramText = 'Baseline model'
        else: paramText = texEscapeText("+".join(parambatch[0].param_set))
        section = '\\newpage\\section{ ' + paramText + '}'
    else: section = ''
    if not args.compare is None:
        compares = Opts.filterForDataCompare(parambatch, args.compare)
        if len(compares) == len(args.compare):
            lines.append(section)
            lines += compareTable(compares, args.titles)
        else: print 'no matches for compare: ' + paramtag
    else:
        lines.append(section)
        for jobItem in parambatch:
            if (os.path.exists(jobItem.distPath) or args.bestfitonly) and (args.converge == 0 or jobItem.hasConvergeBetterThan(args.converge)):
                if not args.forpaper: lines.append('\\subsection{ ' + texEscapeText(jobItem.name) + '}')
                tableLines = paramResultTable(jobItem)
                if args.separate_tex: ResultObjs.textFile(tableLines).write(jobItem.distRoot + '.tex')
                lines += tableLines

if not args.forpaper: lines.append('\\end{document}')

if outfile.find('.') < 0: outfile += '.tex'
(outdir, outname) = os.path.split(outfile)
if len(outdir) > 0 and not os.path.exists(outdir): os.makedirs(outdir)
ResultObjs.textFile(lines).write(outfile)
(root, _) = os.path.splitext(outfile)

if not args.forpaper:
    print 'Now converting to PDF...'
    delext = ['aux', 'log', 'out', 'toc']
    if len(outdir)>0: dodir= 'cd ' + outdir + '; '
    else: dodir='';
    os.system(dodir+'pdflatex ' + outname)
    # #again to get table of contents
    os.system(dodir + 'pdflatex ' + outname)
    # and again to get page numbers
    os.system(dodir +'pdflatex ' + outname)
    for ext in delext:
        if os.path.exists(root + '.' + ext):
            os.remove(root + '.' + ext)


