import os, sys, batchJob, ResultObjs


if len(sys.argv) < 3:
    print 'Usage: python/makeTables.py directory_with_outputs outfile [-bestfitonly, -bestfit]'
    sys.exit()

outfile = sys.argv[2]
if outfile.find('.') < 0: outfile += '.tex'

opt = ''
if len(sys.argv) > 3: opt = sys.argv[3]

batchPath = os.path.abspath(sys.argv[1]) + os.sep
batch = batchJob.readobject(batchPath + 'batch.pyobj')

# this is just for the latex labelsm set None to use those in chain .paramnames
paramNameFile = 'clik_latex.paramnames'

lines = []
lines.append('\\documentclass[11pt]{article}')
lines.append('\\usepackage{fullpage}')
lines.append('\\usepackage[landscape]{geometry}')
lines.append('\\renewcommand{\\arraystretch}{1.5}')
lines.append('\\begin{document}')

def addResultTable(caption, bf_file, marge_file):
    if bf_file != '' and os.path.exists(bf_file):
        bf = ResultObjs.bestFit(bf_file, paramNameFile)
        caption += ' (Best-fit $\\chi^2_{\\rm eff} = ' + ('%.2f' % bf.logLike) + '$)'
    else: bf = None
    if opt == '-bestfitonly':
        lines.append(bf.resultTable(3, caption).tableTex())
    else:
        if os.path.exists(marge_file):
            marge = ResultObjs.margeStats(marge_file, paramNameFile)
            if not bf is None and opt == '-bestfit': marge.addBestFit(bf)
            lines.append(marge.resultTable(3, caption).tableTex())
        else: print 'missing: ' + marge_file
    lines.append('\\newpage')



for jobItem in batch.items():
    caption = jobItem.name.replace('_', '-')
    bf_file = jobItem.chainRoot + '.minimum'
    marge_file = jobItem.distRoot() + '.margestats'
    addResultTable(caption, bf_file, marge_file)
    for imp in batch.importanceRuns:
        tag = '_post_' + imp
        bf_file = ''
        marge_File = jobItem.distRoot() + tag + '.margestats'
        caption = (jobItem.name + tag).replace('_', '-')
        addResultTable(caption, bf_file, marge_file)


lines.append('\\end{document}')

textFileHandle = open(outfile, 'w')
textFileHandle.write("\n".join(lines))
textFileHandle.close()

print 'Now converting to PDF...'
os.system('pdflatex ' + outfile)
