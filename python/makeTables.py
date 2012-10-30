import os, sys, batchJob, ResultObjs

def checkDir(fname):
    if not os.path.exists(fname): os.makedirs(fname)


if len(sys.argv) < 3:
    print 'Usage: python/makeTables.py directory_with_outputs outfile [-bestfit]'
    sys.exit()

outfile = sys.argv[2]

batchPath = os.path.abspath(sys.argv[1]) + os.sep
batch = batchJob.readobject(batchPath + 'batch.pyobj')

lines = []
lines.append('\\documentclass[11pt]{article}')
lines.append('\\begin{document}')
lines.append('\\usepackage{fullpage}')
lines.append('\\usepackage[landscape]{geometry}')
lines.append('\\renewcommand{\arraystretch}{1.5}')


for jobItem in batch.items():
    caption = jobItem.name.replace('_', '-')
    if len(sys.argv) > 3 and sys.argv[3] == '-bestfit':
        bf = ResultObjs.bestFit(jobItem.chainRoot + '.minimum')
        caption = '$\\chi^2_{\\rm eff} = ' + ('%.2f' % bf.logLike) + '$'
        lines.append(bf.resultTable(3, caption).tableTex())
    else:
        fname = jobItem.distPath + jobItem.name + '.margestats'
        if os.path.exists(fname):
            marge = ResultObjs.margeStats(fname)
            lines.append(marge.resultTable(3, caption).tableTex())
        else: print 'missing: ' + fname
    lines.append('\\newpage')

lines.append('\\end{document}')

textFileHandle = open(outfile, 'w')
textFileHandle.write("\n".join(lines))
textFileHandle.close()

