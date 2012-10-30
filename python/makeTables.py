import os, sys, batchJob, ResultObjs

def checkDir(fname):
    if not os.path.exists(fname): os.makedirs(fname)


if len(sys.argv) < 3:
    print 'Usage: python/makeTables.py directory_with_outputs outfile'
    sys.exit()

outfile = sys.argv[2]

batchPath = os.path.abspath(sys.argv[1]) + os.sep
batch = batchJob.readobject(batchPath + 'batch.pyobj')

lines = []
lines.append('\\documentclass[11pt]{article}')
lines.append('\\begin{document}')

for jobItem in batch.items():
    caption = jobItem.name.replace('_', '-')
    bf = ResultObjs.bestFit(jobItem.chainRoot + '.minimum')
    lines.append(bf.resultTable(3, caption).tableTex())

lines.append('\\end{document}')

textFileHandle = open(outfile, 'w')
textFileHandle.write("\n".join(lines))
textFileHandle.close()

# m = margeStats('z:/base_nrun_r_planck_CAMspec_post_BAO.margestats')
