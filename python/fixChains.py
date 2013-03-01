import os, batchJobArgs

Opts = batchJobArgs.batchArgs('remove bad lines (old bug)')

(batch, args) = Opts.parseForBatch()

def writeChain(outfile, lines):
        textFileHandle = open(outfile, 'w')
        textFileHandle.write("".join(lines))
        textFileHandle.close()


def fixChain(fileName):

        textFileHandle = open(fileName)
        outline = []
        skip = 0
        hasSkip = False
        for i, line in enumerate(textFileHandle):
            if i == 0: corr = len(line)
            if skip != 0:skip -= 1
            if skip == 0:
                if len(line) != corr:
                    hasSkip = True
                    skip = 10
                else: outline.append(line)
        textFileHandle.close()
        if hasSkip: return outline
        else: return None


for jobItem in Opts.filteredBatchItems():
    if not jobItem.isImportanceJob and jobItem.chainExists():
        print jobItem.name
        for i in range(1, 8 + 1):
            name = jobItem.chainRoot + '_' + str(i) + '.txt'
            if os.path.exists(name):
                fix = fixChain(name)
                if fix is not None:
                    print 'fixing: ' + name
                    writeChain(name, fix)
