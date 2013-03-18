# checking .minimum files, also the code I used to extract data for PLA

import shutil, os, sys, glob, batchJobArgs, paramNames, GetDistPlots, ResultObjs


renames = [['CamSpec CovSpec_v1_cuts', 'CamSpec'], ['CamSpec CovSpec_v1', 'CamSpec'], ['actspt_2013_01.clik', 'highL'], ['lowlike_v22.cldf.clik', 'lowLike'],
	 ['commander_v4.1_lm49.cldf.clik', 'lowl'], ['actspt_2013_01.cldf.clik', 'highL'],
        ['commander_v4.1_lm49.clik', 'lowl'], ['lensing_likelihood_v4_ref.clik_lensing', 'lensing'], ['lowlike_v222.clik', 'lowLike'], ['SNLS April_2012', 'SNLS']]

Opts = batchJobArgs.batchArgs('Make plots from getdist outputs', importance=True, converge=True)

(batch, args) = Opts.parseForBatch()

newdir = ''


def writeChain(outfile, lines):
        if not Debug:
            textFileHandle = open(outfile, 'w')
            textFileHandle.write("".join(lines))
            textFileHandle.close()


def UnBurnChain(fileName, frac=0.3):
        print 'Unburn: ' + fileName
        textFileHandle = open(fileName)
        outline = []
        for i, line in enumerate(textFileHandle):
            if i == 0:
                corr = len(line)
                burn = os.path.getsize(fileName) * frac / corr
            elif i > burn: outline.append(line)
        textFileHandle.close()
        return outline

Debug = False

def fixOtherFiles():
    testmin = 'testmin'
    testmin2 = 'testmin2'
    testmin3 = 'testmin3'

    for jobItem in Opts.filteredBatchItems(chainExist=True):
        s1 = jobItem.chainRoot + '.minimum'
        s2 = s1.replace('planck_', 'testmin_').replace('main', testmin)
        s3 = s1.replace('planck_', 'testmin_').replace('main', testmin2)
        s4 = s1.replace('planck_', 'testmin_').replace('main', testmin3)
        bfcontent = open(s1, 'r').read()[0:2500]
        tests = [s2, s3, s4]
        for i, test in enumerate(tests):
            if os.path.exists(test):
                bfcontent2 = open(test, 'r').read()[0:2500]
                if (bfcontent2 == bfcontent):
                    print jobItem.name, i
                    (sout, _) = s1.split('.')
                    (sin, _) = test.split('.')
                    shutil.copy(sin + '.bestfit_cl', sout + '.bestfit_cl')
                    shutil.copy(sin + '.minimum.inputparams', sout + '.minimum.inputparams')

def getBestMins():
    base = []
    testmin = 'testmin'
    testmin2 = 'testmin2'
    testmin3 = 'testmin3'

    for jobItem in Opts.filteredBatchItems(chainExist=True):
        s1 = jobItem.chainRoot + '.minimum'
        s2 = s1.replace('planck_', 'testmin_').replace('main', testmin)
        s3 = s1.replace('planck_', 'testmin_').replace('main', testmin2)
        s4 = s1.replace('planck_', 'testmin_').replace('main', testmin3)
        bf = ResultObjs.bestFit(s1)
        best = bf
        tests = [s2, s3, s4]
        bestf = s1
        for i, test in enumerate(tests):
            if os.path.exists(test):
                bf2 = ResultObjs.bestFit(test)
    #        print s2
    #        if bf.logLike-bf2.logLike>0.2: print jobItem.name, bf.logLike*2, bf2.logLike*2
                if best.logLike - bf2.logLike > 0.1:
                    if jobItem.paramtag == 'base': base.append(jobItem.name)
                    print i, jobItem.name, best.logLike * 2, bf2.logLike * 2
                    if jobItem.name != 'base_planck_lowl_lowLike_highL':
                        best = bf2
                        bestf = test
        if best is not bf:
            print jobItem.name + ' -> ' + bestf
        content = open(bestf, 'r').read()
        for rename in renames:
            content = content.replace(rename[0], rename[1])
        if not Debug: open(s1, 'w').write(content)
        print s1


def copyCleanedChains(newdir):
    for jobItem in Opts.filteredBatchItems(chainExist=True):
        besf_fname = jobItem.chainRoot.replace(jobItem.batchPath, newdir)
        (outdir, _) = os.path.split(besf_fname)
        if len(outdir) > 0 and not os.path.exists(outdir): os.makedirs(outdir)
        if not os.path.exists(outdir + os.sep + 'dist/'): os.makedirs(outdir + os.sep + 'dist/')
        for i in range(1, 9):
            if jobItem.chainExists(i):
                fname = jobItem.chainName(i)
                newName = fname.replace(jobItem.batchPath, newdir)
                if not os.path.exists(newName):
                    print newName
                    if not jobItem.isImportanceJob:
                        writeChain(newName, UnBurnChain(fname))
                    else:
                        if not Debug: shutil.copy(fname, newName)
        files = glob.iglob(jobItem.chainRoot + ".*")
        for afile in files:
                if os.path.isfile(afile):
                        print afile.replace(jobItem.batchPath, newdir)
                        if not Debug: shutil.copy(afile, afile.replace(jobItem.batchPath, newdir))

        files = glob.iglob(jobItem.distRoot + ".*")
        for afile in files:
            if os.path.isfile(afile):
                print afile.replace(jobItem.batchPath, newdir)
                if not Debug: shutil.copy(afile, afile.replace(jobItem.batchPath, newdir))


Debug = False
# fixOtherFiles()
newdir = '/scratch/aml1005/tmp/' + args.name[0] + os.sep
print newdir
if not os.path.exists(newdir): os.makedirs(newdir)
copyCleanedChains(newdir)


#        else: print 'OK: ' + jobItem.name
#        print s
