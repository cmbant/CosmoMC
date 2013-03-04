
import os, sys, pickle, ResultObjs, time


def readobject(directory=None):
    if directory == None:
        directory = sys.argv[1]
    with open(os.path.abspath(directory) + os.sep + 'batch.pyobj', 'rb') as inp:
        return pickle.load(inp)

def saveobject(obj, filename):
        with open(filename, 'wb') as output:
            pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

class jobItem:

    def __init__(self, path, param_set, data_set, base='base'):
        self.param_set = param_set
        self.data_set = data_set
        self.base = base
        paramtag = self.base
        for param in param_set:
            paramtag = paramtag + '_' + param
        self.dataname_set = data_set[0]
        datatag = "_".join(self.dataname_set)
        self.datatag = datatag
        self.paramtag = paramtag
        self.name = paramtag + '_' + datatag
        self.batchPath = path
        self.chainPath = path + paramtag + os.sep + datatag + os.sep
        self.chainRoot = self.chainPath + self.name
        self.distPath = self.chainPath + 'dist' + os.sep
        self.distRoot = self.distPath + self.name
        self.isImportanceJob = False
        self.importanceItems = []
        self.result_converge = None
        self.group = None
        self.makeIDs()

    def iniFile(self, variant=''):
        if not self.isImportanceJob:
            return self.batchPath + 'iniFiles' + os.sep + self.name + variant + '.ini'
        else: return self.batchPath + 'postIniFiles' + os.sep + self.name + variant + '.ini'

    def makeImportance(self, importanceRuns):
        self.importanceItems = []
        for (imp, ini, arr) in [(x[0], x[1], x) for x in importanceRuns]:
            if len(arr) > 2 and not arr[2].wantImportance(self): continue
            if len(set(imp).intersection(self.dataname_set)) > 0:
                print 'importance job duplicating parent data set:' + self.name
                continue
            job = jobItem(self.batchPath, self.param_set, self.data_set)
            job.importanceTag = "_".join(imp)
            job.importanceSettings = ini
            tag = '_post_' + job.importanceTag
            job.name = self.name + tag
            job.chainRoot = self.chainRoot + tag
            job.distRoot = self.distRoot + tag
            job.datatag = self.datatag + tag
            job.isImportanceJob = True
            job.parent = self
            job.dataname_set = imp + job.dataname_set
            job.group = self.group
            job.makeIDs()
            self.importanceItems.append(job)

    def makeIDs(self):
        self.normed_params = "_".join(sorted(self.param_set))
        self.normed_data = "_".join(sorted(self.dataname_set))
        self.normed_name = self.base
        if len(self.normed_params) > 0: self.normed_name += '_' + self.normed_params
        self.normed_name += '_' + self.normed_data


    def matchesDatatag(self, tagList):
        if self.datatag in tagList or self.normed_data in tagList: return True
        return self.datatag.replace('_post', '') in  [tag.replace('_post', '') for tag in tagList]

    def importanceJobs(self):
        return self.importanceItems

    def makeChainPath(self):
        if not os.path.exists(self.chainPath): os.makedirs(self.chainPath)
        return self.chainPath

    def writeIniLines(self, f):
        outfile = open(self.iniFile(), 'w')
        outfile.write("\n".join(f))
        outfile.close()

    def chainName(self, chain=1):
        return self.chainRoot + '_' + str(chain) + '.txt'

    def chainExists(self):
        fname = self.chainName()
        return os.path.exists(fname) and os.path.getsize(fname) > 0

    def chainFileDate(self, name, chain=1):
        return os.path.getmtime(self.chainName(chain))

    def chainsDodgy(self, interval=600):
        dates = []
        i = 1
        while os.path.exists(self.chainName(i)):
            dates.append(os.path.getmtime(self.chainName(i)))
            i += 1
        return os.path.exists(self.chainName(i + 1)) or max(dates) - min(dates) > interval

    def notRunning(self):
        if not self.chainExists(): return False  # might be in queue
        lastWrite = self.chainFileDate(self.chainName())
        return lastWrite < time.time() - 5 * 60

    def chainMinimumExists(self):
        fname = self.chainRoot + '.minimum'
        return os.path.exists(fname) and os.path.getsize(fname) > 0

    def chainMinimumConverged(self):
        fname = self.chainRoot + '.minimum'
        if not os.path.exists(fname) or os.path.getsize(fname) == 0: return False
        bf = ResultObjs.bestFit(fname)
        return bf.logLike < 1e29

    def convergeStat(self):
        fname = self.chainRoot + '.converge_stat'
        if not os.path.exists(fname): return None, None
        textFileHandle = open(fname)
        textFileLines = textFileHandle.readlines()
        textFileHandle.close()
        return float(textFileLines[0].strip()), len(textFileLines) > 1 and textFileLines[1].strip() == 'Done'

    def chainFinished(self):
        done = self.convergeStat()[1]
        if done is None: return False
        return done

    def wantCheckpointContinue(self):
        R, done = self.convergeStat()
        if R is None: return False
        if not os.path.exists(self.chainRoot + '_1.chk'): return False
        return not done

    def getDistExists(self):
        return os.path.exists(self.distRoot + '.margestats')

    def R(self):
        if self.result_converge is None:
            fname = self.distRoot + '.converge'
            if not os.path.exists(fname) or os.path.getsize(fname) == 0: return None
            self.result_converge = ResultObjs.convergeStats(fname)
        return float(self.result_converge.worstR())

    def hasConvergeBetterThan(self, R, returnNotExist=False):
        chainR = self.R()
        if chainR is None: return returnNotExist
        return chainR <= R

    def loadJobItemResults(self, paramNameFile=None, bestfit=True, bestfitonly=False, noconverge=False, silent=False):
        marge_root = self.distRoot
        self.result_converge = None
        self.result_marge = None
        self.result_bestfit = None
        bf_file = self.chainRoot + '.minimum'
        if os.path.exists(bf_file):
            self.result_bestfit = ResultObjs.bestFit(bf_file, paramNameFile)
        if not bestfitonly:
            if os.path.exists(marge_root + '.margestats'):
                if not noconverge: self.result_converge = ResultObjs.convergeStats(marge_root + '.converge')
                self.result_marge = ResultObjs.margeStats(marge_root + '.margestats', paramNameFile)
                if not self.result_bestfit is None and bestfit: self.result_marge.addBestFit(self.result_bestfit)
            elif not silent: print 'missing: ' + marge_root


class batchJob:

    def __init__(self, path):
        self.batchPath = path
        self.skip = []
        self.basePath = os.path.dirname(sys.path[0]) + os.sep
        self.commonPath = self.basePath + 'batch1/'
        self.subBatches = []
        self.jobItems = None

    def makeItems(self, dataAndParams):
            self.jobItems = []
            for group in dataAndParams:
                for data_set in group.datasets:
                    for param_set in group.params:
                        item = jobItem(self.batchPath, param_set, data_set)
                        if hasattr(group, 'groupName'): item.group = group.groupName
                        if not item.name in self.skip:
                            item.makeImportance(group.importanceRuns)
                            self.jobItems.append(item)
            for item in self.items():
                for x in [imp for imp in item.importanceJobs()]:
                    if self.has_normed_name(x.normed_name):
                        print 'replacing importance sampling run with full run: ' + x.name
                        item.importanceItems.remove(x)
            for item in self.items():
                for x in [imp for imp in item.importanceJobs()]:
                    if self.has_normed_name(x.normed_name, wantImportance=True, exclude=x):
                        print 'removing duplicate importance sampling run: ' + x.name
                        item.importanceItems.remove(x)


    def items(self, wantSubItems=True, wantImportance=False):
        for item in self.jobItems:
            yield(item)
            if wantImportance:
                for imp in item.importanceJobs():
                    if not imp.name in self.skip: yield(imp)

        if wantSubItems:
            for subBatch in self.subBatches:
                for item in subBatch.items(wantSubItems, wantImportance): yield(item)

    def hasName(self, name, wantSubItems=True):
        for jobItem in self.items(wantSubItems):
            if jobItem.name == name: return True
        return False

    def has_normed_name(self, name, wantSubItems=True, wantImportance=False, exclude=None):
        return self.normed_name_item(name, wantSubItems, wantImportance, exclude) is not None

    def normed_name_item(self, name, wantSubItems=True, wantImportance=False, exclude=None):
        for jobItem in self.items(wantSubItems, wantImportance):
            if jobItem.normed_name == name and not jobItem is exclude: return jobItem
        return None

    def normalizeDataTag(self, tag):
        return "_".join(sorted(tag.replace('_post', '').split('_')))

    def save(self, filename=''):
        saveobject(self, (self.batchPath + 'batch.pyobj', filename)[filename != ''])


    def makeDirectories(self):
            if not os.path.exists(self.batchPath):
                os.makedirs(self.batchPath)

            if not os.path.exists(self.batchPath + 'iniFiles'):
                os.makedirs(self.batchPath + 'iniFiles')

            if not os.path.exists(self.batchPath + 'postIniFiles'):
                os.makedirs(self.batchPath + 'postIniFiles')

