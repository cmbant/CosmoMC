
import os, sys, pickle


def readobject(directory=None):
    if directory == None:
        directory = sys.argv[1]
    with open(os.path.abspath(directory) + os.sep + 'batch.pyobj', 'rb') as inp:
        return pickle.load(inp)

def saveobject(obj, filename):
        with open(filename, 'wb') as output:
            pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

class jobItem:

    def __init__(self, path, param_set, data_set):
        self.param_set = param_set
        self.data_set = data_set
        paramtag = 'base'
        for param in param_set:
            paramtag = paramtag + '_' + param
        datatag = "_".join(data_set[0])
        self.datatag = datatag
        self.paramtag = paramtag
        self.name = paramtag + '_' + datatag
        self.batchPath = path
        self.chainPath = path + paramtag + '/' + datatag + '/'
        self.chainRoot = self.chainPath + self.name
        self.distPath = self.chainPath + 'dist/'
        self.distRoot = self.distPath + self.name
        self.isImportanceJob = False
        self.importanceItems = []

    def iniFile(self, variant=''):
        if not self.isImportanceJob:
            return self.batchPath + 'iniFiles/' + self.name + variant + '.ini'
        else: return self.batchPath + 'postIniFiles/' + self.name + variant + '.ini'

    def makeImportance(self, importanceRuns):
        self.importanceItems = []
        for (imp, ini) in [(x[0], x[1]) for x in importanceRuns]:
            job = jobItem(self.batchPath, self.param_set, self.data_set)
            job.importanceTag = imp
            job.importaceSettings = ini
            tag = '_post_' + imp
            job.name = self.name + tag
            job.chainRoot = self.chainRoot + tag
            job.distRoot = self.distRoot + tag
            job.datatag = self.datatag + tag
            job.isImportanceJob = True
            job.parent = self
            self.importanceItems.append(job)

    def importanceJobs(self):
        return self.importanceItems

    def makeChainPath(self):
        if not os.path.exists(self.chainPath): os.makedirs(self.chainPath)
        return self.chainPath

    def writeIniLines(self, f):
        outfile = open(self.iniFile(), 'w')
        outfile.write("\n".join(f))
        outfile.close()

    def chainExists(self):
        fname = self.chainRoot + '_1.txt'
        return os.path.exists(fname) and os.path.getsize(fname) > 0

    def getDistExists(self):
        return os.path.exists(self.distRoot + '.margestats')


class batchJob:

    def __init__(self, path):
        self.batchPath = path
        self.extparams = []
        self.datasets = []
        self.skip = []
        self.basePath = os.path.dirname(sys.path[0]) + os.sep
        self.commonPath = self.basePath + 'batch1/'
        self.subBatches = []
        self.jobItems = None

    def makeItems(self, importanceRuns=[]):
            self.jobItems = []
            for data_set in self.datasets:
                for param_set in self.extparams:
                    item = jobItem(self.batchPath, param_set, data_set)
                    if not item.name in self.skip:
                        item.makeImportance(importanceRuns)
                        self.jobItems.append(item)


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

    def save(self, filename=''):
        saveobject(self, (self.batchPath + 'batch.pyobj', filename)[filename != ''])


    def makeDirectories(self):
            if not os.path.exists(self.batchPath):
                os.makedirs(self.batchPath)

            if not os.path.exists(self.batchPath + 'iniFiles'):
                os.makedirs(self.batchPath + 'iniFiles')

            if not os.path.exists(self.batchPath + 'postIniFiles'):
                os.makedirs(self.batchPath + 'postIniFiles')

