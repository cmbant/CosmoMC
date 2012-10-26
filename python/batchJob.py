
import os, sys, pickle

def readobject(filename):
    with open(filename, 'rb') as inp:
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

    def iniFile(self):
        return self.batchPath + 'iniFiles/' + self.name + '.ini'

    def makeChainPath(self):
        if not os.path.exists(self.chainPath): os.makedirs(self.chainPath)
        return self.chainPath

    def writeIniLines(self, f):
        outfile = open(self.iniFile(), 'w')
        outfile.write("\n".join(f))
        outfile.close()



class batchJob:

    def __init__(self, path):
        self.batchPath = path
        self.extparams = []
        self.datasets = []
        self.basePath = os.path.dirname(sys.path[0]) + os.sep
        self.commonPath = self.basePath + 'batch1/'

    def items(self):
        for data_set in self.datasets:
            for param_set in self.extparams:
                    yield(jobItem(self.batchPath, param_set, data_set))

    def save(self, filename=''):
        saveobject(self, (self.batchPath + 'batch.pyobj', filename)[filename != ''])


    def makeDirectories(self):
            if not os.path.exists(self.batchPath):
                os.makedirs(self.batchPath)

            if not os.path.exists(self.batchPath + 'iniFiles'):
                os.makedirs(self.batchPath + 'iniFiles')


