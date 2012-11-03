import sys, batchJob
try: import argparse
except:
    print 'use "module load" to load python 2.7'
    sys.exit()

class batchArgs():

        def __init__(self, desc='', importance=False, noBatchPath=False):
            self.parser = argparse.ArgumentParser(description=desc)
            if not noBatchPath: self.parser.add_argument('batchPath')
            self.importanceParameter = importance;

        def parseForBatch(self):
            if self.importanceParameter: self.parser.add_argument('--importance', nargs='*', default=None)
            self.args = self.parser.parse_args()
            return (batchJob.readobject(self.args.batchPath), self.args)

        def wantImportance(self, importanceTag):
            return self.args.importance is None or len(self.args.importance) == 0 or importanceTag in self.args.importance

        def jobItemWanted(self, jobItem):
            return not jobItem.isImportanceJob and (self.args.importance is None) or jobItem.isImportanceJob and self.wantImportance(jobItem.importanceTag)




