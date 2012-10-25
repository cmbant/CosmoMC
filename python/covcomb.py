# usage:
# python covcmb.py out.covmat in1.covmat in2.covmat
# Nb. in1 values take priority over in2

import sys
import numpy
import covMat

if len(sys.argv) < 3:
    print 'Usage: python covcmb.py out.covmat in1.covmat in2.covmat'

foutname = sys.argv[1]
finname1 = sys.argv[2]
finname2 = sys.argv[3]
    
cov1 = covMat()    
cov2 = covMat()    
cov1.loadFromFile(finname1)
cov2.loadFromFile(finname2) 
cov = covMat()

cov.paramNames.add(cov1.paramNames)   
  
   
for param in cov2.paramNames:
        if param not in cov.paramNames:
            cov.paramNames.append(param)

l1 = len(cov1.paramNames)
l2 = len(cov2.paramNames)
l = len(cov.paramNames)

params1=cov1.paramNames
params2=cov2.paramNames

map1 = dict(zip(params2,range(0,l1)))
map2 = dict(zip(params2,range(0,l2)))
covmap = dict(zip(range(0,l),cov.paramNames))

cov.matrix = numpy.zeros((l,l))

for i in range(0, l):
    for j in range(0, l):
        if cov.paramNames[i] in params1 and cov.paramNames[j] in params1:
            cov.matrix[i][j] = cov1.matrix[map1[covmap[i]]][map1[covmap[j]]]
        elif cov.paramNames[i] in params2 and cov.paramNames[j] in params2:
            cov[i][j] = cov2[map2[map[i]]][map2[map[j]]]

cov.saveToFile(foutname)

