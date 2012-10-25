# usage:
# python covcmb.py out.covmat in1.covmat in2.covmat
# Nb. in1 values take priority over in2

import sys
import numpy
import covMat

#if len(sys.argv) < 3:
#    print 'Usage: python covcmb.py out.covmat in1.covmat in2.covmat'

foutname = sys.argv[1]
finname1 = sys.argv[2]
finname2 = sys.argv[3]

  
cov1 = covMat.covMat(finname1)    
cov2 = covMat.covMat(finname2)    
C = covMat.covMat()

C.paramNames.extend(cov1.paramNames)   
     
for param in cov2.paramNames:
        if param not in C.paramNames:
            C.paramNames.append(param)

params1=cov1.paramNames
params2=cov2.paramNames

l1 = len(params2)
l2 = len(params2)
l = len(C.paramNames)

map1 = dict(zip(params2,range(0,l1)))
map2 = dict(zip(params2,range(0,l2)))
covmap = dict(zip(range(0,l),C.paramNames))

C.matrix = numpy.zeros((l,l))

for i in range(0, l):
    for j in range(0, l):
        if C.paramNames[i] in params1 and C.paramNames[j] in params1:
            C.matrix[i,j] = cov1.matrix[map1[covmap[i]],map1[covmap[j]]]
        elif C.paramNames[i] in params2 and C.paramNames[j] in params2:
            C.matrix[i,j] = cov2.matrix[map2[covmap[i]],map2[covmap[j]]]

C.saveToFile(foutname)

