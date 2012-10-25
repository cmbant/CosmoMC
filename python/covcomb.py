# usage:
# python covcmb.py out.covmat in1.covmat in2.covmat
# Nb. in1 values take priority over in2

import sys
import numpy

foutname = sys.argv[1]
finname1 = sys.argv[2]
finname2 = sys.argv[3]
    
fin1 = open(finname1, 'r')
fin2 = open(finname2, 'r')

params = []
params1 = []
params2 = []

line = fin1.readline()
words = line.split()
for word in words:
    if word !='#':
        params1.append(word)
        if word not in params:
            params.append(word)
            
line = fin2.readline()
words = line.split()
for word in words:
    if word !='#':
        params2.append(word)
        if word not in params:
            params.append(word)


#print params1
#print params2
#print params

l1 = len(params1)
l2 = len(params2)
l = len(params)

map1 = dict(zip(params1,range(0,l1)))
map2 = dict(zip(params2,range(0,l2)))
map = dict(zip(range(0,l),params))

cov1 = numpy.empty((l1,l1))
cov2 = numpy.empty((l2,l2))
cov = numpy.empty((l,l))

for i in range(0,l1):
    data = fin1.readline().split()
    for j in range(0,l1):
        cov1[i][j]=data[j]

for i in range(0,l2):
    data = fin2.readline().split()
    for j in range(0,l2):
        cov2[i][j]=data[j]

fin1.close()
fin2.close()

for i in range(0, l):
    for j in range(0, l):
        if params[i] in params1 and params[j] in params1:
            cov[i][j] = cov1[map1[map[i]]][map1[map[j]]]
        elif params[i] in params2 and params[j] in params2:
            cov[i][j] = cov2[map2[map[i]]][map2[map[j]]]
        else:
            cov[i][j] = 0.

fout = open(foutname, 'w')
delimiter = ' '
fout.write('# '+delimiter.join(params)+'\n')
for i in range(0,l):
    nums = ["%E" % number for number in cov[i]]
    fout.write(delimiter.join(nums)+'\n')
    
fout.close
    

