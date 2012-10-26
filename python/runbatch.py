#!/usr/bin/env python
import sys
import os


if len(sys.argv) < 2:
    print 'Usage: python/runbatch.py directory_containing_inifiles [num_nodes]'
    sys.exit()

subScript = 'runMPI_HPCS.pl'

noOfMpiNodes = '2';
if len(sys.argv) > 2: noOfMpiNodes = sys.argv[2]

iniDir = os.path.abspath(sys.argv[1]) + os.sep

dirList = os.listdir(iniDir)
for fname in dirList:
    command = 'perl ' + subScript + ' ' + iniDir + fname + ' ' + noOfMpiNodes
    print 'Submitting...' + command
    os.system(command)
