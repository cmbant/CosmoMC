#!/usr/bin/env python
import sys
import os


if len(sys.argv) < 2:
	print 'Usage: test/runbatch.py directory_containing_inifiles'

subScript = 'runMPI_HPCS.pl'

noOfMpiNodes = 2;

iniDir=os.path.abspath(sys.argv[1])

dirList=os.listdir(iniDir)
for fname in dirList:
	print '***** submitting...' + fname
	print  'test:'+ 'perl '+subScript +' ' + fname + ' '+noOfMpiNodes
#	os.system('perl '+subScript +' ' + fname + ' '+noOfMpiNodes)

