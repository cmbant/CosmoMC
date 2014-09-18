import numpy as np
import subprocess
import fileinput

root = '/home/adammoss/work/code/PLANCK/cosmomcplanck/'

angles = np.linspace(1,301,30)
NL = [['F','1'],['T','1'],['T','2'],['T','3'],['T','4']]
filename = root+'data/CFHTLENS/2D/cut_values_test.dat'
program = [root+'/cosmomc',root+'ancillary/test_2D.ini']

for i in angles:
	for j in angles:
		result = ''
		for k in NL:
			with open(root+'batch2/WL_test_2D.ini') as a:
				for x in a:
					if 'use_non_linear' in x:
						for line in fileinput.input(root+'batch2/WL_test_2D.ini',inplace=True):
							print line.replace(x,'use_non_linear = '+k[0]+'\n'),
			with open(root+'batch2/WL_test_2D.ini') as a:
				for x in a:
					if 'halofit_version' in x:
						for line in fileinput.input(root+'batch2/WL_test_2D.ini',inplace=True):
							print line.replace(x,'halofit_version = '+k[1]+'\n'),
			file = open(filename,'w')
			file.write(str(i)+'\t'+str(j))
			file.close
			file = open(root+'ancillary/output_2D.dat','w')
			p = subprocess.call(program,stdout=file)
			file.close
			with open(root+'ancillary/output_2D.dat') as file:
				for l in file:
					if 'WL: CFHTLENS_1bin' in l:
						out = l.split()[1];
						print k[0],k[1],i,j,out
						result += ', '+str(out)
		outfile = open(root+'ancillary/result_2D.dat','a')
		outfile.write(str(i)+', '+str(j)+result+'\n')
		outfile.close()
		


