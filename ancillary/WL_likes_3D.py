import numpy as np
import subprocess
import fileinput

root = '/home/adammoss/work/code/PLANCK/cosmomcplanck/'

theta1 = np.linspace(1,31,5)
theta2= np.linspace(1,31,5)
theta3 = np.linspace(1,31,5)
theta4 = np.linspace(1,31,5)
NL = [['F','1'],['T','1'],['T','2'],['T','3'],['T','4']]
filename = root+'data/CFHTLENS/blu_sample/cut_values_test.dat'
program = [root+'/cosmomc',root+'ancillary/WL_test_3D.ini']

for t1 in theta1:
	for t2 in theta2:
		for t3 in theta3:
			for t4 in theta4:
				result = ''
				for k in NL:
					with open(root+'batch2/WL_test_3D.ini') as a:
						for x in a:
							if 'use_non_linear' in x:
								for line in fileinput.input(root+'batch2/WL_test_3D.ini',inplace=True):
									print line.replace(x,'use_non_linear = '+k[0]+'\n'),
					with open(root+'batch2/WL_test_3D.ini') as a:
						for x in a:
							if 'halofit_version' in x:
								for line in fileinput.input(root+'batch2/WL_test_3D.ini',inplace=True):
									print line.replace(x,'halofit_version = '+k[1]+'\n'),
					file = open(filename,'w')
					file.write(str(t1)+'\t'+str(t2)+'\n')
					file.write(str(t1)+'\t'+str(t2)+'\n')
					file.write(str(t3)+'\t'+str(t2)+'\n')
					file.write(str(t3)+'\t'+str(t2)+'\n')
					file.write(str(t3)+'\t'+str(t4)+'\n')
					file.write(str(t3)+'\t'+str(t4)+'\n')
					file.close
					file = open(root+'ancillary/output_3D.dat','w')
					p = subprocess.call(program,stdout=file)
					file.close
					with open(root+'ancillary/output_3D.dat') as file:
						for l in file:
							if 'WL: CFHTLENS_6bin' in l:
								out = l.split()[1];
								print k[0],k[1],t1,t2,t3,t4,out
								result += ', '+str(out)
				outfile = open(root+'ancillary/result_3D.dat','a')
				outfile.write(str(t1)+', '+str(t2)+', '+str(t3)+', '+str(t4)+result+'\n')
				outfile.close()
		


