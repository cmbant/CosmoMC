import os, sys

#path to put input files and outputchains, make if not exist
batchPath=  'z://test1/'
#'/scratch/aml1005/gridtest/'

#sets of parameters to vary in addition to baseline
extparams=[[],['mnu'],['omegak']]

#dataset names
planck='planck_CAMspec'
highL='highL'

datasets=[]
#lists of dataset names to combine, with corresponding sets of inis to include
datasets.append([[planck],['CAMspec_defaults.ini']])
datasets.append([[planck,highL],['CAMspec_ACTSPT_defaults.ini']])

#priors and widths for parameters which are varied
params=dict()
params['mnu']='0 0 5 0.1 0.03'
params['omegak']='0 -0.3 0.3 0.001 0.001'
params['w']='-1 -3 -0.3 0.02 0.02'
params['nnu']='3.046 0 10 0.1 0.1'
params['nrun']='0 -1 1 0.001 0.001'
params['r']='0 0 2 0.03 0.03'
params['Alens']='0 0 10 0.05 0.05'
params['yhe']='0.245 0.1 0.5 0.006 0.006'
params['alpha1']='0 -1 1 0.0003 0.0003'
params['deltazrei']='0.5 0.1 3 0.3 0.3'
params['wa']='0 -2 2 0.3 0.3'


#ini files you want to base each set of runs on
defaults=['common_batch1.ini']


basePath = os.path.dirname(sys.path[0]) +'/../'

commonPath = basePath + 'batch1/'

if not os.path.exists(batchPath):
    os.makedirs(batchPath)
    
if not os.path.exists(batchPath+'iniFiles'):
    os.makedirs(batchPath+'iniFiles')



for dataset in datasets:
    for param_set in extparams:
        paramtag = 'base';
        f = []
    
        for param in param_set:
            f.append('param['+param+']=' + params[param])
            paramtag= paramtag + '_' + param
        for deffile in defaults:
            f.append('DEFAULT(' + commonPath + deffile + ')')
        datatag="_".join(dataset[0])
        for ini in dataset[1]:
            f.append('INCLUDE('+commonPath+ ini+')')
                
        if 'mnu' in param_set:
            f.append('num_massive_neutrinos=3')        
                
        chainPath = batchPath +paramtag+'/' +datatag  +'/'
        os.makedirs(chainPath)
        fulltag=paramtag+'_'+datatag 
        f.append('file_root = ' + chainPath + fulltag)
        
        f.append('propose_matrix =' + datatag+'.covmat') 
    
        outfile = open(batchPath+'iniFiles/' +fulltag+ '.ini', 'w')
        outfile.write("\n".join(f))
        outfile.close()

print 'Done... to run do: python python/runbatch.py '+batchPath+'iniFiles'
