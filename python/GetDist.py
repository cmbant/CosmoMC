# GetDist.py

import os
import sys
import numpy as np

import iniFile

import paramNames
from MCSamples import MCSamples

# 




# Initialization
# --------------

# use batchJobArgs.batchArgs

if (len(sys.argv)<2):
    print 'Usage: python/GetDist.py ini_file [chain_root]'
    sys.exit()
    #sys.exit('No parameter input file')

# Parameter file
ini_file = os.path.abspath(sys.argv[1])
if not os.path.isfile(ini_file):
    print 'Parameter file do not exist ', ini_file
    sys.exit()

ini = iniFile.iniFile()
ini.readFile(ini_file)
dataset = ini.params

# file_root
if (len(sys.argv)>2):
    file_root = sys.argv[2]
else:
    file_root = ini.params['file_root']
if (file_root=="") or (not os.path.isfile(file_root)):
     print 'Root file do not exist ', file_root
     sys.exit()


parameter_names_labels = ini.string('parameter_names_labels')

nparams = ini.int('nparams')
columnnum = ini.int('columnnum',0)
# TODO ncols = 

single_column_chain_files = ini.bool('single_column_chain_files', False)


num_bins = ini.int('num_bins')
num_bins_2D = ini.int('num_bins_2D', num_bins)
smooth_scale_1D = ini.float('smooth_scale_1D', -1.)
smooth_scale_2D = ini.float('smooth_scale_2D', -1.)
# TODO: ???
#    if (smooth_scale_1D>0 .and. smooth_scale_1D>1) write(*,*) 'WARNING: smooth_scale_1D>1 is oversmoothed'
#    if (smooth_scale_1D>0 .and. smooth_scale_1D>1.9) stop 'smooth_scale_1D>1 is now in stdev units'


credible_interval_threshold = ini.float('credible_interval_threshold', 0.05)

ignorerows = ini.float('ignore_rows',0.)

adjust_priors = ini.bool('adjust_priors', False)

plot_ext = ini.string('plot_ext','py')

plot_output = ini.string('plot_output','pdf')

#subplot_size_inch = ini.float('subplot_size_inch')
#subplot_size_inch2 = ini.float('subplot_size_inch2', subplot_size_inch)
#subplot_size_inch3 = ini.float('subplot_size_inch3', subplot_size_inch)

font_scale = ini.float('font_scale',1.)
finish_run_command = ini.string('finish_run_command')

auto_label = ini.bool('auto_label',False)

prob_label = ini.bool('prob_label',False)

samples_are_chains = ini.bool('samples_are_chains',True)

no_plots = ini.bool('no_plots',False)
plots_only = ini.bool('plots_only',False)
no_tests = plots_only or ini.bool('no_tests',False)
line_labels = ini.bool('line_labels',False)

thin_factor = ini.int('thin_factor',0)
thin_cool = ini.float('thin_cool',1.)


make_single_samples = ini.bool('make_single_samples', False)
single_thin = ini.int('single_thin',1)
cool = ini.float('cool',1.)


# plotparams_num deprecated message
# plot_params 


plot_2D_param = 0 # todo 


triangle_plot = ini.bool('triangle_plot',False)
if (triangle_plot):
    no_triangle_axis_labels = ini.bool('no_triangle_axis_labels',False)


exclude_chain = ini.bool('exclude_chain') 
# ...


map_params = ini.bool('map_params',False)
if (map_params):
    print 'WARNING: Mapping params - .covmat file is new params.'

shade_meanlikes = ini.bool('shade_meanlikes',False)



out_dir = ini.bool('out_dir')

out_root = ini.bool('out_root')
if (out_root!=''):
    rootname = out_root
    print 'producing files with with root ', out_root


plot_data_dir = ini.bool('plot_data_dir')
if (plot_data_dir==''):
    plot_data_dir = 'plot_data/'


rootdirname = os.path.join(out_dir, rootname)

num_contours = ini.int('num_contours',2)


if (not no_tests):
    converge_test_limit = ini.float('converge_test_limit',contours(num_contours))
    corr_length_thin = ini.int('corr_length_thin',corr_length_thin)
    corr_length_steps = ini.int('corr_length_steps',corr_length_steps)
    

force_twotail = ini.bool('force_twotail',False)
if (force_twotail):
    print 'Computing two tail limits'



plot_meanlikes = ini.bool('plot_meanlikes',False)



# PCA 

PCA_num = ini.int('PCA_num',0)
if (PCA_num!=0):
    if (PCA_num<2):
        print 'Can only do PCA for 2 or more parameters'
        sys.exit(1)
        
    # ...



num_3D_plots = ini.int('num_3D_plots',0)
# ...


make_scatter_samples = ini.bool('make_scatter_samples',False)
max_scatter_points = ini.int('max_scatter_points',2000)





def getLastChainIndex(file_root):
    # TODO: get max index of chain
    return 0

first_chain = ini.int('first_chain',1)
last_chain = ini.int('chain_num',-1) 
# -1 means keep reading until one not found
if(last_chain==-1):
    last_chain = getLastChainIndex(file_root)



# Read in the chains
chains = loadChains(file_root, range(first_chain, last_chain+1), ignore_rows=ignorerows)









self.nrows = 0
self.num_chains_used = 0





if (self.num_chains_used > self.max_chains):
    sys.exit('Increase max_chains in GetDist')

    #chain_indices(num_chains_used) = nrows
    #chain_numbers(num_chains_used) = chain_ix


if (single_column_chain_files):
    # Use used for WMAP 5-year chains suppled on LAMBDA; code from Mike Nolta
    # Standard CosmoMC case below

    first_haschain = 0
    for ip in range(self.ncols):
            #infile = concat(File%CheckTrailingSlash(concat(in_root,chain_ix)), pname(ip))
        if not os.path.exists(infile):
            print 'skipping missing ', infile 
                #coldata(ip,:) = 0
                #nrows2(ip) = -1
        else:
            print 'reading ', infile 
                #call ChainFile%Open(infile)
            if (first_haschain==0): first_haschain = ip
                # ...
            idx = 0
                # ...
                # TODO: read data from file


        if (first_haschain==0):
            sys.exit('no chain parameter files read!')

        for ip in range(2, self.ncols+1): # ?
            # ...
            pass

        #nrows = nrows2(first_haschain)
        print 'all columns match, nrows = ', self.nrows

    else:
        # Not single column chain files (usual cosmomc format)
        # This increments nrows by number read in

        # ...
        pass


    if (map_params):
        # ...
        pass


    if (ignorerows<1) and (ignorerows!=0):
        # ... 
        pass


if (self.nrows==0):
    sys.exit('No un-ignored rows! (check number of chains/burn in)')

if (cool!=1):
    self.CoolChain(cool)

# Adjust weights if requested
if (adjust_priors):
    self.AdjustPriors()  


# See which parameters are fixed
self.GetUsedCols()


# ...


if (not no_tests):
    self.DoConvergeTests(converge_test_limit)

if (adjust_priors):
    self.DeleteZeros()

print 'mean input multiplicity = ', mean_mult


# Output thinned data if requested
# Must do this with unsorted output
if (thin_factor!=0):
    self.ThinData(thin_factor)
    self.WriteThinData(self.rootdirname+'_thin.txt', thin_cool) # ?


# Produce file of weight-1 samples if requested

# ...


# Only use variables whose labels are not empty (and in list of plotparams if plotparams_num /= 0)
num_vars = 0

if (plotparams_num!=0):

    # ...
    pass

else:

    # ...
    pass


for j in range(num_vars):
    #mean(j) = sum(coldata(1,0:nrows-1)*coldata( colix(j),0:nrows-1))/numsamp
    #sddev(j)  = sqrt(sum(coldata(1,0:nrows-1)*(coldata(colix(j),0:nrows-1) -mean(j))**2)/numsamp)
    pass


if (make_single_samples):
    self.MakeSingleSamples(single_thin)


# IO_WriteBounds

# Sort data in order of likelihood of points
#call SortColData(2)

#numsamp = sum(coldata(1,0:nrows-1))


# Get ND confidence region (index into sorted coldata)
counts = 0
ND_cont1, ND_cont2 = -1, -1

for j in range(0, self.nrows-1):
    # ...
    pass


triangle_plot = triangle_plot and (num_vars>1)
if (triangle_plot):
    # ...
    pass

print 'using ',nrows,' rows, processing ',num_vars,' parameters'
if (indep_thin!=0):
    print 'Approx indep samples: ', int(numsamp/indep_thin) # equiv. to nint ?
else:
    print  'effective number of samples (assuming indep): ', int(numsamp/max_mult)  # equiv. to nint ?

# Get covariance matrix and correlation matrix
self.GetCovMatrix()

if (PCA_num>0) and not plots_only:
    self.PCA(PCA_params,PCA_num,PCA_func, PCA_NormParam)


# Find best fit, and mean likelihood
self.GetChainLikeSummary(stdout)

if (not no_plots):
    # Matlab only ? 
    pass


LowerUpperLimits = 0


# 1D 
# --


# Do 1D bins
for j in range(num_vars):
    # ... 

    self.Get1DDensity(j)







# 2D 
# --



# 3D 
# --



# Statistics
# ----------




import pdb
pdb.set_trace()




