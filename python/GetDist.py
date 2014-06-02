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
if (out_root<>''):
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
if (PCA_num<>0):
    if (PCA_num<2):
        print 'Can only do PCA for 2 or more parameters'
        sys.exit(1)
        
    # ...



num_3D_plots = ini.int('num_3D_plots',0)
# ...


make_scatter_samples = ini.bool('make_scatter_samples',False)
max_scatter_points = ini.int('max_scatter_points',2000)

BW = ini.bool('B&W', False)
do_shading = ini.bool('do_shading',True)


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
chains.getChainsStats()

# FIXME?: chains to exclude 
# FIXME?: distinction single column vs. multiple columns
# FIXME?: ignorerows 
# FIXME?: check if no chains were read


#TODO: CoolChain(cool)

# Adjust weights if requesteda
if (adjust_priors):
    # TODO: AdjustPriors()
    pass

# See which parameters are fixed (see Chains.loadWMAPChain)
#TODO: GetUsedCols()


# ...


if (not no_tests):
    #TODO: DoConvergeTests(converge_test_limit)
    pass

if (adjust_priors):
    #TODO: DeleteZeros
    pass


#print 'mean input multiplicity = ',mean_mult


# Output thinned data if requested
# Must do this with unsorted output
if (thin_factor<>0):
    #TODO ThinData(thin_factor)
    # Loop over chain + call chain.thin(thin_factor)
    filename = rootdirname + '_thin.txt'
    #TODO WriteThinData(trim(rootdirname)//'_thin.txt',thin_cool)



# Produce file of weight-1 samples if requested

# ...


# Only use variables whose labels are not empty (and in list of plotparams if plotparams_num /= 0)
num_vars = 0

if (plotparams_num<>0):

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
#TODO call SortColData(2)

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
if (indep_thin<>0):
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



# Do 1D bins

for j in range(num_vars):
    # ... 

    #self.Get1DDensity(j)
    
    #call to twoTailLimits in chains.py
    pass

if (not no_plots):
    if (plot_ext=='py'):
        text = "g.plots_1d(roots)"

if (not no_plots):
    # Output files for 1D plots
    filename = rootdirname + plot_ext
    textFileHandle = open(filename, 'w')
    if (plot_ext=='py'):
        text = "g.plots_1d(roots)"
        # ...

    
    if (triangle_plot):
        filename = rootdirname + '_tri.' + plot_ext
        textFileHandle = open(filename, 'w')
        if (plot_ext=='py'):
            text = "g.triangle_plot(roots, %s)"%data # todo 


# Do 2D bins
if (plot_2D_param==0) and (num_cust2D_plots==0) and (not no_plots):
    # In this case output the most correlated variable combinations
    print "doing 2D plots for most correlated variables"
    
    # ...

    
if (num_cust2D_plots==0):
    num_2D_plots = 0

    # ...
else:
    num_2D_plots = num_cust2D_plots
      

if (num_2D_plots>0) and (not no_plots):
    print 'Producing ', num_2D_plots,' 2D plots'
    filename = rootdirname + '_2D.' + plot_ext
    textFileHandle = open(filename, 'w')
    if (plot_ext=='py'):
        text = 'pairs=[]'

    # ... 
            
if (triangle_plot) and (not no_plots):
    # Add the off-diagonal 2D plots
    for i in range(triangle_num):
        # Matlab only ???
        pass



# Do 3D plots (i.e. 2D scatter plots with coloured points)

if (num_3D_plots<>0 and not no_plots):
    print 'producing ',num_3D_plots, '2D colored scatter plots'
    filename = rootdirname + '_3D.' + plot_ext
    textFileHandle = open(filename, 'w')
    # call to WritePlotFileInit
    text  = ""
    text += "sets=[]"
    for j in range(num_3D_plots):
        text += "sets.append([%s,%s,%s])"%(v1, v2, v3) # todo
    text += "g.plots_3d(roots,sets)"
    textFileHandle.close()




def WritePlotFileInit(handle):
    text = """
import GetDistPlots, os
g=GetDistPlots.GetDistPlotter('%s')
g.settings.setWithSubplotSize(%f)
outdir='%s'
roots=['%s', '%s']
"""
    return text




# Write out stats marginalized
# TODO: margeStats

# TODO: write file paramnames in plot_data/ dir

# Limits from global likelihood
# TODO: LikeFile 


# System command
if (finish_run_command):
    finish_run_command = finish_run_command.replace('%ROOTNAME%',rootname)
    finish_run_command = finish_run_command.replace('%PLOTDIR%',plot_data_dir)
    finish_run_command = finish_run_command.replace('%PLOTROOT%', os.path.join(plot_data_dir, rootname))
    os.system(finish_run_command)
    






