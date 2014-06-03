# GetDist.py

import os
import sys
import glob
import numpy as np

import chains as ch
import iniFile
import paramNames
from MCSamples import MCSamples

#MT FIXME?: use batchJobArgs.batchArgs
if (len(sys.argv)<2):
    print 'Usage: python/GetDist.py ini_file [chain_root]'
    sys.exit()

# Parameter file
ini_file = os.path.abspath(sys.argv[1])
if not os.path.isfile(ini_file):
    print 'Parameter file do not exist ', ini_file
    sys.exit()

# Input parameters
ini = iniFile.iniFile()
ini.readFile(ini_file)
#dataset = ini.params

parameter_names_file = ini.string('parameter_names')
if (parameter_names_file):
    p = paramNames.paramNames(parameter_names_file)
parameter_names_labels = ini.string('parameter_names_labels', False) # MT: False added here

if (len(sys.argv)>2):
    in_root = sys.argv[2]
else:
    in_root = ini.params['file_root']
if (in_root==""): # or (not os.path.isfile(file_root)): # MT FIXME test file existence? 
     print 'Root file do not exist ', in_root
     sys.exit()
rootname = os.path.basename(in_root)

if (ini.params.has_key('nparams')):
    ncols = ini.int('nparams') + 2 
    if (ini.params.has_key('columnnum')):
        print  'specify only one of nparams or columnnum'
        sys.exit(1)
else:
    ncols = ini.int('columnnum',0)

single_column_chain_files = ini.bool('single_column_chain_files', False)

num_bins = ini.int('num_bins')
num_bins_2D = ini.int('num_bins_2D', num_bins)
smooth_scale_1D = ini.float('smooth_scale_1D', -1.)
smooth_scale_2D = ini.float('smooth_scale_2D', -1.)
#MT TODO: ???
#    if (smooth_scale_1D>0 .and. smooth_scale_1D>1) write(*,*) 'WARNING: smooth_scale_1D>1 is oversmoothed'
#    if (smooth_scale_1D>0 .and. smooth_scale_1D>1.9) stop 'smooth_scale_1D>1 is now in stdev units'

credible_interval_threshold = ini.float('credible_interval_threshold', 0.05)

ignorerows = ini.float('ignore_rows',0.)

adjust_priors = ini.bool('adjust_priors', False)

plot_ext = ini.string('plot_ext','py')

plot_output = ini.string('plot_output','pdf')

subplot_size_inch = ini.float('subplot_size_inch', 3.0)
subplot_size_inch2 = ini.float('subplot_size_inch2', subplot_size_inch)
subplot_size_inch3 = ini.float('subplot_size_inch3', subplot_size_inch)

font_scale = ini.float('font_scale',1.)
finish_run_command = ini.string('finish_run_command', '') # MT added '' here

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

#MT TODO
#    do ix = 1, Ini%Read_Int('compare_num',0)
#        call ComparePlots%Add(File%ExtractName(Ini%Read_String(numcat('compare',ix))))


has_limits_top = False
has_limits_bot = False

bin_limits = ini.string('all_limits')
markers = 0
has_markers = False

#MT TODO
# ...

if (ini.params.has_key('plotparams_num')):
    print 'plotparams_num deprectated; just use plot_params'
    sys.exit(1)

line = ini.string('plot_params', 0) # MT added 0 here
if (line<>""):
    plotparams_num = -1

line = ini.string('plot_2D_param')
if (line==""):
    plot_2D_param = 0
else:
    #MT TODO
    #call NameMapping%ReadIndices(InLine, tmp_params, 1)
    plot_2D_param = 123456789 # todo tmp_params(1)
    #if (plot_2D_param/=0 .and. plotparams_num/=0 .and. &
    #count(plotparams(1:plotparams_num)==plot_2D_param)==0) &
    #stop 'plot_2D_param not in plotparams'
    
if (plot_2D_param<>0):
    plot_2D_param = plot_2D_param + 2
    num_cust2D_plots = 0
else:
    # Use custom array of specific plots
    num_cust2D_plots = ini.int('plot_2D_num',0)
    for i in range(1, num_cust2D_plots+1):
        line = ini.string('plot'+str(i))
        #MT todo
        #call NameMapping%ReadIndices(InLine, tmp_params, 2)
        #if (plotparams_num/=0 .and. (count(plotparams(1:plotparams_num)==tmp_params(1))==0 .or. &
        #count(plotparams(1:plotparams_num)==tmp_params(2))==0)) then
        #   write(*,*) trim(numcat('plot',ix)) //': parameter not in plotparams'
        #   stop
        #cust2DPLots(ix) = tmp_params(1)+2 + (tmp_params(2)+2)*1000

triangle_plot = ini.bool('triangle_plot',False)
if (triangle_plot):
    no_triangle_axis_labels = ini.bool('no_triangle_axis_labels',False)
    line = ini.string('triangle_params')
    triangle_num = -1
    if (line<>""):
        #call NameMapping%ReadIndices(InLine, triangle_params, triangle_num, unknown_value=-1)
        pass

exclude_chain = ini.string('exclude_chain') 
# ...

map_params = ini.bool('map_params',False)
if (map_params):
    print 'WARNING: Mapping params - .covmat file is new params.'

shade_meanlikes = ini.bool('shade_meanlikes',False)

out_dir = ini.string('out_dir')

out_root = ini.string('out_root')
if (out_root<>''):
    rootname = out_root
    print 'producing files with with root ', out_root

plot_data_dir = ini.string('plot_data_dir')
if (plot_data_dir==''):
    plot_data_dir = 'plot_data/'

rootdirname = os.path.join(out_dir, rootname)

num_contours = ini.int('num_contours',2)
contours = []
for i in range(1, num_contours+1):
    contours.append(ini.float('contour'+str(i)))
    #MT TODO
    #max_frac_twotail(i) = Ini%Read_Double(numcat('max_frac_twotail',i), exp(-dinvnorm((1-contours(i))/2)**2/2))
contours_str = "; ".join([ str(c) for c in contours ]) 

if (not no_tests):
    converge_test_limit = ini.float('converge_test_limit',contours[num_contours-1])
    corr_length_thin = ini.int('corr_length_thin', 0)
    corr_length_steps = ini.int('corr_length_steps', 15)
    
force_twotail = ini.bool('force_twotail',False)
if (force_twotail): print 'Computing two tail limits'

#MT TODO
#if (Ini%Read_String('cov_matrix_dimension')=='') then
#   if (NameMapping%nnames/=0) covmat_dimension = NameMapping%num_MCMC
#else
#   covmat_dimension = Ini%Read_Int('cov_matrix_dimension',0)
#   if (covmat_dimension == -1) covmat_dimension = ncols-2
#end if

plot_meanlikes = ini.bool('plot_meanlikes',False)

if (ini.params.has_key('do_minimal_1d_intervals')):
    print 'do_minimal_1d_intervals no longer used; set credible_interval_threshold instead'
    sys.exit(1)

PCA_num = ini.int('PCA_num',0)
if (PCA_num<>0):
    if (PCA_num<2):
        print 'Can only do PCA for 2 or more parameters'
        sys.exit(1)
    line = ini.string('PCA_params')
    PCA_func = ini.string('PCA_func')
    # Characters representing functional mapping
    if (PCA_func==''):
        PCA_func = ['N'] * PCA_num  # No mapping
    if (line.lower()=="all"):
        PCA_params = range(1, PCA_num+1)
    else:
        PCA_params = [] #MT todo find indexes
    line = ini.string('PCA_normparam')
    if (line==""):
        PCA_NormParam = 0
    else:
        PCA_NormParam = 123456789 # todo tmp_params(1)

num_3D_plots = ini.int('num_3D_plots',0)
#plot_3D = [ ini.string('3D_plot'+str(ix)) of ix in range(1, num_3D_plots+1) ]
plot_3D = []
for ix in range(1, num_3D_plots+1):
    plot_3D.append(ini.string('3D_plot'+str(ix)))

make_scatter_samples = ini.bool('make_scatter_samples',False)
max_scatter_points = ini.int('max_scatter_points',2000)

BW = ini.bool('B&W', False)
do_shading = ini.bool('do_shading',True)


def getLastChainIndex(in_root):
    chain_files = glob.glob(in_root+"_*.txt")
    if not chain_files: return 0
    names_files = [ os.path.basename(f) for f in chain_files ]
    basename = os.path.basename(in_root)
    indexes = [ int(f.replace(basename+'_', '').replace('.txt', '')) for f in names_files ]
    return max(indexes)

first_chain = ini.int('first_chain',1) 
last_chain = ini.int('chain_num',-1) 
# -1 means keep reading until one not found
if(last_chain==-1):
    last_chain = getLastChainIndex(in_root)

# Read in the chains
chains = ch.loadChains(in_root, range(first_chain, last_chain+1), ignore_frac=ignorerows)
#chains.getChainsStats()

# FIXME?: chains to exclude 
# FIXME?: distinction single column vs. multiple columns
# FIXME?: ignorerows 
# FIXME?: check if no chains were read



import pdb; pdb.set_trace()

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
    






