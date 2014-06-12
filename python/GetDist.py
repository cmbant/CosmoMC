# GetDist.py

import os
import sys
import glob
import math
import numpy as np
from scipy.stats import norm

import chains as ch
import iniFile
import paramNames
import MCSamples

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

# File root
if (len(sys.argv)>2):
    in_root = sys.argv[2]
else:
    in_root = ini.params['file_root']
if (in_root==''): # or (not os.path.isfile(file_root)): 
     print 'Root file do not exist ', in_root
     sys.exit()
rootname = os.path.basename(in_root)

# Create instance of MCSamples
mc = MCSamples.MCSamples(in_root)

if (ini.params.has_key('nparams')):
    ncols = ini.int('nparams') + 2 
    if (ini.params.has_key('columnnum')):
        print  'specify only one of nparams or columnnum'
        sys.exit()
else:
    ncols = ini.int('columnnum', 0)

single_column_chain_files = ini.bool('single_column_chain_files', False)

num_bins = ini.int('num_bins')
num_bins_2D = ini.int('num_bins_2D', num_bins)
smooth_scale_1D = ini.float('smooth_scale_1D', -1.0)
# Smoothing scale in terms of bin scale
smooth_scale_2D = ini.float('smooth_scale_2D', -1.0)
if (smooth_scale_1D>0) and (smooth_scale_1D>1): 
    print 'WARNING: smooth_scale_1D>1 is oversmoothed'
if (smooth_scale_1D>0) and (smooth_scale_1D>1.9):
    print 'smooth_scale_1D>1 is now in stdev units'
    sys.exit()

credible_interval_threshold = ini.float('credible_interval_threshold', 0.05)

ignorerows = ini.float('ignore_rows', 0.0)

adjust_priors = ini.bool('adjust_priors', False)

plot_ext = ini.string('plot_ext', 'py')

plot_output = ini.string('plot_output', 'pdf')

subplot_size_inch  = ini.float('subplot_size_inch' , 3.0)
subplot_size_inch2 = ini.float('subplot_size_inch2', subplot_size_inch)
subplot_size_inch3 = ini.float('subplot_size_inch3', subplot_size_inch)

font_scale = ini.float('font_scale', 1.0)
finish_run_command = ini.string('finish_run_command', '') 

auto_label = ini.bool('auto_label', False)

prob_label = ini.bool('prob_label', False)

samples_are_chains = ini.bool('samples_are_chains', True)

no_plots = ini.bool('no_plots', False)
plots_only = ini.bool('plots_only', False)
no_tests = plots_only or ini.bool('no_tests', False)
line_labels = ini.bool('line_labels', False)

thin_factor = ini.int('thin_factor', 0)
thin_cool = ini.float('thin_cool', 1.0)

make_single_samples = ini.bool('make_single_samples', False)
single_thin = ini.int('single_thin', 1)
cool = ini.float('cool', 1.0)

bin_limits = ini.string('all_limits')

indexes = mc.index2name.keys()
indexes.sort()
for ix in indexes:
    name = mc.index2name[ix]
    mini = mc.ranges.min(name)
    maxi = mc.ranges.max(name)
    if (mini and maxi and mini<>maxi):
        mc.limmin[name] = mini
        mc.limmax[name] = maxi
    if (bin_limits<>''):
        line = bin_limits
    else:
        line = ''
        if ini.params.has_key('limits[%s]'%name.strip()):
            line = ini.string('limits[%s]'%name.strip())
    if (line<>''):
        limits = [ s for s in line.split(' ') if s<>'' ]
        if len(limits)==2:
            if limits[0]<>'N': mc.limmin[name] = float(limits[0])
            if limits[1]<>'N': mc.limmax[name] = float(limits[1])
    if ini.params.has_key('marker[%s]'%name.strip()):
        line = ini.string('marker[%s]'%name.strip())
        if (line<>''):
            mc.markers[name] = float(line)

if (ini.params.has_key('plotparams_num')):
    print 'plotparams_num deprectated; just use plot_params'
    sys.exit()

line = ini.string('plot_params', 0)
if (line not in ['', 0]):
    plotparams_num = -1
    plotparams = [ s for s in line.split(' ') if s<>'' ] 
    plotparams = [ mc.index[name] for name in plotparams]
else:
    plotparams_num = 0

line = ini.string('plot_2D_param')
if (line==''):
    plot_2D_param = 0
else:
    tmp_params = [ s for s in line.split(' ') if s<>'' ]
    plot_2D_param = int(tmp_params[0])
    if (plot_2D_param<>0 and plotparams_num<>0 and plotparams.count(plot_2D_param)==0):
        print 'plot_2D_param not in plotparams'
        sys.exit()
    
if (plot_2D_param<>0):
    plot_2D_param = plot_2D_param + 2
    num_cust2D_plots = 0
else:
    # Use custom array of specific plots
    num_cust2D_plots = ini.int('plot_2D_num', 0)
    cust2DPlots = []
    for i in range(1, num_cust2D_plots+1):
        line = ini.string('plot'+str(i))
        tmp_params = [ s for s in line.split(' ') if s<>'' ]
        tmp_params = [ mc.index[name] for name in tmp_params]
        if (plotparams_num<>0 and plotparams.count(tmp_params[0])==0):
            print 'plot'+str(i), ': parameter not in plotparams'
            sys.exit()
        cust2DPlots.append(tmp_params[0]+2 + (tmp_params[1]+2)*1000)

triangle_plot = ini.bool('triangle_plot', False)
if (triangle_plot):
    no_triangle_axis_labels = ini.bool('no_triangle_axis_labels', False)
    line = ini.string('triangle_params')
    triangle_num = -1
    if (line<>''):
        triangle_params = [ s for s in line.split(' ') if s<>'' ]
        triangle_params = [ mc.index[name] for name in triangle_params if mc.index.has_key(name) ]

exclude_chain = ini.string('exclude_chain') 
chain_exclude = [ int(s) for s in exclude_chain.split(' ') if s<>'' ]
num_exclude = len(chain_exclude) - chain_exclude.count(0)

map_params = ini.bool('map_params', False)
if (map_params):
    print 'WARNING: Mapping params - .covmat file is new params.'

shade_meanlikes = ini.bool('shade_meanlikes', False)

out_dir = ini.string('out_dir')
if (out_dir<>''):
    print 'producing files in directory ', out_dir

out_root = ini.string('out_root')
if (out_root<>''):
    rootname = out_root
    print 'producing files with with root ', out_root

plot_data_dir = ini.string('plot_data_dir')
if (plot_data_dir==''):
    plot_data_dir = 'plot_data/'

#abs_plot_data_dir = os.path.join(out_dir, plot_data_dir)
abs_plot_data_dir = plot_data_dir
if not os.path.isdir(abs_plot_data_dir):
    os.mkdir(abs_plot_data_dir)

rootdirname = os.path.join(out_dir, rootname)

num_contours = ini.int('num_contours', 2)
contours = []
max_frac_twotail = []
for i in range(1, num_contours+1):
    contours.append(ini.float('contour'+str(i)))
    max_frac = ini.float('max_frac_twotail'+str(i), math.exp(math.pow(-1.0*norm.ppf((1-contours[i-1])/2), 2)/2))
    max_frac_twotail.append(max_frac)
contours_str = '; '.join([ str(c) for c in contours ]) 

if (not no_tests):
    converge_test_limit = ini.float('converge_test_limit', contours[num_contours-1])
    corr_length_thin = ini.int('corr_length_thin', 0)
    corr_length_steps = ini.int('corr_length_steps', 15)
    
force_twotail = ini.bool('force_twotail', False)
if (force_twotail): print 'Computing two tail limits'

if (ini.params.has_key('cov_matrix_dimension')):
    covmat_dimension = len(mc.paramNames.list())
else:
    covmat_dimension = ini.int('cov_matrix_dimension', 0)
    if (covmat_dimension==-1):
        covmat_dimension = ncols - 2

plot_meanlikes = ini.bool('plot_meanlikes', False)

if (ini.params.has_key('do_minimal_1d_intervals')):
    print 'do_minimal_1d_intervals no longer used; set credible_interval_threshold instead'
    sys.exit()

PCA_num = ini.int('PCA_num',0)
if (PCA_num<>0):
    if (PCA_num<2):
        print 'Can only do PCA for 2 or more parameters'
        sys.exit()
    line = ini.string('PCA_params')
    PCA_func = ini.string('PCA_func')
    # Characters representing functional mapping
    if (PCA_func==''):
        PCA_func = ['N'] * PCA_num  # No mapping
    if (line.lower()=='all'):
        PCA_params = range(1, PCA_num+1)
    else:
        PCA_params = mc.index.keys() 
    line = ini.string('PCA_normparam')
    if (line==''):
        PCA_NormParam = 0
    else:
        PCA_NormParam = 123456789 # todo tmp_params(1)

num_3D_plots = ini.int('num_3D_plots',0)
#plot_3D = [ ini.string('3D_plot'+str(ix)) of ix in range(1, num_3D_plots+1) ]
plot_3D = []
for ix in range(1, num_3D_plots+1):
    plot_3D.append(ini.string('3D_plot'+str(ix)))

make_scatter_samples = ini.bool('make_scatter_samples', False)
max_scatter_points = ini.int('max_scatter_points', 2000)

BW = ini.bool('B&W', False)
do_shading = ini.bool('do_shading', True)

# Chain files
chain_files = glob.glob(in_root+'_*.txt')

def getLastChainIndex(in_root):
    if not chain_files: return 0
    names_files = [ os.path.basename(f) for f in chain_files ]
    basename = os.path.basename(in_root)
    indexes = [ int(f.replace(basename+'_', '').replace('.txt', '')) for f in names_files ]
    return max(indexes)

first_chain = ini.int('first_chain',1) 
last_chain = ini.int('chain_num',-1) 
# -1 means keep reading until one not found
if(last_chain==-1): last_chain = getLastChainIndex(in_root)

# Read in the chains
ok = mc.loadChains(in_root, chain_files)
#if (not ok): print ''
if (mc.numrows==0):
    print 'No un-ignored rows! (check number of chains/burn in)'
    sys.exit()

mc.makeSingle()

if (cool<>1):
    mc.CoolChain(cool)

# Adjust weights if requesteda
if (adjust_priors):
    mc.AdjustPriors()

mc.GetUsedCols()

mean_mult = mc.norm/mc.numrows
max_mult = (mean_mult*mc.numrows)/min(mc.numrows/2,500)
outliers = len(mc.weights[np.where(mc.weights>max_mult)])
if (outliers<>0):
    print 'outlier fraction ', float(outliers)/mc.numrows
max_mult = np.max(mc.weights)
numsamp = np.sum(mc.weights)

indep_thin = 0
if (not no_tests):
    #indep_thin modified here ...
    mc.DoConvergeTests(converge_test_limit)

adjust_priors = True
if (adjust_priors):
    mc.DeleteZeros()

print 'mean input multiplicity = ', mean_mult

# Output thinned data if requested
# Must do this with unsorted output
if (thin_factor<>0):
    mc.thin_indices(thin_factor) # ThinData => thin_indices
    filename = rootdirname + '_thin.txt'
    mc.WriteThinData(filename, thin_cool)

# Produce file of weight-1 samples if requested
if ((num_3D_plots<>0 and not make_single_samples or make_scatter_samples) and not no_plots):
    make_single_samples = True
    #single_thin = max(1, round(numsamp/max_mult)/max_scatter_points)

# Only use variables whose labels are not empty (and in list of plotparams if plotparams_num /= 0)
num_vars = 0

if (plotparams_num<>0):

    # ...
    pass

else:

    # ...
    pass


#
#mc.getChainsStats()

if (make_single_samples):
    filename = os.path.join(plot_data_dir, rootname.strip()+'_single.txt')
    mc.WriteSingleSamples(filename, single_thin)

# IO_WriteBounds
filename = os.path.join(plot_data_dir, rootname.strip()+'.bounds')
mc.WriteBounds(filename)

# Sort data in order of likelihood of points
mc.SortColData(1)

#numsamp = sum(coldata(1,0:nrows-1))

# Get ND confidence region (index into sorted coldata)
counts = 0
ND_cont1, ND_cont2 = mc.GetConfidenceRegion()

triangle_plot = triangle_plot and (num_vars>1)
if (triangle_plot):
    # ...
    pass

#print 'using ', nrows,' rows, processing ', num_vars,' parameters'
if (indep_thin<>0):
    print 'Approx indep samples: ', round(numsamp/indep_thin)
else:
    print  'effective number of samples (assuming indep): ', round(numsamp/max_mult)

# Get covariance matrix and correlation matrix
mc.GetCovMatrix()

if (PCA_num>0) and not plots_only:
    mc.PCA(PCA_params, PCA_num, PCA_func, PCA_NormParam)

# Find best fit, and mean likelihood
mc.GetChainLikeSummary(toStdOut=True)

LowerUpperLimits = 0

# Do 1D bins
for j in range(num_vars):

    #mc.Get1DDensity(j)
    
    #call to twoTailLimits in chains.py
    pass


if (not no_plots):
    # Output files for 1D plots
    filename = rootdirname + '.' + plot_ext
    textFileHandle = open(filename, 'w')
    if (plot_ext=='py'):
        text = 'g.plots_1d(roots)'
        # ...

    if (triangle_plot):
        filename = rootdirname + '_tri.' + plot_ext
        textFileHandle = open(filename, 'w')
        textInit = MCSamples.WritePlotFileInit()
        textFileHandle.write(textInit%(plot_data_dir, subplot_size_inch, out_dir, rootname))
        if (plot_ext=='py'):
            #todo: make string for tuple of names
            text = 'g.triangle_plot(roots, %s)'%data 
        textExport = MCSamples.WritePlotFileExport()
        fname = rootname + tag + '.ext???'
        textFileHandle.write(textExport%(fname))
    textFileHandle.close()
  



# Do 2D bins
if (plot_2D_param==0) and (num_cust2D_plots==0) and (not no_plots):
    # In this case output the most correlated variable combinations
    print 'doing 2D plots for most correlated variables'
    
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
    textInit = MCSamples.WritePlotFileInit()
    textFileHandle.write(textInit%(plot_data_dir, subplot_size_inch2, out_dir, rootname))
    if (plot_ext=='py'):
        textFileHandle.write('pairs=[]\n')
        for i in range(num_vars):
            textFileHandle.write("pairs.append(['%s','%s'])\n"%(name1, name2))
        textFileHandle.write('g.plots_2d(roots,param_pairs=pairs)\n')
        textExport = MCSamples.WritePlotFileExport()
        fname = rootname + '_2D.' + plot_ext
        textFileHandle.write(textExport%(fname))
    textFileHandle.close()
        

            
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
    textInit = MCSamples.WritePlotFileInit()
    textFileHandle.write(textInit%(plot_data_dir, subplot_size_inch3, out_dir, rootname))
    textFileHandle.write('sets=[]\n')
    for j in range(num_3D_plots):
        v1, v2, v3 = 0, 0, 0
        text += 'sets.append([%s,%s,%s])'%(v1, v2, v3) # todo
    text += 'g.plots_3d(roots,sets)'
    fname = rootname + '_3D.' + plot_ext
    textExport = MCSamples.WritePlotFileExport()
    textFileHandle.write(textExport%(fname))
    textFileHandle.close()


# Write out stats marginalized
if (not plots_only):
    mc.OutputMargeStats()

# Write paramNames file
filename = os.path.join(plot_data_dir, rootname + '.paramnames')
mc.WriteParamNames(filename)

# Limits from global likelihood
if (not plots_only):
    filename = rootdirname + '.likestats'
    textFileHandle = open(filename, 'w')
    textInit = mc.GetChainLikeSummary(toStdOut=False)
    textFileHandle.write(textInit)
    textFileHandle.write('param  bestfit        lower1         upper1         lower2         upper2')
    for j in range(num_vars):
        #todo: values
        textFileHandle.write('%5i%15.7E%15.7E%15.7E%15.7E%15.7E'%(j, v1, v2, v3, v4, v5))
    textFileHandle.close()

# System command
if (finish_run_command):
    finish_run_command = finish_run_command.replace('%ROOTNAME%',rootname)
    finish_run_command = finish_run_command.replace('%PLOTDIR%',plot_data_dir)
    finish_run_command = finish_run_command.replace('%PLOTROOT%', os.path.join(plot_data_dir, rootname))
    os.system(finish_run_command)
    
