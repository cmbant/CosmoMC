# GetDist.py

import os
import sys
import iniFile
import chains
import MCSamples

# ==============================================================================

if (len(sys.argv) < 2):
    print 'Usage: python/GetDist.py ini_file [chain_root]'
    sys.exit()

# Parameter file
ini_file = os.path.abspath(sys.argv[1])
if not os.path.isfile(ini_file):
    print 'Parameter file does not exist ', ini_file
    sys.exit()

# Input parameters
ini = iniFile.iniFile()
ini.readFile(ini_file)

# File root
if (len(sys.argv) > 2):
    in_root = sys.argv[2]
else:
    in_root = ini.params['file_root']
if (in_root == ''):  # or (not os.path.isfile(file_root)):
    print 'Root file does not exist ', in_root
    sys.exit()
rootname = os.path.basename(in_root)

ignorerows = ini.float('ignore_rows', 0.0)

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


mc.initParameters(ini)


adjust_priors = ini.bool('adjust_priors', False)

plot_ext = ini.string('plot_ext', 'py')

plot_output = ini.string('plot_output', 'pdf')

font_scale = ini.float('font_scale', 1.0)
finish_run_command = ini.string('finish_run_command', '')

auto_label = ini.bool('auto_label', False)

prob_label = ini.bool('prob_label', False)

samples_are_chains = ini.bool('samples_are_chains', True)

no_plots = ini.bool('no_plots', False); mc.no_plots = no_plots
plots_only = ini.bool('plots_only', False)
no_tests = plots_only or ini.bool('no_tests', False)
line_labels = ini.bool('line_labels', False)

thin_factor = ini.int('thin_factor', 0)
thin_cool = ini.float('thin_cool', 1.0)

make_single_samples = ini.bool('make_single_samples', False)
single_thin = ini.int('single_thin', 1)
cool = ini.float('cool', 1.0)

# Compute limits
mc.initLimits(ini)

if (ini.params.has_key('plotparams_num')):
    print 'plotparams_num deprectated; just use plot_params'
    sys.exit()

plotparams = []
line = ini.string('plot_params', 0)
if (line not in ['', 0]):
    plotparams_num = -1
    plotparams = [ s for s in line.split(' ') if s <> '' ]
    plotparams = [ mc.index[name] for name in plotparams ]
else:
    plotparams_num = 0

line = ini.string('plot_2D_param')
if (line == ''):
    plot_2D_param = 0
else:
    tmp_params = [ s for s in line.split(' ') if s <> '' ]
    plot_2D_param = int(tmp_params[0])
    if (plot_2D_param <> 0 and plotparams_num <> 0 and plotparams.count(plot_2D_param) == 0):
        print 'plot_2D_param not in plotparams'
        sys.exit()

if (plot_2D_param <> 0):
    plot_2D_param = plot_2D_param + 2
    num_cust2D_plots = 0
else:
    # Use custom array of specific plots
    num_cust2D_plots = ini.int('plot_2D_num', 0)
    cust2DPlots = []
    for i in range(1, num_cust2D_plots + 1):
        line = ini.string('plot' + str(i))
        tmp_params = [ s for s in line.split(' ') if s <> '' ]
        tmp_params = [ mc.index[name] for name in tmp_params ]
        if (plotparams_num <> 0 and plotparams.count(tmp_params[0]) == 0):
            print 'plot' + str(i), ': parameter not in plotparams'
            sys.exit()
        cust2DPlots.append(tmp_params[0] + 2 + (tmp_params[1] + 2) * 1000)

triangle_params = []
triangle_plot = ini.bool('triangle_plot', False)
if (triangle_plot):
    no_triangle_axis_labels = ini.bool('no_triangle_axis_labels', False)
    line = ini.string('triangle_params')
    triangle_num = -1
    if (line <> ''):
        triangle_params = [ s for s in line.split(' ') if s <> '' ]
        triangle_params = [ mc.index[name] for name in triangle_params if mc.index.has_key(name) ]
        triangle_num = len(triangle_params)

exclude_chain = ini.string('exclude_chain')
chain_exclude = [ int(s) for s in exclude_chain.split(' ') if s <> '' ]
num_exclude = len(chain_exclude) - chain_exclude.count(0)

map_params = ini.bool('map_params', False)
if (map_params):
    print 'WARNING: Mapping params - .covmat file is new params.'

shade_meanlikes = ini.bool('shade_meanlikes', False)
mc.shade_meanlikes = shade_meanlikes

out_dir = ini.string('out_dir')
if (out_dir <> ''):
    print 'producing files in directory ', out_dir
mc.out_dir = out_dir

out_root = ini.string('out_root')
if (out_root <> ''):
    rootname = out_root
    print 'producing files with with root ', out_root
mc.rootname = rootname

plot_data_dir = ini.string('plot_data_dir')
if (plot_data_dir == ''):
    plot_data_dir = 'plot_data/'

abs_plot_data_dir = plot_data_dir
if not os.path.isdir(abs_plot_data_dir):
    os.mkdir(abs_plot_data_dir)
mc.plot_data_dir = plot_data_dir

rootdirname = os.path.join(out_dir, rootname); mc.rootdirname = rootdirname

mc.initContours(ini)

if (not no_tests):
    converge_test_limit = ini.float('converge_test_limit', mc.contours[mc.num_contours - 1])
    corr_length_thin = ini.int('corr_length_thin', 0)
    corr_length_steps = ini.int('corr_length_steps', 15)
    mc.corr_length_thin = corr_length_thin
    mc.corr_length_steps = corr_length_steps

if (ini.params.has_key('cov_matrix_dimension')):
    covmat_dimension = len(mc.paramNames.list())
else:
    covmat_dimension = ini.int('cov_matrix_dimension', 0)
    if (covmat_dimension == -1):
        covmat_dimension = ncols - 2
mc.covmat_dimension = covmat_dimension

if (ini.params.has_key('do_minimal_1d_intervals')):
    print 'do_minimal_1d_intervals no longer used; set credible_interval_threshold instead'
    sys.exit()

PCA_num = ini.int('PCA_num', 0)
if (PCA_num <> 0):
    if (PCA_num < 2):
        print 'Can only do PCA for 2 or more parameters'
        sys.exit()
    line = ini.string('PCA_params')
    PCA_func = ini.string('PCA_func')
    # Characters representing functional mapping
    if (PCA_func == ''):
        PCA_func = ['N'] * PCA_num  # No mapping
    if (line.lower() == 'all'):
        PCA_params = range(1, PCA_num + 1)
    else:
        names = [ s for s in line.split(' ') if s <> '' ]
        PCA_params = [ mc.index[name] for name in names ]
    line = ini.string('PCA_normparam')
    if (line == ''):
        PCA_NormParam = 0
    else:
        tmp_params = [ s for s in line.split(' ') if s <> '' ]
        tmp_params = [ mc.index[name] for name in tmp_params]
        PCA_NormParam = int(tmp_params[0])

num_3D_plots = ini.int('num_3D_plots', 0)
plot_3D = []
for ix in range(1, num_3D_plots + 1):
    line = ini.string('3D_plot' + str(ix))
    plot_3D.append([ s for s in line.split(' ') if s <> '' ])

make_scatter_samples = ini.bool('make_scatter_samples', False)

BW = ini.bool('B&W', False)
do_shading = ini.bool('do_shading', True)

# ==============================================================================

# Chain files
chain_files = chains.chainFiles(in_root)

def getLastChainIndex(in_root):
    if not chain_files: return 0
    names_files = [ os.path.basename(f) for f in chain_files ]
    basename = os.path.basename(in_root)
    indexes = [ int(f.replace(basename + '_', '').replace('.txt', '')) for f in names_files ]
    return max(indexes)

first_chain = ini.int('first_chain', 1)
last_chain = ini.int('chain_num', -1)
# -1 means keep reading until one not found
if(last_chain == -1): last_chain = getLastChainIndex(in_root)

# Read in the chains
ok = mc.loadChains(in_root, chain_files)
# if (not ok): print ''
# if (mc.numrows==0):
#    print 'No un-ignored rows! (check number of chains/burn in)'
#    sys.exit()

mc.removeBurnFraction(ignorerows)
if (not no_tests):
    mc.DoConvergeTests(converge_test_limit)

#
mc.makeSingle()

if (cool <> 1):
    mc.CoolChain(cool)

# Adjust weights if requested
if (adjust_priors):
    mc.AdjustPriors()

# See which parameters are fixed
mc.GetUsedCols()

mc.ComputeMultiplicators()

if (adjust_priors):
    mc.DeleteZeros()

print 'mean input multiplicity = ', mc.mean_mult

# Output thinned data if requested
# Must do this with unsorted output
if (thin_factor <> 0):
    thin_ix = mc.thin_indices(thin_factor)
    filename = rootdirname + '_thin.txt'
    mc.WriteThinData(filename, thin_ix, thin_cool)

# Produce file of weight-1 samples if requested
if ((num_3D_plots <> 0 and not make_single_samples or make_scatter_samples) and not no_plots):
    make_single_samples = True
    single_thin = max(1, int(round(mc.numsamp / mc.max_mult)) / mc.max_scatter_points)

# Compute means and std dev.
mc.ComputeStats()

if (make_single_samples):
    filename = os.path.join(plot_data_dir, rootname.strip() + '_single.txt')
    mc.MakeSingleSamples(filename, single_thin)

# IO_WriteBounds
filename = os.path.join(plot_data_dir, rootname.strip() + '.bounds')
mc.WriteBounds(filename)

# Sort data in order of likelihood of points
mc.SortColData(1)

mc.ComputeNumSamp()

# Get ND confidence region (index into sorted coldata)
counts = 0
mc.GetConfidenceRegion()

triangle_plot = triangle_plot and (mc.num_vars > 1)
if (triangle_plot):
    if (triangle_num == -1):
        triangle_num = mc.num_vars
        triangle_params = mc.paramNames.list()
    else:
        ix = triangle_num

        triangle_plot = triangle_num > 1

num_parameters = mc.isused[mc.isused == True].shape[0]
print 'using ', mc.numrows, ' rows, processing ', num_parameters, ' parameters'
if (mc.indep_thin <> 0):
    print 'Approx indep samples: ', round(mc.numsamp / mc.indep_thin)
else:
    print  'effective number of samples (assuming indep): ', round(mc.numsamp / mc.max_mult)

# Get covariance matrix and correlation matrix
mc.GetCovMatrix()

if (PCA_num > 0) and not plots_only:
    mc.PCA(PCA_params, PCA_func, PCA_NormParam)

# Find best fit, and mean likelihood
mc.GetChainLikeSummary(toStdOut=True)

# Initialize variables for 1D bins
mc.Init1DDensity()

# Do 1D bins
mc.Do1DBins(mc.max_frac_twotail)

if (not no_plots):
    # Output files for 1D plots
    filename = rootdirname + '.' + plot_ext
    mc.WriteScriptPlots1D(filename)

    if (triangle_plot):
        filename = rootdirname + '_tri.' + plot_ext
        mc.WriteScriptPlotsTri(filename, triangle_params)

# Do 2D bins
if (plot_2D_param == 0) and (num_cust2D_plots == 0) and (not no_plots):
    # In this case output the most correlated variable combinations
    print 'doing 2D plots for most correlated variables'
    num_cust2D_plots_0 = 12
    cust2DPlots, num_cust2D_plots = mc.GetCust2DPlots(num_cust2D_plots_0)

if (num_cust2D_plots == 0):
    num_2D_plots = 0
    for j in range(mc.num_vars):
        if (mc.ix_min[j] <> mc.ix_max[j]):
            for j2 in range(j + 1, mc.num_vars):
                if (mc.ix_min[j2] <> mc.ix_max[j2]):
                    if (plot_2D_param in [0, j, j2]):
                        num_2D_plots += 1
else:
    num_2D_plots = num_cust2D_plots

if ((num_2D_plots > 0) and (not no_plots)):
    print 'Producing ', num_2D_plots, ' 2D plots'
    filename = rootdirname + '_2D.' + plot_ext
    mc.WriteScriptPlots2D(filename, plot_2D_param, num_cust2D_plots, cust2DPlots, plots_only)

if (triangle_plot and not no_plots):
    # Add the off-diagonal 2D plots
    for i in range(triangle_num):
        for i2 in range(i + 1, triangle_num):
            j = triangle_params[i]
            j2 = triangle_params[i2]
            if (not mc.isused[j]) or (not mc.isused[j2]): continue
            if (not mc.done2D[j2][j] and not plots_only): mc.Get2DPlotData(j2, j)

# Do 3D plots (i.e. 2D scatter plots with coloured points)
if (num_3D_plots <> 0 and not no_plots):
    print 'producing ', num_3D_plots, '2D colored scatter plots'
    filename = rootdirname + '_3D.' + plot_ext
    mc.WriteScriptPlots3D(filename, num_3D_plots, plot_3D)

# Write out stats marginalized
if (not plots_only):
    mc.OutputMargeStats()

# Write paramNames file
filename = os.path.join(plot_data_dir, rootname + '.paramnames')
mc.WriteParamNames(filename)

# Limits from global likelihood
if (not plots_only):
    filename = rootdirname + '.likestats'
    mc.WriteGlobalLikelihood(filename)

# System command
if (finish_run_command):
    finish_run_command = finish_run_command.replace('%ROOTNAME%', rootname)
    finish_run_command = finish_run_command.replace('%PLOTDIR%', plot_data_dir)
    finish_run_command = finish_run_command.replace('%PLOTROOT%', os.path.join(plot_data_dir, rootname))
    os.system(finish_run_command)

