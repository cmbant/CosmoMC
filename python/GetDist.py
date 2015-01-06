# GetDist.py

import os
import sys
import iniFile
import chains
import MCSamples
import numpy as np
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

single_column_chain_files = ini.bool('single_column_chain_files', False)

mc.initParameters(ini)


adjust_priors = ini.bool('adjust_priors', False)

plot_ext = ini.string('plot_ext', 'py')

plot_output = ini.string('plot_output', 'pdf')

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
        PCA_params = mc.paramNames.list()
    else:
        PCA_params = line.split()
    line = ini.string('PCA_normparam')
    PCA_NormParam = line or None

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
mc.deleteFixedParams()

if (not no_tests):
    mc.DoConvergeTests(converge_test_limit)

mc.makeSingle()

def filterPars(names):
    return [ name for name in names if mc.paramNames.parWithName(name) ]

if (cool <> 1):
    mc.CoolChain(cool)

# Adjust weights if requested
if (adjust_priors):
    mc.AdjustPriors()

plotparams = []
line = ini.string('plot_params', '')
if (line not in ['', 0]):
    plotparams = filterPars(line.split())

line = ini.string('plot_2D_param', '')
plot_2D_param = None
if line.strip() and line <> '0':
    plot_2D_param = line.strip()

cust2DPlots = []
if not plot_2D_param:
    # Use custom array of specific plots
    num_cust2D_plots = ini.int('plot_2D_num', 0)
    for i in range(1, num_cust2D_plots + 1):
        line = ini.string('plot' + str(i))
        pars = filterPars(line.split())
        if len(pars) <> 2: raise Exception('plot_2D_num parameter not found, not varied, or not wrong number of parameters')
        cust2DPlots.append(pars)

triangle_params = []
triangle_plot = ini.bool('triangle_plot', False)
if (triangle_plot):
    no_triangle_axis_labels = ini.bool('no_triangle_axis_labels', False)
    line = ini.string('triangle_params')
    if line: triangle_params = filterPars(line.split())
    triangle_num = len(triangle_params)
    triangle_plot = triangle_num > 1

num_3D_plots = ini.int('num_3D_plots', 0)
plot_3D = []
for ix in range(1, num_3D_plots + 1):
    line = ini.string('3D_plot' + str(ix))
    pars = filterPars(line.split())
    if len(pars) <> 3: raise Exception('3D_plot parameter not found, not varied, or not wrong number of parameters')
    plot_3D.append(line.split)


if (adjust_priors):
    mc.DeleteZeros()

mc.updateChainBaseStatistics(ini)

mc.writeCovMatrix()
mc.writeCorrMatrix()

# Output thinned data if requested
# Must do this with unsorted output
if (thin_factor <> 0):
    thin_ix = mc.thin_indices(thin_factor)
    filename = rootdirname + '_thin.txt'
    mc.WriteThinData(filename, thin_ix, thin_cool)

# Produce file of weight-1 samples if requested
if ((num_3D_plots and not make_single_samples or make_scatter_samples) and not no_plots):
    make_single_samples = True
    single_thin = max(1, int(round(mc.numsamp / mc.max_mult)) / mc.max_scatter_points)

if (make_single_samples):
    filename = os.path.join(plot_data_dir, rootname.strip() + '_single.txt')
    mc.MakeSingleSamples(filename, single_thin)

print 'mean input multiplicity = ', mc.mean_mult

num_parameters = mc.paramNames.numParams()
print 'using ', mc.numrows, ' rows, processing ', num_parameters, ' parameters'
if (mc.indep_thin <> 0):
    print 'Approx indep samples: ', round(mc.numsamp / mc.indep_thin)
else:
    print  'effective number of samples (assuming indep): ', round(mc.numsamp / mc.max_mult)


filename = os.path.join(plot_data_dir, rootname.strip() + '.bounds')
mc.writeBounds(filename)

if (PCA_num > 0) and not plots_only:
    mc.PCA(PCA_params, PCA_func, PCA_NormParam, writeDataToFile=True)

# Do 1D bins
mc.Do1DBins(mc.max_frac_twotail, writeDataToFile=True)

if not no_plots:
    # Output files for 1D plots
    print 'Calculating plot data...'
    filename = rootdirname + '.' + plot_ext
    mc.WriteScriptPlots1D(filename, plotparams, ext=plot_output)

# Do 2D bins
if plot_2D_param == 'corr' and not no_plots:
    # In this case output the most correlated variable combinations
    print 'doing 2D plots for most correlated variables'
    num_cust2D_plots_0 = 12
    cust2DPlots = mc.GetCust2DPlots(num_cust2D_plots_0)
    plot_2D_param = None
elif plot_2D_param:
    mc.paramNames.parWithName(plot_2D_param, error=True)  # just check

if (cust2DPlots or plot_2D_param) and  not no_plots:
    filename = rootdirname + '_2D.' + plot_ext
    mc.WriteScriptPlots2D(filename, plot_2D_param, cust2DPlots, plots_only, ext=plot_output)

if (triangle_plot and not no_plots):
    # Add the off-diagonal 2D plots
    mc.WriteScriptPlotsTri(rootdirname + '_tri.' + plot_ext, triangle_params, ext=plot_output)
    for i in range(triangle_num):
        for i2 in range(i + 1, triangle_num):
            j = mc.index[triangle_params[i]]
            j2 = mc.index[triangle_params[i2]]
            if (mc.done2D is None or not mc.done2D[j2][j]) and not plots_only: mc.Get2DPlotData(j2, j, writeDataToFile=True)

# Do 3D plots (i.e. 2D scatter plots with coloured points)
if (num_3D_plots <> 0 and not no_plots):
    print 'producing ', num_3D_plots, '2D colored scatter plots'
    filename = rootdirname + '_3D.' + plot_ext
    mc.WriteScriptPlots3D(filename, plot_3D, ext=plot_output)

# Write out stats marginalized
if not plots_only: mc.saveMargeStats()

# Write paramNames file
mc.paramNames.saveAsText(os.path.join(plot_data_dir, rootname + '.paramnames'))

# Limits from global likelihood
if not plots_only:
    mc.WriteGlobalLikelihood(rootdirname + '.likestats')

# System command
if (finish_run_command):
    finish_run_command = finish_run_command.replace('%ROOTNAME%', rootname)
    finish_run_command = finish_run_command.replace('%PLOTDIR%', plot_data_dir)
    finish_run_command = finish_run_command.replace('%PLOTROOT%', os.path.join(plot_data_dir, rootname))
    os.system(finish_run_command)

