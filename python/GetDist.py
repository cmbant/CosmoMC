#!/usr/bin/env python
import os
import sys
import iniFile
import subprocess
from getdist import MCSamples, chains
# ==============================================================================

if len(sys.argv) < 2:
    print 'Usage: python/GetDist.py ini_file [chain_root]'
    sys.exit()

# Parameter file
ini_file = os.path.abspath(sys.argv[1])
if not os.path.isfile(ini_file):
    print 'Parameter file does not exist ', ini_file
    sys.exit()

# Input parameters
ini = iniFile.iniFile(ini_file)

# File root
if len(sys.argv) > 2:
    in_root = sys.argv[2]
else:
    in_root = ini.params['file_root']
if in_root == '':
    print 'Root file does not exist ', in_root
    sys.exit()
rootname = os.path.basename(in_root)

ignorerows = ini.float('ignore_rows', 0.0)

# Create instance of MCSamples
mc = MCSamples.MCSamples(in_root)

mc.initParameters(ini)


if ini.bool('adjust_priors', False) or ini.bool('map_params', False):
    print 'To adjust priors or define new parameters, use a separate python script; see the python getdist docs for examples'
    sys.exit()

plot_ext = ini.string('plot_ext', 'py')
finish_run_command = ini.string('finish_run_command', '')

auto_label = ini.bool('auto_label', False)

prob_label = ini.bool('prob_label', False)

samples_are_chains = ini.bool('samples_are_chains', True)

no_plots = ini.bool('no_plots', False);
plots_only = ini.bool('plots_only', False)
no_tests = plots_only or ini.bool('no_tests', False)
make_plots = ini.bool('make_plots', False)

thin_factor = ini.int('thin_factor', 0)
thin_cool = ini.float('thin_cool', 1.0)

make_single_samples = ini.bool('make_single_samples', False)
single_thin = ini.int('single_thin', 1)
cool = ini.float('cool', 1.0)

chain_exclude = ini.int_list('exclude_chain')
num_exclude = len(chain_exclude) - chain_exclude.count(0)

shade_meanlikes = ini.bool('shade_meanlikes', False)
mc.shade_meanlikes = shade_meanlikes

out_dir = ini.string('out_dir')
if out_dir:
    print 'producing files in directory ', out_dir
mc.out_dir = out_dir

out_root = ini.string('out_root')
if out_root:
    rootname = out_root
    print 'producing files with with root ', out_root
mc.rootname = rootname

plot_data_dir = ini.string('plot_data_dir') or 'plot_data/'

abs_plot_data_dir = plot_data_dir
if not os.path.isdir(abs_plot_data_dir):
    os.mkdir(abs_plot_data_dir)
mc.plot_data_dir = plot_data_dir

rootdirname = os.path.join(out_dir, rootname); mc.rootdirname = rootdirname

if ini.params.has_key('do_minimal_1d_intervals'):
    print 'do_minimal_1d_intervals no longer used; set credible_interval_threshold instead'
    sys.exit()

line = ini.string('PCA_params', '')
if line.lower() == 'all':
    PCA_params = mc.paramNames.list()
else:
    PCA_params = line.split()
PCA_num = ini.int('PCA_num', len(PCA_params))
if PCA_num <> 0:
    if PCA_num < 2:
        print 'Can only do PCA for 2 or more parameters'
        sys.exit()
    PCA_func = ini.string('PCA_func', '')
    # Characters representing functional mapping
    if PCA_func == '':
        PCA_func = ['N'] * PCA_num  # No mapping
    PCA_NormParam = ini.string('PCA_normparam', '') or None

make_scatter_samples = ini.bool('make_scatter_samples', False)

# ==============================================================================

# Chain files
chain_files = chains.chainFiles(in_root)

def getLastChainIndex(in_root):
    if not chain_files: return 0
    names_files = [ os.path.basename(f) for f in chain_files ]
    basename = os.path.basename(in_root)
    indexes = [ int(f.replace(basename + '_', '').replace('.txt', '')) for f in names_files ]
    return max(indexes)

def runScript(fname):
    subprocess.Popen(['python', fname])

first_chain = ini.int('first_chain', 1)
last_chain = ini.int('chain_num', -1)
# -1 means keep reading until one not found
if last_chain == -1: last_chain = getLastChainIndex(in_root)

mc.loadChains(in_root, chain_files)

mc.removeBurnFraction(ignorerows)
mc.deleteFixedParams()
mc.makeSingle()


def filterPars(names):
    return [ name for name in names if mc.paramNames.parWithName(name) ]

if cool <> 1:
    print 'Cooling chains by ', cool
    mc.cool(cool)

plotparams = []
line = ini.string('plot_params', '')
if line not in ['', '0']:
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
if triangle_plot:
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
    plot_3D.append(pars)

mc.updateBaseStatistics()

if not no_tests:
    mc.getConvergeTests(mc.converge_test_limit, writeDataToFile=True, feedback=True)

mc.writeCovMatrix()
mc.writeCorrelationMatrix()

# Output thinned data if requested
# Must do this with unsorted output
if thin_factor <> 0:
    thin_ix = mc.thin_indices(thin_factor)
    filename = rootdirname + '_thin.txt'
    mc.WriteThinData(filename, thin_ix, thin_cool)

# Produce file of weight-1 samples if requested
if (num_3D_plots and not make_single_samples or make_scatter_samples) and not no_plots:
    make_single_samples = True
    single_thin = max(1, int(round(mc.norm / mc.max_mult)) / mc.max_scatter_points)

if make_single_samples:
    filename = os.path.join(plot_data_dir, rootname.strip() + '_single.txt')
    mc.makeSingleSamples(filename, single_thin)

print mc.getNumSampleSummaryText().strip()
print mc.likeStats.likeSummary().strip()

if PCA_num > 0 and not plots_only:
    mc.PCA(PCA_params, PCA_func, PCA_NormParam, writeDataToFile=True)

# Do 1D bins
mc.setDensitiesandMarge1D(writeDataToFile=not no_plots)

if not no_plots:
    # Output files for 1D plots
    print 'Calculating plot data...'

    mc.getBounds().saveToFile(os.path.join(plot_data_dir, rootname.strip() + '.bounds'))

    filename = rootdirname + '.' + plot_ext
    mc.WriteScriptPlots1D(filename, plotparams)
    if make_plots: runScript(filename)

    # Do 2D bins
    if plot_2D_param == 'corr':
        # In this case output the most correlated variable combinations
        print '...doing 2D plots for most correlated variables'
        cust2DPlots = mc.getCorrelatedVariable2DPlots()
        plot_2D_param = None
    elif plot_2D_param:
        mc.paramNames.parWithName(plot_2D_param, error=True)  # just check

    if cust2DPlots or plot_2D_param:
        print '...producing 2D plots'
        filename = rootdirname + '_2D.' + plot_ext
        mc.WriteScriptPlots2D(filename, plot_2D_param, cust2DPlots, plots_only)
        if make_plots: runScript(filename)

    if triangle_plot:
        # Add the off-diagonal 2D plots
        print '...producing triangle plot'
        filename = rootdirname + '_tri.' + plot_ext
        mc.WriteScriptPlotsTri(filename, triangle_params)
        for i in range(triangle_num):
            for i2 in range(i + 1, triangle_num):
                j = mc.index[triangle_params[i]]
                j2 = mc.index[triangle_params[i2]]
                if mc.done2D is None or not mc.done2D[j2][j] and not plots_only:
                    mc.get2DPlotData(j2, j, writeDataToFile=True)
        if make_plots: runScript(filename)

    # Do 3D plots (i.e. 2D scatter plots with coloured points)
    if num_3D_plots <> 0:
        print '...producing ', num_3D_plots, '2D colored scatter plots'
        filename = rootdirname + '_3D.' + plot_ext
        mc.WriteScriptPlots3D(filename, plot_3D)
        if make_plots: runScript(filename)

# Write paramNames file
    mc.paramNames.saveAsText(os.path.join(plot_data_dir, rootname + '.paramnames'))

if not plots_only:
# Write out stats marginalized
    mc.getMargeStats().saveAsText(rootdirname + '.margestats')

# Limits from global likelihood
    mc.getLikeStats().saveAsText(rootdirname + '.likestats')

# System command
if finish_run_command:
    finish_run_command = finish_run_command.replace('%ROOTNAME%', rootname)
    finish_run_command = finish_run_command.replace('%PLOTDIR%', plot_data_dir)
    finish_run_command = finish_run_command.replace('%PLOTROOT%', os.path.join(plot_data_dir, rootname))
    os.system(finish_run_command)

