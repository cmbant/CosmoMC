from __future__ import absolute_import
from __future__ import print_function
import os
import subprocess
import getdist
from getdist import IniFile
import time
from paramgrid import batchjob_args


def checkDir(fname):
    if not os.path.exists(fname): os.makedirs(fname)


Opts = batchjob_args.batchArgs('Run getdist over the grid of models', notExist=True)
Opts.parser.add_argument('--update_only', action='store_true')
Opts.parser.add_argument('--make_plots', action='store_true', help='run generated script plot files to make PDFs')
Opts.parser.add_argument('--norun', action='store_true')
Opts.parser.add_argument('--plot_data', default=None,
                         help="directory to store the plot_data in for each chain. Default None to generate on the fly.")
Opts.parser.add_argument('--burn_removed', action='store_true', help="if burn in has already been removed from chains")
Opts.parser.add_argument('--burn_remove', type=float,
                         help="fraction of chain to remove as burn in (if not importance sampled or already done)")

Opts.parser.add_argument('--no_plots', action='store_true',
                         help="just make non-plot outputs (faster if using old plot_data)")
Opts.parser.add_argument('--delay', type=int, help="run after delay of some number of seconds")
Opts.parser.add_argument('--procs', type=int, default=1, help="number of getdist instances to run in parallel")
Opts.parser.add_argument('--base_ini', default=getdist.default_getdist_settings, help="default getdist settings")
Opts.parser.add_argument('--command', default='python', help="program to run")
Opts.parser.add_argument('--command_params', nargs='*',
                         default=['python/GetDist.py'], help="arguments program to run (excl. ini name)")
Opts.parser.add_argument('--exist', action='store_true', help="Silently skip all chains that don't exist")

(batch, args) = Opts.parseForBatch()

plot_ext = 'py'
plot_cmd = 'python'

plot_types = ['.', '_2D.']
# you don't need these for python plots generated separately
# '_tri' is very slow for so many

if args.plot_data is None and getdist.use_plot_data:
    data_dir = os.path.abspath(args.plot_data) + os.sep
elif args.plot_data is not None:
    data_dir = os.path.abspath(args.plot_data) + os.sep
else:
    data_dir = None

if data_dir: checkDir(data_dir)

ini_dir = batch.batchPath + 'getdist' + os.sep

checkDir(ini_dir)

if args.delay: time.sleep(args.delay)
processes = set()

for jobItem in Opts.filteredBatchItems():
    ini = IniFile()
    ini.params['file_root'] = jobItem.chainRoot
    ini.params['batch_path'] = jobItem.batchPath
    checkDir(jobItem.distPath)
    ini.params['out_dir'] = jobItem.distPath
    if data_dir: ini.params['plot_data_dir'] = data_dir
    custom_plot = batch.commonPath + 'plots' + os.sep + jobItem.paramtag + '.ini'
    custom_plot2 = batch.commonPath + 'plots' + os.sep + jobItem.name + '.ini'
    if os.path.exists(custom_plot2):
        ini.includes.append(custom_plot2)
    elif os.path.exists(custom_plot):
        ini.includes.append(custom_plot)
    if os.path.exists(args.base_ini):
        ini.defaults.append(args.base_ini)
    elif os.path.exists(batch.commonPath + args.base_ini):
        ini.defaults.append(batch.commonPath + args.base_ini)
    else:
        raise ValueError("base_ini file not found")
    if hasattr(batch, 'getdist_options'):
        ini.params.update(batch.getdist_options)
    tag = ''
    if jobItem.isImportanceJob or args.burn_removed or jobItem.isBurnRemoved():
        ini.params['ignore_rows'] = 0
    elif args.burn_remove is not None:
        ini.params['ignore_rows'] = args.burn_remove

    if jobItem.isImportanceJob:
        ini.params['compare_num'] = 1
        ini.params['compare1'] = jobItem.parent.chainRoot
    if args.no_plots: ini.params['no_plots'] = True
    if args.make_plots: ini.params['make_plots'] = True
    fname = ini_dir + jobItem.name + tag + '.ini'
    ini.params.update(jobItem.dist_settings)
    ini.saveFile(fname)
    if not args.norun and (not args.notexist or not jobItem.getDistExists()) and (
                not args.update_only or jobItem.getDistNeedsUpdate()):
        if jobItem.chainExists():
            print("running: " + fname)
            processes.add(subprocess.Popen([args.command] + args.command_params + [fname]))
            while len(processes) >= args.procs:
                time.sleep(.1)
                processes.difference_update([p for p in processes if p.poll() is not None])
        else:
            if not args.exist: print("Chains do not exist yet: " + jobItem.chainRoot)
