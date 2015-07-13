from __future__ import absolute_import
from __future__ import print_function
import os
import fnmatch
import shutil
from paramgrid import batchjob_args

Opts = batchjob_args.batchArgs('copy all files of a given type from all getdist output directories in the batch',
                               importance=True, converge=True)

Opts.parser.add_argument('target_dir')
Opts.parser.add_argument('file_extension', nargs='+')
Opts.parser.add_argument('--normalize_names', action='store_true', help='replace actual name tags with normalized names')
Opts.parser.add_argument('--tag_replacements', nargs='+', help="XX YY XX2 YY2 replaces name XX with YY, XX2 with YY2 etc.")


(batch, args) = Opts.parseForBatch()

target_dir = os.path.abspath(args.target_dir) + os.sep
if not os.path.exists(target_dir): os.makedirs(target_dir)

if args.tag_replacements is not None:
    replacements = dict()
    for i, val in enumerate(args.tag_replacements[::2]):
        replacements[val] = args.tag_replacements[i * 2 + 1]
else: replacements = None

for ext in args.file_extension:
    if not '.' in ext: pattern = '.' + ext
    else: pattern = ext
    for jobItem in Opts.filteredBatchItems():
        if os.path.exists(jobItem.distPath) and (args.converge == 0 or jobItem.hasConvergeBetterThan(args.converge)):
            for f in os.listdir(jobItem.distPath):
                if fnmatch.fnmatch(f, jobItem.name + pattern):
                    print(jobItem.distPath + f)
                    if args.normalize_names:
                        fout = jobItem.makeNormedName(replacements)[0] + os.path.splitext(f)[1]
                    else:
                        fout = f
                    shutil.copyfile(jobItem.distPath + f, target_dir + fout)
