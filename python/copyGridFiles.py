from __future__ import absolute_import
from __future__ import print_function
import os
import fnmatch
import shutil
import zipfile
from datetime import datetime, timedelta

from paramgrid import batchjob, batchjob_args

Opts = batchjob_args.batchArgs('copy or zip chains and optionally other files', importance=True, converge=True)

Opts.parser.add_argument('target_dir', help="output root directory or zip file name")

Opts.parser.add_argument('--dist', action='store_true', help="include getdist outputs")
Opts.parser.add_argument('--chains', action='store_true', help="include chain files")
Opts.parser.add_argument('--sym_link', action='store_true', help="just make symbolic links to source directories")
Opts.parser.add_argument('--no_config', action='store_true', help="don't copy grid config info")

Opts.parser.add_argument('--remove_burn_fraction', default=0.0, type=float,
                         help="fraction at start of chain to remove as burn in")

Opts.parser.add_argument('--file_extensions', nargs='+', default=['.*'], help='extensions to include')
Opts.parser.add_argument('--skip_extensions', nargs='+',
                         default=['.data', '.chk', '.chk_tmp', '.log', '.corr', '.py', '.m', '.py_mcsamples',
                                  '.pysamples'])
Opts.parser.add_argument('--max_age_days', default=0.0, type=float,
                         help="only include files with date stamp at most max_age_days old")
Opts.parser.add_argument('--dryrun', action='store_true')
Opts.parser.add_argument('--verbose', action='store_true')
Opts.parser.add_argument('--zip', action='store_true',
                         help='make a zip file. Not needed if target_dir is a filename ending in .zip')

(batch, args) = Opts.parseForBatch()

if '.zip' in args.target_dir: args.zip = True
if args.max_age_days:
    max_age = datetime.now() - timedelta(days=args.max_age_days)

sizeMB = 0

if args.zip:
    zipper = zipfile.ZipFile(args.target_dir, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True)
    target_dir = None
else:
    zipper = None
    target_dir = os.path.abspath(args.target_dir) + os.sep
    batchjob.makePath(target_dir)

if args.sym_link and (args.remove_burn_fraction or args.zip): raise Exception('option not compatible with --sym_link')


def fileMatches(f, name):
    for ext in args.file_extensions:
        if fnmatch.fnmatch(f, name + ext):
            for ext2 in args.skip_extensions:
                if fnmatch.fnmatch(f, name + ext2): return False
            return True
    return False


def doCopy(source, dest, f, hasBurn=False):
    global sizeMB
    if args.verbose: print(source + f)
    frac = 1
    if not args.dryrun:
        if args.remove_burn_fraction and hasBurn:
            lines = open(source + f).readlines()
            lines = lines[int(len(lines) * args.remove_burn_fraction):]
            frac = 1 - args.remove_burn_fraction
        else:
            lines = None
        destf = dest + f
        if args.zip:
            if lines:
                zipper.writestr(destf, "".join(lines))

            else:
                zipper.write(source + f, destf)
        else:
            if lines:
                open(target_dir + destf, 'w').writelines(lines)
            else:
                if args.sym_link:
                    if os.path.islink(target_dir + destf): os.unlink(target_dir + destf)
                    os.symlink(os.path.realpath(source + f), target_dir + destf)
                else:
                    shutil.copyfile(source + f, target_dir + destf)
    elif args.remove_burn_fraction and hasBurn:
        frac = 1 - args.remove_burn_fraction
    sizeMB += os.path.getsize(source + f) / 1024. ** 2 * frac


def writeIni(iniName, props):
    if args.dryrun: return
    if args.zip:
        zipper.writestr(outdir + iniName, "\n".join(props.fileLines()))
    else:
        props.saveFile(target_dir + outdir + iniName)


if not args.no_config:
    config_path = os.path.join(batch.batchPath, 'config/')
    if not args.dryrun and not args.zip: batchjob.makePath(target_dir + 'config')
    for f in os.listdir(config_path):
        doCopy(config_path, 'config/', f)

for jobItem in Opts.filteredBatchItems():
    if args.converge == 0 or jobItem.hasConvergeBetterThan(args.converge):
        print(jobItem.name)
        chainfiles = 0
        infofiles = 0
        distfiles = 0
        doneProperties = False
        outdir = jobItem.relativePath
        if not args.zip: batchjob.makePath(target_dir + outdir)
        if args.chains and jobItem.chainExists() and (not args.max_age_days or datetime.fromtimestamp(jobItem.chainFileDate()) > max_age):
            i = 1
            while os.path.exists(jobItem.chainRoot + '_%d.txt' % i):
                f = jobItem.name + '_%d.txt' % i
                chainfiles += 1
                doCopy(jobItem.chainPath, outdir, f, not jobItem.isImportanceJob)
                i += 1
            if not jobItem.isImportanceJob and args.remove_burn_fraction:
                props = jobItem.propertiesIni()
                props.params['burn_removed'] = True
                writeIni(jobItem.name + '.properties.ini', props)
                doneProperties = True

        for f in os.listdir(jobItem.chainPath):
            if fileMatches(f, jobItem.name):
                if doneProperties and '.properties.ini' in f: continue
                if not args.max_age_days or datetime.fromtimestamp(os.path.getmtime(jobItem.chainPath + f)) > max_age:
                    infofiles += 1
                    if args.verbose: print(jobItem.chainPath + f)
                    doCopy(jobItem.chainPath, outdir, f)
        if args.dist and os.path.exists(jobItem.distPath):
            outdir += 'dist' + os.sep
            if not args.zip: batchjob.makePath(target_dir + outdir)
            for f in os.listdir(jobItem.distPath):
                if fileMatches(f, jobItem.name) and (not args.max_age_days or
                                                     datetime.fromtimestamp(
                                                         os.path.getmtime(jobItem.distPath + f)) > max_age):
                    distfiles += 1
                    doCopy(jobItem.distPath, outdir, f)
        print('... %d chain files, %d other files and %d dist files' % (chainfiles, infofiles, distfiles))

if zipper: zipper.close()

print('Total size: %u MB' % sizeMB)
