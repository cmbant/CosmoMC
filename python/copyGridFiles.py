import os, fnmatch, shutil, batchJobArgs, batchJob, zipfile

Opts = batchJobArgs.batchArgs('copy chains and optionally other files', importance=True, converge=True)

Opts.parser.add_argument('target_dir', help="output root directory or zip file name")

Opts.parser.add_argument('--dist', action='store_true')
Opts.parser.add_argument('--chains', action='store_true', default=True)
Opts.parser.add_argument('--file_extensions', nargs='+', default=['.*'])
Opts.parser.add_argument('--skip_extensions', nargs='+', default=['.data', '.chk', '.log', '.corr', '.py', '.m'])
Opts.parser.add_argument('--dryrun', action='store_true')
Opts.parser.add_argument('--verbose', action='store_true')
Opts.parser.add_argument('--zip', action='store_true')

(batch, args) = Opts.parseForBatch()

if not args.chains and not args.dist: print 'use --chains or --dist!'

if '.zip' in args.target_dir: args.zip = True

if args.zip:
    zipper = zipfile.ZipFile(args.target_dir, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True)
else:
    target_dir = os.path.abspath(args.target_dir) + os.sep
    batchJob.makePath(target_dir)


def fileMatches(f, name):
    for ext in args.file_extensions:
        if fnmatch.fnmatch(f, name + ext):
            for ext2 in args.skip_extensions:
                if fnmatch.fnmatch(f, name + ext2): return False
            return True
    return False

def doCopy(source, dest, name):
    if args.verbose: print source + f
    if not args.dryrun:
        if args.zip:
            zipper.write(source + f, dest + f)
        else:
            shutil.copyfile(source + f, target_dir + dest + f)




for jobItem in Opts.filteredBatchItems():
    if (args.converge == 0 or jobItem.hasConvergeBetterThan(args.converge)):
        print jobItem.name
        chainfiles = 0
        infofiles = 0
        distfiles = 0
        outdir = jobItem.relativePath
        if args.chains:
            if not args.zip: batchJob.makePath(target_dir + outdir)
            i = 1
            while os.path.exists(jobItem.chainRoot + ('_%d.txt') % i):
                f = jobItem.name + ('_%d.txt') % i
                chainfiles += 1
                doCopy(jobItem.chainPath, outdir, f)
                i += 1
            for f in os.listdir(jobItem.chainPath):
                if fileMatches(f, jobItem.name):
                    infofiles += 1
                    if args.verbose: print jobItem.chainPath + f
                    doCopy(jobItem.chainPath, outdir, f)
        if args.dist:
            outdir += 'dist' + os.sep
            if not args.zip: batchJob.makePath(target_dir + outdir)
            for f in os.listdir(jobItem.distPath):
                if fileMatches(f, jobItem.name):
                    distfiles += 1
                    doCopy(jobItem.distPath, outdir, f)
        print '... %d chain files, %d other files and %d dist files' % (chainfiles, infofiles, distfiles)

if args.zip: zipper.close()
