#!/usr/bin/env python
import os, jobQueue, batchJobArgs

parser = batchJobArgs.argParser("Submit a single job to queue")

parser.add_argument('iniFile', nargs='+')

jobQueue.addArguments(parser)

args = parser.parse_args()

omp = jobQueue.getArgsOmp(args, msg=True, runsPerJob=len(args.iniFile))
ini = [ini.replace('.ini', '') for ini in args.iniFile]

if not args.dryrun:
    jobQueue.submitJob(os.path.basename(ini[0]), ini, omp=omp, **args.__dict__)
