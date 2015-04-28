#!/usr/bin/env python
import os

from paramgrid import batchJobArgs, jobQueue


parser = batchJobArgs.argParser("Submit a single job to queue")

parser.add_argument('iniFile', nargs='+')

jobQueue.addArguments(parser)

args = parser.parse_args()

ini = [ini.replace('.ini', '') for ini in args.iniFile]

jobQueue.submitJob(os.path.basename(ini[0]), ini, msg=True, **args.__dict__)
