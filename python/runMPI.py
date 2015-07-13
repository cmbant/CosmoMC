#!/usr/bin/env python

import os
from paramgrid import batchjob_args, jobqueue


parser = batchjob_args.argParser("Submit a single job to queue")

parser.add_argument('iniFile', nargs='+')

jobqueue.addArguments(parser)

args = parser.parse_args()

ini = [ini.replace('.ini', '') for ini in args.iniFile]

jobqueue.submitJob(os.path.basename(ini[0]), ini, msg=True, **args.__dict__)
