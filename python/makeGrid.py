#!/usr/bin/env python

from paramgrid import gridconfig

args = gridconfig.getArgs()
args.interactive = True
gridconfig.makeGrid(**args.__dict__)
