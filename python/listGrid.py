from __future__ import absolute_import
from __future__ import print_function
from paramgrid import batchjob_args


Opts = batchjob_args.batchArgs('List items in a grid', importance=True, converge=True, notExist=True)
Opts.parser.add_argument('--exists', action='store_true', help='chain must exist')
Opts.parser.add_argument('--normed', action='store_true', help='Output normed names')

(batch, args) = Opts.parseForBatch()
items = Opts.sortedParamtagDict(chainExist=args.exists)

for paramtag, parambatch in items:
    for jobItem in parambatch:
        if hasattr(jobItem, 'group'):
            tag = '(%s)' % jobItem.group
        else:
            tag = ''
        if args.normed:
            print(jobItem.normed_name, tag)
        else:
            print(jobItem.name, tag)
