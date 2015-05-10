from __future__ import absolute_import
from __future__ import print_function
import os
import re
from paramgrid import batchjob_args
from getdist import paramnames


Opts = batchjob_args.batchArgs('rename parameter in all .paramnames files in grid', importance=True)
Opts.parser.add_argument('--old_new', nargs='+', help="list of oldname newname oldname2 newname2...")
Opts.parser.add_argument('--labelNames', default=None, help=".paramnames file for new param labels")
Opts.parser.add_argument('--map_file', help="file with rows of oldname newname label")
Opts.parser.add_argument('--confirm', action='store_true', help="true to replace .paramnames files")


(batch, args) = Opts.parseForBatch()

if args.old_new and len(args.old_new) < 2: raise Exception('Must have at least one pair of parameters to rename')

if args.labelNames:
    labels = paramnames.ParamNames(args.labelNames)
else:
    labels = None

mapper = dict()
if args.map_file:
    with open(args.map_file) as f:
        for line in f:
            if line.strip():
                old, new, label = [s.strip() for s in line.split(None, 2)]
                mapper[old] = (new, label)
if args.old_new:
    for old, new in zip(args.old_new[::2], args.old_new[1::2]):
        mapper[old] = (new, None)


for jobItem in Opts.filteredBatchItems():
    name = jobItem.chainRoot + '.paramnames'
    if os.path.exists(name):
        names = paramnames.ParamNames(name)
        has = False
        if mapper:
            for p in names.names:
                new = mapper.get(p.name, None)
                if new:
                    p.name, lab = new
                    if lab is not None: p.label = lab
                    has = True
        if labels:
            for p in names.names:
                plab = labels.parWithName(p.name)
                if plab:
                    has = True
                    p.label = plab.label
        if has:
            print(jobItem.chainRoot)
            if args.confirm: names.saveAsText(name)
        bestfit = jobItem.chainRoot + '.minimum'
        if os.path.exists(bestfit):
            lines = open(bestfit).readlines()
            has = False
            for i, line in enumerate(lines):
                line = line.rstrip()
                if i < 3: continue
                items = re.split(r'(\S+)', line, maxsplit=4)
                if len(items) > 1 and items[1] == '-log(Like)': break
                if len(items) < 6: continue
                name = items[5].strip()
                new = mapper.get(name, None)
                lab = None
                if new:
                    newname, lab = new
                    items[5] = newname
                    items[6] = ' ' * (len(items[6]) + len(name) - len(newname))
                    has = has or newname != name
                if labels:
                    plab = labels.parWithName(name)
                    if plab: lab = plab.label
                if lab is not None:
                    if len(items) < 7:
                        items.append(lab)
                        has = True
                    else:
                        has = has or "".join(items[7:]) != lab
                        items[7] = lab
                        del items[8:]
                lines[i] = "".join(items) + "\n"
            if has:
                print(bestfit)
                if args.confirm: open(bestfit).write("".join(lines))

if args.confirm:
    print('Done. Re-run getdist to update getdist outputs.')
else:
    print('... run with --confirm to actually replace .paramname files')
