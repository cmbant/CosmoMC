from __future__ import absolute_import
from __future__ import print_function
import os
import sys
from paramgrid import batchjob

if len(sys.argv) < 3:
    print('Usage: python/addGridBatch.py directory_with_outputs directory_with_output_to_add [and_another..]')
    sys.exit()

batch = batchjob.readobject()

for subBatch in sys.argv[2:]:
    batchPath2 = os.path.abspath(subBatch) + os.sep
    batch2 = batchjob.readobject(batchPath2)
    batch.subBatches.append(batch2)

for jobItem in list(batch.items()):
    for x in [imp for imp in jobItem.importanceJobsRecursive()]:
        if batch.hasName(x.name.replace('_post', '')):
            print('replacing importance sampling run (not deleting files): ' + x.name)
            jobItem.removeImportance(x)

batch.save()
