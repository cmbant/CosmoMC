import os, sys, batchJob, fnmatch, shutil

if len(sys.argv) < 4:
    print 'Usage: python/extractFiles.py directory_with_outputs target_dir [pdf] [ext2] [ext3..]'
    sys.exit()


batch = batchJob.readobject()

target_dir = os.path.abspath(sys.argv[2]) + os.sep
if not os.path.exists(target_dir): os.makedirs(target_dir)

for ext in sys.argv[3:]:
    pattern = '*.' + ext
    for jobItem in batch.items(wantImportance=True):
        for f in os.listdir(jobItem.distPath):
            if fnmatch.fnmatch(f, jobItem.name + pattern):
                print jobItem.distPath + f
                shutil.copyfile(jobItem.distPath + f, target_dir + f)
