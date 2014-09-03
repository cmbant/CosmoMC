
import os
import sys
import glob

import numpy as np

# ==============================================================================
 
plot_data_a = "/home/mtartar/freelance/cosmomc/outputs_f903"
plot_data_b = "/home/mtartar/freelance/cosmomc/python/outputs_python3"

# ==============================================================================
# .corr



# ==============================================================================
# .covmat 

skiprows_covmat = 1


# ==============================================================================


# .likestats 

skiprows_likestats = 6
usecols_likestats = (1, 2, 3, 4, 5) # for 2 contours  




# ==============================================================================
# .margestats

skiprows_margestats = 3 
usecols_margestats = (1, 2, 3, 4, 6, 7) 

a = None
files_a = glob.glob(os.path.join(plot_data_a, "*.margestats"))
file_a = ""
if files_a: 
    file_a = files_a[0]
    a = np.loadtxt(file_a, skiprows=skiprows_margestats, usecols=usecols_margestats)

b = None
files_b = glob.glob(os.path.join(plot_data_b, "*.margestats"))
file_b = ""
if files_b: 
    file_b = files_b[0]
    b = np.loadtxt(file_b, skiprows=skiprows_margestats, usecols=usecols_margestats)

if (a is not None) and (b is not None):
    diff = np.abs(a - b)
    mini = np.min(diff)
    maxi = np.max(diff)
    print "Files .margestats, min, max = ", mini, maxi
else:
    print "Can not compare .margestats"


# ==============================================================================



