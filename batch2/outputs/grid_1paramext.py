import planckStyle as s
from pylab import *
import sys

g = s.getSubplotPlotter()
g.settings.legend_frac_subplot_margin = 0.1
# g.settings.solid_colors = [g.settings.solid_colors[0]]+g.settings.solid_colors

# g.settings.solid_colors = ['#009966', '#000866', '#000866']

bases = ['omegabh2', 'omegach2', 'ns', 'H0', 'sigma8']
params = ['omegak', 'mnu', 'nnu', 'yhe', 'nrun', 'r02']
paramchains = ['omegak', 'mnu', 'nnu', 'yhe', 'nrun', 'r']
defs = [0, 0.06, 3.046, 0.2449, 0, 0, -1]

# get means for LCDM
refJobItem = g.getJobItem('', s.defdata_all)
bfs = []
for base in bases:
    bfs.append(refJobItem.result_marge.parWithName(base).mean)

roots = [[g.getRoot(param, s.defdata_TT), g.getRoot(param, s.defdata_all), g.getRoot(param, s.defdata_all + '_BAO')] for param in paramchains]

labels = [s.datalabel[s.defdata_TT], s.datalabel[s.defdata_all], s.datalabel[s.defdata_all] + '+BAO']

g.rectangle_plot(bases, params, yroots=roots, ymarkers=defs, xmarkers=bfs, filled=True, legend_labels=labels)
g.export()

