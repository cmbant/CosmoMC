import planckStyle as s
from pylab import *
g = s.getSinglePlotter()

labels = [s.planckTT, s.planckall]
roots = [s.defdata_TT, s.defdata_all]
roots = ['base_nrun_' + root for root in roots]

g.plot_3d(roots, ['ns02', 'nrun', 'ns'])
# g.add_legend(labels, legend_loc='upper right')
g.add_y_marker(0)
g.export()
