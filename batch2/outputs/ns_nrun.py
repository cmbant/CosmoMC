import planckStyle as s
from pylab import *
g = s.getSinglePlotter()

dataroots = [s.defdata_TT, s.defdata_all]
roots = ['base_nrun_' + root for root in dataroots]

g.plot_3d(roots, ['ns02', 'nrun', 'ns'])

g.add_y_marker(0)
g.export()
