import planckStyle as s
from pylab import *
g = s.getSinglePlotter()

labels = [s.planckTT, '+lensing', s.planckall, '+lensing', '+BAO+H0+JLA' ]
roots = [s.defdata_TT, s.defdata_TT + '_lensing', s.defdata_all, s.defdata_all + '_lensing', s.defdata_all + '_BAO_H070p6_JLA_lensing' ]
roots = [g.getRoot('mnu', root) for root in roots]

g.plot_1d(roots, 'mnu')
g.add_legend(labels, legend_loc='upper right')

g.export()

