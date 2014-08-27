import planckStyle as s
from pylab import *
g = s.getSinglePlotter()

labels = [s.planckTT, '+lensing', s.planckall, '+lensing', '+BAO+HST+JLA' ]
roots = [s.defdata_TT, s.defdata_TT + '_lensing', s.defdata_all, s.defdata_all + '_lensing', s.defdata_all + '_BAO_HST70p6_JLA_post_lensing' ]
roots = ['base_mnu_' + root for root in roots]

g.plot_1d(roots, 'mnu')
g.add_legend(labels, legend_loc='upper right')

g.export()

