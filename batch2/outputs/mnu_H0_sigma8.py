import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.settings.param_names_for_labels = 'clik_Hunits.paramnames'


base = 'base_mnu_'
roots = [base + s.defdata]


mnu = [0, 2.3]
H, sigma = s.H0_gpe
g.add_y_bands(H, sigma, xlim=mnu)


g.plot_3d(roots, ['mnu', 'H0', 'sigma8'])

root = g.getRoot('mnu', s.defdata + '_lensing')
g.add_2d_contours(root, 'mnu', 'H0', plotno=0)

root = g.getRoot('mnu', s.defdata + '_BAO_lensing')
g.add_2d_contours(root, 'mnu', 'H0', filled=True, zorder=2, alpha=0.6)

# gca().set_xticks([2, 2.5,3.0,3.5,4])
ylim([51, 79])
gca().set_yticks([55, 60, 65, 70, 75])

g.export()
