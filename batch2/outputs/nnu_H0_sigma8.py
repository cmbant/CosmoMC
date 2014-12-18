import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.settings.param_names_for_labels = 'clik_Hunits.paramnames'

base = 'base_nnu_'
roots = [base + s.defdata]

neff = [1, 5]
H = 70.6
sigma = 3.3
g.add_y_bands(H, sigma, xlim=neff)


g.plot_3d(roots, ['nnu', 'H0', 'sigma8'])
norm = 3.046
g.add_x_marker(norm, ls='-')
g.add_x_marker(norm + 0.39)
g.add_x_marker(norm + 0.57)
g.add_x_marker(norm + 1)

g.add_2d_contours(g.getRoot('nnu', s.defdata_all + '_BAO'), 'nnu', 'H0', color='black')
# g.add_2d_contours(g.getRoot('nnu', s.defdata_all + '_abundances'), 'nnu', 'H0', color='red', ls='--')
# g.add_2d_contours(g.getRoot('nnu', s.defdata + '_H073p9'), 'nnu', 'H0', color='magenta', ls='-')

# g.add_legend([s.planckall + '+BAO'], legend_loc='upper left', colored_text=True)

gca().set_xticks([2, 2.5, 3.0, 3.5, 4, 4.5])
ylim([57, 82])

g.export()
