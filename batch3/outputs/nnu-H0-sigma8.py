import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.settings.param_names_for_labels = 'clik_Hunits.paramnames'

base = 'base_nnu_'
roots = [base + s.defdata_all]

neff = [1, 5]

updated = True

if updated:
    g.add_y_bands(74.03, 1.42)
    g.add_text_left('Riess et al. (2019)', 0.03, 0.8, color='k', fontsize=7)
else:
    g.add_y_bands(73.45, 1.66)
    g.add_text_left('Riess et al. (2018)', 0.03, 0.76, color='k', fontsize=7)

g.plot_3d(roots, ['nnu', 'H0', 'sigma8'])
norm = 3.046
g.add_x_marker(norm, ls='-')
# g.add_x_marker(norm + 0.39)
# g.add_x_marker(norm + 0.57)
g  # .add_x_marker(norm + 1)

g.add_2d_contours(g.getRoot('nnu', s.defdata_all_lensing + '_BAO'), 'nnu', 'H0', color='black')
# g.add_2d_contours(g.getRoot('nnu', s.defdata_all + '_abundances'), 'nnu', 'H0', color='red', ls='--')
# g.add_2d_contours(g.getRoot('nnu', s.defdata + '_H073p9'), 'nnu', 'H0', color='magenta', ls='-')

g.add_2d_contours(g.getRoot('nnu', s.defdata_all + '_Riess18_post_BAO_lensing'), 'nnu', 'H0', color='darkblue', ls='--',
                  alpha=0.7)

# g.add_2d_contours(g.getRoot('nnu', s.defdata_all + '_BAO_Cooke17_Aver15'), 'nnu', 'H0', color='green', filled=False,
#                  alpha=1)

# g.add_legend([s.planckall + '+BAO'], legend_loc='upper left', colored_text=True)

gca().set_xticks([2, 2.5, 3.0, 3.5, 4])
ylim([57, 78])
xlim([2, 4.2])

g.export()
