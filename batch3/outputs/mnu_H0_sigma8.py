import planckStyle as s
import pylab as plt

g = s.getSinglePlotter()

g.settings.param_names_for_labels = 'clik_Hunits.paramnames'

base = 'base_mnu_'
roots = [base + s.defdata_all]

H = 73.45
sigma = 1.66
# g.add_y_bands(H, sigma)
g.add_x_marker(0.0982, ls=':')

plt.gca().axvspan(0, 0.056, color='lightgray', zorder=-10, alpha=0.3)

g.plot_3d(roots, ['mnu', 'H0', 'sigma8'])

# root = g.getRoot('mnu', s.defdata_all)
# g.plot_2d(root, 'mnu', 'H0', plotno=0)

root = g.getRoot('mnu', s.defdata_all_lensing)
g.add_2d_contours(root, 'mnu', 'H0', plotno=0)

root = g.getRoot('mnu', s.defdata_all + '_BAO_lensing')
g.add_2d_contours(root, 'mnu', 'H0', filled=False, zorder=2, ls='--', color='darkblue')

g.add_2d_contours('base_nnu_mnu_' + s.defdata_all_lensing + '_BAO', 'mnu', 'H0', ls='--', color='green', alpha=0.5)

# g.add_text('Riess 2018', 0.97, 0.85)

plt.gca().text(0.06, 60.7, 'NH', horizontalalignment='left', fontsize=7, color='darkgray')
plt.gca().text(0.105, 60.7, 'NH or IH', horizontalalignment='left', fontsize=7, color='darkgray')

# gca().set_xticks([2, 2.5,3.0,3.5,4])
plt.ylim([60, 73])
plt.xlim([0, 0.44])

plt.gca().set_yticks([60, 62, 64, 66, 68, 70])

g.add_legend([s.shortall + '+lensing', '+BAO', r'$+N_{\rm eff}$'], legend_loc='upper right', fontsize=7,
             align_right=True)

g.export()
