import planckStyle as s
import pylab as plt

g = s.getSinglePlotter()
g.settings.param_names_for_labels = 'clik_Hunits.paramnames'

updated = True
legacy = False

if updated:
    g.add_y_bands(74.03, 1.42)
    g.add_text_left('Riess et al. (2019)', 0.03, 0.8, color='k', fontsize=7)
else:
    g.add_y_bands(73.45, 1.66)
    g.add_text_left('Riess et al. (2018)', 0.03, 0.76, color='k', fontsize=7)


roots = []
roots.append('base_BAO_Cooke17_Pantheon18')
roots.append('base_lensing_lenspriors_BAO_post_Pantheon18')
roots.append('base_BAO_Cooke17_Pantheon18_theta')
roots.append('base_plikHM_TTTEEE_lowl_lowE')
g.plot_2d(roots, [u'omegam', u'H0'], filled=True, shaded=False)

if not legacy:
    g.add_2d_contours('base_BAO_Cooke17_JLA',u'omegam', u'H0', color='darkgreen', ls='--', alpha=0.3)

#g.add_2d_contours('base_mnu_lensing_lenspriors_BAO_post_Pantheon','omegam','H0', ls='--', filled=False)
g.add_legend(['BAO+Pantheon+D/H BBN', 'BAO+Pantheon+D/H BBN+lensing', r'BAO+Pantheon+D/H BBN+$\theta_{\rm MC}$', s.planckall],
             colored_text=True, fontsize=6.5, legend_loc='lower left')
plt.ylim(58, 78)
plt.xlim(0.23, 0.4)
g.export()
