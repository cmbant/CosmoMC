import planckStyle as s

g = s.getSubplotPlotter(subplot_size=2)

g.settings.solid_colors = [('#8CD3F5', '#006FED'), ('#F7BAA6', '#E03424'), ('#B1B1B1', '#515151'), 'g']

roots = []
roots.append('base_DESlens_lenspriors_BAO')
roots.append('base_lensing_lenspriors_BAO')
roots.append('base_DESlens_lenspriors_lensing_BAO')
roots.append('base_plikHM_TTTEEE_lowl_lowE')
params = [u'H0', u'omegam', u'sigma8']
g.triangle_plot(roots, params, filled=True,
                legend_labels=['DES lensing+BAO', r'$\textit{Planck}$ lensing+BAO', r'(DES+$\textit{Planck}$) lensing + BAO', s.planckall])
g.export()
