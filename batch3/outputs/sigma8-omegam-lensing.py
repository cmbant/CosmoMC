import planckStyle as s
import pylab as plt

g = s.getSinglePlotter(ratio=1)
roots = []
roots.append('base_plikHM_TTTEEE_lowl_lowE')
roots.append('base_plikHM_TTTEEE_lowl_lowE_lensing')
g.plot_2d(roots, [u'omegam', u'sigma8'], filled=True, shaded=False)
g.add_2d_contours('base_plikHM_TTTEEE_lowl_lowE_lensing_post_zre6p5', u'omegam', u'sigma8', ls='--', color='k')
# roots.append('base_plikHM_TTTEEE_lowl_lowE_lensing_post_zre6p5')
g.add_legend([s.planckall, s.planckall + '+lensing', s.planckall + r'+lensing, $z_{\rm re} > 6.5$'], colored_text=True)

plt.gca().set_yticks([0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85])
g.export()
