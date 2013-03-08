import planckStyle as s

g=s.plotter

roots = ['base_omegak_planck_lowl_lowLike_highL','base_omegak_planck_lowl_lowLike_highL_lensing','base_omegak_planck_lowl_lowLike_highL_BAO_post_BAO']
params = g.get_param_array(roots[0], ['omegam', 'omegal', 'H0'])
g.settings.setWithSubplotSize(4)
g.settings.lw_contour = 0.2
g.setAxes(params, lims=[0, 1, 0, 1])
g.add_3d_scatter(roots[0], params)
g.add_2d_contours(roots[1], params[0], params[1], filled=False, alpha =0.9)
g.add_2d_contours(roots[2], params[0], params[1], filled=True, color='#991100')
g.add_line([1, 0], [0, 1], zorder=1,ls='--',color='gray')
g.export('Omegam-Omegal-H0')

g.newPlot()
g.setAxes(params, lims=[0.2, 0.5, 0.5, 0.8])
g.add_3d_scatter(roots[0], params)
g.add_2d_contours(roots[1], params[0], params[1], filled=False, alpha =0.9)
g.add_2d_contours(roots[2], params[0], params[1], filled=True, color='#991100')
g.add_line([1, 0], [0, 1], zorder=1,ls='--',color='gray')
g.export('Omegam-Omegal-H0_zoom')
