import GetDistPlots
g=GetDistPlots.GetDistPlotter('main/plot_data')

roots = ['base_planck_CAMspec_lowl','base_planck_CAMspec_lowl_lowLike_post_lensing']
params = g.get_param_array(roots[0], ['omegam', 'H0', 'ns'])
g.settings.setWithSubplotSize(4)
g.settings.lw_contour = 0.2
g.setAxes(params)
g.add_3d_scatter(roots[0], params)
g.add_2d_contours(roots[1], params[0], params[1], filled=False)
g.export('outputs/Omegam-H0-ns.eps')


