import planckStyle as s

#import GetDistPlots
#g=GetDistPlots.GetDistPlotter('main/plot_data')

g=s.plotter

roots = ['base_nnu_meffsterile_planck_lowl_lowLike_highL','base_nnu_meffsterile_planck_lowl_lowLike_highL_post_BAO']

params = g.get_param_array(roots[0], ['meffsterile', 'nnu', 'omegach2'])
g.settings.setWithSubplotSize(4)
g.settings.lw_contour = 0.2
g.setAxes(params, lims=[0, 3, 3.046, 4.9])
g.add_3d_scatter(roots[0], params)
#g.add_2d_contours(roots[1],params[0],params[1])

g.export('msterile-nnu')


