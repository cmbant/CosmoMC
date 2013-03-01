import planckStyle as s

#import GetDistPlots
#g=GetDistPlots.GetDistPlotter('main/plot_data')

g=s.plotter

roots = ['base_nnu_meffsterile_planck_lowl_lowLike']

params = g.get_param_array(roots[0], ['meffsterile', 'nnu', 'omegach2'])
g.settings.setWithSubplotSize(4)
g.settings.lw_contour = 0.2
g.setAxes(params, lims=[0, 5, 3.046, 5.5])
g.add_3d_scatter(roots[0], params)

g.export('msterile-nnu')


