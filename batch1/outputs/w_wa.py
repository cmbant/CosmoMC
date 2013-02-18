import planckStyle as s

g=s.plotter

roots = ['base_w_wa_planck_CAMspec_lowl_lowLike_BAO','base_w_wa_planck_CAMspec_lowl_lowLike_SNLS','base_w_wa_planck_CAMspec_lowl_lowLike_Union2']
params = g.get_param_array(roots[0], ['w','wa'])
g.setAxes(params,lims=[-2, -0.3 , -1.6, 2])
g.settings.setWithSubplotSize(4)
g.settings.lw_contour = 0.2
g.add_2d_contours(roots[0], params[0], params[1], filled=False, color='#ff0000')
g.add_2d_contours(roots[2], params[0], params[1], plotno=1,filled=True,color='#00ff00' )
g.add_2d_contours(roots[1], params[0], params[1], plotno=2,filled=True,color='#0000ff' )
g.add_line([-1., -1.], [-2., 2.], zorder=3)
g.add_line([-2., -.3], [0.,0.], zorder=3)
g.export('w_wa_test')

