import planckStyle as s

g=s.getSinglePlotter()

g.settings.param_names_for_labels = 'clik_units_for_w_wa.paramnames'

roots = ['base_w_wa_planck_lowl_lowLike_BAO','base_w_wa_planck_lowl_lowLike_SNLS','base_w_wa_planck_lowl_lowLike_highL_BAO_post_HST']
params = g.get_param_array(roots[0], ['w','wa','H0'])
#g.settings.setWithSubplotSize(4)

g.setAxes(params,lims=[-2, -0.3 , -1.6, 2])
g.settings.lw_contour = 0.8
g.add_2d_contours(roots[0], params[0], params[1], filled=False, color='#000000')
#g.add_2d_contours(roots[1], params[0], params[1], filled=False,color='#0000ff' )
g.add_3d_scatter(roots[0],params)
#g.add_2d_contours(roots[2], params[0], params[1], filled=False,color='#000000' )

#g.add_line([-1., -1.], [-2., 2.], zorder=3)
#g.add_line([-2., -.3], [0.,0.], zorder=3)
g.add_x_marker(-1)
g.add_y_marker(0)

#savefig('w_wa_H0_test.eps',transparent=True)
g.export('w_wa_H0')

