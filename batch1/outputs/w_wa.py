import planckStyle as s

g=s.getSinglePlotter()

g.settings.param_names_for_labels = 'clik_units_for_w_wa.paramnames'

roots = ['base_w_wa_planck_lowl_lowLike_BAO','base_w_wa_planck_lowl_lowLike_Union2','base_w_wa_planck_lowl_lowLike_SNLS']

params = g.get_param_array(roots[0], ['w','wa'])



g.plot_2d(roots, param_pair=params, filled=True,lims=[-2, -0.3 , -1.6, 2])

#g.add_2d_contours(roots[0], params[0], params[1], filled=False, color='#ff0000')
#g.add_2d_contours(roots[2], params[0], params[1], plotno=1,filled=True,color='#00ff00' )
#g.add_2d_contours(roots[1], params[0], params[1], plotno=2,filled=True,color='#0000ff' )

g.add_x_marker(-1)
g.add_y_marker(0)
g.export('w_wa')

