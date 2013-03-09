import planckStyle as s

g=s.getSinglePlotter()

roots = ['base_planck_lowl','base_planck_lowl_lowLike']
params = g.get_param_array(roots[0], ['omegam', 'H0', 'ns'])
g.setAxes(params)
g.add_3d_scatter(roots[0], params)
g.add_2d_contours(roots[1], params[0], params[1], filled=False)
g.export('Omegam-H0-ns')


