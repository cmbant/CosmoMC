import planckStyle as s
g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)

roots = [g.getRoot('', s.defdata_all), 'base_lensonly', g.getRoot('', s.defdata_all + '_lensing')]

params = ['s8omegamp25', 'rmsdeflect', 'H0']

g.plot_3d(roots, params)

# g.add_2d_contours('base_lensonly',param_pair=['s8omegamp25','rmsdeflect'])
g.export()
