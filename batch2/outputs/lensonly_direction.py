import planckStyle as s
g = s.getSinglePlotter()

roots = [g.getRoot('', s.defdata), 'base_lensonly', g.getRoot('', s.defdata + '_lensing')]

params = ['s8omegamp25', 'rmsdeflect', 'H0']

g.plot_3d(roots, params)

# g.add_2d_contours('base_lensonly',param_pair=['s8omegamp25','rmsdeflect'])
g.export()
