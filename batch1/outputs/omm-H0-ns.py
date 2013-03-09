import planckStyle as s
from pylab import *
g=s.getSinglePlotter()

roots = ['base_planck_lowl','base_planck_lowl_lowLike']
params = g.get_param_array(roots[0], ['omegam', 'H0', 'ns'])
g.add_3d_scatter(roots[0], params)
g.add_2d_contours(roots[1], params[0], params[1], filled=False)
g.setAxes(params)
yticks([64, 66, 68, 70,72])
xticks([0.26, 0.3, 0.34, 0.38])

g.export('Omegam-H0-ns')


