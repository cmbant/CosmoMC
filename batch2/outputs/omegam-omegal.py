import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

roots = [g.getRoot('omegak', s.defdata_all),
       g.getRoot('omegak', s.defdata_all + '_lensing'),
       g.getRoot('omegak', s.defdata_all + '_lensing_BAO')]

params = g.get_param_array(roots[0], ['omegam', 'omegal', 'H0'])

g.newPlot()
# g.setAxes(params, lims=[0.2, 0.5, 0.5, 0.8])
g.setAxes(params, lims=[0.2, 0.55, 0.45, 0.8])

g.add_3d_scatter(roots[0], params)
g.add_2d_contours(roots[1], params[0], params[1], filled=False, zorder=1)
g.add_2d_contours(roots[2], params[0], params[1], filled=True, zorder=2, alpha=0.85)

g.add_line([1, 0], [0, 1], zorder=0)

g.add_legend(['+lensing', '+lensing+BAO'])
g.export()
