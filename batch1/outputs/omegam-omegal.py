import planckStyle as s
from pylab import *

g=s.getSinglePlotter()

roots = ['base_omegak_planck_lowl_lowLike_highL','base_omegak_planck_lowl_lowLike_highL_lensing','base_omegak_planck_lowl_lowLike_highL_lensing_post_BAO']

params = g.get_param_array(roots[0], ['omegam', 'omegal', 'H0'])

g.setAxes(params, lims=[0, 1, 0, 1])
g.add_3d_scatter(roots[0], params)
g.add_2d_contours(roots[1], params[0], params[1], filled=False)
#g.add_2d_contours(roots[2], params[0], params[1], filled=True)
g.add_line([1, 0], [0, 1], zorder=1)
g.export('Omegam-Omegal-H0')

g.newPlot()
g.setAxes(params, lims=[0.2, 0.5, 0.5, 0.8])

g.add_3d_scatter(roots[0], params)
g.add_2d_contours(roots[1], params[0], params[1], filled=False, zorder=1)
g.add_2d_contours(roots[2], params[0], params[1], filled=True, zorder=2, alpha=0.85)

g.add_line([1, 0], [0, 1], zorder=0)

g.add_legend(['+lensing','+lensing+BAO'])
g.export('Omegam-Omegal-H0_zoom')
