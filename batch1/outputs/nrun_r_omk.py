import planckStyle as s
from pylab import *


g=s.getSinglePlotter()
g.settings.colorbar_rotation =-90
g.settings.lab_fontsize=10

roots = ['base_nrun_r_omegak_planck_lowl_lowLike_highL','base_nrun_r_planck_lowl_lowLike_highL','base_nrun_r_omegak_planck_lowl_lowLike_highL_BAO' ]
params = g.get_param_array(roots[0], ['nrun','r02','omegak','ns'])

g.setAxes(params)
g.settings.lw_contour = 0.2
g.add_2d_contours(roots[0], params[0], params[1], filled=False, color='#ff0000')
#g.add_2d_contours(roots[1], params[0], params[1], filled=False, color='#00ff00')
g.add_2d_contours(roots[2], params[0], params[1], filled=False, color='#0000ff')
#g.add_2d_contours(roots[2], params[0], params[1], filled=False,color='#0000ff' )
g.add_3d_scatter(roots[0],params)
#g.add_2d_contours(roots[2], params[0], params[1], filled=False,color='#000000' )
g.add_x_marker(0)
ylim([0, 1.4])
ylabel(r'Tensor-to-Scalar Ratio ($r$)')
xlabel(r'Running spectral index ($dn_{\rm s}/d\ln k$)')
g.rotate_yticklabels()
for ticklabel in gca().yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")
#g.add_line([-1., 1.], [0., 0.], zorder=3)
#g.add_line([0., 0.], [-1.,2.], zorder=3)
#savefig('w_wa_H0_test.eps',transparent=True)
g.export('nrun_r_omk')

