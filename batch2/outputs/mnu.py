import planckStyle as s
from pylab import *

g=s.plotter
g.settings.line_labels = False
g.settings.setWithSubplotSize(4)

g.settings.tight_layout=True
g.settings.lab_fontsize=32
g.settings.axes_fontsize = 24
g.settings.lw1 = 2

g.setAxes(lims=[0.0, 2.0, 0.0, 1.1])

roots=['base_mnu_planck_lowl_lowLike_highL','base_mnu_planck_lowl_lowLike_highL_lensing', 'base_mnu_Alens_planck_lowl_lowLike_highL','base_mnu_Alens_planck_lowl_lowLike_highL_post_lensing']
g.add_1d(roots[0],'mnu');
g.add_1d(roots[1],'mnu',color='b');
g.add_1d(roots[2],'mnu',color='r');
g.add_1d(roots[3],'mnu',color='g');

text(0.9,1.02, s.WPhighL, fontsize=18)
text(0.9,0.94, s.WPhighLlensing, color='b', fontsize=18)
text(0.9,0.85, s.WPhighL+'($A_{\mathrm L}$)', color='r', fontsize=18)
text(0.9,0.77, s.WPhighLlensing+'($A_{\mathrm L}$)', color='g', fontsize=18)

xlabel(r'$\Sigma m_\nu\,[\mathrm{eV}]$',fontsize=32)
ylabel(r'$P/P_{\rm max}$',fontsize=32)

g.export('mnu')
