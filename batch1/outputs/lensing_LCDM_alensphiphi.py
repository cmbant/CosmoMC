import planckStyle as s

from pylab import *


g=s.plotter
g.settings.line_labels = False
g.settings.setWithSubplotSize(6)

roots=['base_Aphiphi_planck_lowl_lowLike_highL','base_Alens_planck_lowl_lowLike_highL_post_lensing']

g.plots_1d(roots,['Alens'])
g.add_1d_marker(1,color='k',ls='--')
text(1.3,1,s.WPhighLlensing,color='r')
text(1.3,0.9, s.WPhighL,color='k')
gca().set_xticks([0.9,1,1.1,1.2, 1.3, 1.4, 1.5, 1.6])

g.export('lensing_LCDM_alensphiphi')


