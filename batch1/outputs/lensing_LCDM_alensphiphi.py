import planckStyle as s

from pylab import *


g=s.plotter
g.settings.line_labels = False
g.settings.setWithSubplotSize(4)

roots=['base_Aphiphi_planck_lowl_lowLike_highL_lensing','base_Alens_planck_lowl_lowLike_highL_post_lensing','base_Alens_planck_lowl_lowLike_highL']
#g.setAxes(params,lims=[, -0.3 , -1.6, 2])

g.add_1d(roots[0],'Aphiphi',color='b');
g.add_1d(roots[1],'Alens',color='r');
g.add_1d(roots[2],'Alens',color='c',ls='--');

#g.plots_1d(roots,['Alens'])
g.add_1d_marker(1,color='k',ls='--')
text(1.3,1,s.WPhighLlensing,color='r')
text(1.3,0.9, s.WPhighL,color='b')
#gca().set_xticks([0.9,1,1.1,1.2, 1.3, 1.4, 1.5, 1.6])
gca().set_yticks([])

g.export('lensing_LCDM_alensphiphi')


