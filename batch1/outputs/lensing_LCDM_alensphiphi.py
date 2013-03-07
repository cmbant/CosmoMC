import planckStyle as s

from pylab import *


g=s.plotter
g.settings.line_labels = False
g.settings.setWithSubplotSize(4)

roots=['base_Aphiphi_planck_lowl_lowLike_highL_lensing','base_Alens_planck_lowl_lowLike_highL_post_lensing','base_Alens_planck_lowl_lowLike_highL']
g.setAxes(lims=[0.7,1.8 , 0, 1.1])

g.add_1d(roots[0],'Aphiphi',color='b');
g.add_1d(roots[1],'Alens',color='r');
g.add_1d(roots[2],'Alens',color='c',ls='--');

Aphiphi = g.param_latex_label(roots[0], 'Aphiphi')
Alens = g.param_latex_label(roots[1], 'Alens')

#g.plots_1d(roots,['Alens'])
g.add_x_marker(1,color='k',ls='--')

text(1.35,1, Aphiphi+' ('+s.WPhighLlensing+')' ,color='b')
text(1.35,0.9, Alens +' ('+s.WPhighLlensing+')' ,color='r')
text(1.35,0.8, Alens +' ('+s.WPhighL+')' ,color='c')
#text(1.3,1,s.WPhighLlensing,color='r')
#gca().set_xticks([0.9,1,1.1,1.2, 1.3, 1.4, 1.5, 1.6])
xlabel('Amplitude relative to physical',fontsize=14)

g.export('lensing_LCDM_alensphiphi')


