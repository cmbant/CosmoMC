import GetDistPlots
from pylab import *
g=GetDistPlots.GetDistPlotter('main/plot_data')
g.settings.line_labels = False
g.settings.setWithSubplotSize(6)

labels=['Planck TT','Planck+WP','Planck+WP+highL','Planck+Lensing']
roots=['base_Alens_planck_CAMspec_lowl_lowLike_highL','base_Alens_planck_CAMspec_lowl_lowLike_highL_post_lensing']

g.plots_1d(roots,['Alens'],legend_labels=labels)
g.add_1d_marker(1,color='k',ls='--')
text(1.3,1,'Planck+WP+highL+lensing',color='r')
text(1.3,0.9, 'Planck+WP+highL',color='k')
gca().set_xticks([0.9,1,1.1,1.2, 1.3, 1.4, 1.5, 1.6])

g.export('outputs/lensing_LCDM_alensphiphi.eps')


