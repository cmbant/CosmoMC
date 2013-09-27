import planckStyle as s
from pylab import *

g=s.plotter

labels=[s.WPhighL, s.WP]
roots=['base_planck_lowl_lowLike_highL','base_planck_lowl_lowLike']
#g.plots_1d(roots,paramList='batch1/outputs/nuisance.paramnames',nx=5, legend_labels=labels)
#g.export('nuisance_1d_LCDM')


#g.newPlot()
#g.plots_1d(roots,paramList='batch1/outputs/foregrounds.paramnames',nx=3, legend_labels=labels)
#g.exportExtra('foregrounds_1d_LCDM')

#g.newPlot()
#g.plots_1d(roots,paramList='batch1/outputs/camspec.paramnames',nx=5, legend_labels=labels)
#g.exportExtra('camspec_1d_LCDM')

g.newPlot()
g.settings.line_labels = False
g.settings.figure_legend_frame =False
g.settings.axes_fontsize+=3
g.settings.lab_fontsize+=2

g.plots_1d(roots,paramList='batch1/outputs/camspec_foregrounds.paramnames',nx=3)

g.settings.lab_fontsize-=4
subplot(3,3,9)
gca().axis('off')
g.add_legend(labels,figure=False,legend_loc='right')

g.export('camspec_foregrounds_1d_LCDM')
