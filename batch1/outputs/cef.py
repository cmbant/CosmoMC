import planckStyle as s
from pylab import *

g=s.plotter
g.settings.line_labels = False
g.settings.setWithSubplotSize(4)

g.settings.tight_layout=True
g.settings.lab_fontsize=32
g.settings.axes_fontsize = 24

g.setAxes(lims=[0.94,1.06, 0.0, 1.1])

roots = ['test_H0_alpha_cosmoplanckmass_WMAP9', 'test_H0_alpha_cosmoplanckmass_lowllikeCAMspec_new', 'test_H0_alpha_cosmoplanckmass_lowllikeCAMspec_new_BAO']
g.add_1d(roots[0],'cef',color='r');
g.add_1d(roots[1],'cef',color='b');
g.add_1d(roots[2],'cef',color='g');

text(0.942,1.02, r'\textit{WMAP}9', fontsize=24)
text(0.942,0.92, s.WP, color='b', fontsize=24)
text(0.942,0.82, s.WP+'+BAO', color='r', fontsize=24)

xlabel(r'$\alpha/\alpha_0$',fontsize=32)
ylabel(r'$P/P_{\rm max}$',fontsize=32)

g.export('cef')
