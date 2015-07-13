import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.settings.lineM = ['-b', '-r', '-g', '--r', ':r', '--g', '-c', '-y']

labels = [s.planckTT, '+lensing', '+lensing+BAO', s.defplanck, s.defplanck + '+WP', s.defplanck + '+BAO']
roots = [s.defdata_TTonly,
         s.defdata_TTonly + '_lensing',
         s.defdata_TTonly + '_lensing_BAO',
         s.defdata,
         'plikHM_TT_WMAPTEB',
         s.defdata + '_BAO',
         ]

roots = [g.getRoot('', root) for root in roots]

g.plot_1d(roots, 'tau', normalized=True)
g.add_legend(labels, legend_loc='upper right')
# xlim([0.01, 0.3])

g.export()

