import planckStyle as s
from pylab import *
g = s.getSinglePlotter()

g.settings.lineM = ['-b', '-r', '-g', '--b', '--r', '--g', '-c', '-y']

labels = [s.planckTT, '+lensing', '+lensing+BAO' , s.planck + '\\ TT +WP', s.defplanck, '+lensing+BAO' ]
roots = [s.defdata_TTonly,
         s.defdata_TTonly + '_lensing',
         s.defdata_TTonly + '_lensing_BAO',
         'plikHM_TT_WMAPTEB',
         s.defdata,
         s.defdata + '_lensing_BAO',
]

roots = [g.getRoot('', root) for root in roots]

g.plot_1d(roots, 'tau', normalized=True)
g.add_legend(labels, legend_loc='upper right')


g.export()

