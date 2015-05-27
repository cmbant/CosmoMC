import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.settings.colorbar_rotation = -90
dataroots = [s.defdata_TT, s.defdata_TT + '_lensing_BAO', s.defdata_TT + '_BKP_lensing_BAO']

roots = [g.getRoot('nrun_r', s.defdata_TT),
         g.getRoot('nrun_r', s.defdata_TT + '_lensing_BAO'),
         g.getRoot('r', s.defdata_TT + '_lensing_BAO'),
         g.getRoot('nrun_r', s.defdata_TT + '_BKP_lensing_BAO'),
         g.getRoot('r', s.defdata_TT + '_BKP_lensing_BAO'),
         ]

nrun = g.param_latex_label(roots[0], 'nrun')
labels = [s.LCDM + '+running+tensors', s.LCDM + '+tensors']

# roots = ['base_nrun_r_' + s.defdata_all, 'base_r_' + s.defdata_all]

g.plot_3d(roots,
          params_for_plots=[['ns02', 'r02', 'nrun'], ['ns02', 'r02'], ['ns', 'r02'], ['ns02', 'r02'], ['ns', 'r02']],
          colors=['k', 'b', 'k', 'b'], ls=['-', '-', '--', '--'])

g.add_legend(labels, legend_loc='upper left', colored_text=True)

gca().set_xticks([0.9, 0.95, 1, 1.05, 1.1])

g.export()
