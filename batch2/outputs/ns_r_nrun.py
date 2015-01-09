import planckStyle as s
from pylab import *
g = s.getSinglePlotter()

g.settings.colorbar_rotation = -90
roots = [s.defdata_TT, s.defdata_all]
roots = ['base_nrun_r_' + root for root in roots] + ['base_r_' + s.defdata_all]

nrun = g.param_latex_label(roots[0], 'nrun')
labels = [s.LCDM + '+running+tensors', s.LCDM + '+tensors' ]

# roots = ['base_nrun_r_' + s.defdata_all, 'base_r_' + s.defdata_all]

g.plot_3d(roots, params_for_plots=[['ns02', 'r02', 'nrun'], ['ns02', 'r02'], ['ns', 'r02']], colors=['k', 'b'])

g.add_legend(labels, legend_loc='upper left', colored_text=True)
g.export()
