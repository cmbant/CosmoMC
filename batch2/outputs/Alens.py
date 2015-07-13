import planckStyle as s
from matplotlib.pyplot import *

g = s.getSinglePlotter()
g.settings.lineM = ['-b', '-.b', '--r', '--g', '--k', '-m', '-y', '--y']

roots = [s.defdata,
         s.defdata_TE,
         s.defdata_EE,
         s.defdata_all,
         ]

roots = [g.getRoot('Alens', root) for root in roots]
roots = roots[:1] + [g.getRoot('Aphiphi', s.defdata + '_lensing')] + roots[1:]

g.plot_1d(roots, 'Alens', normalized=True, param_renames={'Alensf': 'Alens', 'Aphiphi': 'Alens'},
          colors=['b', 'navy', 'r', 'g', 'k'])

Aphiphi = g.param_latex_label(roots[1], 'Aphiphi', labelParams='clik_latex.paramnames')

g.add_legend([s.deflabel, '+lensing (' + Aphiphi + ')', s.datalabel[s.defdata_TE], s.datalabel[s.defdata_EE],
              s.datalabel[s.defdata_all]],
             legend_loc='upper right', colored_text=True)

ext = []
for lab in ['TT', 'TE', 'EE', 'TTTEEE']:
    ext.append('base_Alens_CamSpecHM_' + lab + '_lowEB')
g.plot_1d(ext, 'Alens', normalized=True, colors=['b', 'r', 'g', 'k'], ls=[':', ':', ':', ':'])

g.add_x_marker(1, ls='-')
gca().set_xticks([0.6, 1, 1.4, 1.8, 2.2])
xlim([0.3, 2.5])
g.export()

g.newPlot()

g.settings.lineM = ['-b', '--r', '--g', '-k', '-c', '-m', ':c', '-.c']

roots = ['base_' + s.defdata_TTonly,
         'base_' + s.defdata_TE,
         'base_' + s.defdata_EE,
         'base_lensonly',
         'base_' + s.defdata,
         'base_' + s.defdata + '_lensing',
         # 'base_' + s.defdata_all + '_lensing',
         # 'base_lensonly_BAO_theta',
         'base_Alensf_' + s.defdata,
         'base_Alens_' + s.defdata,
         ]

g.plot_1d(roots, 'rmsdeflect', normalized=True, param_renames={'Alensf': 'Alens'})
g.add_legend(['TT', 'TE+lowP', 'EE+lowP', 'lensing', 'TT+lowP', 'TT+lowP+lensing'],
             legend_loc='upper right', colored_text=True)
xlim([1.9, 3.1])
g.export(tag='deflect')
