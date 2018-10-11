import planckStyle as s
from matplotlib.pyplot import *

g = s.getSinglePlotter()

roots = [s.defdata,
         s.defdata_TE,
         s.defdata_EE,
         s.defdata_all,
         s.defdata_all_lensing
         ]

roots = [g.getRoot('Alens', root) for root in roots]

g.plot_1d(roots, 'Alens', normalized=True, colors=['C0', 'C1', 'C2', 'k', 'C3'], ls=['-'] * 4 + ['--'])

roots =[x.replace('plikHM','CamSpecHM') for x in roots[:-1]]
g.plot_1d(roots, 'Alens', normalized=True, colors=['C0', 'C1', 'C2', 'k', 'C3'], ls=[':']*4)


g.add_legend([s.shortTT, s.shortTE, s.shortEE,
              s.shortall, s.shortall + '+lensing'],
             legend_loc='upper right', colored_text=True, fontsize=8)

g.add_x_marker(1, ls='-')
gca().set_xticks([0.2, 0.6, 1, 1.4, 1.8, 2.2, 2.6])
# xlim([0.3, 2.5])
g.export()
