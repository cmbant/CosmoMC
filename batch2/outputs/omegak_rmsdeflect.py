import planckStyle as s
from pylab import *
g = s.getSinglePlotter()

roots = []
roots.append('base_omegak_plikHM_TT_lowTEB')
roots.append('base_omegak_plikHM_TT_lowTEB_lensing')
roots.append('base_omegak_plikHM_TT_lowTEB_BAO')
roots.append('base_omegak_plikHM_TT_lowTEB_BAO_post_lensing')
g.plot_3d(roots, ['omegak', 'rmsdeflect', 'H0'], filled=[False, False, True], colors=['k', 'r', g.settings.solid_colors[0]])
g.add_legend([s.defplanck + '+lensing', s.defplanck + '+BAO', s.defplanck + '+lensing+BAO' ], legend_loc='lower left', colored_text=True)
g.add_x_marker(0)
xlim(-0.18, 0.025)
ylim(2.2, 2.99)

g.export()
