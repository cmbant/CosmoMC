import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)
base = 'base_nnu_meffsterile_'
roots = [ g.getRoot('nnu_meffsterile', s.defdata),
          g.getRoot('nnu_meffsterile', s.defdata + '_BA0'),
          g.getRoot('nnu_mmu', s.defdata + '_BAO'),
          g.getRoot('', s.defdata + '_BAO')
                  ]

omm = np.arange(0.2, 0.45, 0.01)


# plotBounds(galaxygalaxy,'green')
s.plotBounds(omm, s.CFTHlens_mnu, 'burlywood')
s.plotBounds(omm, s.PLSZ, 'gray')
# s.plotBounds(s8, s.planck_lensing, 'cyan')


g.plot_3d(roots, ['omegam', 'sigma8', 'H0'])

# g.add_2d_contours(base + s.defdata_all + '_post_BAO', 'sigma8', 'omegam')
# g.add_2d_contours(base+s.defdata_all+'_post_HST70p6', 'sigma8', 'omegam', plotno=1)


mnu = g.param_latex_label(roots[2], 'mnu', labelParams='params_CMB.paramnames')
nnu = g.param_latex_label(roots[2], 'nnu')
meff = g.param_latex_label(roots[1], 'meffsterile', labelParams='clik_latex.paramnames')

legends = [meff + '+' + nnu, mnu + '+' + nnu, s.LCDM]

g.add_legend(legends, legend_loc='upper right', colored_text=True);


g.export()
