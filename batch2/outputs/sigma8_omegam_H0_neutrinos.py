import planckStyle as s
from pylab import *

g = s.getSubplotPlotter()

g.settings.legend_frac_subplot_margin = 0.1

for TT in [False, True]:
    g.newPlot()
    if not TT:
        basedat = s.defdata_all
        fname = 'Planckall'
    else:
        basedat = s.defdata_TT
        fname = 'PlanckTT'
    basedatname = s.datalabel[basedat]

    ref = 'base_' + basedat + '_lensing_post_BAO'
    vars = ['nnu', 'mnu', 'nnu_mnu', 'nnu_meffsterile']
    roots = [[g.getRoot(var, basedat), g.getRoot(var, basedat + '_lensing'), g.getRoot(var, basedat + '_lensing_BAO')] for var in vars]
    labels = [basedatname, '+lensing', '+lensing+BAO', '$\Lambda$CDM']

    mnu = r'$\Sigma m_\nu$'
    nnu = g.param_latex_label(roots[2][0], 'nnu', labelParams='clik_latex.paramnames')
    meff = g.param_latex_label(roots[3][0], 'meffsterile', labelParams='clik_latex.paramnames')

    legends = [s.LCDM , '+' + nnu, '+' + mnu, '+' + nnu + '+' + mnu, '+' + nnu + '+' + meff]

    plotroots = [ [root + [ref] for root in roots[0:2]], [root + [ref] for root in roots[2:]]] * 2
    plottexts = [ [legend for legend in legends[1:3]], [legend for legend in legends[3:]] ] * 2

    lims = {'sigma8':[0.61, 0.97], 'H0':[52, 79], 'omegam':[0.21, 0.55]}
    g.rectangle_plot(['omegam', 'omegam', 'H0', 'H0'], ['sigma8'] * 2, plot_roots=plotroots, plot_texts=plottexts,
                      param_limits=lims, legend_labels=labels, legend_ncol=4)

    g.export(tag=fname)
    g.newPlot()
