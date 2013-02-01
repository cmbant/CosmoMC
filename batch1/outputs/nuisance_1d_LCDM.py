import GetDistPlots
g=GetDistPlots.GetDistPlotter('main/plot_data')
g.settings.plot_meanlikes=False

labels=['Planck+WP+highL','Planck+WP']
roots=['base_planck_CAMspec_lowl_lowLike_highL','base_planck_CAMspec_lowl_lowLike']
g.plots_1d(roots,paramList='batch1/outputs/nuisance.paramnames',nx=5, legend_labels=labels)
g.export('nuisance_1d_LCDM.eps')


g.newPlot()
g.plots_1d(roots,paramList='batch1/outputs/foregrounds.paramnames',nx=3, legend_labels=labels)
g.export('foregrounds_1d_LCDM.eps')

g.newPlot()
g.plots_1d(roots,paramList='batch1/outputs/camspec.paramnames',nx=5, legend_labels=labels)
g.export('camspec_1d_LCDM.eps')


g.newPlot()
g.plots_1d(roots,paramList='batch1/outputs/camspec_foregrounds.paramnames',nx=3, legend_labels=labels)
g.export('camspec_foregrounds_1d_LCDM.eps')
