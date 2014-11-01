import planckStyle as s
g = s.plotter


g.settings.lineM = ['-g', '-r', '-b', '-k', '--r', '--b']

pol = ['TT', 'TE', 'EE', 'TTTEEE']
dataroots = [ getattr(s, 'defdata_' + p) for p in pol]
labels = [s.datalabel[r] for r in dataroots]
labels += [r'\textit{Planck} TE$+$ lowEB', r'\textit{Planck} EE$+$ lowEB']

for par in ['', 'nnu', 'mnu', 'Alens', 'r', 'yhe', 'nrun']:
    g.newPlot()
    base = 'base_'
    if par: base += par + '_'
    roots = [base + dat for dat in dataroots]

    roots += [roots[1].replace('TEB', 'EB'), roots[2].replace('TEB', 'EB')]

    g.settings.legend_frac_subplot_margin = 0.15

    plotpars = [ 'zrei', 'H0', 'omegabh2', 'thetastar', 'A', 'tau', 'omegam', 'omegach2', 'ns', 'sigma8']
    if par: plotpars[0] = par
    g.plots_1d(roots, plotpars, nx=5, legend_ncol=len(roots), legend_labels=labels, share_y=True)

    g.export(tag=par)


