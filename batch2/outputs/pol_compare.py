import planckStyle as s
g = s.getSubplotPlotter()

g.settings.legend_fontsize -= 3.5

g.settings.lineM = ['-g', '-r', '-b', '-k', '--r', '--b']

pol = ['TT', 'TE', 'EE', 'TTTEEE']
dataroots = [ getattr(s, 'defdata_' + p) for p in pol]
dataroots += [dataroots[1].replace('lowEB', 'lowTEB'), dataroots[2].replace('lowEB', 'lowTEB')]

for par, marker in zip(['', 'nnu', 'mnu', 'Alens', 'r', 'yhe', 'nrun'], [None, 3.046, 0.06, 1, None, 0.2449, 0]):
    g.newPlot()
    base = 'base_'
    if par: base += par + '_'
    roots = [base + dat for dat in dataroots]
    labels = [s.datalabel[r] for r in dataroots]

    g.settings.legend_frac_subplot_margin = 0.15

    plotpars = [ 'zrei', 'H0', 'omegabh2', 'thetastar', 'A', 'tau', 'omegam', 'omegach2', 'ns', 'sigma8']
    if par: plotpars[0] = par
    g.plots_1d(roots, plotpars, nx=5, legend_ncol=len(roots), legend_labels=labels, share_y=True, markers=[marker])

    g.export(tag=par)


