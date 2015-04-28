
import planckStyle as s
g = s.getSubplotPlotter()


params = ['theta', 'omegabh2', 'omegach2', 'logA', 'ns', 'tau']

for par in ['']:
    g.newPlot()
    dataroots = [s.defdata_root + '_EE_lowEB', s.defdata_root + '_TE_lowEB', s.defdata_root + '_TT_lowTEB', s.defdata_root + '_TTTEEE_lowTEB']
    roots = [ g.getRoot(par, root) for root in dataroots]

    labs = [s.datalabel[t] for t in dataroots]
    if par: params = [par] + params

    g.triangle_plot(roots, params, filled_compare=True, legend_labels=labs)

    g.export(tag=par)


