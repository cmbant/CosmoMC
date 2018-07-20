
import planckStyle as s
g = s.getSubplotPlotter()


params = ['theta', 'omegabh2', 'omegach2', 'logA', 'ns', 'tau']

for camspec in [True, False]:
    for par in ['']:
        g.newPlot()
        dataroots = [s.defdata_root + '_EE_lowE', s.defdata_root + '_TE_lowE', s.defdata_root + '_TT_lowl_lowE', s.defdata_root + '_TTTEEE_lowl_lowE']
        labs = [s.datalabel[t] for t in dataroots]

        if camspec:
            dataroots = [x.replace('plikHM','CamSpecHM') for x in dataroots]
        print(dataroots)
        roots = [ g.getRoot(par, root) for root in dataroots]

        if par: params = [par] + params

        g.triangle_plot(roots, params, filled_compare=True, legend_labels=labs)

        g.export(tag=par + ('_CamSpec' if camspec else ''))


