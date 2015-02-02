
import planckStyle as s
g = s.getSubplotPlotter()


params = ['theta', 'omegabh2', 'omegach2', 'logA', 'ns', 'tau']


# for par in ['', 'nnu']:
for par in ['']:
    g.newPlot()
    dataroots = [s.defdata_root + '_EE_lowEB', s.defdata_root + '_TE_lowEB', s.defdata_root + '_TT_lowTEB', s.defdata_root + '_TTTEEE_lowTEB']
    roots = [ g.getRoot(par, root) for root in dataroots]

    labs = [s.datalabel[t] for t in dataroots]
    params = ['theta', 'omegabh2', 'omegach2', 'logA', 'ns', 'tau']
    if par: params = [par] + params
    # g.settings.plot_args = [{'color':'olive'}, {'color':'red'}, {'color': 'black', 'alpha':1}]
    g.triangle_plot(roots, params, filled_compare=True,
    #                contour_args=[{'filled':False}, {'filled':False}, {'filled':False}, {'filled':True}],
                        legend_labels=labs)

    g.export(tag=par)


