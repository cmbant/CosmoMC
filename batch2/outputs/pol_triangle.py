
import planckStyle as s
g = s.getSubplotPlotter()


params = ['theta', 'omegabh2', 'omegach2', 'logA', 'ns', 'tau']


for par in ['', 'nnu']:
    g.newPlot()
    roots = [ g.getRoot(par, s.defdata_root + '_EE_lowEB'),
         g.getRoot(par, s.defdata_root + '_TE_lowEB'),
         g.getRoot(par, s.defdata_root + '_TT_lowTEB'),
         g.getRoot(par, s.defdata_root + '_TTTEEE_lowTEB'),
       ]


    labs = ['EE + lowEB', 'TE + lowEB', 'TT + lowTEB', 'TTTEEE + lowTEB']
    labs = [s.planck + ' ' + t for t in labs]
    params = ['theta', 'omegabh2', 'omegach2', 'logA', 'ns', 'tau']
    if par: params = [par] + params
    # g.settings.plot_args = [{'color':'olive'}, {'color':'red'}, {'color': 'black', 'alpha':1}]
    g.triangle_plot(roots, params, filled_compare=True,
    #                contour_args=[{'filled':False}, {'filled':False}, {'filled':False}, {'filled':True}],
                        legend_labels=labs)
    g.export(tag=par)


