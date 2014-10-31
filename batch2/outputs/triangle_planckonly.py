
import planckStyle as s
g = s.getSubplotPlotter()



roots = [g.getRoot('', s.defdata_TTonly),
       g.getRoot('', s.defdata_allNoLowE),
       g.getRoot('', s.defdata_all)]

labs = [s.planckTT, s.NoLowLE, s.planckall]

g.settings.plot_args = [{'color':'olive'}, {'color':'red'}, {'color': 'black', 'alpha':1}]
g.triangle_plot(roots, ['theta', 'omegabh2', 'omegach2', 'logA', 'ns', 'tau'], plot_3d_with_param='H0', filled_compare=False,
                    legend_labels=labs)
g.export()


