import GetDistPlots
g=GetDistPlots.GetDistPlotter('main/plot_data')
g.settings.setWithSubplotSize(5)
g.settings.plot_meanlikes=False
g.settings.tight_layout=True
g.settings.prob_label='Probability'
g.settings.lab_fontsize=12
#g.settings.tight_gap_fraction = .05
#g.settings.line_labels=False

labels=['Planck+WP','Planck+WP+BAO','Planck+WP+SNLS','Planck+WP+Union2.1']
roots=['base_w_planck_CAMspec_lowl_lowLike','base_w_planck_CAMspec_lowl_lowLike_BAO', 'base_w_planck_CAMspec_lowl_lowLike_SNLS','base_w_planck_CAMspec_lowl_lowLike_Union2']
#roots=['base_w_planck_CAMspec_lowl_lowLike','base_w_planck_CAMspec_lowl_lowLike_BAO' ]
#g.plots_1d(roots, legend_labels=labels)
#g.export('planck_datasets_1d.pdf')
#g.newPlot()
g.plots_1d(roots,['w'],legend_labels=labels)
g.add_line([-1., -1.], [-2., 2.], zorder=1)
g.export('w_test2.pdf')
g.export('w_test2.eps')
