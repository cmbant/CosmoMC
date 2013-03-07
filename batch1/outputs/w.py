import planckStyle as s

g=s.plotter
g.settings.setWithSubplotSize(4)

g.settings.tight_layout=True
g.settings.lab_fontsize=10
g.settings.legend_frac_subplot_margin=-0.02

labels=[s.WP,s.WP+'+BAO',s.WP+'+SNLS',s.WP+'+Union2.1']
roots=['base_w_planck_lowl_lowLike','base_w_planck_lowl_lowLike_BAO', 'base_w_planck_lowl_lowLike_SNLS','base_w_planck_lowl_lowLike_Union2']
#roots=['base_w_planck_lowl_lowLike','base_w_planck_lowl_lowLike_BAO' ]

g.plots_1d(roots,['w'],legend_labels=labels,legend_ncol=2)
g.add_x_marker(-1)

g.export('w_test2')


