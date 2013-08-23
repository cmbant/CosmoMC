import planckStyle as s

g=s.plotter
g.settings.setWithSubplotSize(3.5)

g.settings.lab_fontsize=10
g.settings.legend_frac_subplot_margin=-0.02
g.settings.legend_fontsize = 9

labels=[s.WP+'+BAO',s.WP+'+Union2.1',s.WP+'+SNLS',s.WP]
roots=['base_w_planck_lowl_lowLike_BAO', 'base_w_planck_lowl_lowLike_Union2','base_w_planck_lowl_lowLike_SNLS','base_w_planck_lowl_lowLike']
#roots=['base_w_planck_lowl_lowLike','base_w_planck_lowl_lowLike_BAO' ]

g.plots_1d(roots,['w'],legend_labels=labels,legend_ncol=2)

g.add_x_marker(-1)

g.export('w')


