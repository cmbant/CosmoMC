import planckStyle as s

g=s.plotter

roots = ['base_planck_lowl_lowLike', 'base_planck_lowl_lowLike_post_BAO','base_WMAP']

g.settings.lab_fontsize += 4
g.settings.axes_fontsize+=1
g.settings.plot_args= [None,None, {'color': 'gray'}]
g.triangle_plot(roots, ['omegabh2', 'omegach2', 'ns', 'tau', 'omegal', ], plot_3d_with_param='H0', filled_compare=False,
                    legend_labels=[s.WP, s.WP+'+BAO','WMAP9'])
g.export('triangle_plot_vs_WMAP')

#import GetDistPlots
#roots = ['base_Alens_planck_lowl_lowLike', 'base_Alens_planck_lowl_lowLike_highL','base_Alens_planck_lowl_lowLike_#post_lensing']
#g.triangle_plot(roots, ['omegabh2', 'omegach2', 'Alens', 'omegal'], legend_labels=['Planck', 'Planck+highL','Planck_lensing'])
#g.export('Alens-highL-compare.pdf')
