import GetDistPlots
g=GetDistPlots.GetDistPlotter('main/plot_data')
roots = ['base_planck_CAMspec_lowl_lowLike', 'base_planck_CAMspec_lowl_lowLike_post_BAO','base_WMAP']

g.settings.lab_fontsize += 4
g.settings.axes_fontsize+=1
g.settings.plot_meanlikes = False
g.settings.plot_args= [None,None, {'color': 'gray'}]
g.triangle_plot(roots, ['omegabh2', 'omegach2', 'ns', 'tau', 'omegal', ], plot_3d_with_param='H0', filled_compare=False,
                    legend_labels=['Planck+WP', 'Planck+WP+BAO','WMAP9'])
g.export('triangle_plot_vs_WMAP.eps')

#import GetDistPlots
#roots = ['base_Alens_planck_CAMspec_lowl_lowLike', 'base_Alens_planck_CAMspec_lowl_lowLike_highL','base_Alens_planck_CAMspec_lowl_lowLike_#post_lensing']
#g.triangle_plot(roots, ['omegabh2', 'omegach2', 'Alens', 'omegal'], legend_labels=['Planck', 'Planck+highL','Planck_lensing'])
#g.export('Alens-highL-compare.pdf')
