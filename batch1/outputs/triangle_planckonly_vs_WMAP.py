#import GetDistPlots
#g=GetDistPlots.GetDistPlotter('main/plot_data')

import planckStyle as s
g=s.plotter

roots = ['base_planck_CAMspec_lowl_post_lensing', 'base_planck_CAMspec_lowl_lowLike','base_WMAP']
labs=[s.lensing, s.WP,'WMAP9']
#labs=['Planck+lensing','Planck+WP','WMAP9']

g.settings.plot_args= [None,None, {'color': 'gray', 'alpha':0.6}]
g.triangle_plot(roots, ['omegabh2', 'omegach2', 'ns', 'tau', 'omegal', ], plot_3d_with_param='H0', filled_compare=False,
                    legend_labels=labs)
g.export('triangle_planckonly_vs_WMAP')
#g.export('triangle_planckonly_vs_WMAP.eps')


