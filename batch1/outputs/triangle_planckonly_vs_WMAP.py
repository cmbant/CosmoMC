#import GetDistPlots
#g=GetDistPlots.GetDistPlotter('main/plot_data')

import planckStyle as s
g=s.getSubplotPlotter()


roots = ['base_planck_lowl_post_lensing', 'base_planck_lowl_lowLike','base_WMAP']
labs=[s.lensing, s.WP,'WMAP-9']

g.settings.plot_args= [None,None, {'color': 'gray', 'alpha':0.8}]
g.triangle_plot(roots, ['omegabh2', 'omegach2', 'ns', 'tau', 'omegal', ], plot_3d_with_param='H0', filled_compare=False,
                    legend_labels=labs)
g.export('triangle_planckonly_vs_WMAP')


