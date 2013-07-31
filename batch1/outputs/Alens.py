import planckStyle as s
g=s.plotter

labels=['Planck+WP ($A_L$ free)','Planck+WP+highL ($A_L$ free)','Planck+Lensing ($A_L$ free)','Planck+WP+lensing ($A_L$ free)','Planck+WP ($A_L$=1)']
roots=['base_Alens_planck_lowl_lowLike','base_Alens_planck_lowl_lowLike_highL','base_Alens_planck_lowl_post_lensing','base_Alens_planck_lowl_lowLike_post_lensing','base_planck_lowl_lowLike']
g.plots_1d(roots, legend_labels=labels)

g.exportExtra('Alens_compare')

