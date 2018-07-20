import planckStyle as s

g=s.getSinglePlotter()
roots = ['base_nrun_r_plikHM_TTTEEE_lowl_lowE_lensing', 'base_nrun_plikHM_TTTEEE_lowl_lowE_post_lensing']
g.plot_3d(roots, ['ns02', 'nrun', 'ns'])
g.add_y_marker(0)
g.export()
