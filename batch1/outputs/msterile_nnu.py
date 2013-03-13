import planckStyle as s

g=s.getSinglePlotter()

g.plot_3d('base_nnu_meffsterile_planck_lowl_lowLike_highL', ['meffsterile', 'nnu', 'omegach2'], lims=[0, 3, 3.046, 4.9])

g.export('msterile-nnu')

