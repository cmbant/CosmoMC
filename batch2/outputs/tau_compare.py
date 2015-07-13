import planckStyle as s

g = s.getPlotter()

g.settings.lineM = ['-k', '-r', '-b', '-g', '-m', '--c']

labels = [s.planckTT, s.NoLowLE, s.planckTT + '+lensing', s.NoLowLE + '+lensing', s.NoLowLE + '+lensing+BAO',
          s.planckTTlowTEB]
roots = [s.defdata_TTonly,
         s.defdata_allNoLowE,
         s.defdata_TTonly + '_lensing',
         s.defdata_allNoLowE + '_lensing',
         s.defdata_allNoLowE + '_lensing_post_BAO',
         s.defdata,
         ]
roots = ['base_' + root for root in roots]

g.settings.legend_frac_subplot_margin = 0.05
g.plots_1d(roots, ['omegabh2', 'thetastar', 'A', 'tau', 'omegam', 'omegach2', 'ns', 'sigma8', 'zrei', 'H0'], nx=5,
           legend_ncol=len(roots), legend_labels=labels, share_y=True)

g.export()


