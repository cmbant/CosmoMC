import planckStyle as s

g = s.getSubplotPlotter()
roots = ['base_Alens_plikHM_TT_lowl_lowE', 'base_Alens_plikHM_TTTEEE_lowl_lowE', 'base_plikHM_TTTEEE_lowl_lowE',
         'base_Alens_CamSpecHM_TTTEEE_lowl_lowE']
for i, root in enumerate(roots):
    samples = g.getSamples(root)
    p = samples.getParams()
    samples.addDerived(p.rmsdeflect ** 2, 'vardeflect', label=r'$\langle |\nabla\phi|^2\rangle\,[{\rm arcmin}^2]$')
    roots[i] = samples

yparams = [u'Alens', u'vardeflect']
xparams = [u'omegabh2', u'omegach2', 'ns', u'H0', u'omegam', u'sigma8']
g.rectangle_plot(xparams, yparams, roots=roots, ymarkers=[1, None], filled=[True] * 3 + [False],
                 colors=g.settings.solid_colors[:3] + ['k'], ls=['-'] * 3 + ['--'],
                 legend_labels=[s.planckTT + r' ($\Lambda{\rm CDM}+A_L$)', s.planckall + r' ($\Lambda{\rm CDM}+A_L$)',
                                s.planckall + r' ($\Lambda{\rm CDM}$)'])
g.export()
