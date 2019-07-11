import planckStyle as s
import pylab as plt

for with_KIDS in [True, False]:

    g = s.getSinglePlotter(chain_dir=[r'C:\Tmp\Planck\KiDs', r'C:\Tmp\Planck\2017\Dec17'])
    roots = []
    roots.append('base_DESlens_DESpriors')
    roots.append('base_lensing_DESpriors')
    roots.append('base_DESlens_DESpriors_lensing')
    # roots.append('base_DES')
    roots.append('base_' + s.defdata_all)
    g.plot_2d(roots, [u'omegam', u'sigma8'], filled=True, shaded=False)
    g.add_2d_contours('base_DES_DESpriors', u'omegam', u'sigma8', ls='--')

    if with_KIDS:
        g.add_2d_contours('KiDS_lcdm_DESpriors', u'omegam', u'sigma8', ls=':', color='black', alpha=0.5)
        g.add_legend(['DES lensing', r'$\textit{Planck}$ lensing', r'DES lensing+$\textit{Planck}$ lensing',
                      s.planckall, r'DES joint', 'KiDS-450'], align_right=True)
        plt.ylim(None, 1.29)
        g.export(tag='with_KIDS')
    else:
        g.add_legend(['DES lensing', r'$\textit{Planck}$ lensing', r'(DES+$\textit{Planck}$) lensing', s.planckall,
                      r'DES joint'])
        g.export()
