import planckStyle as s
from pylab import *

if True:
    g = s.getSubplotPlotter()
    samps = g.sampleAnalyser.samplesForRoot('base_plikHM_TTTEEE_lowl_lowE_lensing')
    p = samps.getParams()
    for i in range(100):
        power = i/100.
        newp=p.sigma8*p.omegam**power
        newp2 = p.sigma8 * (p.omegam/0.3) ** power
        print(r'\sigma_8 \Omm^{%.2f} = %.4f \pm %.4f, \sigma_8 (\Omm/0.3)^{%.2f} = %.3f \pm %.3f'%(power,samps.mean(newp),samps.std(newp),power, samps.mean(newp2),samps.std(newp2)))

    import sys
    sys.exit()

    samps = g.sampleAnalyser.samplesForRoot('base_mnu_plikHM_TTTEEE_lowl_lowE_lensing_BAO')
    p = samps.getParams()
    samps.filter(p.mnu > 0.06)
    samps.ranges.setRange('mnu', (0.06, None))
    samps.updateBaseStatistics()
    print(samps.getInlineLatex('mnu', limit=2))

g = s.getSinglePlotter()
g.settings.norm_prob_label = r'Probability density [${\rm eV}^{-1}$]'
g.settings.lineM = ['-b', '-r', '-g', '--b', '--r', '--g', '-c', '-y']

labels = [s.planckall, '+lensing', '+BAO', ]
roots = [s.defdata_all,
         s.defdata_all + '_lensing',
         s.defdata_all + '_lensing_BAO']
roots = [g.getRoot('mnu', root) for root in roots]

g.plot_1d(roots, 'mnu', normalized=True)
g.add_legend(labels, legend_loc='upper right')

xlim([0, 0.6])

g.export()
