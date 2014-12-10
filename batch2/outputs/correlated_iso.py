import planckStyle as s
from pylab import *
g = s.getSinglePlotter()

dataroots = [s.defdata_TT, s.defdata_all]
roots = ['base_alpha1_' + root for root in dataroots]

if False:
    samples = g.sampleAnalyser.samplesForRoot(roots[1])
    p = samples.getParams()
    samples.filter(p.alpha1 < 0)
    samples.updateChainBaseStatistics(ini_settings={'limits[alpha1]':'N 0'})
    p = samples.getParams()
    print samples.confidence(p.alpha1, 0.05, upper=False)


g.plot_2d(roots, ['ns', 'alpha1'], filled=True)
g.add_legend([s.datalabel[x] for x in dataroots], legend_loc='lower right')
g.add_y_marker(0)
g.export()
