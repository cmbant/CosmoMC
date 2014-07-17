import planckStyle as s
from pylab import *

g=s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)
outdir=''
roots=['base_nnu_meffsterile_planck_lowl_lowLike']

g.plot_3d(roots, ['nnu', 'ns', 'H0'])

g.export(os.path.join(outdir,'meffsterile_ns_H0.pdf'))
