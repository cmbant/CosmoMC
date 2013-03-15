import planckStyle as s
from pylab import *

g=s.getSinglePlotter()

g.plot_3d('base_nnu_meffsterile_planck_lowl_lowLike_highL', ['meffsterile', 'nnu', 'omegach2'], lims=[0, 3, 3.046, 4.9])

m = arange(0.0, 3, 0.1)
dN = (m / 6) ** 4.3
plot(m, dN, '--k')

g.export('msterile-nnu')

