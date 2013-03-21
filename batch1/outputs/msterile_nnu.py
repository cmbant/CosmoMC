import planckStyle as s
from pylab import *

g=s.getSinglePlotter()

root='base_nnu_meffsterile_planck_lowl_lowLike_highL'

g.plot_3d(root, ['meffsterile', 'nnu', 'omegach2'], lims=[0, 3, 3.046, 4.9])

meff = g.param_latex_label(root, 'meffsterile')

m = np.arange(0, 3.0, 0.01)
N = np.arange(3.047, 4.9, 0.01)
N, m = np.meshgrid(N, m)
z = m / (N - 3.046) ** 0.75
CS = contour(m, N, z, origin='lower', levels=[0.5, 1, 2, 5, 10], colors='k', linestyles='--',linewidths=0.3, extent=[0,3,3.046,4.9])
clabel(CS, CS.levels, inline=True,fmt='%1.1f',fontsize=7)

z = m / (N - 3.046)
CS = contour(m, N, z, origin='lower', levels=[0.5, 1, 2, 5, 10], colors='gray', linestyles=':',linewidths=0.3, extent=[0,3,3.046,4.9])
#clabel(CS, CS.levels, inline=True,fmt='%1.1f',fontsize=7)

xlim([0, 3])
ylim([3.046, 4.9])
#xlabel(meff+r'$(=\Omega_\nu h^2/94.1)$')

g.export('msterile-nnu')

