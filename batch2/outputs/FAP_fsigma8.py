
import planckStyle as s
from paramgrid import batchjob
from pylab import *
from getdist.densities import Density2D


roots = ['base_' + s.defdata + '_lensing']
g = s.getSinglePlotter(ratio=1)

pars = g.get_param_array(roots[0], ['FAP057', 'fsigma8z057'])

def RSDdensity(FAPbar, f8bar, covfile):
    incov = loadtxt(covfile)
    invcov = inv(inv(incov)[1:, 1:])

    FAPv = np.arange(0.56, 0.78, 0.003)
    f8v = np.arange(0.28, 0.63, 0.003)

    FAP, f8 = np.meshgrid(FAPv, f8v)
    like = (FAP - FAPbar) ** 2 * invcov[0, 0] + 2 * (FAP - FAPbar) * (f8 - f8bar) * invcov[0, 1] + (f8 - f8bar) ** 2 * invcov[1, 1]

    density = Density2D(FAPv, f8v, exp(-like / 2))
    density.contours = exp(-np.array([1.509, 2.4477]) ** 2 / 2)
    return density



FAPbar = 0.6725
f8bar = 0.4412
density = RSDdensity(FAPbar, f8bar, batchjob.getCodeRootPath() + 'data/sdss_DR11CMASS_RSD_bao_invcov_Samushia.txt')
g.add_2d_contours(roots[0], 'FAP057', 'fsigma8z057', filled=True, density=density)


# CS = contourf(FAP, f8, like, origin='lower', levels=[2.279, 5.991], colors='r')


FAPbar = .683
f8bar = 0.422
density = RSDdensity(FAPbar, f8bar, batchjob.getCodeRootPath() + 'data/sdss_DR11CMASS_RSD_bao_invcov_Beutler.txt')
g.add_2d_contours(roots[0], 'FAP057', 'fsigma8z057', filled=False, density=density, ls=':', alpha=0.5)

g.add_2d_contours(roots[0], 'FAP057', 'fsigma8z057', filled=True, plotno=3)


g.add_legend(['BOSS CMASS (Samushia et al.)', 'BOSS CMASS (Beutler et al.)', s.defplanck + '+lensing'], legend_loc='upper left')
g.setAxes(params=pars)
g.export()






