from __future__ import print_function
import planckStyle as s
import pylab as plt
import numpy as np
from numpy.linalg import inv
from scipy import stats

g = s.getSinglePlotter()


def quantify_shift(pars, single, joint):
    print(single.name_tag, joint.name_tag, pars)
    Csingle = single.cov(pars)
    Cjoint = joint.cov(pars)
    p1 = single.mean(pars)
    p2 = joint.mean(pars)
    delta = p1 - p2
    diffcov = Csingle - Cjoint
    chi2 = inv(diffcov).dot(delta).dot(delta)
    print('Parameter shift chi2:', chi2, 1 - stats.chi2.cdf(chi2, len(pars)))
    print('d chi^2 single expected', np.trace(inv(Csingle).dot(diffcov)))
    print('\n')


if False:
    planck = g.sampleAnalyser.samplesForRoot('base_plikHM_TTTEEE_lowl_lowE_lensing')
    DESplanck = g.sampleAnalyser.samplesForRoot('base_plikHM_TTTEEE_lowl_lowE_DES_post_lensing')
    plancknolens = g.sampleAnalyser.samplesForRoot('base_plikHM_TTTEEE_lowl_lowE')
    DESplancknolens = g.sampleAnalyser.samplesForRoot('base_plikHM_TTTEEE_lowl_lowE_DES')

    DES = g.sampleAnalyser.samplesForRoot('base_DES_DESpriors')
    DESlens = g.sampleAnalyser.samplesForRoot('base_DESlens_DESpriors')

    lowl = g.sampleAnalyser.samplesForRoot('base_plikHM_TTTEEE_lowl')
    lowl_lowE = g.sampleAnalyser.samplesForRoot('base_plikHM_TTTEEE_lowl_lowE')

    quantify_shift(['S8', 'omegam'], planck, DESplanck)
    quantify_shift(['S8', 'omegam'], plancknolens, DESplancknolens)

    quantify_shift(['S8', 'omegam', 'logA', 'omegabh2', 'ns'], planck, DESplanck)
    quantify_shift(['S8', 'omegam', 'logA', 'omegabh2', 'ns'], plancknolens, DESplancknolens)

    quantify_shift(['S8', 'omegam'], DES, DESplanck)
    quantify_shift(['S8', 'omegam', 'logA', 'omegabh2', 'ns'], DES, DESplanck)

    pars = ['S8', 'omegam']
    delta = planck.mean(pars) - DES.mean(pars)
    chi2 = inv(planck.cov(pars) + DES.cov(pars)).dot(delta).dot(delta)
    print('Independent shift chi2(2)', chi2, 1 - stats.chi2.cdf(chi2, len(pars)))

    pars = ['S8', 'omegam', 'logA', 'omegabh2', 'ns']
    delta = planck.mean(pars) - DES.mean(pars)
    chi2 = inv(planck.cov(pars) + DES.cov(pars)).dot(delta).dot(delta)
    print('Independent shift chi2(5)', chi2, 1 - stats.chi2.cdf(chi2, len(pars)))

if True:
    roots = []
    roots.append('base_DESwt_DESpriors')
#    roots.append('base_DESlens_DESpriors')
    roots.append('base_DES_DESpriors')

    # roots.append(g.getRoot('','DES_DESpriors_CookeDH_lensing_BAO'))
    p2='sigma8'

    roots.append('base_plikHM_TT_lowl_lowE')
    roots.append('base_plikHM_TTTEEE_lowl_lowE')
    #roots.append('base_plikHM_TTTEEE_lowl_lowE_DES_post_lensing')

    #roots.append('base_plikHM_TTTEEE_lowl_lowE_post_BAO_zre6p5')
    g.plot_2d(roots, [u'omegam', p2], filled=True, shaded=False)

    g.add_2d_contours('base_plikHM_TTTEEE_lowl_lowE_DES_post_lensing', u'omegam', p2,
                      ls='-')

    g.add_2d_contours('base_CamSpecHM_TTTEEE_lowl_lowE','omegam',p2, ls=':',color='navy')

    g.add_legend(
        ['DES $\gamma_t+w$', 'DES (joint)', s.planckTT, s.planckall,
         '+DES+lensing'], colored_text=True, align_right=True, label_order=[2, 3, 4, 0, 1])
    plt.xlim(0.17, 0.39)
    plt.ylim(0.61, 1.07)
    g.export()
