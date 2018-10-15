import planckStyle as s
import numpy as np
import pylab as plt
import os

g = s.getSinglePlotter()


def fit_file(root):
    for chain_dir in g.sampleAnalyser.chain_dirs:
        jobItem = chain_dir.resolveRoot(root)
        if jobItem:
            return np.loadtxt(jobItem.chainRoot + '.minimum.theory_cl')
            break


bands = np.loadtxt(os.path.join(os.path.dirname(__file__),
                                '../../data/planck_lensing_2017/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8_bandpowers.dat'))
agr = np.loadtxt(os.path.join(os.path.dirname(__file__),
                              '../../data/planck_lensing_2017/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_agr2_bandpowers.dat'))

bands[:, 4:6] *= 1e7
agr[:, 4:6] *= 1e7

ALens_root = 'base_Alens_plikHM_TTTEEE_lowl_lowE'
fit = fit_file('base_plikHM_TTTEEE_lowl_lowE')
fit_lens = fit_file('base_lensing_lenspriors')

fit_ALens = fit_file(ALens_root)

plt.semilogx(fit_lens[:, 0], fit_lens[:, 5] * 1e7, lw=0.7, color='C2', ls='-')
plt.semilogx(fit[:, 0], fit[:, 5] * 1e7, lw=0.7, color='C0', ls='-.')

plt.semilogx(fit_ALens[:, 0], fit_ALens[:, 5] * 1e7, ls='--', lw=0.7, color='C0', alpha=0.4)

plt.scatter(bands[:, 3], bands[:, 4], marker='.', s=3, color='k', zorder=100)

for lmin, lmax, PP, error in zip(bands[:, 1].astype(np.int), bands[:, 2].astype(np.int),
                                 bands[:, 4], bands[:, 5]):
    plt.fill_between(np.arange(lmin, lmax + 1), PP - error, PP + error, color='gray',
                     alpha=0.3, linewidth=0)

plt.errorbar(agr[:, 3], agr[:, 4], agr[:, 5], fmt='.', markersize=3, alpha=0.85, color='C1', zorder=100)

ALens = g.param_latex_label(ALens_root, 'Alens', labelParams='clik_latex.paramnames')
legend = plt.legend(['lensing+priors',s.planckall, s.planckall +' ('+ ALens+')'], fontsize=6)
legend.get_frame().set_edgecolor('none')

plt.xlim(8, 2048)
#plt.axhline(0, color='k', lw=0.5)
plt.ylim(0,None)
plt.xlabel('$L$')
plt.ylabel(r'$10^7 [L(L+1)]^2 C_L^{\phi\phi}/2\pi$')

g.export()
