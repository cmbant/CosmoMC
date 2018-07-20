import planckStyle as s
import pylab as plt

g = s.getSinglePlotter()

g.plot_3d('base_omegak_plikHM_TTTEEE_lowl_lowE', ['omegak', 'omegam', 'H0'], alpha_samples=True)

plt.axvline(0, c='k', ls='--', color='gray', alpha=0.5, lw=0.7)

g.add_2d_contours('base_omegak_plikHM_TTTEEE_lowl_lowE', 'omegak', 'omegam', filled=False, ls='--', color='k')

g.add_2d_contours('base_omegak_plikHM_TTTEEE_lowl_lowE_lensing', 'omegak', 'omegam', filled=False, ls='-', color='g')

# g.add_2d_contours('base_omegak_plikHM_TTTEEE_lowl_lowE_BAO_post_lensing','omegak','omegam',filled=True, alpha=0.85)

g.add_2d_contours('base_omegak_plikHM_TTTEEE_lowl_lowE_BAO_post_lensing', 'omegak', 'omegam',
                  filled=True, color='purple', alpha=0.85)

g.add_2d_contours('base_omegak_CamSpecHM_TTTEEE_lowl_lowE', 'omegak', 'omegam', filled=False, ls=':', color='k')

g.add_legend([s.shortlabel[s.defdata_all], '+lensing', '+BAO'])
plt.gca().set_xticks([-0.1, -0.05, 0])

g.export()
