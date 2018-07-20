import planckStyle as s

g = s.getSinglePlotter()

g.settings.solid_colors = [('#8CD3F5', '#006FED'), ('#F7BAA6', '#E03424'), 'g', 'cadetblue', 'olive',
                           'darkcyan']

roots = []
roots.append('base_lensing_lenspriors')
roots.append('base_plikHM_TE_lowE')
roots.append('base_plikHM_TT_lowl_lowE')
roots.append('base_plikHM_TTTEEE_lowl_lowE_lensing')
g.plot_2d(roots, [u'omegamh2', u's8omegamp25'], filled=True, shaded=False)
g.add_legend(
    [s.planck + ' ' + s.lensonly, s.datalabel[s.defdata_TE], s.planckTT, s.planckall + '+lensing'],
    colored_text=True, legend_loc='lower right', align_right=True)
g.export()
