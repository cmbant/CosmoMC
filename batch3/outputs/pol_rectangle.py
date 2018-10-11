import planckStyle as s

for dat in ['plikHM', 'CamSpecHM']:

    g = s.getSubplotPlotter(subplot_size=2)

    roots = []
    roots.append('base_%s_EE_lowE_BAO' % dat)
    roots.append('base_%s_TE_lowE' % dat)
    roots.append('base_%s_TT_lowl_lowE' % dat)
    roots.append('base_%s_TTTEEE_lowl_lowE' % dat)
    xparams = [u'omegabh2', u'omegach2', u'theta', u'tau', u'ns', 'logA']
    yparams = [u'H0', u'omegam', u'sigma8']
    labels = [s.datalabel[s.defdata_EE] + '+BAO', s.datalabel[s.defdata_TE], s.datalabel[s.defdata_TT], s.planckall]
    filled = True
    # ranges = {'omegabh2':}
    g.rectangle_plot(xparams, yparams, roots=roots, filled=filled, legend_labels=labels)
    if dat == 'CamSpecHM':
        for int, ax1, ax in zip(ints, axs, g.subplots.reshape(-1)):
            ax.set_xlim(int[0])
            ax.set_ylim(int[1])
            ax.set_yticks(ax1.get_yticks())
            ax.set_xticks(ax1.get_xticks())

    axs = g.subplots.reshape(-1)
    ints = [[ax.xaxis.get_view_interval(), ax.yaxis.get_view_interval()] for ax in g.subplots.reshape(-1)]
    ax = g.get_axes_for_params('omegabh2', 'omegam')
    print ax,g.get_axes_for_params('omegam', 'omegab2')
    if ax:
        ax.set_yticks([0.28, 0.30, 0.32, 0.34, 0.36])
    g.export(tag='' if dat == 'plikHM' else dat)
