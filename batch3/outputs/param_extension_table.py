import planckStyle as s

g = s.getSinglePlotter()

pars = [[], ['r'], ['nrun'], ['nrun', 'r'], ['nrunrun', 'nrun'], ['nnu'], ['nnu', 'nrun'], ['mnu'], ['mnu', 'nnu'],
        ['meffsterile', 'nnu'], ['alpha1'], ['w'], ['omegak'], ['yhe'], ['yhe', 'nnu'], ['Alens']]

toppars = ['omegabh2', 'omegach2', 'theta', 'H0', 'ns', 'logA']

datatag = 'plikHM_TTTEEE_lowl_lowE_lensing'

lines = []
heading = ''
for i, par in enumerate(pars):
    print(par)
    samples = g.getSamples('_'.join(par), datatag)
    samples.paramNames.setLabelsAndDerivedFromParamNames(g.settings.param_names_for_labels)
    varied = [samples.paramNames.parWithName(p).getLabel() for p in par]
    line = ', '.join(varied).replace(r'[\mathrm{eV}]', '')
    latex = samples.getLatex(toppars, limit=1, err_sig_figs=2)
    if i == 0:
        heading += '& ' + ' & '.join(latex[0]) + ''

    line += r' & ' + ' & '.join(latex[1])
    lines.append(line)
print(heading + '\\cr\n\\cr\n' + '\\cr\n'.join(lines) + '\\cr\n')
