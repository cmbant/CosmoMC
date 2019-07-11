from __future__ import absolute_import
from __future__ import print_function
import planckStyle as s

g = s.getSinglePlotter()


class comb(object):
    def __init__(self, varname, prior, title, blockstart=False):
        self.varname = varname
        self.prior = prior
        self.title = title
        self.blockstart = blockstart


items = []

items += [comb('', 'lenspriors', r'MV conservative $8\le L \le 400$')]
items += [comb('DESlens', 'lenspriors', r'DES lensing joint')]
items += [comb('DES', 'lenspriors', r'DES combined joint')]
items += [comb('theta', 'lenspriors', '$100\\thetaMC = 1.0409 \pm 0.0006$ joint')]
items += [comb('', 'plikHM_TT_lowl_lowE_lensing', r'\planckTT\ joint')]
items += [comb('', 'plikHM_TTTEEE_lowl_lowE_lensing', r'\planckall\ joint')]
items += [comb('conslmin40', 'lenspriors', r'MV conservative $40\le L \le 400$ ', blockstart=True)]
items += [comb('agrlmax425', 'lenspriors', r'MV aggressive $8\le L \le 425$')]
items += [comb('agr2', 'lenspriors', r'MV aggressive $8\le L \le 2048$')]
items += [comb('ptt', 'lenspriors', r'TT conservative $8\le L \le 400$')]
items += [comb('pttagr2', 'lenspriors', r'TT aggressive $8\le L \le 2048$')]
items += [comb('Apr6', 'lenspriors', r"CompSep mask $8\le L \le 400$")]
items += [comb('', 'DESpriors', r'DES priors', blockstart=True)]
items += [comb('', 'DESpriors_CookeDH', r"$'$$'$ + ($\Omega_{\rm b}h^2=0.0222\pm 0.0005$)")]
items += [comb('bfcl', 'lenspriors', r'Best-fit $C^{\rm CMB}_\ell$')]
items += [comb('agr2bfcl', 'lenspriors', r"$'$$'$ (MV aggressive $8\le L \le 2048$)")]
items += [comb('takahashi', 'lenspriors', r"Takahashi {\HALOFIT}")]
items += [comb('agr2takahashi', 'lenspriors', r"$'$$'$ (MV aggressive $8\le L \le 2048$)")]
items += [comb('linear', 'lenspriors', r'Linear theory')]
items += [comb('acc', 'lenspriors', r'Higher accuracy')]
items += [comb('agr2acc', 'lenspriors', r"$'$$'$ (MV aggressive $8\le L \le 2048$)")]

lines = []
heading = ''
for i, item in enumerate(items):
    line = item.title
    if item.blockstart:
        line = '\\hline\n' + line
    roots = ['base_lensing_%s', 'base_lensing_%s_BAO',
             'base_mnu_lensing_%s', 'base_mnu_lensing_%s_BAO']
    for root, pars in zip(roots, [['s8omegamp25'], ['sigma8', 'H0', 'omegam'], ['s8omegamp25'], ['sigma8', 'mnu']]):
        if 'plikHM' in item.prior:
            root = root.replace('lensing_%s', item.prior)
        else:
            root = root % item.prior
            if item.varname: root += '_' + item.varname
        print(root)
        try:
            if 'base_mnu' in root:
                paramtag = 'mnu'
                datatag = root.replace('base_mnu_', '')
            else:
                paramtag = ''
                datatag = root.replace('base_', '')
            root = g.getRoot(paramtag, datatag)
            samples = g.sampleAnalyser.samplesForRoot(root)
            samples.paramNames.setLabelsAndDerivedFromParamNames(g.settings.param_names_for_labels)
            if 'planck' in item.title:
                latex = samples.getLatex(params=pars, limit=1, err_sig_figs=1)
            else:
                latex = samples.getLatex(params=pars, limit=1)
                if 'DESpriors' in item.prior and 'mnu' in pars:
                    latex[1][1] = '\\text{--}'
        except Exception as e:
            print(e)
            print('Missing root:' + root)
            latex = [[''] * len(pars), [''] * len(pars)]
        if i == 0: heading += '& $' + '$ & $'.join(latex[0]) + '$'
        line += r' & $' + '$ & $'.join(latex[1]) + '$'
    lines.append(line)
print(heading + '\\\\\n\\hline\n' + '\\\\\n'.join(lines) + '\\\\\n')
