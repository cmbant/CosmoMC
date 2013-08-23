import planckStyle as s
import sys

g=s.getSubplotPlotter()

g.settings.lineM = ['-k', '-r', '-b', '-g', '--m', '-y']

markers=[None, 1, 3.046,0,0.06,0.247816]

pars=['','Alens','nnu','nrun','mnu','yhe']

labels=[r'{\tt CamSpec}',r'100 + 143 ($\ell_{\rm max}=2000$)',r'100 + 217 ($\ell_{\rm max}=2000$)',r'100 + 217 ($\ell_{\rm max}=2500$)']

exts=['']

comparePars=['ns']+pars[1:]
for ext in exts:
    parroots=[]
    for ix,par in enumerate(comparePars):
        base = 'base_'
        if len(par)>0 and ix!=0: base+=par+'_'
        roots=[base+'planck_lowl_lowLike'+ext,base+'planck_143only_lowl_lowLike'+ext,base+'planck_217onlylmax2000_lowl_lowLike'+ext,base+'planck_217only_lowl_lowLike'+ext]
        parroots.append(roots)

    g.plots_1d(parroots, comparePars,legend_labels=labels, nx=3,legend_ncol=2,roots_per_param=True, markers=markers)
    g.exportExtra('param_extensions_camspec_freq_compare'+ext)
    g.export('param_extensions_camspec_freq_compare'+ext+'.png')
    g.newPlot()

