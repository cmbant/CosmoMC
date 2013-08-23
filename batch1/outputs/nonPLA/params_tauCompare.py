import planckStyle as s
import sys

g=s.getSubplotPlotter()
g.settings.lineM = ['-k', '--k', '-r', '--r', '--m', '-y']

markers=[None, 1, 3.046,0,0.06,0.247816]


pars=['','Alens','nnu','nrun','mnu','yhe']

labels=[r'\textit{Planck}+WP',r'\textit{Planck}+$(\tau=0.07\pm 0.013)$',r'\textit{Planck}+WP+highL',r'\textit{Planck}+$(\tau=0.07\pm 0.013)$+highL']

exts=['']

comparePars=['ns']+pars[1:]
for ext in exts:
    parroots=[]
    for ix,par in enumerate(comparePars):
        base = 'base_'
        if len(par)>0 and ix!=0: base+=par+'_'
        roots=[base+'planck_lowl_lowLike'+ext,base+'planck_lowl_tau07'+ext,base+'planck_lowl_lowLike_highL'+ext,base+'planck_lowl_tau07_highL'+ext]
        parroots.append(roots)

    g.plots_1d(parroots, comparePars,legend_labels=labels, nx=3,legend_ncol=2,roots_per_param=True, markers=markers)
    g.exportExtra('param_extensions_camspec_tau_compare'+ext)
    g.export('param_extensions_camspec_tau_compare'+ext+'.png')
    g.newPlot()
