import planckStyle as s
from pylab import *
import sys

g=s.getSubplotPlotter()

g.settings.lineM = ['-k', '-r', '-b', '-g', '--m', '-y']

markers=[None, 1, 3.046,0,0.06,0.247816]


pars=['','Alens','nnu','nrun','mnu','yhe']

labels=[r'{\tt CamSpec}',r'$\ell_{\rm min}=1200$ for 143, 217',r'$\ell_{\rm max}=2000$',r'No $217\times217$',r'{\tt Plik}']

exts=['','_highL']

comparePars=['ns']+pars[1:]
xlims=[(0.93,0.99),(0.81,1.9),(2,5.3),(-0.059,0.027),(0,1.7),(0.16,0.39) ]
for i,ext in enumerate(exts):
    parroots=[]
    for ix,par in enumerate(comparePars):
        base = 'base_'
        if len(par)>0 and ix!=0: base+=par+'_'
        roots=[base+'planck_lowl_lowLike'+ext,base+'planck_lmin1200_lowl_lowLike'+ext,base+'planck_lmax2000_lowl_lowLike'+ext,base+'planck_no217auto_lowl_lowLike'+ext,base+'plik_lowl_lowLike'+
ext]
        parroots.append(roots)

    g.plots_1d(parroots, comparePars,legend_labels=labels, nx=3,legend_ncol=3,xlims=xlims,roots_per_param=True, markers=markers)
    g.exportExtra('param_extensions_camspec_cuts_compare'+ext)
    g.export('param_extensions_camspec_cuts_compare'+ext+'.png')
    g.newPlot()
#     

sys.exit()
for ext in exts:
    for par in pars:
        base = 'base_'
        if len(par)>0: base+=par+'_'
        roots=[base+'planck_lowl_lowLike'+ext,base+'planck_lmin1200_lowl_lowLike'+ext,base+'planck_lmax2000_lowl_lowLike'+ext,base+'planck_no217auto_lowl_lowLike'+ext,base+'plik_lowl_lowLike'+ext]

        g.plots_1d(roots, legend_labels=labels, nx=5,legend_ncol=5)
        g.exportExtra(base+'camspec_cuts_compare'+ext)
        g.newPlot()
#        if par=='': 
#            g.plots_1d(roots, legend_labels=labels,nx=4, legend_ncol=3)
#            g.exportExtra(base+'camspec_cuts_compare_thin'+ext)
#            g.newPlot()
