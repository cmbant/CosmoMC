import planckStyle as s
import sys

g=s.getSubplotPlotter()

g.settings.lineM = ['-k', '-r', '-b', '-g', '--m', '-y']
g.settings.legend_frac_subplot_margin=0.3
g.settings.param_names_for_labels = clik_cib.paramnames

markers=[None, 1, 3.046,0,0.06,0.247816]


pars=['','Alens','nnu','nrun','mnu','yhe'] 

labels=[r'$\gamma^{\rm CIB}=0.7\pm 0.2$, $d\gamma^{\rm CIB}/d\ln \ell=0$',r'$\gamma^{\rm CIB}=0.7\pm 0.2$, free running',r'no $\gamma^{\rm CIB}$ prior and free running']

exts=['']

comparePars=['ns']+pars[1:] +  ['acib143','aps217','acib217','ncib','nruncib','asz143']
for ext in exts:
    parroots=[]
    p=[]
    for ix,par in enumerate(comparePars):
        base = 'base_'
        if len(par)>0 and ix!=0 and ix<6: base+=par+'_'

        roots=[base+'planck_lowl_lowLike_highL'+ext,base+'planck_cibrun_lowl_lowLike_highL'+ext,base+'planck_cibrunnoprior_lowl_lowLike_highL'+ext]
        parroots.append(roots)
        print par, roots[0]
        p.append(g.check_param(roots[1],par))

    g.plots_1d(parroots, p,legend_labels=labels, nx=3,legend_ncol=2,roots_per_param=True, markers=markers)
    g.exportExtra('param_extensions_camspec_cib_compare'+ext)
    g.export('param_extensions_camspec_cib_compare'+ext+'.png')
    g.newPlot()


