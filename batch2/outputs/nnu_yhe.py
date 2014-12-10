import planckStyle as s
from pylab import *

g = s.getSinglePlotter()

roots = [g.getRoot('nnu_yhe', s.defdata),
         g.getRoot('nnu_yhe', s.defdata_all),
         g.getRoot('nnu_yhe', s.defdata_all + '_BAO'),
         ]


g.plot_2d(roots, param_pair=['nnu', 'yhe'], filled=True)

g.add_2d_contours(g.getRoot('nnu', s.defdata), 'nnu', 'yheused')

nnu = g.param_latex_label(roots[0], 'nnu')
# yhe = g.param_latex_label(roots[0], 'yheused')

bbn = s.defplanck + ' (' + s.LCDM + '+' + nnu + '+BBN'
labels = [s.defplanck, s.planckall, '+BAO' ]

yhe = 0.2465 * 0.2499 / 0.24732
sigma = 0.0097
g.add_y_bands(yhe, sigma, color='darkblue', alpha1=0.25)


yhe = 0.2551 * 0.2499 / 0.24732
sigma = 0.0022
g.add_y_bands(yhe, sigma, color='darkgreen', alpha1=0.5, alpha2=0.2)


g.add_x_marker(3.046)
g.add_legend(labels, legend_loc='upper right', colored_text=True);

g.add_text_left(s.LCDM + '+' + nnu + '+BBN', y=0.15, color='black')
g.add_text_left('Izotov et al. (2014)', y=0.1, color='darkgreen')
g.add_text_left('Aver et al. (2012)', y=0.05, color='darkblue')

ylim([0.15, 0.359])
xlim([1.5, 4.99])
g.export()
