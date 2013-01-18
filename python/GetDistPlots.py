import os, paramNames
import matplotlib
matplotlib.use('Agg')
from pylab import *

class GetDistPlotSettings:

    def __init__(self, subplot_size_inch=2):

        self.plot_meanlikes = True
        self.shade_meanlikes = True
        self.prob_label = None
        # self.prob_label = 'Probability'
        if self.plot_meanlikes:
            self.lineM = ['-k', '-r', '-b', '-g', '-m', '-y']
        else:
            self.lineM = ['-k', '--r', '-.b', ':g', '--m', '-.y']
        self.plot_args = None
        self.lineL = [':k', ':r', ':b', ':m', ':g', ':y']
        self.solid_colors = ['#009966', '#66CC99', '#336600', '#006633' , 'g', 'm', 'r']
        self.printex = 'pdf'
        self.line_labels = True
        self.x_label_rotation = 0
        self.num_shades = 80
        self.setWithSubplotSize(subplot_size_inch)
        self.progress = False
        self.tight_layout = True
        self.no_triangle_axis_labels = True
        self.solid_contour_palefactor = 0.5
# see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
        self.colormap = cm.Blues
        self.colormap_scatter = cm.jet
        self.param_names_for_labels = None  # 'clik_latex.paramnames'

    def setWithSubplotSize(self, size_inch):
        self.subplot_size_inch = size_inch
        self.lab_fontsize = 3 + 3 * self.subplot_size_inch
        self.axes_fontsize = 2 + 2 * self.subplot_size_inch
        self.lw1 = self.subplot_size_inch / 3.0
        self.lw_contour = self.lw1 * 0.4
        self.lw_likes = self.subplot_size_inch / 6.0
        self.scatter_size = 1 + self.subplot_size_inch


class Density1D():
    def bounds(self): return min(self.x), max(self.x)

class Density2D():
    def bounds(self):
        return (self.x1.min(), self.x2.min()), (self.x1.max(), self.x2.max())

class SampleAnalysisGetDist():

    def __init__(self, plot_data):
        self.plot_data = plot_data
        self.newPlot()

    def newPlot(self):
        self.single_samples = dict()

    def get_density_grid(self, root, param1, param2, conts=True, likes=False):
        result = Density2D()
        if likes:  result.pts = self.load_2d(root, param1, param2, '_likes')
        else: result.pts = self.load_2d(root, param1, param2)
        if result.pts is None: return None
        result.x1 = self.load_1d(root, param1)[:, 0]
        result.x2 = self.load_1d(root, param2)[:, 0]
        if (conts): result.contours = self.load_2d(root, param1, param2, '_cont')
        return result

    def get_density(self, root, param, likes=False):
        result = Density1D()
        pts = self.load_1d(root, param)
        if pts is None: return None
        result.x = pts[:, 0]
        result.pts = pts[:, 1]
        if (likes): result.likes = self.load_1d(root, param, '.likes')[:, 1]
        return result

    def load_single_samples(self, root):
        if not root in self.single_samples: self.single_samples[root] = loadtxt(self.plot_data_file(root) + '_single.txt')
        return self.single_samples[root]

    def paramsForRoot(self, root, labelParams=None):
        names = paramNames.paramNames(self.plot_data_file(root) + '.paramnames')
        if labelParams is not None: names.setLabelsAndDerivedFromParamNames(labelParams)
        return names

    def plot_data_file(self, root):
        return self.plot_data + os.sep + root

    def plot_data_file_1D(self, root, name):
        return self.plot_data_file(root) + '_p_' + name

    def plot_data_file_2D(self, root, name1, name2):
        fname = self.plot_data_file(root) + '_2D_' + name2 + '_' + name1
        if not os.path.exists(fname):
            return self.plot_data_file(root) + '_2D_' + name1 + '_' + name2, True
        else: return fname, False

    def load_1d(self, root, param, ext='.dat'):
        fname = self.plot_data_file_1D(root, param.name) + ext
        if not hasattr(param, 'plot_data'): param.plot_data = dict()
        if not fname in param.plot_data:
            if not os.path.exists(fname): param.plot_data[fname] = None
            else: param.plot_data[fname] = loadtxt(fname)
        return param.plot_data[fname]

    def load_2d(self, root, param1, param2, ext=''):
        fname, transpose = self.plot_data_file_2D(root, param1.name, param2.name)
        if not os.path.exists(fname + ext): return None
        pts = loadtxt(fname + ext)
        if transpose: pts = pts.transpose()
        return pts



class GetDistPlotter():

    def __init__(self, plot_data, settings=None,):
        if settings is None: self.settings = GetDistPlotSettings()
        else: self.settings = settings
        self.plot_data = plot_data
        self.sampleAnalyser = SampleAnalysisGetDist(plot_data)
        self.newPlot()

    def newPlot(self):
        clf()
        self.extra_artists = []
        self.param_name_sets = dict()
        self.sampleAnalyser.newPlot()

    def get_plot_args(self, plotno, **kwargs):
        if not self.settings.plot_args is None and len(self.settings.plot_args) > plotno:
            args = self.settings.plot_args[plotno]  #
            if args is None: args = dict()
        else: args = dict()
        args.update(kwargs)
        return args

    def get_color(self, plotno, **kwargs):
        args = self.get_plot_args(plotno, **kwargs)
        return args.get('color', self.settings.lineM[plotno][-1])

    def paramNamesForRoot(self, root):
        if not root in self.param_name_sets: self.param_name_sets[root] = self.sampleAnalyser.paramsForRoot(root, labelParams=self.settings.param_names_for_labels)
        return self.param_name_sets[root]

    def add_1d(self, root, param, plotno=0, **kwargs):
        param = self.check_param(roots, param)
        density = self.sampleAnalyser.get_density(root, param, likes=self.settings.plot_meanlikes)
        if density is None: return None;

        args = self.get_plot_args(plotno)
        args.update(kwargs)
        kwargs = args
        plot(density.x, density.pts, self.settings.lineM[plotno], linewidth=self.settings.lw1, **kwargs)
        if self.settings.plot_meanlikes:
            plot(density.x, density.likes, self.settings.lineL[plotno], linewidth=self.settings.lw_likes, **kwargs)

        return density.bounds()

    def add_2d_contours(self, root, param1, param2, plotno=0, filled=False, color=None, cols=None, alpha=0.5, **kwargs):
        param1, param2 = self.get_param_array(root, [param1, param2])

        density = self.sampleAnalyser.get_density_grid(root, param1, param2, conts=True, likes=False)
        if density is None: return None

        if filled:
            if cols is None:
                if color is None: color = self.settings.solid_colors[plotno]
                cols = [matplotlib.colors.colorConverter.to_rgb(color)]
                for _ in range(1, len(density.contours)):
                    cols = [[c * (1 - self.settings.solid_contour_palefactor) + self.settings.solid_contour_palefactor for c in cols[0]]] + cols
            contourf(density.x1, density.x2, density.pts, sorted(np.append([density.pts.max() + 1], density.contours)), colors=cols, alpha=alpha)
        else:
            if color is None: color = self.get_color(plotno, **kwargs)
            cols = [color]
        kwargs = self.get_plot_args(plotno, **kwargs)
        contour(density.x1, density.x2, density.pts, density.contours, colors=cols , linewidth=self.settings.lw_contour, alpha=alpha, **kwargs)

        return density.bounds()

    def add_2d_shading(self, root, param1, param2):
        param1, param2 = self.get_param_array(root, [param1, param2])
        density = self.sampleAnalyser.get_density_grid(root, param1, param2, conts=False, likes=self.settings.shade_meanlikes)

#        pcolor(x1, x2, shade_dat, cmap=self.settings.colormap, vmax=shade_dat.max(), vmin=shade_dat.min())
        contourf(density.x1, density.x2, density.pts, self.settings.num_shades, cmap=self.settings.colormap)
# doing contour gets rid of annoying wehite lines
        contour(density.x1, density.x2, density.pts, self.settings.num_shades, cmap=self.settings.colormap)
#        contour(cs, hold='on')

    def plot_2d(self, roots, param_pair, shaded=True, filled=False, **ax_args):
        param_pair = self.get_param_array(roots[0], param_pair)
        if self.settings.progress: print 'plotting: ', [param.name for param in param_pair]
        if shaded and not filled: self.add_2d_shading(roots[0], param_pair[0], param_pair[1])
        for i, root in enumerate(roots):
            res = self.add_2d_contours(root, param_pair[0], param_pair[1], i, filled=filled)
            if res is None: continue
            if i == 0: mins, maxs = res
            elif not shaded:
                mins = min(res[0], mins)
                maxs = max(res[1], maxs)
        if not 'lims' in ax_args: ax_args['lims'] = [mins[0], maxs[0], mins[1], maxs[1]]
        self.setAxes(param_pair, **ax_args)

    def add_1d_marker(self, marker, marker_color='k'):
        gca().add_line(Line2D([marker, marker], [0, 2], color=marker_color))

    def setAxisProperties(self, axis, x):
        formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        axis.set_major_formatter(formatter)
        tick_params(axis='both', which='major', labelsize=self.settings.axes_fontsize)
        if x and self.settings.x_label_rotation != 0:setp(xticks()[1], rotation=self.settings.x_label_rotation)
        else:
            xmin, xmax = xlim()
            if (abs(xmax - xmin) < 0.01 or max(abs(xmin), abs(xmax)) >= 1000):
                axis.set_major_locator(MaxNLocator(self.settings.subplot_size_inch / 2 + 3))
            else:
                axis.set_major_locator(MaxNLocator(4 + self.settings.subplot_size_inch / 2))

    def setAxes(self, params, lims=None, do_xlabel=True, do_ylabel=True, no_label_no_numbers=False):
        if lims is not None: axis(lims)
        self.setAxisProperties(gca().xaxis, True)
        if do_xlabel:self.set_xlabel(params[0])
        elif no_label_no_numbers: gca().set_xticklabels([])
        if len(params) > 1:
            self.setAxisProperties(gca().yaxis, False)
            if do_ylabel:self.set_ylabel(params[1])
            elif no_label_no_numbers: gca().set_yticklabels([])

    def set_xlabel(self, param):
        xlabel(r'$' + param.label + '$', fontsize=self.settings.lab_fontsize)

    def set_ylabel(self, param):
        ylabel(r'$' + param.label + '$', fontsize=self.settings.lab_fontsize)

    def plot_1d(self, roots, param, marker=None, marker_color='k', **ax_args):
        param = self.check_param(roots[0], param)
        xmin, xmax = self.add_1d(roots[0], param, 0)
        for i, root in enumerate(roots[1:]):
            bounds = self.add_1d(root, param, i + 1)
            if bounds is not None:
                xmin = min(xmin, bounds[0])
                xmax = max(xmax, bounds[1])

        if marker is not None: self.add_1d_marker(marker, marker_color)
        if not 'lims' in ax_args:ax_args['lims'] = [xmin, xmax, 0, 1.1]
        self.setAxes([param], **ax_args)

        if self.settings.prob_label is not None: ylabel(self.settings.prob_label)
        gca().set_yticks([])

    def make_figure(self, nplot=1, nx=None, ny=None, xstretch=1):
        self.newPlot()
        if nx is None: plot_col = round(sqrt(nplot / 1.4))
        else: plot_col = nx
        if ny is None: plot_row = (nplot + plot_col - 1) / plot_col
        else: plot_row = ny
        self.fig = figure(figsize=(self.settings.subplot_size_inch * plot_col * xstretch, self.settings.subplot_size_inch * plot_row))
        return plot_col, plot_row

    def get_param_array(self, root, in_params):
        if in_params is None or len(in_params) == 0: return self.paramNamesForRoot(root).names
        else:
            if not isinstance(in_params[0], paramNames.paramInfo): return self.paramNamesForRoot(root).parsWithNames(in_params)
        return in_params

    def check_param(self, root, param):
        if not isinstance(param, paramNames.paramInfo): return self.paramNamesForRoot(root).parWithName(param)
        return param

    def finish_plot(self, legend_labels=[], legend_loc='upper center', line_offset=0, no_gap=False):
        if self.settings.tight_layout:
            if no_gap: tight_layout(h_pad=0, w_pad=0)
            else: tight_layout()

        if self.settings.line_labels and len(legend_labels) > 1:
            lines = []
            for i in enumerate(legend_labels):
                color = self.get_color(i[0] + line_offset)
                col = self.settings.lineM[i[0] + line_offset]

                lines.append(Line2D([0, 1], [0, 1], color=color, ls=col[0:-1]))
            self.extra_artists = [self.fig.legend(lines, legend_labels, legend_loc, prop={'size':self.settings.lab_fontsize})]


    def plots_1d(self, roots, params=None, legend_labels=None, nx=None):
        params = self.get_param_array(roots[0], params)

        nparam = len(params)
        plot_col, plot_row = self.make_figure(nparam, nx=nx)

        for i, param in enumerate(params):
            subplot(plot_row, plot_col, i + 1)
            self.plot_1d(roots, param)

        self.finish_plot([legend_labels, roots][legend_labels is None])

        return plot_col, plot_row

    def plots_2d(self, roots, param1=None, params2=None, param_pairs=None, nx=None, legend_labels=None, filled=False):
        pairs = []
        if param_pairs is None:
            if param1 is not None:
                param1 = self.check_param(roots[0], param1)
                params2 = self.get_param_array(roots[0], params2)
                for param in params2:
                    if param.name != param1.name: pairs.append((param1, param))
            else: raise Exception('No parameter or parameter pairs for 2D plot')
        else:
            for pair in param_pairs:
                pairs.append((self.check_param(roots[0], pair[0]), self.check_param(roots[0], pair[1])))

        plot_col, plot_row = self.make_figure(len(pairs), nx=nx)

        for i, pair in enumerate(pairs):
            subplot(plot_row, plot_col, i + 1)
            self.plot_2d(roots, param_pair=pair, filled=filled)

        self.finish_plot([legend_labels, roots][legend_labels is None])

        return plot_col, plot_row

    def triangle_plot(self, roots, in_params=None, legend_labels=None, plot_3d_with_param=None, filled_compare=False, shaded=True):
        params = self.get_param_array(roots[0], in_params)
        plot_col = len(params)
        if plot_3d_with_param is not None: col_param = self.check_param(roots[0], plot_3d_with_param)
        self.make_figure(nx=plot_col, ny=plot_col)
        lims = dict()
        for i, param in enumerate(params):
            subplot(plot_col, plot_col, i * plot_col + i + 1)
            self.plot_1d(roots, param, do_xlabel=i == plot_col - 1, no_label_no_numbers=self.settings.no_triangle_axis_labels)
            if not shaded:
                xmin, xmax = xlim()
                r = xmax - xmin
                xlim([xmin - r / 10, xmax + r / 10])
            lims[i] = xlim()
        for i, param in enumerate(params):
            for i2 in range(i + 1, len(params)):
                param2 = params[i2]
                subplot(plot_col, plot_col, i2 * plot_col + i + 1)
                if plot_3d_with_param is not None:
                    self.plot_3d(roots, [param, param2, col_param], color_bar=False, line_offset=1,
                      do_xlabel=i2 == plot_col - 1, do_ylabel=i == 0, filled=filled_compare, no_label_no_numbers=self.settings.no_triangle_axis_labels)
                else:
                    self.plot_2d(roots, param_pair=[param, param2], do_xlabel=i2 == plot_col - 1, do_ylabel=i == 0,
                                        no_label_no_numbers=self.settings.no_triangle_axis_labels, shaded=shaded)
                xlim(lims[i])
                ylim(lims[i2])

        if self.settings.no_triangle_axis_labels:subplots_adjust(wspace=0, hspace=0)
        if plot_3d_with_param is not None:
            cb = self.fig.colorbar(self.last_scatter, cax=self.fig.add_axes([0.9, 0.5, 0.03, 0.35]))
            self.add_colorbar_label(cb, col_param)

        self.finish_plot([legend_labels, roots][legend_labels is None], legend_loc=None, no_gap=self.settings.no_triangle_axis_labels)

    def add_colorbar(self, param):
        cb = colorbar()
        self.add_colorbar_label(cb, param)
        return cb

    def add_line(self, P1, P2, zorder=0, color='k', **kwargs):
            gca().add_line(Line2D(P1, P2, color=color, zorder=zorder, **kwargs))

    def add_colorbar_label(self, cb, param):
        cb.set_label(r'$' + param.label + '$', fontsize=self.settings.lab_fontsize)
        setp(getp(cb.ax, 'ymajorticklabels'), fontsize=self.settings.axes_fontsize)

    def add_3d_scatter(self, root, in_params, color_bar=True):
        params = self.get_param_array(root, in_params)
        pts = self.sampleAnalyser.load_single_samples(root)
        names = self.paramNamesForRoot(root)
        cols = []
        for param in params:
            cols.append(names.numberOfName(param.name) + 2)
        self.last_scatter = scatter(pts[:, cols[0]], pts[:, cols[1]], edgecolors='none',
                s=self.settings.scatter_size, c=pts[:, cols[2]], cmap=self.settings.colormap_scatter)
        if color_bar: self.add_colorbar(params[2])

    def plot_3d(self, roots, in_params, color_bar=True, line_offset=0, filled=False, **ax_args):
        params = self.get_param_array(roots[0], in_params)
        self.add_3d_scatter(roots[0], params, color_bar=color_bar)
        for i, root in enumerate(roots[1:]):
            self.add_2d_contours(root, params[0], params[1], i + line_offset, filled=filled)
        self.setAxes(params, **ax_args)


    def plots_3d(self, roots, param_sets, nx=None, filled_compare=False, legend_labels=None):
        sets = [[self.check_param(roots[0], param) for param in param_group] for param_group in param_sets]
        plot_col, plot_row = self.make_figure(len(sets), nx=nx, xstretch=1.3)

        for i, triplet in enumerate(sets):
            subplot(plot_row, plot_col, i + 1)
            self.plot_3d(roots, triplet, filled=filled_compare)

        self.finish_plot([legend_labels, roots[1:]][legend_labels is None])
        return plot_col, plot_row

    def export(self, fname):
        savefig(fname, bbox_extra_artists=self.extra_artists, bbox_inches='tight')



def sample_plots():
    g = GetDistPlotter('main/plot_data')
    g.settings.setWithSubplotSize(3)
    g.settings.param_names_for_labels = 'clik_latex.paramnames'

    roots = ['base_omegak_planck_CAMspec_lowl_lowLike', 'base_omegak_planck_CAMspec_lowl_lowLike_post_lensing', 'base_omegak_planck_CAMspec_lowl_lowLike_BAO']
    params = g.get_param_array(roots[0], ['omegam', 'omegal', 'H0'])
    g.add_3d_scatter(roots[0], params)
    g.add_line([1, 0], [0, 1])
    g.add_2d_contours(roots[1], params[0], params[1], filled=False, color='g')
    g.add_2d_contours(roots[2], params[0], params[1], filled=True, color='#CC1100')
    g.setAxes(params, lims=[0, 1, 0, 1])
    g.export('omegam-omegal-H0_planck.pdf')

    g.triangle_plot(roots, ['omegabh2', 'omegach2', 'ns', 'omegal', 'omegak', ], plot_3d_with_param='H0', filled_compare=True,
                    legend_labels=['Planck', 'Planck+BAO'])
    g.export(roots[0] + '.pdf')

    roots = ['base_omegak_planck_CAMspec_lowl_lowLike', 'base_omegak_planck_CAMspec_lowl_lowLike_post_lensing', 'base_omegak_planck_CAMspec_lowl_lowLike_BAO', 'base_planck_CAMspec_lowl_lowLike']
    g.plots_1d(roots)
    g.export(roots[0] + '_1d.pdf')

    roots = ['base_omegak_planck_CAMspec_lowl_lowLike', 'base_omegak_planck_CAMspec_lowl_lowLike_BAO']
    g.plots_2d(roots, param_pairs=[('omegabh2', 'ns'), ('logA', 'tau')], nx=2)
    g.export(roots[0] + '_2d.pdf')

    roots = ['base_omegak_planck_CAMspec_lowl_lowLike', 'base_omegak_planck_CAMspec_lowl_lowLike_post_lensing', 'base_omegak_planck_CAMspec_lowl_lowLike_BAO', 'base_planck_CAMspec_lowl_lowLike']
    g.triangle_plot(roots, ['omegabh2', 'omegach2', 'ns', 'omegal', 'omegak', ], plot_3d_with_param='H0',
                    legend_labels=['Planck', 'Planck+lensing', 'Planck+BAO', 'Planck (flat)'])
    g.export(roots[0] + '_tri.pdf')

    roots = ['base_planck_CAMspec_lowl_lowLike', 'base_planck_CAMspec_lowl_lowLike_post_BAO']
    g.plots_3d(roots, [['omegabh2', 'omegach2', 'ns'], ['omegach2', 'logA', 'tau']])
    g.export(roots[0] + '_3d.pdf')


# sample_plots()
