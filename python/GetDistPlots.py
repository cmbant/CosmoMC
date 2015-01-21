import os
import sys
import iniFile
import matplotlib
import copy as cp
matplotlib.use('Agg')
from pylab import *
import batchJob
from getdist import MCSamples, paramNames

"""Plotting scripts for GetDist outputs"""


class GetDistPlotSettings(object):
    """Default sizes, font, styles etc settings for use by plots"""


    def __init__(self, subplot_size_inch=2, fig_width_inch=None):
        # if fig_width_inch set, forces fixed size, subplot_size_inch then just determines font sizes etc
        # otherwise width as wide as neccessary to show all subplots at specified size

        self.plot_meanlikes = False
        self.shade_meanlikes = False
        self.prob_label = None
        self.norm_prob_label = 'P'
        self.prob_y_ticks = False
        # self.prob_label = 'Probability'
        self.lineM = ['-k', '-r', '-b', '-g', '-m', '-c', '-y']
        self.plot_args = None
        self.solid_colors = ['#006FED', '#E03424', 'gray', '#009966', '#000866', '#336600', '#006633' , 'm', 'r']  # '#66CC99'
        self.default_dash_styles = {'--':(3, 2), '-.':(4, 1, 1, 1)}
        self.line_labels = True
        self.x_label_rotation = 0
        self.num_shades = 80
        self.fig_width_inch = fig_width_inch  # if you want to force specific fixed width
        self.progress = False
        self.tight_layout = True
        self.no_triangle_axis_labels = True
# see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
        self.colormap = cm.Blues
        self.colormap_scatter = cm.jet
        self.colorbar_rotation = None  # e.g. -90
        self.colorbar_label_pad = 0
        self.colorbar_label_rotation = -90  # seems to cause problems with some versions, can set to zero

        self.setWithSubplotSize(subplot_size_inch)

        self.param_names_for_labels = 'clik_latex.paramnames'
        self.tick_prune = None  # 'lower' or 'upper'
        self.tight_gap_fraction = 0.13  # space between ticks and the edge

        self.legend_loc = 'best'
        self.figure_legend_loc = 'upper center'
        self.legend_frame = True
        self.figure_legend_frame = True
        self.figure_legend_ncol = 1

        self.legend_rect_border = False
        self.legend_frac_subplot_margin = 0.2
        self.legend_frac_subplot_line = 0.1
        self.legend_fontsize = None

        self.num_contours = 2
        self.solid_contour_palefactor = 0.6
        self.alpha_filled_add = 0.85
        self.alpha_factor_contour_lines = 0.5

        self.axis_marker_color = 'gray'
        self.axis_marker_ls = '--'
        self.axis_marker_lw = 0.5

    def setWithSubplotSize(self, size_inch=3.5, size_mm=None):
        if size_mm: size_inch = size_mm * 0.0393700787
        self.subplot_size_inch = size_inch
        self.lab_fontsize = 7 + 2 * self.subplot_size_inch
        self.axes_fontsize = 4 + 2 * self.subplot_size_inch
        self.legend_fontsize = self.axes_fontsize
        self.font_size = self.lab_fontsize
        self.lw1 = self.subplot_size_inch / 3.0
        self.lw_contour = self.lw1 * 0.6
        self.lw_likes = self.subplot_size_inch / 6.0
        self.scatter_size = 3
        if size_inch > 4: self.scatter_size = size_inch * 2
        self.colorbar_axes_fontsize = self.axes_fontsize
        if  self.colorbar_label_rotation:  self.colorbar_label_pad = size_inch * 3

    def rcSizes(self, axes_fontsize=None, lab_fontsize=None, legend_fontsize=None):
        self.font_size = rcParams['font.size']
        self.legend_fontsize = legend_fontsize or rcParams['legend.fontsize']
        self.lab_fontsize = lab_fontsize  or rcParams['axes.labelsize']
        self.axes_fontsize = axes_fontsize or rcParams['xtick.labelsize']
        if isinstance(self.axes_fontsize, int): self.colorbar_axes_fontsize = self.axes_fontsize - 1
        else: self.colorbar_axes_fontsize = 'smaller'


defaultSettings = GetDistPlotSettings()

def getSinglePlotter(ratio=3 / 4., plot_data=None, chain_dir=None, width_inch=3.464, settings=None, mcsamples=False, analysis_settings=None):
    """ 
    Wrapper functions to get plotter to make single plot of fixed width (default: one page column width)
    """
    plotter = GetDistPlotter(plot_data=plot_data, chain_dir=chain_dir, settings=settings or GetDistPlotSettings(), mcsamples=mcsamples, analysis_settings=analysis_settings)
    plotter.settings.setWithSubplotSize(width_inch)
    plotter.settings.fig_width_inch = width_inch
#    if settings is None: plotter.settings.rcSizes()
    plotter.make_figure(1, xstretch=1 / ratio)
    return plotter

def getSubplotPlotter(plot_data=None, chain_dir=None, subplot_size=2, settings=None, mcsamples=False, analysis_settings=None, width_inch=None):
    """ 
    Wrapper functions to get plotter to make array of subplots 
    if width_inch is None, just makes plot as big as needed
    """
    plotter = GetDistPlotter(plot_data=plot_data, chain_dir=chain_dir, settings=settings or GetDistPlotSettings(), mcsamples=mcsamples, analysis_settings=analysis_settings)
    plotter.settings.setWithSubplotSize(subplot_size)
    if width_inch:
        plotter.settings.fig_width_inch = width_inch
        if not settings: plotter.settings.rcSizes()
    if subplot_size < 3 and settings is None:
        plotter.settings.axes_fontsize += 2
        plotter.settings.colorbar_axes_fontsize += 2
        plotter.settings.legend_fontsize = plotter.settings.lab_fontsize + 1
    return plotter


class Density1D(object):
    def bounds(self): return min(self.x), max(self.x)

class Density2D(object):
    def bounds(self):
        return (self.x1.min(), self.x2.min()), (self.x1.max(), self.x2.max())

    def xy_bounds(self):
        return (self.x1.min(), self.x1.max()), (self.x2.min(), self.x2.max())

class paramBounds(object):
    def __init__(self, fileName):
        self.upper = dict()
        self.lower = dict()
        if os.path.exists(fileName):
            textFileHandle = open(fileName)
            for line in textFileHandle:
                name, low, high = [text.strip() for text in  line.split()]
                if low != 'N': self.lower[name] = float(low)
                if high != 'N': self.upper[name] = float(high)
            textFileHandle.close()

    def getUpper(self, name):
        if name in self.upper: return self.upper[name]
        else: return None

    def getLower(self, name):
        if name in self.lower: return self.lower[name]
        else: return None


class SampleAnalysisGetDist(object):

    def __init__(self, plot_data):
        self.plot_data = plot_data
        self.newPlot()
        self.paths = dict()

    def newPlot(self):
        self.single_samples = dict()

    def get_density_grid(self, root, param1, param2, conts=2, likes=False):
        if likes:  res = self.load_2d(root, param1, param2, '_likes')
        else: res = self.load_2d(root, param1, param2)
        if res is None: return None
        result = Density2D()
        (result.pts, result.x1, result.x2) = res
        if conts > 0: result.contours = self.load_2d(root, param1, param2, '_cont', no_axes=True)[0:conts]
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
        if not root in self.single_samples: self.single_samples[root] = loadtxt(self.plot_data_file(root) + '_single.txt')[:, 2:]
        return self.single_samples[root]

    def paramsForRoot(self, root, labelParams=None):
        names = paramNames.paramNames(self.plot_data_file(root) + '.paramnames')
        if labelParams is not None: names.setLabelsAndDerivedFromParamNames(labelParams)
        return names

    def boundsForRoot(self, root):
        return paramBounds(self.plot_data_file(root) + '.bounds')

    def plot_data_file(self, root):
        # find first match to roots that exist in list of plot_data paths
        if os.sep in root: return root
        path = self.paths.get(root, None)
        if path is not None: return path
        for datadir in self.plot_data:
            path = datadir + os.sep + root
            if os.path.exists(path + '.paramnames'):
                self.paths[root] = path
                return path
        return self.plot_data[0] + os.sep + root

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

    def load_2d(self, root, param1, param2, ext='', no_axes=False):
        fname, transpose = self.plot_data_file_2D(root, param1.name, param2.name)
        if not os.path.exists(fname + ext): return None
        pts = loadtxt(fname + ext)
        if transpose: pts = pts.transpose()
        if no_axes: return pts
        x = loadtxt(fname + '_x')
        y = loadtxt(fname + '_y')
        if transpose: return (pts, y, x)
        else: return (pts, x, y)


class MCSampleAnalysis(object):

    def __init__(self, chain_dir=None, settings=None):
        self.chain_dir = chain_dir
        self.batch = None
        ini = None
        if chain_dir and isinstance(chain_dir, basestring):
            import makeGrid
            if makeGrid.pathIsGrid(chain_dir):
                self.batch = batchJob.readobject(chain_dir)
                ini = iniFile.iniFile(self.batch.commonPath + 'getdist_common.ini')
            else:
                self.chain_dir = [chain_dir]
        if settings and isinstance(settings, basestring):
            ini = iniFile.iniFile(settings)

        if not ini: ini = iniFile.iniFile()
        if isinstance(settings, dict): ini.params.update(settings)
        self.ini = ini
        self.reset()

    def reset(self):
        self.mcsamples = {}
        # Dicts. 1st key is root; 2nd key is param
        self.densities_1D = dict()
        self.densities_2D = dict()
        self.single_samples = dict()

    def samplesForRoot(self, root, file_root=None, cache=True):
        if os.path.isabs(root):
            root = os.path.basename(root)
        if root in self.mcsamples and cache: return self.mcsamples[root]
        jobItem = None
        dist_setings = {}
        if not file_root:
            if self.batch:
                jobItem = self.batch.resolveRoot(root)
                if not jobItem: raise Exception('chain not found: ' + root)
                file_root = jobItem.chainRoot
                dist_setings = jobItem.dist_settings
            else:
                for directory in self.chain_dir:
                    name = os.path.join(directory, root)
                    if os.path.exists(name + '_1.txt'):
                        file_root = name
                        break
                if not file_root:  raise Exception('chain not found: ' + root)
        self.mcsamples[root] = MCSamples.loadMCSamples(file_root, self.ini, jobItem, dist_settings=dist_setings)
        return self.mcsamples[root]

    def addRootGrid(self, root):
        if self.batch is None: return None
        return self.samplesForRoot(root)

    def addRoots(self, roots):
        for root in roots:
            self.addRoot(root)

    def addRoot(self, file_root):
        return self.samplesForRoot(os.path.basename(file_root), file_root)

    def removeRoot(self, file_root):
        root = os.path.basename(file_root)
        print "remove root for %s" % root
        if self.mcsamples.has_key(root):
            del self.mcsamples[root]

    def newPlot(self):
        pass

    def get_1d(self, root, param, likes=False):
        rootdata = self.densities_1D.get(root)
        if rootdata is None:
            rootdata = {}
            self.densities_1D[root] = rootdata

        name = param.name
        samples = self.samplesForRoot(root)
        density = rootdata.get(name)
        if density is None:
            index = samples.index.get(name)
            if index is None: return None
            density = samples.Get1DDensity(index)

        if density is None: return None
        dat, likedata = density
        if likes:
            return likedata
        else:
            return dat

    def get_density_grid(self, root, param1, param2, conts=2, likes=False):
        rootdata = self.densities_2D.get(root)
        if not rootdata:
            rootdata = {}
            self.densities_2D[root] = rootdata
        key = (param1.name, param2.name)
        density = rootdata.get(key)
        if not density:
            samples = self.samplesForRoot(root)
            index1 = samples.index.get(param1.name)
            index2 = samples.index.get(param2.name)
            if index1 is None or index2 is None: return None
            samples.initParamRanges(index1)
            samples.initParamRanges(index2)
            density = samples.Get2DPlotData(index2, index1)
            if density is None: return None
            rootdata[key] = density
        result = Density2D()
        dat, likes, cont, result.x1, result.x2 = density
        if likes:
            result.pts = likes
        else:
            result.pts = dat
        if conts > 0: result.contours = cont[0:conts]
        return result

    def get_density(self, root, param, likes=False):
        result = Density1D()
        pts = self.get_1d(root, param)
        if pts is None: return None
        result.x = pts[:, 0]
        result.pts = pts[:, 1]
        if (likes): result.likes = self.get_1d(root, param, True)[:, 1]
        return result

    def load_single_samples(self, root):
        if not root in self.single_samples:
            self.single_samples[root] = self.samplesForRoot(root).MakeSingleSamples()
        return self.single_samples[root]


    def paramsForRoot(self, root, labelParams=None):
        samples = self.samplesForRoot(root)
        names = samples.paramNames
        if labelParams is not None:
            names.setLabelsAndDerivedFromParamNames(os.path.join(batchJob.getCodeRootPath(), labelParams))
        return names

    def boundsForRoot(self, root):
        lower, upper = self.samplesForRoot(root).getBounds()
        bounds = paramBounds("")
        bounds.lower = lower
        bounds.upper = upper
        return bounds


class GetDistPlotter(object):

    def __init__(self, plot_data=None, settings=None, mcsamples=False, chain_dir=None, analysis_settings=None):
        """
        Set plot_data to directory name if you have pre-compputed plot_data/ directory from GetDist
        Set chain_dir to directly to use chains in the given directory (can also be a list of directories to search)
        """
        if settings is None: self.settings = defaultSettings
        else: self.settings = settings
        if isinstance(plot_data, basestring): self.plot_data = [plot_data]
        else: self.plot_data = plot_data
        if chain_dir is not None or mcsamples:
            self.sampleAnalyser = MCSampleAnalysis(chain_dir, analysis_settings)
        else:
            self.sampleAnalyser = SampleAnalysisGetDist(self.plot_data)
        self.newPlot()

    def newPlot(self):
        clf()
        self.extra_artists = []
        self.contours_added = []
        self.lines_added = dict()
        self.param_name_sets = dict()
        self.param_bounds_sets = dict()
        self.sampleAnalyser.newPlot()
        self.fig = None
        self.subplots = None

    def show_all_settings(self):
        print 'Python version:', sys.version
        print '\nMatplotlib version:', matplotlib.__version__
        print '\nGetDist Plot Settings:'
        sets = self.settings.__dict__
        for key, value in sets.items():
            print key, ':', value
        print '\nRC params:'
        for key, value in matplotlib.rcParams.items():
            print key, ':', value

    def get_plot_args(self, plotno, **kwargs):
        if not self.settings.plot_args is None and len(self.settings.plot_args) > plotno:
            args = self.settings.plot_args[plotno]  #
            if args is None: args = dict()
        else: args = dict()
        args.update(kwargs)
        return args

    def get_dashes_for_ls(self, ls):
        return self.settings.default_dash_styles.get(ls, None)

    def get_default_ls(self, plotno=0):
        try:
            return self.settings.lineM[plotno]
        except IndexError:
            print 'Error adding line ' + plotno + ': Add more default line stype entries to settings.lineM'
            raise

    def get_line_styles(self, plotno, **kwargs):
        args = self.get_plot_args(plotno, **kwargs)
        if not 'ls' in args: args['ls'] = self.get_default_ls(plotno)[:-1]
        if not 'dashes' in args:
            dashes = self.get_dashes_for_ls(args['ls'])
            if dashes is not None: args['dashes'] = dashes
        if not 'color' in args:
            args['color'] = self.get_default_ls(plotno)[-1]
        if not 'lw' in args: args['lw'] = self.settings.lw1
        return args

    def get_color(self, plotno, **kwargs):
        return self.get_line_styles(plotno, **kwargs)['color']

    def get_linestyle(self, plotno, **kwargs):
        return self.get_line_styles(plotno, **kwargs)['ls']

    def get_alpha2D(self, plotno, **kwargs):
        args = self.get_plot_args(plotno, **kwargs)
        if kwargs.get('filled') and plotno > 0: default = self.settings.alpha_filled_add
        else: default = 1
        return args.get('alpha', default)

    def paramNamesForRoot(self, root):
        if not root in self.param_name_sets: self.param_name_sets[root] = self.sampleAnalyser.paramsForRoot(root, labelParams=self.settings.param_names_for_labels)
        return self.param_name_sets[root]

    def paramBoundsForRoot(self, root):
        if not root in self.param_bounds_sets: self.param_bounds_sets[root] = self.sampleAnalyser.boundsForRoot(root)
        return self.param_bounds_sets[root]

    def checkBounds(self, root, name, xmin, xmax):
        d = self.paramBoundsForRoot(root)
        low = d.getLower(name)
        if low is not None: xmin = max(xmin, low)
        up = d.getUpper(name)
        if up is not None: xmax = min(xmax, up)
        return xmin, xmax

    def add_1d(self, root, param, plotno=0, normalized=False, **kwargs):
        param = self.check_param(root, param)
        density = self.sampleAnalyser.get_density(root, param, likes=self.settings.plot_meanlikes)
        if density is None: return None;
        pts = density.pts
        if normalized:
            norm = (pts[0] + pts[-1]) / 2 + sum(pts[1:-1])
            norm *= density.x[1] - density.x[0]
            pts /= norm

        kwargs = self.get_line_styles(plotno, **kwargs)
        self.lines_added[plotno] = kwargs
        l, = plot(density.x, pts, **kwargs)
        if kwargs.get('dashes'):
            l.set_dashes(kwargs['dashes'])
        if self.settings.plot_meanlikes:
            kwargs['lw'] = self.settings.lw_likes
            plot(density.x, density.likes, **kwargs)

        return density.bounds()

    def add_2d_contours(self, root, param1=None, param2=None, plotno=0, of=None, cols=None,
                            add_legend_proxy=True, param_pair=None, density=None, alpha=None, **kwargs):
        param1, param2 = self.get_param_array(root, param_pair or [param1, param2])

        if not density: density = self.sampleAnalyser.get_density_grid(root, param1, param2, conts=self.settings.num_contours, likes=False)
        if density is None:
            if add_legend_proxy: self.contours_added.append(None)
            return None
        if alpha is None: alpha = self.get_alpha2D(plotno, **kwargs)

        if add_legend_proxy:
            proxyIx = len(self.contours_added)
            self.contours_added.append(None)
        elif None in self.contours_added and self.contours_added.index(None) == plotno:
            proxyIx = plotno
        else: proxyIx = -1

        if kwargs.get('filled'):
            linestyles = ['-']
            if cols is None:
                color = kwargs.get('color')
                if color is None:
                    if of is not None:color = self.settings.solid_colors[of - plotno - 1]
                    else: color = self.settings.solid_colors[plotno]
                if isinstance(color, basestring):
                    cols = [matplotlib.colors.colorConverter.to_rgb(color)]
                    for _ in range(1, len(density.contours)):
                        cols = [[c * (1 - self.settings.solid_contour_palefactor) + self.settings.solid_contour_palefactor for c in cols[0]]] + cols
                else: cols = color
            levels = sorted(np.append([density.pts.max() + 1], density.contours))
            CS = contourf(density.x1, density.x2, density.pts, levels, colors=cols, alpha=alpha, **kwargs)
            if proxyIx >= 0: self.contours_added[proxyIx] = (Rectangle((0, 0), 1, 1, fc=CS.tcolors[1][0]))
            contour(density.x1, density.x2, density.pts, levels[:1], colors=CS.tcolors[1],
                    linewidths=self.settings.lw_contour, alpha=alpha * self.settings.alpha_factor_contour_lines, **kwargs)
        else:
            args = self.get_line_styles(plotno, **kwargs)
#            if color is None: color = self.get_color(plotno, **kwargs)
#            cols = [color]
#            if ls is None: ls = self.get_linestyle(plotno, **kwargs)
            linestyles = [args['ls']]
            cols = [args['color']]
            kwargs = self.get_plot_args(plotno, **kwargs)
            kwargs['alpha'] = alpha
            CS = contour(density.x1, density.x2, density.pts, density.contours, colors=cols , linestyles=linestyles, linewidths=self.settings.lw_contour, **kwargs)
            dashes = args.get('dashes')
            if dashes:
                for c in CS.collections:
                    c.set_dashes([(0, dashes)])
            if proxyIx >= 0:
                line = Line2D([0, 1], [0, 1], ls=linestyles[0], lw=self.settings.lw_contour, color=cols[0], alpha=args.get('alpha'))
                if dashes: line.set_dashes(dashes)
                self.contours_added[proxyIx] = line

        return density.xy_bounds()

    def add_2d_shading(self, root, param1, param2):
        param1, param2 = self.get_param_array(root, [param1, param2])
        density = self.sampleAnalyser.get_density_grid(root, param1, param2, conts=0, likes=self.settings.shade_meanlikes)
        if density is None: return

#        pcolor(x1, x2, shade_dat, cmap=self.settings.colormap, vmax=shade_dat.max(), vmin=shade_dat.min())
        contourf(density.x1, density.x2, density.pts, self.settings.num_shades, cmap=self.settings.colormap)
# doing contour gets rid of annoying wehite lines
        contour(density.x1, density.x2, density.pts, self.settings.num_shades, cmap=self.settings.colormap)
#        contour(cs, hold='on')

    def updateLimit(self, bounds, curbounds):
        if not bounds: return curbounds
        if curbounds is None or curbounds[0] is None: return bounds
        return min(curbounds[0], bounds[0]), max(curbounds[1], bounds[1])

    def updateLimits(self, res, xlims, ylims, doResize=True):
        if res is None: return xlims, ylims
        if xlims is None and ylims is None: return res
        if not doResize: return xlims, ylims
        else: return self.updateLimit(res[0], xlims), self.updateLimit(res[1], ylims)


    def _make_line_args(self, nroots, **kwargs):
        line_args = kwargs.get('line_args')
        if line_args is None: line_args = kwargs.get('contour_args')
        if line_args is None: line_args = [{}] * nroots
        elif isinstance(line_args, dict): line_args = [line_args] * nroots
        if len(line_args) < nroots: line_args += [{}] * (nroots - len(line_args))
        colors = kwargs.get('colors')
        lws = kwargs.get('lws')
        alphas = kwargs.get('alphas')
        ls = kwargs.get('ls')
        for i, args in enumerate(line_args):
            c = cp.copy(args)  # careful to copy before modifying any
            line_args[i] = c
            if colors and i < len(colors) and colors[i]:
                if isinstance(colors[i], basestring):
#                    c['color'] = matplotlib.colors.colorConverter.to_rgb(colors[i])
                    c['color'] = colors[i]
                else:
                    c['color'] = colors[i]
            if ls and i < len(ls) and ls[i]: c['ls'] = ls[i]
            if alphas and i < len(alphas) and alphas[i]: c['alpha'] = alphas[i]
            if lws and i < len(lws) and lws[i]: c['lw'] = lws[i]
        return line_args

    def _make_contour_args(self, nroots, **kwargs):
        contour_args = self._make_line_args(nroots, **kwargs)
        filled = kwargs.get('filled')
        if filled and not isinstance(filled, bool):
            for cont, fill in zip(contour_args, filled):
                cont['filled'] = fill
        for cont in contour_args:
            if cont.get('filled') is None: cont['filled'] = filled or False
        return contour_args

    def plot_2d(self, roots, param1=None, param2=None, param_pair=None, shaded=False, add_legend_proxy=True, **kwargs):
        if self.fig is None: self.make_figure()
        if isinstance(roots, basestring):roots = [roots]
        if isinstance(param1, list):
            param_pair = param1
            param1 = None
        param_pair = self.get_param_array(roots[0], param_pair or [param1, param2])
        if self.settings.progress: print 'plotting: ', [param.name for param in param_pair]
        if shaded and not kwargs.get('filled'): self.add_2d_shading(roots[0], param_pair[0], param_pair[1])
        xbounds, ybounds = None, None
        contour_args = self._make_contour_args(len(roots), **kwargs)
        for i, root in enumerate(roots):
            res = self.add_2d_contours(root, param_pair[0], param_pair[1], i, of=len(roots),
                                       add_legend_proxy=add_legend_proxy, **contour_args[i])
            xbounds, ybounds = self.updateLimits(res, xbounds, ybounds, doResize=not shaded)
        if xbounds is None: return
        if not 'lims' in kwargs:
            lim1 = self.checkBounds(roots[0], param_pair[0].name , xbounds[0], xbounds[1])
            lim2 = self.checkBounds(roots[0], param_pair[1].name , ybounds[0], ybounds[1])
            kwargs['lims'] = [lim1[0], lim1[1], lim2[0], lim2[1]]
        self.setAxes(param_pair, **kwargs)
        return xbounds, ybounds

    def add_1d_marker(self, marker, color=None, ls=None):
        self.add_x_marker(marker, color, ls)

    def add_x_marker(self, marker, color=None, ls=None, lw=None):
        if color is None: color = self.settings.axis_marker_color
        if ls is None: ls = self.settings.axis_marker_ls
        if lw is None: lw = self.settings.axis_marker_lw
        axvline(marker, ls=ls, color=color, lw=lw)

    def add_y_marker(self, marker, color=None, ls=None, lw=None):
        if color is None: color = self.settings.axis_marker_color
        if ls is None: ls = self.settings.axis_marker_ls
        if lw is None: lw = self.settings.axis_marker_lw
        axhline(marker, ls=ls, color=color, lw=lw)

    def add_y_bands(self, y, sigma, xlim=None, color='gray', ax=None, alpha1=0.15, alpha2=0.1):
        ax = ax or gca()
        if xlim is None: xlim = ax.xaxis.get_view_interval()
        one = array([1, 1])
        c = color
        if alpha2 > 0: ax.fill_between(xlim, one * (y - sigma * 2), one * (y + sigma * 2), facecolor=c, alpha=alpha2, edgecolor=c, lw=0)
        if alpha1 > 0: ax.fill_between(xlim, one * (y - sigma), one * (y + sigma), facecolor=c, alpha=alpha1, edgecolor=c, lw=0)

    def set_locator(self, axis, x=False, prune=None):
        if x: xmin, xmax = axis.get_view_interval()
        if (x and (abs(xmax - xmin) < 0.01 or max(abs(xmin), abs(xmax)) >= 1000)):
            axis.set_major_locator(MaxNLocator(self.settings.subplot_size_inch / 2 + 3, prune=prune))
        else:
            axis.set_major_locator(MaxNLocator(self.settings.subplot_size_inch / 2 + 4, prune=prune))


    def setAxisProperties(self, axis, x, prune):
        formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        axis.set_major_formatter(formatter)
        tick_params(axis='both', which='major', labelsize=self.settings.axes_fontsize)
        if x and self.settings.x_label_rotation != 0:setp(xticks()[1], rotation=self.settings.x_label_rotation)
        self.set_locator(axis, x, prune=prune)

    def setAxes(self, params=[], lims=None, do_xlabel=True, do_ylabel=True, no_label_no_numbers=False, pos=None, prune=None,
                color_label_in_axes=False, ax=None, **other_args):
        if lims is not None: axis(lims)
        if prune is None: prune = self.settings.tick_prune
        ax = ax or gca()
        self.setAxisProperties(ax.xaxis, True, prune)
        if pos is not None: ax.set_position(pos)  # # set [left, bottom, width, height] for the figure
        if do_xlabel and len(params) > 0:self.set_xlabel(params[0])
        elif no_label_no_numbers: ax.set_xticklabels([])
        if len(params) > 1:
            self.setAxisProperties(ax.yaxis, False, prune)
            if do_ylabel:self.set_ylabel(params[1])
            elif no_label_no_numbers: ax.set_yticklabels([])
        if color_label_in_axes and len(params) > 2: self.add_text(r'$' + params[2].label + '$')
        return ax

    def set_xlabel(self, param):
        xlabel(r'$' + param.label + '$', fontsize=self.settings.lab_fontsize,
                verticalalignment='baseline', labelpad=4 + self.settings.font_size)  # test_size because need a number not e.g. 'medium'

    def set_ylabel(self, param):
        ylabel(r'$' + param.label + '$', fontsize=self.settings.lab_fontsize)

    def plot_1d(self, roots, param, marker=None, marker_color=None, label_right=False,
                no_ylabel=False, no_ytick=False, no_zero=False, normalized=False, param_renames={}, **kwargs):
        if isinstance(roots, basestring): roots = [roots]
        if self.fig is None: self.make_figure()
        plotparam = None
        plotroot = None
        line_args = self._make_line_args(len(roots), **kwargs)
        xmin, xmax = None, None
        for i, root in enumerate(roots):
            root_param = self.check_param(root, param, param_renames)
            if not root_param: continue
            bounds = self.add_1d(root, root_param, i, normalized=normalized, **line_args[i])
            xmin, xmax = self.updateLimit(bounds, (xmin, xmax))
            if bounds is not None and not plotparam:
                    plotparam = root_param
                    plotroot = root
        if plotparam is None: raise Exception('No roots have parameter: ' + str(param))
        if marker is not None: self.add_x_marker(marker, marker_color)
        if not 'lims' in kwargs:
            xmin, xmax = self.checkBounds(plotroot, plotparam.name, xmin, xmax)
            if normalized: mx = gca().yaxis.get_view_interval()[-1]
            else: mx = 1.099
            kwargs['lims'] = [xmin, xmax, 0, mx]
        ax = self.setAxes([plotparam], **kwargs)

        if normalized: lab = self.settings.norm_prob_label
        else: lab = self.settings.prob_label
        if lab and not no_ylabel:
            if label_right:
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
                ylabel(lab)
            else: ylabel(lab)
        if no_ytick or not self.settings.prob_y_ticks: ax.set_yticks([])
        elif no_ylabel: ax.set_yticklabels([])
        elif no_zero and not normalized:
            ticks = ax.get_yticks()
            if ticks[-1] > 1: ticks = ticks[:-1]
            ax.set_yticks(ticks[1:])

    def make_figure(self, nplot=1, nx=None, ny=None, xstretch=1, ystretch=1):
        self.newPlot()
        if nx is None: self.plot_col = int(round(sqrt(nplot / 1.4)))
        else: self.plot_col = nx
        if ny is None: self.plot_row = (nplot + self.plot_col - 1) / self.plot_col
        else: self.plot_row = ny
        if self.settings.fig_width_inch is not None:
            self.fig = figure(figsize=(self.settings.fig_width_inch, (self.settings.fig_width_inch * self.plot_row * ystretch) / (self.plot_col * xstretch)))
        else:
            self.fig = figure(figsize=(self.settings.subplot_size_inch * self.plot_col * xstretch, self.settings.subplot_size_inch * self.plot_row * ystretch))
        self.subplots = np.ndarray((self.plot_row, self.plot_col), dtype=object)
        self.subplots[:, :] = None
        return self.plot_col, self.plot_row

    def get_param_array(self, root, in_params=None, renames={}):
        if in_params is None or len(in_params) == 0: return self.paramNamesForRoot(root).names
        else:
            if not all([isinstance(param, paramNames.paramInfo) for param in in_params]):
                return self.paramNamesForRoot(root).parsWithNames(in_params, error=True, renames=renames)
        return in_params

    def check_param(self, root, param, renames={}):
        if not isinstance(param, paramNames.paramInfo):
            return self.paramNamesForRoot(root).parWithName(param, error=True, renames=renames)
        elif renames: return self.paramNamesForRoot(root).parWithName(param.name, error=False, renames=renames)
        return param

    def param_latex_label(self, root, param, labelParams=None):
        if labelParams is not None:
            p = self.sampleAnalyser.paramsForRoot(root, labelParams=labelParams).parWithName(param)
        else:
            p = self.check_param(root, param)
        if not p: raise Exception('Parameter not found: ' + param)
        return r'$' + p.label + r'$'

    def add_legend(self, legend_labels, legend_loc=None, line_offset=0, legend_ncol=None, colored_text=False,
                   figure=False, ax=None, label_order=None, align_right=False, fontsize=None):
            if legend_loc is None:
                if figure: legend_loc = self.settings.figure_legend_loc
                else: legend_loc = self.settings.legend_loc
            if legend_ncol is None: legend_ncol = self.settings.figure_legend_ncol
            lines = []
            if len(self.contours_added) == 0:
                for i in enumerate(legend_labels):
                    args = self.lines_added.get(i[0]) or self.get_line_styles(i[0] + line_offset)
                    args.pop('filled', None)
                    lines.append(Line2D([0, 1], [0, 1], **args))
            else: lines = self.contours_added
            args = {'ncol':legend_ncol}
            if fontsize or self.settings.legend_fontsize: args['prop'] = {'size':fontsize or self.settings.legend_fontsize}
            if colored_text:
                args['handlelength'] = 0
                args['handletextpad'] = 0
            if label_order is not None:
                if label_order == '-1': label_order = range(len(lines)).reverse()
                lines = [lines[i] for i in label_order]
                legend_labels = [legend_labels[i] for i in label_order]
            if figure:
#                args['frameon'] = self.settings.figure_legend_frame
                self.legend = self.fig.legend(lines, legend_labels, legend_loc, **args)
                if not self.settings.figure_legend_frame:
                    # this works with tight_layout
                    self.legend.get_frame().set_edgecolor('none')
            else:
                args['frameon'] = self.settings.legend_frame and not colored_text
                self.legend = (ax or gca()).legend(lines, legend_labels, legend_loc, **args)
            if align_right:
                vp = self.legend ._legend_box._children[-1]._children[0]
                for c in vp._children:
                    c._children.reverse()
                vp.align = "right"
            if not self.settings.legend_rect_border:
                for rect in self.legend.get_patches():
                    rect.set_edgecolor(rect.get_facecolor())
            if colored_text:
                for h, text in zip(self.legend.legendHandles, self.legend.get_texts()):
                    h.set_visible(False)
                    if isinstance(h, Line2D):
                        c = h.get_color()
                    elif isinstance(h, matplotlib.patches.Patch):
                        c = h.get_facecolor()
                    else: continue
                    text.set_color(c)
            return self.legend

    def finish_plot(self, legend_labels=None, legend_loc=None, line_offset=0, legend_ncol=None, label_order=None, no_gap=False, no_extra_legend_space=False, no_tight=False):
        has_legend = self.settings.line_labels and legend_labels and len(legend_labels) > 1
        if self.settings.tight_layout and not no_tight:
            if no_gap: tight_layout(h_pad=0, w_pad=0)
            else: tight_layout()

        if has_legend:
            if legend_ncol is None: legend_ncol = self.settings.figure_legend_ncol
            if legend_loc is None: legend_loc = self.settings.figure_legend_loc
            self.extra_artists = [self.add_legend(legend_labels, legend_loc, line_offset, legend_ncol, label_order=label_order, figure=True)]
            if self.settings.tight_layout and not no_extra_legend_space:
                frac = self.settings.legend_frac_subplot_margin + (len(legend_labels) / legend_ncol) * self.settings.legend_frac_subplot_line
                if 'upper' in legend_loc: subplots_adjust(top=1 - frac / self.plot_row)
                elif 'lower' in legend_loc: subplots_adjust(bottom=frac / self.plot_row)

    def _escapeLatex(self, text):
        if matplotlib.rcParams['text.usetex']:
            return text.replace('_', '{\\textunderscore}')
        else:
            return text

    def default_legend_labels(self, legend_labels, roots):
        if legend_labels is None:
            return [self._escapeLatex(root) for root in roots]
        else: return legend_labels

    def plots_1d(self, roots, params=None, legend_labels=None, legend_ncol=None, label_order=None, nx=None,
                 paramList=None, roots_per_param=False, share_y=None, markers=None, xlims=None, param_renames={}):
        if roots_per_param:
            params = [self.check_param(root[0], param, param_renames) for root, param in zip(roots, params)]
        else: params = self.get_param_array(roots[0], params, param_renames)
        if paramList is not None:
            wantedParams = self.paramNameListFromFile(paramList)
            params = [param for param in params if param.name in wantedParams or param_renames.get(param.name, '') in wantedParams]
        nparam = len(params)
        if share_y is None: share_y = self.settings.prob_label is not None and nparam > 1
        plot_col, plot_row = self.make_figure(nparam, nx=nx)
        plot_roots = roots

        for i, param in enumerate(params):
            ax = self.subplot_number(i)
            if roots_per_param: plot_roots = roots[i]
            if markers is not None:
                if isinstance(markers, dict): marker = markers.get(param.name, None)
                elif i < len(markers): marker = markers[i]
            else: marker = None
            self.plot_1d(plot_roots, param, no_ylabel=share_y and  i % self.plot_col > 0, marker=marker, param_renames=param_renames)
            if xlims is not None: ax.set_xlim(xlims[i][0], xlims[i][1])
            if share_y: self.spaceTicks(ax.xaxis, expand=True)

        self.finish_plot(self.default_legend_labels(legend_labels, roots), legend_ncol=legend_ncol, label_order=label_order)
        if share_y: subplots_adjust(wspace=0)
        return plot_col, plot_row

    def plots_2d(self, roots, param1=None, params2=None, param_pairs=None, nx=None, legend_labels=None, legend_ncol=None,
                 label_order=None, filled=False):
        pairs = []
        if isinstance(roots, basestring): roots = [roots]
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
            self.subplot_number(i)
            self.plot_2d(roots, param_pair=pair, filled=filled, add_legend_proxy=i == 0)

        self.finish_plot(self.default_legend_labels(legend_labels, roots), legend_ncol=legend_ncol, label_order=label_order)

        return plot_col, plot_row

    def subplot(self, x, y, **kwargs):
        self.subplots[y, x] = ax = subplot(self.plot_row, self.plot_col, y * self.plot_col + x + 1, **kwargs)
        return ax

    def subplot_number(self, i, **kwargs):
        self.subplots[i / self.plot_col, i % self.plot_col] = ax = subplot(self.plot_row, self.plot_col, i + 1)
        return ax

    def plots_2d_triplets(self, root_params_triplets, nx=None, filled=False, x_lim=None):
        plot_col, plot_row = self.make_figure(len(root_params_triplets), nx=nx)
        for i, (root, param1, param2) in enumerate(root_params_triplets):
            ax = self.subplot_number(i)
            self.plot_2d(root, param_pair=[param1, param2], filled=filled, add_legend_proxy=i == 0)
            if x_lim is not None:ax.set_xlim(x_lim)
        self.finish_plot()
        return plot_col, plot_row

    def spaceTicks(self, axis, expand=True):
            lims = axis.get_view_interval()
            tick = [tick for tick in axis.get_ticklocs() if tick > lims[0] and tick < lims[1]]
            gap_wanted = (lims[1] - lims[0]) * self.settings.tight_gap_fraction
            if expand:
                lims = [min(tick[0] - gap_wanted, lims[0]), max(tick[-1] + gap_wanted, lims[1])]
                axis.set_view_interval(lims[0], lims[1])
            else:
                if tick[0] - lims[0] < gap_wanted: tick = tick[1:]
                if lims[1] - tick[-1] < gap_wanted:tick = tick[:-1]
            axis.set_ticks(tick)
            return tick


    def triangle_plot(self, roots, in_params=None, legend_labels=None, plot_3d_with_param=None, filled=False, filled_compare=False, shaded=False,
                      contour_args=None, contour_colors=None, contour_ls=None, contour_lws=None, line_args=None, label_order=None):
        if isinstance(roots, basestring):roots = [roots]
        params = self.get_param_array(roots[0], in_params)
        plot_col = len(params)
        if plot_3d_with_param is not None: col_param = self.check_param(roots[0], plot_3d_with_param)
        self.make_figure(nx=plot_col, ny=plot_col)
        lims = dict()
        ticks = dict()
        line_args = None
        if filled_compare: filled = filled_compare
        contour_args = self._make_contour_args(len(roots), filled=filled, contour_args=contour_args,
                                               colors=contour_colors, ls=contour_ls, lws=contour_lws)
        if filled and not line_args:
            cols = [self.settings.solid_colors[len(roots) - plotno - 1] for plotno in range(len(params))]
            line_args = []
            for col in cols:
                if isinstance(col, tuple) or isinstance(col, list): col = col[-1]
                line_args += [{'color': col} ]
        for i, param in enumerate(params):
            ax = self.subplot(i, i)
            self.plot_1d(roots, param, do_xlabel=i == plot_col - 1, no_label_no_numbers=self.settings.no_triangle_axis_labels,
                         label_right=True, no_zero=True, no_ylabel=True, no_ytick=True, line_args=line_args)
            # set no_ylabel=True for now, can't see how to not screw up spacing with right-sided y label
            if self.settings.no_triangle_axis_labels: self.spaceTicks(ax.xaxis, expand=not shaded)
            lims[i] = ax.get_xlim()
            ticks[i] = ax.get_xticks()
        for i, param in enumerate(params):
            for i2 in range(i + 1, len(params)):
                param2 = params[i2]
                ax = self.subplot(i, i2)
                if plot_3d_with_param is not None:
                    self.plot_3d(roots, [param, param2, col_param], color_bar=False, line_offset=1, add_legend_proxy=False,
                      do_xlabel=i2 == plot_col - 1, do_ylabel=i == 0, contour_args=contour_args,
                      no_label_no_numbers=self.settings.no_triangle_axis_labels)
                else:
                    self.plot_2d(roots, param_pair=[param, param2], do_xlabel=i2 == plot_col - 1, do_ylabel=i == 0,
                                        no_label_no_numbers=self.settings.no_triangle_axis_labels, shaded=shaded,
                                         add_legend_proxy=i == 1 and i2 == i + 1, contour_args=contour_args)
                ax.set_xticks(ticks[i])
                ax.set_yticks(ticks[i2])
                ax.set_xlim(lims[i])
                ax.set_ylim(lims[i2])

        if self.settings.no_triangle_axis_labels:subplots_adjust(wspace=0, hspace=0)
        if plot_3d_with_param is not None:
            bottom = 0.5
            if len(params) == 2: bottom += 0.1;
            cb = self.fig.colorbar(self.last_scatter, cax=self.fig.add_axes([0.9, bottom, 0.03, 0.35]))
            self.add_colorbar_label(cb, col_param)

        self.finish_plot(self.default_legend_labels(legend_labels, roots), label_order=label_order,
                         legend_loc=None, no_gap=self.settings.no_triangle_axis_labels, no_extra_legend_space=True)

    def rectangle_plot(self, xparams, yparams, yroots=None, roots=None, plot_roots=None, plot_texts=None,
                       ymarkers=None, xmarkers=None, param_limits={}, legend_labels=None, legend_ncol=None,
                       label_order=None, marker_args={}, **kwargs):
            """
                roots uses the same set of roots for every plot in the rectangle
                yroots (list of list of roots) allows use of different set of roots for each row of the plot
                plot_roots allows you to specify (via list of list of list of roots) the set of roots for each individual subplot
            """
            self.make_figure(nx=len(xparams), ny=len(yparams))
#            f, plots = subplots(len(yparams), len(xparams), sharex='col', sharey='row')
            sharey = None
            yshares = []
            xshares = []
            ax_arr = []
            if plot_roots and yroots or roots and yroots or plot_roots and roots:
                raise Exception('rectangle plot: must have one of roots, yroots, plot_roots')
            limits = dict()
            for x, xparam in enumerate(xparams):
                sharex = None
                if plot_roots: yroots = plot_roots[x]
                elif roots: yroots = [roots for _ in yparams]
                axarray = []
                res = None
                for y, (yparam, subplot_roots) in enumerate(zip(yparams, yroots)):
                    if x > 0: sharey = yshares[y]
                    ax = self.subplot(x, y , sharex=sharex, sharey=sharey)
                    if y == 0:
                        sharex = ax
                        xshares.append(ax)
                    res = self.plot_2d(subplot_roots, param_pair=[xparam, yparam], do_xlabel=y == len(yparams) - 1,
                                 do_ylabel=x == 0, add_legend_proxy=x == 0 and y == 0, **kwargs)
                    if ymarkers is not None and ymarkers[y] is not None: self.add_y_marker(ymarkers[y], **marker_args)
                    if xmarkers is not None and xmarkers[x] is not None: self.add_x_marker(xmarkers[x], **marker_args)
                    limits[xparam], limits[yparam] = self.updateLimits(res, limits.get(xparam), limits.get(yparam))
                    if y != len(yparams) - 1: setp(ax.get_xticklabels(), visible=False)
                    if x != 0: setp(ax.get_yticklabels(), visible=False)
                    if x == 0: yshares.append(ax)
                    if plot_texts: self.add_text_left(plot_texts[x][y], y=0.9, ax=ax)
                    axarray.append(ax)
                ax_arr.append(axarray)
            for xparam, ax in zip(xparams, xshares):
                ax.set_xlim(param_limits.get(xparam, limits[xparam]))
                self.spaceTicks(ax.xaxis)
                ax.set_xlim(ax.xaxis.get_view_interval())
            for yparam, ax in zip(yparams, yshares):
                ax.set_ylim(param_limits.get(yparam, limits[yparam]))
                self.spaceTicks(ax.yaxis)
                ax.set_ylim(ax.yaxis.get_view_interval())
            subplots_adjust(wspace=0, hspace=0)
            if roots: legend_labels = self.default_legend_labels(legend_labels, roots)
            self.finish_plot(no_gap=True, legend_labels=legend_labels, label_order=label_order, legend_ncol=legend_ncol or len(legend_labels))
            return ax_arr

    def rotate_yticklabels(self, ax=None, rotation=90):
        if ax is None: ax = gca()
        for ticklabel in ax.get_yticklabels():
            ticklabel.set_rotation(rotation)

    def add_colorbar(self, param, orientation='vertical', **ax_args):
        cb = colorbar(orientation=orientation)
        cb.set_alpha(1)
        cb.draw_all()
        if not ax_args.get('color_label_in_axes'):
            self.add_colorbar_label(cb, param)
            if self.settings.colorbar_rotation is not None:
                self.rotate_yticklabels(cb.ax, self.settings.colorbar_rotation)
                labels = [label.get_text() for label in cb.ax.yaxis.get_ticklabels()[::2]]
                cb.ax.yaxis.set_ticks(cb.ax.yaxis.get_ticklocs()[::2])
                cb.ax.yaxis.set_ticklabels(labels)
        return cb

    def add_line(self, P1, P2, zorder=0, color=None, ls=None, ax=None, **kwargs):
            if color is None: color = self.settings.axis_marker_color
            if ls is None: ls = self.settings.axis_marker_ls
            (ax or gca()).add_line(Line2D(P1, P2, color=color, ls=ls, zorder=zorder, **kwargs))

    def add_colorbar_label(self, cb, param):
        cb.set_label(r'$' + param.label + '$', fontsize=self.settings.lab_fontsize,
                     rotation=self.settings.colorbar_label_rotation, labelpad=self.settings.colorbar_label_pad)
        setp(getp(cb.ax, 'ymajorticklabels'), fontsize=self.settings.colorbar_axes_fontsize)

    def _makeParamObject(self, names, samples):
        class sampleNames(object): pass
        p = sampleNames()
        for i, par in enumerate(names.names):
            setattr(p, par.name, samples[:, i])
        return p

    def add_3d_scatter(self, root, in_params, color_bar=True, alpha=1, extra_thin=1, **ax_args):
        params = self.get_param_array(root, in_params)
        pts = self.sampleAnalyser.load_single_samples(root)
        names = self.paramNamesForRoot(root)
        samples = []
        for param in params:
            if hasattr(param, 'getDerived'):
                samples.append(param.getDerived(self._makeParamObject(names, pts)))
            else:
                samples.append(pts[:, names.numberOfName(param.name)])
        if extra_thin > 1:
            samples = [ pts[::extra_thin] for pts in samples]
        self.last_scatter = scatter(samples[0], samples[1], edgecolors='none',
                s=self.settings.scatter_size, c=samples[2], cmap=self.settings.colormap_scatter, alpha=alpha)
        if color_bar: self.last_colorbar = self.add_colorbar(params[2], **ax_args)
        xbounds = [min(samples[0]), max(samples[0])]
        r = xbounds[1] - xbounds[0]
        xbounds[0] -= r / 20
        xbounds[1] += r / 20
        ybounds = [min(samples[1]), max(samples[1])]
        r = ybounds[1] - ybounds[0]
        ybounds[0] -= r / 20
        ybounds[1] += r / 20
        return [xbounds, ybounds]

    def plot_3d(self, roots, in_params=None, params_for_plots=None, color_bar=True, line_offset=0, add_legend_proxy=True, **kwargs):
        if isinstance(roots, basestring): roots = [roots]
        if params_for_plots:
            params_for_plots = [self.get_param_array(root, p) for p, root in zip(params_for_plots, roots)]
        else:
            if not in_params: raise Exception('No parameters for plot_3d!')
            params = self.get_param_array(roots[0], in_params)
            params_for_plots = [params for root in roots]  # all the same
        if self.fig is None: self.make_figure()
        contour_args = self._make_contour_args(len(roots) - 1, **kwargs)
        xlims, ylims = self.add_3d_scatter(roots[0], params_for_plots[0], color_bar=color_bar, **kwargs)
        for i, root in enumerate(roots[1:]):
            params = params_for_plots[i + 1]
            res = self.add_2d_contours(root, params[0], params[1], i + line_offset, add_legend_proxy=add_legend_proxy, zorder=i + 1, **contour_args[i])
            xlims, ylims = self.updateLimits(res, xlims, ylims)
        if not 'lims' in kwargs:
            params = params_for_plots[0]
            lim1 = self.checkBounds(roots[0], params[0].name , xlims[0], xlims[1])
            lim2 = self.checkBounds(roots[0], params[1].name , ylims[0], ylims[1])
            kwargs['lims'] = [lim1[0], lim1[1], lim2[0], lim2[1]]
        self.setAxes(params, **kwargs)

    def plots_3d(self, roots, param_sets, nx=None, filled_compare=False, legend_labels=None, **kwargs):
        if isinstance(roots, basestring):roots = [roots]
        sets = [[self.check_param(roots[0], param) for param in param_group] for param_group in param_sets]
        plot_col, plot_row = self.make_figure(len(sets), nx=nx, xstretch=1.3)

        for i, triplet in enumerate(sets):
            self.subplot_number(i)
            self.plot_3d(roots, triplet, filled=filled_compare, **kwargs)
        self.finish_plot(self.default_legend_labels(legend_labels, roots[1:]), no_tight=True)
        return plot_col, plot_row

    def plots_3d_z(self, roots, param_x, param_y, param_z=None, max_z=None, **kwargs):
        """Make set of plots of param_x against param_y, each coloured by values of parameters in param_z (all if None)"""
        if isinstance(roots, basestring):roots = [roots]
        param_z = self.get_param_array(roots[0], param_z)
        if max_z is not None and len(param_z) > max_z: param_z = param_z[:max_z]
        param_x, param_y = self.get_param_array(roots[0], [param_x, param_y])
        sets = [[param_x, param_y, z] for z in param_z if z != param_x and z != param_y]
        return self.plots_3d(roots, sets, **kwargs)

    def add_text(self, text_label, x=0.95, y=0.06, ax=None, **kwargs):
        args = {'horizontalalignment':'right', 'verticalalignment':'center'}
        args.update(kwargs)
        if isinstance(ax, int):
            ax = self.fig.axes[ax]
        if isinstance(ax, list):
            ax = self.subplots[ax[0], ax[1]]
        else:
            ax = ax or gca()
        ax.text(x, y, text_label, transform=ax.transAxes, **args)

    def add_text_left(self, text_label, x=0.05, y=0.06, ax=None, **kwargs):
        args = {'horizontalalignment':'left'}
        args.update(kwargs)
        self.add_text(text_label, x, y, ax, **args)

    def export(self, fname=None, adir=None, watermark=None, tag=None):
        if fname is None: fname = os.path.basename(sys.argv[0]).replace('.py', '')
        if tag: fname += '_' + tag
        if not '.' in fname: fname += '.pdf'
        if adir is not None and not os.sep in fname: fname = os.path.join(adir, fname)
        adir = os.path.dirname(fname)
        if adir and not os.path.exists(adir): os.makedirs(adir)
        if watermark:
            gcf().text(0.45, 0.5, self._escapeLatex(watermark), fontsize=30, color='gray', ha='center', va='center', alpha=0.2)

        savefig(fname, bbox_extra_artists=self.extra_artists, bbox_inches='tight')

    def paramNameListFromFile(self, fname):
        p = paramNames.paramNames(fname)
        return [name.name for name in p.names]



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
    g.triangle_plot(roots, ['omegabh2', 'omegach2', 'ns', 'omegal', 'omegak'], plot_3d_with_param='H0',
                    legend_labels=['Planck', 'Planck+lensing', 'Planck+BAO', 'Planck (flat)'])
    g.export(roots[0] + '_tri.pdf')

    roots = ['base_planck_CAMspec_lowl_lowLike', 'base_planck_CAMspec_lowl_lowLike_post_BAO']
    g.plots_3d(roots, [['omegabh2', 'omegach2', 'ns'], ['omegach2', 'logA', 'tau']])
    g.export(roots[0] + '_3d.pdf')


# sample_plots()
