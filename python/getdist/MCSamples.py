# MCSamples.py

import os
import sys
import glob
import math
import logging
import copy
import numpy as np
from scipy.interpolate import splrep, splev
from scipy.stats import norm
from scipy.signal import fftconvolve
import pickle
from iniFile import iniFile
from getdist.chains import chains, chainFiles, lastModified
from getdist import ResultObjs

# =============================================================================

version = 6

default_grid_root = None
output_base_dir = None
cache_dir = None
use_plot_data = False
default_getdist_settings = os.path.join(os.path.dirname(__file__), 'analysis_defaults.ini')
default_plot_output = 'pdf'

config_file = os.environ.get('GETDIST_CONFIG', None)
if not config_file:
    config_file = os.path.join(os.path.dirname(__file__), 'config.ini')
if os.path.exists(config_file):
    config_ini = iniFile(config_file)
    default_grid_root = config_ini.string('default_grid_root', '')
    output_base_dir = config_ini.string('output_base_dir', '')
    cache_dir = config_ini.string('cache_dir', '')
    default_getdist_settings = config_ini.string('default_getdist_settings', default_getdist_settings)
    use_plot_data = config_ini.bool('use_plot_data', use_plot_data)
    default_plot_output = config_ini.string('default_plot_output', default_plot_output)
else:
    config_ini = iniFile()

class MCException(Exception):
    pass

class FileException(MCException):
    pass

class SettingException(MCException):
    pass

class ParamException(MCException):
    pass


def loadMCSamples(file_root, ini=None, jobItem=None, no_cache=False, dist_settings={}):
        files = chainFiles(file_root)
        path, name = os.path.split(file_root)
        path = cache_dir or path
        if not os.path.exists(path): os.mkdir(path)
        cachefile = os.path.join(path, name) + '.py_mcsamples'
        samples = MCSamples(file_root, jobItem=jobItem)
        samples.update_settings(ini, dist_settings)
        allfiles = files + [file_root + '.ranges', file_root + '.paramnames']
        if not no_cache and os.path.exists(cachefile) and lastModified(allfiles) < os.path.getmtime(cachefile):
            try:
                with open(cachefile, 'rb') as inp:
                    cache = pickle.load(inp)
                if cache.version == version and samples.ignore_rows == cache.ignore_rows:
                    changed = len(samples.contours) != len(cache.contours) or \
                                np.any(np.array(samples.contours) != np.array(cache.contours))
                    cache.update_settings(ini, dist_settings, doUpdate=changed)
                    return cache
            except:
                pass
        if not len(files): FileException('No chains found: ' + file_root)
        samples.readChains(files)
        with open(cachefile, 'wb') as output:
                pickle.dump(samples, output, pickle.HIGHEST_PROTOCOL)
        return samples

def get2DContourLevels(bins2D, contours=[0.68, 0.95], missing_norm=0, half_edge=True):
    """
     Get countour levels encolosing "countours" of the probability.
     bins 2D is the density. If half_edge, then edge bins are only half integrated over in each direction.
     missing_norm accounts of any points not included in bin2D (e.g. points in far tails that are not plotted)
    """
    contour_levels = np.zeros(len(contours))
    if half_edge:
        abins = bins2D.copy()
        abins[:, 0] /= 2
        abins[0, :] /= 2
        abins[:, -1] /= 2
        abins[-1, :] /= 2
    else:
        abins = bins2D
    norm = np.sum(abins)
    targets = (1 - np.array(contours)) * norm - missing_norm
    if True:
        bins = abins.reshape(-1)
        indexes = bins2D.reshape(-1).argsort()
        sortgrid = bins[indexes]
        cumsum = np.cumsum(sortgrid)
        ixs = np.searchsorted(cumsum, targets)
        for i, ix in enumerate(ixs):
            if ix == 0:
                raise MCException("Contour level outside plotted ranges")
            h = cumsum[ix] - cumsum[ix - 1]
            d = (cumsum[ix] - targets[i]) / h
            contour_levels[i] = sortgrid[ix] * (1 - d) + d * sortgrid[ix - 1]
    else:
        # old method, twice or so slower
        try_t = np.max(bins2D)
        lastcontour = 0
        for i, (contour, target) in enumerate(zip(contours, targets)):
            if contour < lastcontour: raise SettingException('contour levels must be decreasing')
            lastcontour = contour
            try_b = 0
            lasttry = -1
            while True:
                try_sum = np.sum(abins[bins2D < (try_b + try_t) / 2])
                if try_sum > target:
                    try_t = (try_b + try_t) / 2
                else:
                    try_b = (try_b + try_t) / 2
                if try_sum == lasttry: break
                lasttry = try_sum
            contour_levels[i] = (try_b + try_t) / 2
    return contour_levels


class Ranges(object):

    def __init__(self, fileName=None, setParamNameFile=None):
        self.names = []
        self.mins = {}
        self.maxs = {}
        if fileName is not None: self.loadFromFile(fileName)

    def loadFromFile(self, fileName):
        self.filenameLoadedFrom = os.path.split(fileName)[1]
        f = open(fileName)
        for line in f:
            name, mini, maxi = self.readValues(line.strip())
            if name:
                self.names.append(name)
                self.mins[name] = mini
                self.maxs[name] = maxi

    def readValues(self, line):
        name, mini, maxi = None, None, None
        strings = [ s for s in line.split(" ")  if s <> '' ]
        if len(strings) == 3:
            name = strings[0]
            mini = float(strings[1])
            maxi = float(strings[2])
        return name, mini, maxi

    def listNames(self):
        return self.names

    def min(self, name, error=False):
        if self.mins.has_key(name):
            return self.mins[name]
        if error: ParamException("Name not found:" + name)
        return None

    def max(self, name, error=False):
        if self.maxs.has_key(name):
            return self.maxs[name]
        if error: raise ParamException("Name not found:" + name)
        return None

# =============================================================================

class Kernel1D(object):

    def __init__(self, winw, h):
        self.winw = winw
        self.h = h
        self.x = np.arange(-winw, winw + 1)
        Win = np.exp(-(self.x / h) ** 2 / 2.)
        self.Win = Win / np.sum(Win)


class Density1D(object):

    def __init__(self, xmin, n, spacing, P=None):
        self.n = n
        self.X = xmin + np.arange(n) * float(spacing)
        if P is not None:
            self.P = P
        else:
            self.P = np.zeros(n)
        self.spacing = float(spacing)

    def InitSpline(self):
        self.spl = splrep(self.X, self.P, s=0)

    def Prob(self, x):
        if isinstance(x, np.ndarray):
            return splev(x, self.spl)
        else:
            return splev([x], self.spl)


    def initLimitGrids(self, factor=100):
        class InterpGrid(object): pass

        g = InterpGrid()
        g.factor = factor
        g.bign = (self.n - 1) * factor + 1
        vecx = self.X[0] + np.arange(g.bign) * self.spacing / factor
        g.grid = splev(vecx, self.spl)

        norm = np.sum(g.grid)
        g.norm = norm - (0.5 * self.P[-1]) - (0.5 * self.P[0])

        g.sortgrid = np.sort(g.grid)
        g.cumsum = np.cumsum(g.sortgrid)
        return g

    def Limits(self, p, interpGrid=None):
        g = interpGrid
        if g is None: g = self.initLimitGrids()
        # Return values
        mn, mx = 0., 0.

        target = (1 - p) * g.norm
        ix = np.searchsorted(g.cumsum, target)
        trial = g.sortgrid[ix]
        if ix > 0:
            d = g.cumsum[ix] - g.cumsum[ix - 1]
            frac = (g.cumsum[ix] - target) / d
            trial = (1 - frac) * trial + frac * g.sortgrid[ix + 1]
#         try_t = max(grid)
#         try_b = 0
#         try_last = -1
#         while True:
#             trial = (try_b + try_t) / 2
#             trial_sum = np.sum(grid[grid >= trial])
#             if (trial_sum < p * norm):
#                 try_t = (try_b + try_t) / 2
#             else:
#                 try_b = (try_b + try_t) / 2
#             if (math.fabs(trial_sum / try_last - 1) < 1e-4): break
#             try_last = trial_sum
#         trial = (try_b + try_t) / 2

        lim_bot = (g.grid[0] >= trial)
        if lim_bot:
            mn = self.P[0]
        else:
            for i in range(g.bign):
                if g.grid[i] > trial:
                    mn = self.X[0] + (i - 1) * self.spacing / g.factor
                    break

        lim_top = (g.grid[-1] >= trial)
        if lim_top:
            mx = self.P[-1]
        else:
            indexes = range(g.bign)
            indexes.reverse()
            for i in indexes:
                if (g.grid[i] > trial):
                    mx = self.X[0] + (i - 1) * self.spacing / g.factor
                    break

        return mn, mx, lim_bot, lim_top

# =============================================================================

class MCSamples(chains):

    def __init__(self, root=None, ignore_rows=0, jobItem=None, ini=None):
        chains.__init__(self, root, jobItem=jobItem)

        self.version = version
        self.limmin = []
        self.limmax = []

        self.has_limits = []
        self.has_limits_bot = []
        self.has_limits_top = []

        self.markers = {}

        self.ini = ini

        self.ReadRanges()

        # Other variables
        self.no_plots = False
        self.num_bins = 100
        self.num_bins_2D = 40
        self.smooth_scale_1D = 0.25
        self.smooth_scale_2D = 2
        self.boundary_correction_method = 1
        self.smooth_correct_contours = True
        self.max_corr_2D = 0.95
        self.num_contours = 2
        self.contours = [0.68, 0.95]
        self.max_scatter_points = 2000
        self.credible_interval_threshold = 0.05
        self.plot_meanlikes = False
        self.shade_meanlikes = False

        self.num_vars = 0
        self.numsamp = 0
        self.max_mult = 0
        self.mean_mult = 0
        self.plot_data_dir = ""
        if root:
            self.rootname = os.path.basename(root)
        else:
            self.rootname = ""

        self.rootdirname = ""
        self.meanlike = 0.
        self.mean_loglikes = False
        self.indep_thin = 0
        self.samples = None
        self.ignore_rows = ignore_rows
        self.subplot_size_inch = 4.0
        self.subplot_size_inch2 = self.subplot_size_inch
        self.subplot_size_inch3 = 6.0
        self.plot_output = default_plot_output
        self.out_dir = ""

        self.max_split_tests = 4
        self.force_twotail = False

        self.corr_length_thin = 0
        self.corr_length_steps = 15

        self.done_1Dbins = False
        self.done2D = None

        self.density1D = dict()

        self.update_settings(ini)

    def parName(self, i, starDerived=False):
        return self.paramNames.name(i, starDerived)

    def parLabel(self, i):
        return self.paramNames.names[i].label

    def initParameters(self, ini):

        self.ignore_rows = self.ini.float('ignore_rows', self.ignore_rows)
        self.ignore_lines = int(self.ignore_rows)
        if not self.ignore_lines:
            self.ignore_frac = self.ignore_rows
        else:
            self.ignore_frac = 0

        self.num_bins = ini.int('num_bins', self.num_bins)
        self.num_bins_2D = ini.int('num_bins_2D', self.num_bins_2D)
        self.smooth_scale_1D = ini.float('smooth_scale_1D', self.smooth_scale_1D)
        self.smooth_scale_2D = ini.float('smooth_scale_2D', self.smooth_scale_2D)

        if self.smooth_scale_1D > 0 and self.smooth_scale_1D > 1:
            raise SettingException('mooth_scale_1D>1 is oversmoothed')
        if self.smooth_scale_1D > 0 and self.smooth_scale_1D > 1.9:
            raise SettingException('smooth_scale_1D>1 is now in stdev units')

        self.boundary_correction_method = ini.int('boundary_correction_method',
                                                  getattr(self, 'boundary_correction_method', 1))

        self.smooth_correct_contours = ini.bool('smooth_correct_contours', getattr(self, 'smooth_correct_contours', True))
        self.no_plots = ini.bool('no_plots', False)
        self.shade_meanlikes = ini.bool('shade_meanlikes', False)

        self.force_twotail = ini.bool('force_twotail', False)

        self.plot_meanlikes = ini.bool('plot_meanlikes', self.plot_meanlikes)

        self.max_scatter_points = ini.int('max_scatter_points', self.max_scatter_points)
        self.credible_interval_threshold = ini.float('credible_interval_threshold', self.credible_interval_threshold)

        self.subplot_size_inch = ini.float('subplot_size_inch' , self.subplot_size_inch)
        self.subplot_size_inch2 = ini.float('subplot_size_inch2', self.subplot_size_inch)
        self.subplot_size_inch3 = ini.float('subplot_size_inch3', self.subplot_size_inch)

        self.plot_output = ini.string('plot_output', self.plot_output)

        self.force_twotail = ini.bool('force_twotail', False)
        if self.force_twotail: print 'Computing two tail limits'

        self.max_corr_2D = ini.float('max_corr_2D', self.max_corr_2D)

        if ini and ini.hasKey('num_contours'):
            self.num_contours = ini.int('num_contours', 2)
            self.contours = np.array([ini.float('contour' + str(i + 1)) for i in range(self.num_contours)])
        # how small the end bin must be relative to max to use two tail
        self.max_frac_twotail = []
        for i in range(self.num_contours):
            max_frac = np.exp(-1.0 * math.pow(norm.ppf((1 - self.contours[i]) / 2), 2) / 2)
            if ini:
                max_frac = ini.float('max_frac_twotail' + str(i + 1), max_frac)
            self.max_frac_twotail.append(max_frac)



    def initLimits(self, ini=None):

        bin_limits = ""
        if ini: bin_limits = ini.string('all_limits', '')

        nvars = len(self.paramNames.names)

        self.limmin = np.zeros(nvars)
        self.limmax = np.zeros(nvars)

        self.has_limits = nvars * [False]
        self.has_limits_bot = nvars * [False]
        self.has_limits_top = nvars * [False]

        self.markers = {}

        for ix, name in enumerate(self.paramNames.list()):
            mini = self.ranges.min(name)
            maxi = self.ranges.max(name)
            if mini is not None and maxi is not None and mini <> maxi:
                self.limmin[ix] = mini
                self.limmax[ix] = maxi
                self.has_limits_top[ix] = True
                self.has_limits_bot[ix] = True
            if bin_limits <> '':
                line = bin_limits
            else:
                line = ''
                if ini and ini.params.has_key('limits[%s]' % name):
                    line = ini.string('limits[%s]' % name)
            if line <> '':
                limits = [ s for s in line.split(' ') if s <> '' ]
                if len(limits) == 2:
                    if limits[0] <> 'N':
                        self.limmin[ix] = float(limits[0])
                        self.has_limits_bot[ix] = True
                    if limits[1] <> 'N':
                        self.limmax[ix] = float(limits[1])
                        self.has_limits_top[ix] = True
            if ini and ini.params.has_key('marker[%s]' % name):
                line = ini.string('marker[%s]' % name)
                if line <> '':
                    self.markers[name] = float(line)

    def update_settings(self, ini=None, ini_settings={}, doUpdate=True):
        if not ini:
            ini = self.ini
        elif isinstance(ini, basestring):
            ini = iniFile(ini)
        else:
            ini = copy.deepcopy(ini)
        if not ini: ini = iniFile(default_getdist_settings)
        if ini_settings:
            ini.params.update(ini_settings)
        self.ini = ini
        if ini: self.initParameters(ini)
        if doUpdate and self.samples: self.updateChainBaseStatistics()

    def readChains(self, chain_files):
        # Used for by plotting scripts and gui

        self.loadChains(self.root, chain_files)

        if self.ignore_frac and  (not self.jobItem or
                    (not self.jobItem.isImportanceJob and not self.jobItem.isBurnRemoved())):
            self.removeBurnFraction(self.ignore_frac)

        self.deleteFixedParams()

        # Make a single array for chains
        self.makeSingle()

        self.updateChainBaseStatistics()

        return self

    def updateChainBaseStatistics(self):

        super(MCSamples, self).updateChainBaseStatistics()
        self.initLimits(self.ini)

        self.ComputeMultiplicators()

        # Compute statistics values
        self.ComputeStats()

        self.GetCovMatrix()

        # Find best fit, and mean likelihood
        self.GetChainLikeSummary()

        # Init arrays for 1D densities
        self.Init1DDensity()

        # Get ND confidence region
        self.GetConfidenceRegion()

    def AdjustPriors(self):
        sys.exit('You need to write the AdjustPriors function in MCSamples.py first!')
        # print "Adjusting priors"
        # ombh2 = self.samples[0]  # ombh2 prior
        # chisq = (ombh2 - 0.0213)**2/0.001**2
        # self.weights *= np.exp(-chisq/2)
        # self.loglikes += chisq/2

    def MapParameters(self, invars):
        sys.exit('Need to write MapParameters routine first')

    def CoolChain(self, cool):
        print 'Cooling chains by ', cool
        MaxL = np.max(self.loglikes)
        newL = self.loglikes * cool
        newW = self.weights * np.exp(-(newL - self.loglikes) - (MaxL * (1 - cool)))
        self.weights = np.hstack(newW)
        self.loglikes = np.hstack(newL)

    def DeleteZeros(self):
        indexes = np.where(self.weights == 0)
        self.weights = np.delete(self.weights, indexes)
        self.loglikes = np.delete(self.loglikes, indexes)
        self.samples = np.delete(self.samples, indexes, axis=0)


    def MakeSingleSamples(self, filename="", single_thin=None):
        """
        Make file of weight-1 samples by choosing samples
        with probability given by their weight.
        """
        if single_thin is None:
            single_thin = max(1, self.numsamp / self.max_mult / self.max_scatter_points)
        rand = np.random.random_sample(self.numrows)
        maxmult = np.max(self.weights)

        if filename:
            textFileHandle = open(filename, 'w')
            for i, r in enumerate(rand):
                if (r <= self.weights[i] / maxmult / single_thin):
                    textFileHandle.write("%16.7E" % (1.0))
                    textFileHandle.write("%16.7E" % (self.loglikes[i]))
                    for j in range(self.num_vars):
                        textFileHandle.write("%16.7E" % (self.samples[i][j]))
                    textFileHandle.write("\n")
            textFileHandle.close()
        else:
            # return data
            return self.samples[rand <= self.weights / (maxmult * single_thin)]

    def WriteThinData(self, fname, thin_ix, cool):
        nparams = self.samples.shape[1]
        if(cool <> 1): print 'Cooled thinned output with temp: ', cool
        MaxL = np.max(self.loglikes)
        textFileHandle = open(fname, 'w')
        i = 0
        for thin in thin_ix:
            if (cool <> 1):
                newL = self.loglikes[thin] * cool
                textFileHandle.write("%16.7E" % (
                        np.exp(-(newL - self.loglikes[thin]) - MaxL * (1 - cool))))
                textFileHandle.write("%16.7E" % (newL))
                for j in nparams:
                    textFileHandle.write("%16.7E" % (self.samples[i][j]))
            else:
                textFileHandle.write("%f" % (1.))
                textFileHandle.write("%f" % (self.loglikes[thin]))
                for j in nparams:
                    textFileHandle.write("%16.7E" % (self.samples[i][j]))
            i += 1

        textFileHandle.close()
        print  'Wrote ', len(thin_ix), ' thinned samples'

    # ThinData(self, fac) => chains.thin_indices(factor)

    def GetCovMatrix(self):
        nparam = self.paramNames.numParams()
        paramVecs = [ self.samples[:, i] for i in range(nparam) ]
        self.fullcov = self.cov(paramVecs)
        nparamNonDerived = self.paramNames.numNonDerived()
        self.covmatrix = self.fullcov[:nparamNonDerived, :nparamNonDerived]
        self.corrmatrix = self.corr(paramVecs, cov=self.fullcov)

    def writeCovMatrix(self, filename=None):
        filename = filename or self.rootdirname + ".covmat"
        nparamNonDerived = self.paramNames.numNonDerived()
        with open(filename, "w") as textFileHandle:
            textFileHandle.write("# %s\n" % (" ".join(self.paramNames.list()[:nparamNonDerived])))
            for i in range(nparamNonDerived):
                for j in range(nparamNonDerived):
                    textFileHandle.write("%17.7E" % self.covmatrix[i][j])
                textFileHandle.write("\n")

    def writeCorrMatrix(self, filename=None):
        filename = filename or self.rootdirname + ".corr"
        np.savetxt(filename, self.corrmatrix, fmt="%17.7E")


    def GetFractionIndices(self, weights, n):
        nrows = weights.shape[0]
        numsamp = np.sum(weights)

        tot = 0
        aim = numsamp / n
        num = 0

        fraction_indices = []
        fraction_indices.append(0)

        for i in range(nrows):
            tot = tot + weights[i]
            if (tot > aim):
                num = num + 1
                fraction_indices.append(i)
                if (num == n): break
                aim = aim + numsamp / n

        fraction_indices.append(nrows)
        return fraction_indices


    def PCA(self, params, param_map, normparam=None, writeDataToFile=False, conditional_params=[]):
        """
        Perform principle component analysis. In other words,
        get eigenvectors and eigenvalues for normalized variables
        with optional (log) mapping.
        """

        print 'Doing PCA for ', len(params), ' parameters'
        if len(conditional_params): print 'conditional %u fixed parameters' % len(conditional_params)

        PCAtext = 'PCA for parameters:\n'

        params = [name for name in params if self.paramNames.parWithName(name)]
        nparams = len(params)
        indices = [self.index[param] for param in params]
        conditional_params = [self.index[param] for param in conditional_params]
        indices += conditional_params

        if normparam:
            if normparam in params:
                normparam = params.index(normparam)
            else: normparam = -1
        else: normparam = -1

        n = len(indices)
        corrmatrix = np.zeros((n, n))
        PCdata = self.samples[:, indices]
        PClabs = []

        PCmean = np.zeros(n)
        sd = np.zeros(n)
        newmean = np.zeros(n)
        newsd = np.zeros(n)

        doexp = False
        for i, parix in enumerate(indices):
            if i < nparams:
                label = self.parLabel(parix)
                if (param_map[i] == 'L'):
                    doexp = True
                    PCdata[:, i] = np.log(PCdata[:, i])
                    PClabs.append("ln(" + label + ")")
                elif (param_map[i] == 'M'):
                    doexp = True
                    PCdata[:, i] = np.log(-1.0 * PCdata[:, i])
                    PClabs.append("ln(-" + label + ")")
                else:
                    PClabs.append(label)
                PCAtext += "%10s :%s\n" % (str(parix + 1), str(PClabs[i]))

            PCmean[i] = np.sum(self.weights * PCdata[:, i]) / self.norm
            PCdata[:, i] = PCdata[:, i] - PCmean[i]
            sd[i] = np.sqrt(np.sum(self.weights * np.power(PCdata[:, i], 2)) / self.norm)
            if (sd[i] <> 0): PCdata[:, i] = PCdata[:, i] / sd[i]

        PCAtext += "\n"
        PCAtext += 'Correlation matrix for reduced parameters\n'
        for i, parix in enumerate(indices):
            corrmatrix[i][i] = 1
            for j in range(i):
                corrmatrix[j][i] = np.sum(self.weights * PCdata[:, i] * PCdata[:, j]) / self.norm
                corrmatrix[i][j] = corrmatrix[j][i]
        for i in range(nparams):
            PCAtext += '%12s :' % params[i]
            for j in range(n):
                PCAtext += '%8.4f' % corrmatrix[j][i]
            PCAtext += '\n'

        if len(conditional_params):
            u = np.linalg.inv(corrmatrix)
            u = u[np.ix_(range(len(params)), range(len(params)))]
            u = np.linalg.inv(u)
            n = nparams
            PCdata = PCdata[:, :nparams]
        else:
            u = corrmatrix
        evals, evects = np.linalg.eig(u)
        isorted = evals.argsort()
        u = np.transpose(evects[:, isorted])  # redefining u

        PCAtext += '\n'
        PCAtext += 'e-values of correlation matrix\n'
        for i in range(n):
            isort = isorted[i]
            PCAtext += 'PC%2i: %8.4f\n' % (i + 1, evals[isort])

        PCAtext += '\n'
        PCAtext += 'e-vectors\n'
        for j in range(n):
            PCAtext += '%3i:' % (indices[j] + 1)
            for i in range(n):
                isort = isorted[i]
                '%8.4f' % (evects[j][isort])
            PCAtext += '\n'

        if (normparam <> -1):
            # Set so parameter normparam has exponent 1
            for i in range(n):
                u[i, :] = u[i, :] / u[i, normparam] * sd[normparam]
        else:
            # Normalize so main component has exponent 1
            for i in range(n):
                maxi = np.abs(u[i, :]).argmax()
                u[i, :] = u[i, :] / u[i, maxi] * sd[maxi]

        nrows = PCdata.shape[0]
        for i in range(nrows):
            PCdata[i, :] = np.dot(u, PCdata[i, :])
            if (doexp): PCdata[i, :] = np.exp(PCdata[i, :])

        PCAtext += '\n'
        PCAtext += 'Principle components\n'

        for i in range(n):
            isort = isorted[i]
            PCAtext += 'PC%i (e-value: %f)\n' % (i + 1, evals[isort])
            for j in range(n):
                label = self.parLabel(indices[j])
                if (param_map[j] in ['L', 'M']):
                    expo = "%f" % (1.0 / sd[j] * u[i][j])
                    if (param_map[j] == "M"):
                        div = "%f" % (-np.exp(PCmean[j]))
                    else:
                        div = "%f" % (np.exp(PCmean[j]))
                    PCAtext += '[%f]  (%s/%s)^{%s}\n' % (u[i][j], label, div, expo)
                else:
                    expo = "%f" % (sd[j] / u[i][j])
                    if (doexp):
                        PCAtext += '[%f]   exp((%s-%f)/%s)\n' % (u[i][j], label, PCmean[j], expo)
                    else:
                        PCAtext += '[%f]   (%s-%f)/%s)\n' % (u[i][j], label, PCmean[j], expo)

            newmean[i] = np.sum(self.weights * PCdata[:, i]) / self.norm
            newsd[i] = np.sqrt(np.sum(self.weights * np.power(PCdata[:, i] - newmean[i], 2)) / self.norm)
            PCAtext += '          = %f +- %f\n' % (newmean[i], newsd[i])
            PCAtext += 'ND limits: %9.3f%9.3f%9.3f%9.3f\n' % (
                    np.min(PCdata[0:self.ND_cont1, i]), np.max(PCdata[0:self.ND_cont1, i]),
                    np.min(PCdata[0:self.ND_cont2, i]), np.max(PCdata[0:self.ND_cont2, i]))
            PCAtext += '\n'

        # Find out how correlated these components are with other parameters
        PCAtext += 'Correlations of principle components\n'
        l = [ "%8i" % i for i in range(1, n + 1) ]
        PCAtext += '%s\n' % ("".join(l))

        for i in range(n):
            PCdata[:, i] = (PCdata[:, i] - newmean[i]) / newsd[i]

        for j in range(n):
            PCAtext += 'PC%2i' % (j + 1)
            for i in range(n):
                PCAtext += '%8.3f' % (
                        np.sum(self.weights * PCdata[:, i] * PCdata[:, j]) / self.norm)
            PCAtext += '\n'

        for j in range(self.num_vars):
            PCAtext += '%4i' % (j + 1)
            for i in range(n):
                PCAtext += '%8.3f' % (
                        np.sum(self.weights * PCdata[:, i]
                               * (self.samples[:, j] - self.means[j]) / self.sddev[j]) / self.norm)

            PCAtext += '   (%s)\n' % (self.parLabel(j))

        if writeDataToFile:
            filename = self.rootdirname + ".PCA"
            with open(filename, "w") as f:
                f.write(PCAtext)
        else:
            return PCAtext

    def ComputeMultiplicators(self):
        self.mean_mult = self.norm / self.numrows
        self.max_mult = (self.mean_mult * self.numrows) / min(self.numrows / 2, 500)
        outliers = len(self.weights[self.weights > self.max_mult])
        if (outliers <> 0):
            print 'outlier fraction ', float(outliers) / self.numrows
        self.max_mult = np.max(self.weights)
        self.numsamp = np.sum(self.weights)


    def GetUsedColsChains(self):
        if not hasattr(self, 'chains') or len(self.chains) == 0: return
        nparams = self.chains[0].coldata.shape[1] - 2
        isused = np.zeros(nparams, dtype=bool)
        for ix in range(nparams):
            isused[ix] = not np.all(self.chains[0].coldata[:, ix + 2] == self.chains[0].coldata[0][ix + 2])
        return isused

    def DoConvergeTests(self, limfrac, quick=True):
        """
        Do convergence tests. This is run before combining into one set of samples
        """
        if not hasattr(self, 'chains'): return

        isused = self.GetUsedColsChains()

        # Weights
        weights = np.hstack((chain.coldata[:, 0] for chain in self.chains))

        # Get statistics for individual chains, and do split tests on the samples

        fname = self.rootdirname + '.converge'
        textFileHandle = open(fname, 'w')

        num_chains_used = len(self.chains)
        if (num_chains_used > 1):
            print 'Number of chains used = ', num_chains_used

        nparam = self.paramNames.numParams()
        for chain in self.chains: chain.getCov(nparam)
        mean = np.zeros(nparam)
        norm = np.sum([chain.norm for chain in self.chains])
        nrows = np.sum([ chain.coldata.shape[0] for chain in self.chains ])
        for chain in self.chains:
            mean = mean + chain.means[0:nparam] * chain.norm
        mean /= norm

        fullmean = np.zeros(nparam)
        for j in range(nparam):
            for chain in self.chains:
                fullmean[j] += np.sum(chain.coldata[:, 0] * chain.coldata[:, j + 2]) / norm

        fullvar = np.zeros(nparam)
        for j in range(nparam):
            for chain in self.chains:
                fullvar[j] += np.sum(chain.coldata[:, 0] * np.power(chain.coldata[:, j + 2] - fullmean[j], 2)) / norm

        if (num_chains_used > 1):
            textFileHandle.write(" \n")
            textFileHandle.write(" Variance test convergence stats using remaining chains\n")
            textFileHandle.write(" param var(chain mean)/mean(chain var)\n")
            textFileHandle.write(" \n")

            between_chain_var = np.zeros(nparam)
            in_chain_var = np.zeros(nparam)
            chain_means = np.zeros((num_chains_used, nparam))

            # norm = np.sum([chain.norm for chain in self.chains])
            for j in range(nparam):
                if not isused[j]: continue

                # Get stats for individual chains - the variance of the means over the mean of the variances
                for i in range(num_chains_used):
                    chain = self.chains[i]
                    chain_means[i][j] = np.sum(chain.coldata[:, 0] * chain.coldata[:, j + 2]) / chain.norm
                    between_chain_var[j] += np.power(chain_means[i][j] - mean[j], 2)
                    in_chain_var[j] += np.sum(chain.coldata[:, 0] * np.power(chain.coldata[:, j + 2] - chain_means[i][j], 2))

                between_chain_var[j] /= (num_chains_used - 1)
                in_chain_var[j] /= norm
                label = self.parLabel(j)
                textFileHandle.write("%3i%13.5f  %s\n" % (j + 1, between_chain_var[j] / in_chain_var[j], label))

        nparam = self.paramNames.numNonDerived()
        covmat_dimension = nparam
        if (num_chains_used > 1) and (covmat_dimension > 0):
            # Assess convergence in the var(mean)/mean(var) in the worst eigenvalue
            # c.f. Brooks and Gelman 1997

            meanscov = np.zeros((nparam, nparam))
            cov = np.zeros((nparam, nparam))

            for j in range(nparam):
                for k in range(nparam):
                    meanscov[k, j] = 0.
                    cov[k, j] = 0.
                    for i in range(num_chains_used):
                        chain = self.chains[i]
                        cov[k, j] += np.sum(
                            chain.coldata[:, 0] *
                            (chain.coldata[:, j + 2] - chain_means[i][j]) *
                            (chain.coldata[:, k + 2] - chain_means[i][k]))
                        meanscov[k, j] += (chain_means[i][j] - mean[j]) * (chain_means[i][k] - mean[k])
                    meanscov[j, k] = meanscov[k, j]
                    cov[j, k] = cov[k, j]

            meanscov[:, :] /= (num_chains_used - 1)
            cov[:, :] /= norm

            invertible = np.isfinite(np.linalg.cond(cov))
            if (invertible):
                R = np.linalg.inv(np.linalg.cholesky(cov))
                D = np.linalg.eigvalsh(np.dot(R, meanscov).dot(R.T))
                GelmanRubin = max(np.real(D))

                textFileHandle.write("\n")
                textFileHandle.write("var(mean)/mean(var) for eigenvalues of covariance of means of orthonormalized parameters\n")
                for jj in range(nparam):
                    textFileHandle.write("%3i%13.5f\n" % (jj + 1, D[jj]))
                print " var(mean)/mean(var), remaining chains, worst e-value: R-1 = %13.5F\n" % GelmanRubin
            else:
                print 'WARNING: Gelman-Rubin covariance not invertible'

        textFileHandle.write("\n")

        if quick: return

        # Do tests for robustness under using splits of the samples
        # Return the rms ([change in upper/lower quantile]/[standard deviation])
        # when data split into 2, 3,.. sets
        textFileHandle.write("Split tests: rms_n([delta(upper/lower quantile)]/sd) n={2,3,4}:\n")
        textFileHandle.write("i.e. mean sample splitting change in the quantiles in units of the st. dev.\n")
        textFileHandle.write("\n")

        # Need these values here
        self.weights = np.hstack((chain.coldata[:, 0] for chain in self.chains))
        self.norm = np.sum(self.weights)

        raise MCException('Converge tests not updated yet (and very slow)')

        split_tests = {}
        nparam = self.paramNames.numParams()
        for j in range(nparam):
            if not isused[j]: continue
            coldata = np.hstack((chain.coldata[:, j + 2] for chain in self.chains))
            for endb in [0, 1]:
                for split_n in range(2, self.max_split_tests + 1):
                    frac = self.GetFractionIndices(self.weights, split_n)
                    split_tests[split_n] = 0.
                    confid = self.confidence(coldata, (1 - limfrac) / 2., endb == 0)
                    for i in range(split_n):
                        split_tests[split_n] = split_tests[split_n] + math.pow(self.confidence(coldata[frac[i]:frac[i + 1]], (1 - limfrac) / 2., endb == 0, start=frac[i], end=frac[i + 1]) - confid, 2)

                    split_tests[split_n] = math.sqrt(split_tests[split_n] / split_n / fullvar[j])
                if (endb == 0):
                    typestr = 'upper'
                else:
                    typestr = 'lower'

                textFileHandle.write("%3i" % (j + 1))
                for split_n in range(2, self.max_split_tests + 1):
                    textFileHandle.write("%9.4f" % (split_tests[split_n]))
                label = self.parLabel(j)
                textFileHandle.write("  %s %s\n" % (label, typestr))


        # Now do Raftery and Lewis method
        # See http://www.stat.washington.edu/tech.reports/raftery-lewis2.ps
        # Raw non-importance sampled chains only
        thin_fac = np.zeros(num_chains_used)

        tran = np.zeros((2, 2, 2), dtype=np.int)
        tran2 = np.zeros((2, 2), dtype=np.int)

        epsilon = 0.001

        if (np.all(np.abs(weights - np.round(np.where(weights > 0.6, weights, 0.6))) < 1e-4)):

            nburn = np.zeros(num_chains_used, dtype=np.int)
            markov_thin = np.zeros(num_chains_used, dtype=np.int)
            hardest = -1
            hardestend = 0
            for ix in range(num_chains_used):
                chain = self.chains[ix]
                thin_fac[ix] = int(round(np.max(chain.coldata[:, 0])))

                for j in range(covmat_dimension):
                    coldata = np.hstack((chain.coldata[:, j] for chain in self.chains))
                    if (self.force_twotail or not self.has_limits[j]):
                        for endb in [0, 1]:
                            # Get binary chain depending on whether above or below confidence value
                            u = self.confidence(chain.coldata[:, j], (1 - limfrac) / 2, endb == 0)
                            while(True):
                                thin_ix = self.thin_indices(thin_fac[ix])
                                thin_rows = len(thin_ix)
                                if (thin_rows < 2): break
                                binchain = np.ones(thin_rows)
                                indexes = np.where(coldata[thin_ix] >= u)
                                binchain[indexes] = 0

                                tran[:, :, :] = 0
                                # Estimate transitions probabilities for 2nd order process
                                for i in range(2, thin_rows):
                                    tran[binchain[i - 2]][binchain[i - 1]][binchain[i]] += 1

                                # Test whether 2nd order is better than Markov using BIC statistic
                                g2 = 0
                                for i1 in [0, 1]:
                                    for i2 in [0, 1]:
                                        for i3 in [0, 1]:
                                            if (tran[i1][i2][i3] <> 0):
                                                fitted = float(
                                                    (tran[i1][i2][0] + tran[i1][i2][1]) *
                                                    (tran[0][i2][i3] + tran[1][i2][i3])) \
                                                / float(tran[0][i2][0] + tran[0][i2][1] +
                                                        tran[1][i2][0] + tran[1][i2][1])
                                                focus = float(tran[i1][i2][i3])
                                                g2 += math.log(focus / fitted) * focus
                                g2 = g2 * 2

                                if (g2 - math.log(float(thin_rows - 2)) * 2 < 0): break
                                thin_fac[ix] += 1

                            # Get Markov transition probabilities for binary processes
                            if (np.sum(tran[:, 0, 1]) == 0 or np.sum(tran[:, 1, 0]) == 0):
                                thin_fac[ix] = 0
                                # goto 203

                            alpha = np.sum(tran[:, 0, 1]) / float(np.sum(tran[:, 0, 0]) + np.sum(tran[:, 0, 1]))
                            beta = np.sum(tran[:, 1, 0]) / float(np.sum(tran[:, 1, 0]) + np.sum(tran[:, 1, 1]))
                            probsum = alpha + beta
                            tmp1 = math.log(probsum * epsilon / max(alpha, beta)) / math.log(abs(1.0 - probsum))
                            if (int(tmp1 + 1) * thin_fac[ix] > nburn[ix]):
                                nburn[ix] = int(tmp1 + 1) * thin_fac[ix]
                                hardest = j
                                hardestend = endb

                markov_thin[ix] = thin_fac[ix]

                # Get thin factor to have independent samples rather than Markov
                hardest = max(hardest, 0)
                u = self.confidence(chain.coldata[:, hardest], (1 - limfrac) / 2, hardestend == 0)
                thin_fac[ix] += 1

                while(True):
                    thin_ix = self.thin_indices_chain(chain.coldata[:, 0], thin_fac[ix])
                    thin_rows = len(thin_ix)
                    if (thin_rows < 2): break
                    binchain = np.ones(thin_rows)

                    coldata = np.hstack((chain.coldata[:, hardest] for chain in self.chains))
                    coldata = coldata[thin_ix]
                    indexes = np.where(coldata >= u)
                    binchain[indexes] = 0

                    tran2[:, :] = 0
                    # Estimate transitions probabilities for 2nd order process
                    for i in range(1, thin_rows):
                        tran2[binchain[i - 1]][binchain[i]] += 1

                    # Test whether independence is better than Markov using BIC statistic
                    g2 = 0
                    for i1 in [0, 1]:
                        for i2 in [0, 1]:
                            if (tran2[i1][i2] <> 0):
                                fitted = float(
                                    (tran2[i1][0] + tran2[i1][1]) *
                                    (tran2[0][i2] + tran2[1][i2])) / float(thin_rows - 1)
                                focus = float(tran2[i1][i2])
                                if (fitted <= 0 or focus <= 0):
                                    print 'Raftery and Lewis estimator had problems'
                                    return
                                g2 += math.log(focus / fitted) * focus
                    g2 = g2 * 2

                    if (g2 - math.log(float(thin_rows - 1)) < 0): break

                    thin_fac[ix] += 1
# goto 203
                if (thin_rows < 2): thin_fac[ix] = 0

            textFileHandle.write("\n")
            textFileHandle.write("Raftery&Lewis statistics\n")
            textFileHandle.write("\n")
            textFileHandle.write("chain  markov_thin  indep_thin    nburn\n")

            # Computation of mean_mult
            self.mean_mult = norm / nrows

            for ix in range(num_chains_used):
                if (thin_fac[ix] == 0):
                    textFileHandle.write("%4i      Not enough samples\n" % ix)
                else:
                    textFileHandle.write("%4i%12i%12i%12i" % (
                            ix, markov_thin[ix], thin_fac[ix], nburn[ix]))

            if (not np.all(thin_fac != 0)):
                print 'RL: Not enough samples to estimate convergence stats'
            else:
                print 'RL: Thin for Markov: ', np.max(markov_thin)
                self.indep_thin = np.max(thin_fac)
                print 'RL: Thin for indep samples:  ', str(self.indep_thin)
                print 'RL: Estimated burn in steps: ', np.max(nburn), ' (', int(round(np.max(nburn) / self.mean_mult)), ' rows)'

            # Get correlation lengths
            textFileHandle.write("\n")
            textFileHandle.write("Parameter auto-correlations as function of step separation\n")
            textFileHandle.write("\n")
            if (self.corr_length_thin <> 0):
                autocorr_thin = self.corr_length_thin
            else:
                if (self.indep_thin == 0):
                    autocorr_thin = 20
                elif (self.indep_thin <= 30):
                    autocorr_thin = 5
                else:
                    autocorr_thin = 5 * (self.indep_thin / 30)

            thin_ix = self.thin_indices(autocorr_thin)
            thin_rows = len(thin_ix)
            maxoff = int(min(self.corr_length_steps, thin_rows / (autocorr_thin * num_chains_used)))

            corrs = np.zeros([maxoff, nparam])
            for off in range(maxoff):
                for i in range(off, thin_rows):
                    for j in range(nparam):
                        if not isused[j]: continue
                        coldata = np.hstack((chain.coldata[:, j] for chain in self.chains))
                        corrs[off][j] += (coldata[thin_ix[i]] - fullmean[j]) * (coldata[thin_ix[i - off]] - fullmean[j])
                for j in range(nparam):
                    if not isused[j]: continue
                    corrs[off][j] /= (thin_rows - off) / fullvar[j]

            if (maxoff > 0):
                for i in range(maxoff):
                    textFileHandle.write("%8i" % ((i + 1) * autocorr_thin))
                textFileHandle.write("\n")
                for j in range(nparam):
                    if (isused[j]):
                        label = self.parLabel(j)
                        textFileHandle.write("%3i" % j + 1)
                        for i in range(maxoff):
                            textFileHandle.write("%8.3f" % corrs[i][j])
                        textFileHandle.write("%s\n" % label)
        textFileHandle.close()


    def initParamRanges(self, j, paramConfid=None):

        # Return values
        smooth_1D, end_edge = 0, 0

        paramVec = self.samples[:, j]

        self.param_min[j] = np.min(paramVec)
        self.param_max[j] = np.max(paramVec)
        paramConfid = paramConfid or self.initParamConfidenceData(paramVec)
        self.range_min[j] = min(self.ND_limit_bot[1, j], self.confidence(paramConfid, 0.001, upper=False))
        self.range_max[j] = max(self.ND_limit_top[1, j], self.confidence(paramConfid, 0.001, upper=True))

        width = (self.range_max[j] - self.range_min[j]) / (self.num_bins + 1)
        if width == 0:
            return width, smooth_1D, end_edge

        logging.debug("Smooth scale ... ")
        if self.smooth_scale_1D <= 0:
            # Automatically set smoothing scale from rule of thumb for Gaussian, e.g. see
            # http://en.wikipedia.org/wiki/Kernel_density_estimation
            # 1/5 power is insensitive so just use v crude estimate of effective number
            opt_width = 1.06 / math.pow(max(1.0, self.numsamp / self.max_mult), 0.2) * self.sddev[j]
            smooth_1D = opt_width / width * abs(self.smooth_scale_1D)
            if smooth_1D < 0.5:
                print 'Warning: num_bins not large enough for optimal density - ' + self.parName(j)
            smooth_1D = max(1.0, smooth_1D)
        elif self.smooth_scale_1D < 1.0:
            smooth_1D = self.smooth_scale_1D * self.sddev[j] / width
            if smooth_1D < 1:
                print 'Warning: num_bins not large enough to well sample smoothed density - ' + self.parName(j)
        else:
            smooth_1D = self.smooth_scale_1D

        end_edge = int(round(smooth_1D * 2))

        logging.debug("Limits ... ")
        if self.has_limits_bot[j]:
            if ((self.range_min[j] - self.limmin[j] > (width * end_edge)) and
                 (self.param_min[j] - self.limmin[j] > (width * smooth_1D))):
                # long way from limit
                self.has_limits_bot[j] = False
            else:
                self.range_min[j] = self.limmin[j]

        if self.has_limits_top[j]:
            if ((self.limmax[j] - self.range_max[j] > (width * end_edge)) and
                (self.limmax[j] - self.param_max[j] > (width * smooth_1D))):
                self.has_limits_top[j] = False
            else:
                self.range_max[j] = self.limmax[j]
        self.has_limits[j] = self.has_limits_top[j] or self.has_limits_bot[j]

        if self.has_limits_top[j]:
            self.center[j] = self.range_max[j]
            if self.has_limits_bot[j]:
                width = (self.range_max[j] - self.range_min[j]) / (self.num_bins + 1)  # Feb15
        else:
            self.center[j] = self.range_min[j]

        self.ix_min[j] = int(round((self.range_min[j] - self.center[j]) / width))
        self.ix_max[j] = int(round((self.range_max[j] - self.center[j]) / width))

        if not self.has_limits_bot[j]: self.ix_min[j] -= end_edge
        if not self.has_limits_top[j]: self.ix_max[j] += end_edge

        return width, smooth_1D, end_edge


    def Get1DDensity(self, j, writeDataToFile=False, get_density=False, paramConfid=None):

        logging.debug("1D density for %s (%i)" % (self.parName(j), j))

        paramVec = self.samples[:, j]
        width, smooth_1D, end_edge = self.initParamRanges(j, paramConfid)

        if width == 0: raise MCException("width is 0 in Get1DDensity")

        # In f90, binsraw(ix_min(j):ix_max(j))
        binsraw = np.zeros(self.ix_max[j] - self.ix_min[j] + 1)

        fine_fac = 10
        winw = int(round(2.5 * fine_fac * smooth_1D))
        fine_edge = winw + fine_fac * end_edge
        fine_width = width / fine_fac

        imin = int(round((self.param_min[j] - self.center[j]) / fine_width))
        imax = int(round((self.param_max[j] - self.center[j]) / fine_width))
        imin = min(self.ix_min[j] * fine_fac, imin)
        imax = max(self.ix_max[j] * fine_fac, imax)

        # In f90, finebins(imin-fine_edge:imax+fine_edge)
        fine_min = imin - fine_edge
        finebins = np.zeros((imax + fine_edge) - fine_min + 1)

        if self.plot_meanlikes:
            # In f90, finebinlikes(imin-fine_edge:imax+fine_edge)
            finebinlikes = np.zeros((imax + fine_edge) - fine_min + 1)

        # vectorized code
        ix2 = np.round((paramVec - self.center[j]) / width).astype(np.int)
        minix = np.min(ix2)
        counts = np.bincount(ix2 - minix, weights=self.weights)
        if minix < self.ix_min[j]:
            off = self.ix_min[j] - minix
            alen = min(len(counts) - off, len(binsraw))
            binsraw[:alen] = counts[off:off + alen]
        else:
            off = minix - self.ix_min[j]
            alen = min(len(binsraw) - off, len(counts))
            binsraw[off:off + alen] = counts[:alen]

        ix2 = np.round((paramVec - self.center[j]) / fine_width).astype(np.int)
        minix = np.min(ix2)
        counts = np.bincount(ix2 - minix, weights=self.weights)
        finebins[minix - fine_min:minix - fine_min + len(counts)] = counts
        if self.plot_meanlikes:
            if self.mean_loglikes:
                w = self.weights * self.loglikes
            else:
                w = self.weights * np.exp(self.meanlike - self.loglikes)
            counts = np.bincount(ix2 - minix, weights=w)
            finebinlikes[minix - fine_min:minix - fine_min + len(counts)] = counts

        if self.ix_min[j] <> self.ix_max[j]:
            # account for underweighting near edges
            if (not self.has_limits_bot[j] and binsraw[end_edge - 1] == 0 and
                 binsraw[end_edge] > np.max(binsraw) / 15):
                self.EdgeWarning(j)
            if (not self.has_limits_top[j] and binsraw[self.ix_max[j] - end_edge + 1 - self.ix_min[j]] == 0 and
                 binsraw[self.ix_max[j] - end_edge + 1 - self.ix_min[j]] > np.max(binsraw) / 15):
                self.EdgeWarning(j)

        # High resolution density (sampled many times per smoothing scale)
        if self.has_limits_bot[j]: imin = self.ix_min[j] * fine_fac
        if self.has_limits_top[j]: imax = self.ix_max[j] * fine_fac

        Kernel = Kernel1D(winw, fine_fac * smooth_1D)
        fineused = finebins[imin - winw - fine_min:imax + winw + 1 - fine_min]
        conv = np.convolve(fineused, Kernel.Win, 'valid')

        density1D = Density1D(self.center[j] + imin * fine_width, imax - imin + 1, fine_width, P=conv)
#        density1D.X = self.center[j] + np.arange(imin, imax + 1) * fine_width

        if self.has_limits_bot[j] or self.has_limits_top[j]:
            # correct for cuts allowing for normalization over window
            prior_mask = np.ones(imax - imin + 2 * winw + 1)
            if self.has_limits_bot[j]:
                prior_mask[ winw ] = 0.5
                prior_mask[ : winw ] = 0
            if self.has_limits_top[j]:
                prior_mask[ imax - imin + winw ] = 0.5
                prior_mask[ imax - imin + winw + 1 : ] = 0
            a0 = np.convolve(prior_mask, Kernel.Win, 'valid')
            ix = np.nonzero(a0 * density1D.P)
            a0 = a0[ix]
            normed = density1D.P[ix] / a0
            if self.boundary_correction_method == 1:
                # linear boundary kernel, e.g. Jones 1993, Jones and Foster 1996
                #  www3.stat.sinica.edu.tw/statistica/oldpdf/A6n414.pdf after Eq 1b, expressed for general prior mask
                xWin = Kernel.Win * Kernel.x
                a1 = np.convolve(prior_mask, xWin, 'valid')[ix]
                a2 = np.convolve(prior_mask, xWin * Kernel.x, 'valid')[ix]
                xP = np.convolve(fineused, xWin, 'valid')[ix]
                corrected = (density1D.P[ix] * a2 - xP * a1) / (a0 * a2 - a1 ** 2)
                density1D.P[ix] = normed * np.exp(np.minimum(corrected / normed, 4) - 1)
            elif self.boundary_correction_method == 0:
                density1D.P[ix] = normed
            else: raise SettingException('Unknown boundary_correction_method (expected 0 or 1)')

        maxbin = np.max(density1D.P)
        if maxbin == 0:
            raise MCException('no samples in bin, param: ' + self.parName(j))
        density1D.P /= maxbin

        density1D.InitSpline()
        self.density1D[self.parName(j)] = density1D

        if get_density: return density1D

        logZero = 1e30
        if not self.no_plots:
            bincounts = density1D.P[self.ix_min[j] * fine_fac - imin:self.ix_max[j] * fine_fac - imin + 1:fine_fac]
            if self.plot_meanlikes:
                # In f90, binlikes(ix_min(j):ix_max(j))
                rawbins = conv[self.ix_min[j] * fine_fac - imin:self.ix_max[j] * fine_fac - imin + 1:fine_fac]
                binlikes = np.zeros(self.ix_max[j] - self.ix_min[j] + 1)
                if self.mean_loglikes: binlikes[:] = logZero

                # Output values for plots
                for ix2 in range(self.ix_min[j], self.ix_max[j] + 1):
                    if rawbins[ix2 - self.ix_min[j]] > 0:
                        istart, iend = (ix2 * fine_fac - winw) - fine_min, (ix2 * fine_fac + winw + 1) - fine_min
                        binlikes[ix2 - self.ix_min[j]] = np.dot(Kernel.Win, finebinlikes[istart:iend]) / rawbins[ix2 - self.ix_min[j]]

                if self.mean_loglikes:
                    maxbin = min(binlikes)
                    binlikes = np.where((binlikes - maxbin) < 30, np.exp(-(binlikes - maxbin)), 0)

            x = self.center[j] + np.arange(self.ix_min[j], self.ix_max[j] + 1) * width
            if self.plot_meanlikes:
                maxbin = np.max(binlikes)
                likes = binlikes / maxbin
            else:
                likes = None
            if writeDataToFile:
                logging.debug("Write data to file ...")

                fname = self.rootname + "_p_" + self.parName(j)
                filename = os.path.join(self.plot_data_dir, fname + ".dat")
                with open(filename, 'w') as f:
                    for xval, binval in zip(x, bincounts):
                        f.write("%16.7E%16.7E\n" % (xval, binval))

                if self.plot_meanlikes:
                    filename_like = os.path.join(self.plot_data_dir, fname + ".likes")
                    with open(filename_like, 'w') as f:
                        for xval, binval in zip(x, likes):
                            f.write("%16.7E%16.7E\n" % (xval, binval))
            else:
                return x, bincounts, likes

        return None, None, None


    def Get2DPlotData(self, j, j2, writeDataToFile=False, num_plot_contours=None):
        """
        Get 2D plot data.
        """
        fine_fac_base = 5
        has_prior = self.has_limits[j] or self.has_limits[j2]

        corr = self.corrmatrix[j2][j]
        # keep things simple unless obvious degeneracy
        if abs(corr) < 0.1: corr = 0.
        corr = max(-self.max_corr_2D, corr)
        corr = min(self.max_corr_2D, corr)

        # for tight degeneracies increase bin density
        nbin2D = min(4 * self.num_bins_2D, int(round(self.num_bins_2D / (1 - abs(corr)))))

        widthx = (self.range_max[j] - self.range_min[j]) / (nbin2D + 1)
        widthy = (self.range_max[j2] - self.range_min[j2]) / (nbin2D + 1)
        smooth_scale = (self.smooth_scale_2D * nbin2D) / self.num_bins_2D
        fine_fac = max(2, int(round(fine_fac_base / smooth_scale)))

        ixmin = int(round((self.range_min[j] - self.center[j]) / widthx))
        ixmax = int(round((self.range_max[j] - self.center[j]) / widthx))

        iymin = int(round((self.range_min[j2] - self.center[j2]) / widthy))
        iymax = int(round((self.range_max[j2] - self.center[j2]) / widthy))

        if not self.has_limits_bot[j]: ixmin -= 1
        if not self.has_limits_bot[j2]: iymin -= 1
        if not self.has_limits_top[j]: ixmax += 1
        if not self.has_limits_top[j2]: iymax += 1

        finewidthx = widthx / fine_fac
        finewidthy = widthy / fine_fac
        winw = int(round(2 * fine_fac * smooth_scale))

        ix1s = np.round((self.samples[:, j] - self.center[j]) / finewidthx).astype(np.int)
        ix2s = np.round((self.samples[:, j2] - self.center[j2]) / finewidthy).astype(np.int)
        imin = min(ixmin * fine_fac, np.min(ix1s)) - winw
        imax = max(ixmax * fine_fac, np.max(ix1s)) + winw
        jmin = min(iymin * fine_fac, np.min(ix2s)) - winw
        jmax = max(iymax * fine_fac, np.max(ix2s)) + winw

        flatix = (ix1s - imin) + (ix2s - jmin) * (imax - imin + 1)
        finebins = np.bincount(flatix, weights=self.weights,
                    minlength=(jmax - jmin + 1) * (imax - imin + 1)).reshape((jmax - jmin + 1 , imax - imin + 1))
# equivalent to this slow method
#        finebins = np.zeros((jmax - jmin + 1, imax - imin + 1))
#        for i, (ix1, ix2) in enumerate(zip(ix1s - imin, ix2s - jmin)):
#            finebins[ix2][ix1] += self.weights[i]

        if self.shade_meanlikes:
            likeweights = self.weights * np.exp(self.meanlike - self.loglikes)
            finebinlikes = np.bincount(flatix, weights=likeweights,
                    minlength=(jmax - jmin + 1) * (imax - imin + 1)).reshape((jmax - jmin + 1 , imax - imin + 1))

        # In f90, Win(-winw:winw,-winw:winw)
        Win = np.empty(((2 * winw) + 1, (2 * winw) + 1))
        indexes = np.arange(-winw, winw + 1)
        signorm = 2 * (fine_fac * smooth_scale) ** 2 * (1 - corr ** 2)
        for ix1 in indexes:
            for ix2 in indexes:
                Win[ix2 - (-winw)][ix1 - (-winw)] = np.exp(-(ix1 ** 2 + ix2 ** 2 - 2 * corr * ix1 * ix2) / signorm)
        Win /= np.sum(Win)

        fineused = finebins[iymin * fine_fac - winw - jmin:iymax * fine_fac + winw + 1 - jmin,
                            ixmin * fine_fac - winw - imin:ixmax * fine_fac + winw + 1 - imin]
        bins2D = fftconvolve(fineused, Win, 'valid')

        if self.shade_meanlikes:
            bin2Dlikes = fftconvolve(finebinlikes[iymin * fine_fac - winw - jmin:iymax * fine_fac + winw + 1 - jmin,
                                                  ixmin * fine_fac - winw - imin:ixmax * fine_fac + winw + 1 - imin],
                                                  Win, 'valid')
            del finebinlikes
            mx = 1e-4 * np.max(bins2D)
            bin2Dlikes[bins2D > mx] /= bins2D[bins2D > mx]
            bin2Dlikes[bins2D <= mx] = 0
        else:
            bin2Dlikes = None

        if has_prior:
            # Correct for edge effects
            prior_mask = np.ones((jmax - jmin + 1, imax - imin + 1))
            if self.has_limits_bot[j]:
                prior_mask[:, (ixmin * fine_fac) - imin] /= 2
                prior_mask[:, :(ixmin * fine_fac) - imin] = 0
            if self.has_limits_top[j]:
                prior_mask[:, (ixmax * fine_fac) - imin] /= 2
                prior_mask[:, (ixmax * fine_fac) + 1 - imin:] = 0

            if self.has_limits_bot[j2]:
                prior_mask[(iymin * fine_fac) - jmin, :] /= 2
                prior_mask[:(iymin * fine_fac) - jmin, :] = 0
            if self.has_limits_top[j2]:
                prior_mask[(iymax * fine_fac) - jmin, :] /= 2
                prior_mask[(iymax * fine_fac) + 1 - jmin:, :] = 0

            priorused = prior_mask[iymin * fine_fac - winw - jmin:iymax * fine_fac + winw + 1 - jmin,
                                   ixmin * fine_fac - winw - imin:ixmax * fine_fac + winw + 1 - imin]
            a00 = fftconvolve(priorused, Win, 'valid')
            ix = np.nonzero(a00 * bins2D)
            a00 = a00[ix]
            normed = bins2D[ix] / a00
            if self.boundary_correction_method == 1:
                # linear boundary correction
                y = np.empty(Win.shape)
                for i in range(Win.shape[0]):
                    y[:, i] = indexes
                winx = Win * indexes
                winy = Win * y
                a10 = fftconvolve(priorused, winx, 'valid')[ix]
                a01 = fftconvolve(priorused, winy, 'valid')[ix]
                a20 = fftconvolve(priorused, winx * indexes, 'valid')[ix]
                a02 = fftconvolve(priorused, winy * y, 'valid')[ix]
                a11 = fftconvolve(priorused, winy * indexes, 'valid')[ix]
                xP = fftconvolve(fineused, winx, 'valid')[ix]
                yP = fftconvolve(fineused, winy, 'valid')[ix]
                denom = (a20 * a01 ** 2 + a10 ** 2 * a02 - a00 * a02 * a20 + a11 ** 2 * a00 - 2 * a01 * a10 * a11)
                A = a11 ** 2 - a02 * a20
                Ax = a10 * a02 - a01 * a11
                Ay = a01 * a20 - a10 * a11
                corrected = (bins2D[ix] * A + xP * Ax + yP * Ay) / denom
                bins2D[ix] = normed * np.exp(np.minimum(corrected / normed, 4) - 1)
            elif self.boundary_correction_method == 0:
                # simple boundary correction by normalization
                bins2D[ix] = normed
            else: raise SettingException('unknown boundary_correction_method (expected 0 or 1)')

        mx = np.max(bins2D)
        bins2D = bins2D / mx
        missing_norm = (np.sum(finebins) - np.sum(fineused[winw:-winw, winw:-winw])) / mx

        ncontours = len(self.contours)
        if num_plot_contours: ncontours = min(num_plot_contours, ncontours)
        contours = self.contours[:ncontours]

        # Get contour containing contours(:) of the probability
        contour_levels = get2DContourLevels(bins2D, contours, missing_norm=missing_norm)
        if self.smooth_correct_contours:
            # estimate dependence on smoothing width and roughly correct
            varscale = 2
            Winscale = Win ** varscale
            Winscale /= np.sum(Winscale)
            bins2Droot = fftconvolve(fineused, Winscale, 'valid')
            if has_prior:
                simple_normed = bins2D.copy()
                simple_normed[ix] = normed / mx
                aw00 = fftconvolve(priorused, Winscale, 'valid')
                ix = np.nonzero(aw00 * bins2Droot)
                bins2Droot[ix] = bins2Droot[ix] / aw00[ix]
            else:
                simple_normed = bins2D
            bins2Droot /= np.max(bins2Droot)

            varrat2 = (fine_fac * smooth_scale) ** 2 / (self.fullcov[j2, j2] / finewidthy ** 2)
            varrat = (fine_fac * smooth_scale) ** 2 / (self.fullcov[j, j] / finewidthx ** 2)
            varrat = max(0.02, varrat , varrat2)  # probably better to get from estimate of curvature, but only used for approximate linearization in the smoothing width
            if False:
                # this compres with analytic result for full gaussian assuming varrat is ratio of variances
                p = 1 - np.array(self.contours[:num_plot_contours])
                rat = np.exp(np.log(p) / (1 + varrat)) / p
                contour_levels *= rat
            else:
                for i, contour_level  in enumerate(contour_levels):
                    orig = np.log(np.sum(simple_normed[simple_normed < contour_level]))
                    dpdW = (np.log(np.sum(simple_normed[bins2Droot < contour_level])) - orig) / (1 / (1 + varrat / varscale) - 1 / (1 + varrat))
                    dpdL = (np.log(np.sum(simple_normed[simple_normed < contour_level * 1.1])) - orig) / 0.1
                    if dpdL > 0:
                        fac = dpdW / dpdL * (1 - 1 / (1 + varrat))
                        contour_levels[i] += max(fac * contour_level, 0)

        if True:
            bins2D = bins2D[::fine_fac, ::fine_fac ]
            if self.shade_meanlikes:
                bin2Dlikes = bin2Dlikes[::fine_fac, ::fine_fac ]
            y = self.center[j2] + np.arange(iymin, iymax + 1) * widthy
            x = self.center[j] + np.arange(ixmin, ixmax + 1) * widthx
        else:
            y = self.center[j2] + np.arange(iymin * fine_fac, iymax * fine_fac + 1) * finewidthy
            x = self.center[j] + np.arange(ixmin * fine_fac, ixmax * fine_fac + 1) * finewidthx


        bins2D[bins2D < 1e-30] = 0
        if self.shade_meanlikes:
            bin2Dlikes /= np.max(bin2Dlikes)

#        print 'time 2D:', time.time() - now
        if writeDataToFile:
            # note store things in confusing traspose form
            name = self.parName(j)
            name2 = self.parName(j2)
            plotfile = self.rootname + "_2D_%s_%s" % (name, name2)
            filename = os.path.join(self.plot_data_dir, plotfile)
            np.savetxt(filename, bins2D.T, "%16.7E")
            np.savetxt(filename + "_y", x, "%16.7E")
            np.savetxt(filename + "_x", y, "%16.7E")
            np.savetxt(filename + "_cont-test", np.atleast_2d(contour_levels), "%16.7E")
            if self.shade_meanlikes:
                np.savetxt(filename + "_likes", bin2Dlikes.T , "%16.7E")
        else:
            return bins2D, bin2Dlikes, contour_levels, x, y


    def EdgeWarning(self, i):
        name = self.parName(i)
        print 'Warning: sharp edge in parameter %s - check limits[%s] or limits%i' % (name, name, i + 1)

    def GetChainLikeSummary(self, toStdOut=False):
        text = ""
        maxlike = np.max(self.loglikes)
        text += "Best fit sample -log(Like) = %f\n" % maxlike
        if ((self.loglikes[self.numrows - 1] - maxlike) < 30):
            self.meanlike = np.log(np.sum(np.exp(self.loglikes - maxlike) * self.weights) / self.norm) + maxlike
            text += "Ln(mean 1/like) = %f\n" % (self.meanlike)
        self.meanlike = np.sum(self.loglikes * self.weights) / self.norm
        text += "mean(-Ln(like)) = %f\n" % (self.meanlike)
        self.meanlike = -np.log(np.sum(np.exp(-(self.loglikes - maxlike)) * self.weights) / self.norm) + maxlike
        text += "-Ln(mean like)  = %f\n" % (self.meanlike)
        if toStdOut:
            print text
        else:
            return text


    def ComputeStats(self):
        """
        Compute mean and std dev.
        """
        nparam = self.samples.shape[1]
        self.means = np.zeros(nparam)
        self.sddev = np.zeros(nparam)
        for i in range(nparam): self.means[i] = self.mean(self.samples[:, i])
        for i in range(nparam): self.sddev[i] = self.std(self.samples[:, i])

    def Init1DDensity(self):
        self.done_1Dbins = False
        self.density1D = dict()
        nparam = self.samples.shape[1]
        self.param_min = np.zeros(nparam)
        self.param_max = np.zeros(nparam)
        self.range_min = np.zeros(nparam)
        self.range_max = np.zeros(nparam)
        self.center = np.zeros(nparam)
        self.ix_min = np.zeros(nparam, dtype=np.int)
        self.ix_max = np.zeros(nparam, dtype=np.int)
        #
        self.marge_limits_bot = np.ndarray([self.num_contours, nparam], dtype=bool)
        self.marge_limits_top = np.ndarray([self.num_contours, nparam], dtype=bool)

    def GetConfidenceRegion(self):
                # Sort data in order of likelihood of points
        indexes = self.loglikes.argsort()
        cumsum = np.cumsum(self.weights[indexes])
        self.ND_cont1 = np.searchsorted(cumsum, self.norm * self.contours[0])
        self.ND_cont2 = np.searchsorted(cumsum, self.norm * self.contours[1])

        self.ND_limit_top = np.empty((2, self.num_vars))
        self.ND_limit_bot = np.empty((2, self.num_vars))
        for j in range(self.num_vars):
            region1 = self.samples[indexes[:self.ND_cont1], j]
            region2 = self.samples[indexes[:self.ND_cont2], j]
            self.ND_limit_bot[0, j] = np.min(region1)
            self.ND_limit_bot[1, j] = np.min(region2)
            self.ND_limit_top[0, j] = np.max(region1)
            self.ND_limit_top[1, j] = np.max(region2)

    def ReadRanges(self):
        ranges_file = self.root + '.ranges'
        if os.path.isfile(ranges_file):
            self.ranges = Ranges(ranges_file)
        else:
            self.ranges = Ranges()
#            print "No file %s" % ranges_file

    def getBounds(self):
        upper = dict()
        lower = dict()
        for i, name in enumerate(self.paramNames.list()):
            if self.has_limits_bot[i]:
                lower[name] = self.limmin[i]
            if self.has_limits_top[i]:
                upper[name] = self.limmax[i]
        return lower, upper

    def writeBounds(self, filename):
        lower, upper = self.getBounds()
        with open(filename, 'w') as f:
            for i, name in enumerate(self.paramNames.list()):
                if (self.has_limits_bot[i] or self.has_limits_top[i]):
                    valMin = lower.get(name)
                    if (valMin):
                        lim1 = "%15.7E" % valMin
                    else:
                        lim1 = "    N"
                    valMax = upper.get(name)
                    if (valMax is not None):
                        lim2 = "%15.7E" % valMax
                    else:
                        lim2 = "    N"
                    f.write("%22s%17s%17s\n" % (name, lim1, lim2))


    def getMargeStats(self):
        self.Do1DBins()
        m = ResultObjs.margeStats()
        m.hasBestFit = False
        m.limits = self.contours
        m.names = self.paramNames.names
        return m

    def saveMargeStats(self):
        self.getMargeStats().saveAsText(self.rootdirname + '.margestats')


    # Pass weights as parameter (for chain only)
    def thin_indices_chain(self, weights, factor):
        thin_ix = []
        tot = 0
        i = 0
        numrows = len(weights)
        mult = weights[i]
        if abs(round(mult) - mult) > 1e-4:
                raise SettingException('Can only thin with integer weights')

        while i < numrows:
            if (mult + tot < factor):
                tot += mult
                i += 1
                if i < numrows: mult = weights[i]
            else:
                thin_ix.append(i)
                if mult == factor - tot:
                    i += 1
                    if i < numrows: mult = weights[i]
                else:
                    mult -= (factor - tot)
                tot = 0
        return thin_ix


    def Do1DBins(self, max_frac_twotail=None, writeDataToFile=False):
        if self.done_1Dbins: return
        if max_frac_twotail is None:
            max_frac_twotail = self.max_frac_twotail

        for j in range(self.num_vars):
            paramConfid = self.initParamConfidenceData(self.samples[:, j])
            self.Get1DDensity(j, writeDataToFile, get_density=not writeDataToFile, paramConfid=paramConfid)
            self.setMargeLimits(j, max_frac_twotail, paramConfid=paramConfid)
        self.done_1Dbins = True



    def setMargeLimits(self, j, max_frac_twotail, paramConfid=None):

        # Get limits, one or two tail depending on whether posterior
        # goes to zero at the limits or not
        density1D = self.density1D[self.parName(j)]
        par = self.paramNames.names[j]
        par.mean = self.means[j]
        par.err = self.sddev[j]
        par.limits = []
        interpGrid = None
        paramConfid = paramConfid or self.initParamConfidenceData(self.samples[:, j])
        for ix1 in range(self.num_contours):

            self.marge_limits_bot[ix1][j] = self.has_limits_bot[j] and \
                (not self.force_twotail) and (density1D.P[0] > max_frac_twotail[ix1])
            self.marge_limits_top[ix1][j] = self.has_limits_top[j] and \
                (not self.force_twotail) and (density1D.P[-1] > max_frac_twotail[ix1])

            if not self.marge_limits_bot[ix1][j] or not self.marge_limits_top[ix1][j]:
                # give limit
                if not interpGrid: interpGrid = density1D.initLimitGrids()
                tail_limit_bot, tail_limit_top, marge_bot, marge_top = density1D.Limits(self.contours[ix1], interpGrid)
                self.marge_limits_bot[ix1][j] = marge_bot
                self.marge_limits_top[ix1][j] = marge_top

                limfrac = 1 - self.contours[ix1]

                if (self.marge_limits_bot[ix1][j]):
                    # fix to end of prior range
                    tail_limit_bot = self.range_min[j]
                elif (self.marge_limits_top[ix1][j]):
                    # 1 tail limit
                    tail_limit_bot = self.confidence(paramConfid, limfrac, upper=False)
                else:
                    # 2 tail limit
                    tail_confid_bot = self.confidence(paramConfid, limfrac / 2, upper=False)

                if self.marge_limits_top[ix1][j]:
                    tail_limit_top = self.range_max[j]
                elif self.marge_limits_bot[ix1][j]:
                    tail_limit_top = self.confidence(paramConfid, limfrac, upper=True)
                else:
                    tail_confid_top = self.confidence(paramConfid, limfrac / 2, upper=True)

                if not self.marge_limits_bot[ix1][j] and  not self.marge_limits_top[ix1][j]:
                    # Two tail, check if limits are at very different density
                    if (math.fabs(density1D.Prob(tail_confid_top) -
                                   density1D.Prob(tail_confid_bot))
                        < self.credible_interval_threshold):
                        tail_limit_top = tail_confid_top
                        tail_limit_bot = tail_confid_bot

                lim = [tail_limit_bot, tail_limit_top]
            else:
                # no limit
                lim = [self.range_min[j], self.range_max[j]]

            if self.marge_limits_bot[ix1][j] and self.marge_limits_top[ix1][j]:
                tag = 'none'
            elif self.marge_limits_bot[ix1][j]:
                tag = '>'
            elif self.marge_limits_top[ix1][j]:
                tag = '<'
            else:
                tag = 'two'
            par.limits.append(ResultObjs.paramLimit(lim, tag))


    def GetCust2DPlots(self, num_cust2D_plots):

        try_t = 1e5
        x, y = 0, 0
        cust2DPlots = []
        for j in range(num_cust2D_plots):
            try_b = -1e5
            for ix1 in range(self.num_vars):
                for ix2 in range(ix1 + 1, self.num_vars):
                    if (abs(self.corrmatrix[ix1][ix2]) < try_t) and \
                            (abs(self.corrmatrix[ix1][ix2]) > try_b) :
                        try_b = abs(self.corrmatrix[ix1][ix2])
                        x, y = ix1, ix2
            if (try_b == -1e5):
                num_cust2D_plots = j - 1
                break
            try_t = try_b
            cust2DPlots.append([self.parName(x), self.parName(y)])

        return cust2DPlots


    # Write functions

    def WriteScriptPlots1D(self, filename, plotparams=None, ext=None):
        ext = ext or self.plot_output
        textFileHandle = open(filename, 'w')
        textInit = WritePlotFileInit()
        textFileHandle.write(textInit % (
                self.plot_data_dir, self.subplot_size_inch,
                self.out_dir, self.rootname))
        text = 'markers=' + str(self.markers) + '\n'
        if plotparams:
            text += 'g.plots_1d(roots,[' + ",".join(['\'' + par + '\'' for par in plotparams]) + '], markers=markers)'
        else:
            text += 'g.plots_1d(roots, markers=markers)\n'
        textFileHandle.write(text)
        textExport = WritePlotFileExport()
        fname = self.rootname + '.' + ext
        textFileHandle.write(textExport % (fname))
        textFileHandle.close()


    def WriteScriptPlots2D(self, filename, plot_2D_param, cust2DPlots, plots_only, ext=None):
        ext = ext or self.plot_output
        self.done2D = np.ndarray([self.num_vars, self.num_vars], dtype=bool)
        self.done2D[:, :] = False

        textFileHandle = open(filename, 'w')
        textInit = WritePlotFileInit()
        textFileHandle.write(textInit % (
                self.plot_data_dir, self.subplot_size_inch2,
                self.out_dir, self.rootname))
        textFileHandle.write('pairs=[]\n')
        plot_num = 0
        if cust2DPlots: cuts = [par1 + '__' + par2 for par1, par2 in cust2DPlots]
        for j, par1 in enumerate(self.paramNames.list()):
            if (self.ix_min[j] <> self.ix_max[j]):
                if plot_2D_param or cust2DPlots:
                    if par1 == plot_2D_param: continue
                    j2min = 0
                else:
                    j2min = j + 1

                for j2 in range(j2min, self.num_vars):
                    par2 = self.parName(j2)
                    if (self.ix_min[j2] <> self.ix_max[j2]):
                        if plot_2D_param and par2 <> plot_2D_param: continue
                        if cust2DPlots and (par1 + '__' + par2) not in cuts: continue
                        plot_num += 1
                        self.done2D[j][j2] = True
                        if not plots_only: self.Get2DPlotData(j, j2, writeDataToFile=True)
                        textFileHandle.write("pairs.append(['%s','%s'])\n" % (par1, par2))
        textFileHandle.write('g.plots_2d(roots,param_pairs=pairs)\n')
        textExport = WritePlotFileExport()
        fname = self.rootname + '_2D.' + ext
        textFileHandle.write(textExport % (fname))
        textFileHandle.close()
        print 'Produced ', plot_num, ' 2D plots'


    def WriteScriptPlotsTri(self, filename, triangle_params, ext=None):
        ext = ext or self.plot_output
        textFileHandle = open(filename, 'w')
        textInit = WritePlotFileInit()
        textFileHandle.write(textInit % (
                self.plot_data_dir, self.subplot_size_inch,
                self.out_dir, self.rootname))
        text = 'g.triangle_plot(roots, %s)\n' % triangle_params
        textFileHandle.write(text)
        textExport = WritePlotFileExport()
        fname = self.rootname + '_tri.' + ext
        textFileHandle.write(textExport % (fname))
        textFileHandle.close()


    def WriteScriptPlots3D(self, filename, plot_3D, ext=None):
        ext = ext or self.plot_output
        textFileHandle = open(filename, 'w')
        textInit = WritePlotFileInit()
        textFileHandle.write(textInit % (
                self.plot_data_dir, self.subplot_size_inch3,
                self.out_dir, self.rootname))
        textFileHandle.write('sets=[]\n')
        text = ""
        for v1, v2, v3 in plot_3D:
            text += "sets.append(['%s','%s','%s'])\n" % (v1, v2, v3)
        text += 'g.plots_3d(roots,sets)\n'
        textFileHandle.write(text)
        fname = self.rootname + '_3D.' + ext
        textExport = WritePlotFileExport()
        textFileHandle.write(textExport % (fname))
        textFileHandle.close()


    def WriteGlobalLikelihood(self, filename):
        bestfit_ix = 0  # Since we have sorted the lines
        textFileHandle = open(filename, 'w')
        textInit = self.GetChainLikeSummary(toStdOut=False)
        textFileHandle.write(textInit)
        textFileHandle.write("\n")
        textFileHandle.write('param  bestfit        lower1         upper1         lower2         upper2\n')
        for j in range(self.num_vars):
            best = self.samples[bestfit_ix][j]
            label = self.parLabel(j)
            textFileHandle.write('%5i%15.7E%15.7E%15.7E%15.7E%15.7E   %s\n' % (j + 1, best,
                self.ND_limit_bot[0, j], self.ND_limit_top[0, j], self.ND_limit_bot[1, j], self. ND_limit_top[1, j], label))
        textFileHandle.close()


def WritePlotFileInit():
    text = """import GetDistPlots, os
g=GetDistPlots.GetDistPlotter(plot_data='%s')
g.settings.setWithSubplotSize(%f)
outdir='%s'
roots=['%s']
"""
    return text

def WritePlotFileExport():
    text = "g.export(os.path.join(outdir,'%s'))\n"
    return text


# ==============================================================================

# Usefull functions

def GetChainRootFiles(rootdir):
    pattern = os.path.join(rootdir, '*.paramnames')
    files = [os.path.splitext(f)[0] for f in glob.glob(pattern)]
    files.sort()
    return files

def GetRootFileName(rootdir):
    rootFileName = ""
    pattern = os.path.join(rootdir, '*_*.txt')
    chain_files = glob.glob(pattern)
    chain_files.sort()
    if chain_files:
        chain_file0 = chain_files[0]
        rindex = chain_file0.rindex('_')
        rootFileName = chain_file0[:rindex]
    return rootFileName


# ==============================================================================
