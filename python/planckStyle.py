from __future__ import absolute_import
import os
import copy
import getdist
from getdist import types, plots, IniFile
from getdist.types import BestFit
from matplotlib import rcParams, rc
from matplotlib import pyplot as plt
from paramgrid import batchjob
import numpy as np

# common setup for matplotlib
params = {'backend': 'pdf',
          'axes.labelsize': 9,
          'font.size': 8,
          'legend.fontsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'ytick.major.pad': 4,
          'xtick.major.pad': 4,
          'text.usetex': True,
          'font.family': 'sans-serif',
          # free font similar to Helvetica
          'font.sans-serif': 'FreeSans'}

sfmath = os.path.dirname(os.path.abspath(__file__)) + os.sep + 'sfmath'


# use of Sans Serif also in math mode

def setRc():
    # doesn't seem to work on all systems
    # rc('text.latex', preamble=r'\usepackage{silence}\WarningFilter{latex}{You have requested package}\usepackage{' + sfmath.replace(os.sep, '/') + '}')
    rc('text.latex', preamble=r'\usepackage{' + sfmath.replace(os.sep, '/') + '}')
    rcParams.update(params)


setRc()

if False:
    non_final = True
    version = 'CamSpec v910HM'
    defdata_root = 'CamSpecHM'
else:
    non_final = False
    version = 'clik10.2'
    defdata_root = 'plikHM'

datalabel = dict()
defdata_TT = defdata_root + '_TT_lowl_lowE'
datalabel[defdata_TT] = r'\textit{Planck} TT$+$lowE'
defdata_TE = defdata_root + '_TE_lowE'
datalabel[defdata_TE] = r'\textit{Planck} TE$+$lowE'
defdata_EE = defdata_root + '_EE_lowE'
datalabel[defdata_EE] = r'\textit{Planck} EE$+$lowE'
defdata_TE_TEB = defdata_root + '_TE_lowl_lowE'
datalabel[defdata_TE_TEB] = r'\textit{Planck} TE+$lowT,P'
defdata_EE_TEB = defdata_root + '_EE_lowl_lowE'
datalabel[defdata_EE_TEB] = r'\textit{Planck} EE$+$lowT,P'

defdata_all = defdata_root + '_TTTEEE_lowl_lowE'
datalabel[defdata_all] = r'\textit{Planck} TT,TE,EE$+$lowE'
defdata_TTTEEE = defdata_all
defdata_TTonly = defdata_root + '_TT_lowl'
datalabel[defdata_TTonly] = r'\textit{Planck} TT'
defdata_allNoLowE = defdata_root + '_TTTEEE_lowl'
datalabel[defdata_allNoLowE] = r'\textit{Planck} TT,TE,EE$+$lowT'

defdata = defdata_TT
deflabel = datalabel[defdata_TT]

defdata_lensing = defdata_TT + '_lensing'
datalabel[defdata_lensing] = datalabel[defdata_TT] + '$+$lensing'
defdata_all_lensing = defdata_all + '_lensing'
datalabel[defdata_all_lensing] = datalabel[defdata_all] + '$+$lensing'

planck = r'\textit{Planck}'

planckTT = datalabel[defdata_TT]
planckTTonly = datalabel[defdata_TTonly]

planckTTlowTEB = datalabel[defdata_TT]
planckall = datalabel[defdata_all]
shortall = 'TT,TE,EE$+$lowE'
shortEE = 'EE$+$lowE'
shortTE = 'TE$+$lowE'
shortTT = 'TT$+$lowE'

NoLowLE = datalabel[defdata_allNoLowE]
lensing = datalabel[defdata_lensing]
lensingall = datalabel[defdata_all_lensing]
defplanck = datalabel[defdata]

shortlabel = {}
for key, value in list(datalabel.items()):
    shortlabel[key] = value.replace(planck + ' ', '')

NoLowLhighLtau = r'\textit{Planck}$-$lowL+highL+$\tau$prior'
NoLowLhighL = r'\textit{Planck}$-$lowL+highL'
WPhighLlensing = r'\textit{Planck}+lensing+WP+highL'
WP = r'\textit{Planck}+WP'
WPhighL = r'\textit{Planck}+WP+highL'
NoLowL = r'\textit{Planck}$-$lowL'
lensonly = 'lensing'
HST = r'$H_0$'
BAO = 'BAO'

LCDM = r'$\Lambda$CDM'

s = copy.copy(plots.defaultSettings)
s.legend_frame = False
s.figure_legend_frame = False
s.prob_label = r'$P/P_{\rm max}$'
s.norm_prob_label = 'Probability density'
s.prob_y_ticks = True
s.param_names_for_labels = os.path.join(batchjob.getCodeRootPath(), 'clik_units.paramnames')
s.alpha_filled_add = 0.85
s.solid_contour_palefactor = 0.6

s.solid_colors = [('#8CD3F5', '#006FED'), ('#F7BAA6', '#E03424'), ('#D1D1D1', '#A1A1A1'), 'g', 'cadetblue', 'olive',
                  'darkcyan']
s.axis_marker_lw = 0.6
s.lw_contour = 1

s.param_names_for_labels = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'clik_latex.paramnames'))

use_plot_data = getdist.use_plot_data
rootdir = getdist.default_grid_root or os.path.join(batchjob.getCodeRootPath(), 'main')
output_base_dir = getdist.output_base_dir or batchjob.getCodeRootPath()

H0_gpe = [70.6, 3.3]


# various Omegam sigma8 constraints for plots
def planck_lensing_2015(omm, sigma):
    # g60_full
    return (0.591 + 0.021 * sigma) * omm ** (-0.25)


def planck_lensing(omm, sigma):
    return (0.589 + 0.020 * sigma) * omm ** (-0.25)


def plotBounds(omm, data, c='gray'):
    plt.fill_between(omm, data(omm, -2), data(omm, 2), facecolor=c, alpha=0.15, edgecolor=c, lw=0)
    plt.fill_between(omm, data(omm, -1), data(omm, 1), facecolor=c, alpha=0.25, edgecolor=c, lw=0)


def plotBands(x, means, errors, color='gray', linecolor='k', nbands=2, alphas=[0.5, 0.15, 0.1], lw=0.2):
    if np.isscalar(means): means = np.ones(len(x)) * means
    if nbands > 2: plt.fill_between(x, means - 3 * errors, means + 3 * errors, color=color, alpha=alphas[2], lw=lw)
    if nbands > 1: plt.fill_between(x, means - 2 * errors, means + 2 * errors, color=color, alpha=alphas[1], lw=lw)
    plt.fill_between(x, means - errors, means + errors, color=color, alpha=alphas[0], lw=lw)
    if linecolor: plt.plot(x, means, color=linecolor, lw=0.5)


class planckPlotter(plots.GetDistPlotter):
    def bestfit(self, root, want_fixed=True):
        for chain_dir in self.sampleAnalyser.chain_dirs:
            jobItem = chain_dir.resolveRoot(root)
            if jobItem:
                return BestFit(jobItem.chainRoot + '.minimum', want_fixed=want_fixed)
        return None

    def bestFitCls(self, root):
        from CMBlikes import ClsArray
        for chain_dir in self.sampleAnalyser.chain_dirs:
            jobItem = chain_dir.resolveRoot(root)
            if jobItem:
                return ClsArray(jobItem.chainRoot + '.minimum.theory_cl')
        return None

    def getBatch(self):
        if not hasattr(self, 'batch'): self.batch = batchjob.readobject(rootdir)
        return self.batch

    def doExport(self, fname=None, adir=None, watermark=None, tag=None):
        if watermark is None and non_final:
            watermark = version
        if adir:
            if os.sep not in adir and '/' not in adir: adir = os.path.join(output_base_dir, adir)
        super(planckPlotter, self).export(fname, adir, watermark, tag)

    def export(self, fname=None, tag=None):
        self.doExport(fname, 'outputs', tag=tag)

    def exportExtra(self, fname=None):
        self.doExport(fname, 'plots')

    def getRoot(self, paramtag, datatag, returnJobItem=False):
        return self.getBatch().resolveName(paramtag, datatag, returnJobItem=returnJobItem)

    def getSamples(self, paramtag, datatag=None):
        if datatag is not None:
            return self.sampleAnalyser.samplesForRoot(self.getRoot(paramtag, datatag))
        else:
            return self.sampleAnalyser.samplesForRoot(paramtag)

    def getJobItem(self, paramtag, datatag):
        jobItem = self.getRoot(paramtag, datatag, returnJobItem=True)
        jobItem.loadJobItemResults(paramNameFile=self.settings.param_names_for_labels)
        return jobItem

    def getCacheDir(self):
        cache_dir = os.path.join(self.getBatch().batchPath, 'caches')
        if not os.path.exists(cache_dir):
            os.mkdir(cache_dir)
        return cache_dir

    def getLikelihood(self, paramtag, datatag, want_unbinned=False):
        # get python likelihood, assuming TT or TTTEEE likelihood, plik or Camspec
        jobItem = self.getJobItem(paramtag, datatag)
        ini = IniFile(jobItem.chainRoot + '.inputparams')
        plik_file = ini.string('clik_data_plik', '')
        if plik_file:
            params = ini.string('clik_params_plik', '')
            if not os.path.exists(plik_file):
                plik_file = os.path.join(os.path.dirname(__file__), r'../',
                                         plik_file.replace('%DATASETDIR%', 'data/'))
                params = os.path.join(os.path.dirname(__file__), r'../',
                                      params.replace('%DATASETDIR%', 'data/'))
                if not os.path.exists(plik_file):
                    raise Exception('plik file to found %s' % plik_file)

            from planck.CMB_plik import plik_likelihood
            like = plik_likelihood(plik_file, params)
            if want_unbinned:
                return like, plik_likelihood(plik_file.replace('.clik', '_bin1.clik'), params)
            else:
                return like
        else:
            from planck.CMB_CamSpec import CamSpec_likelihood

            fname = os.path.join(os.path.dirname(__file__), r'../data/planck_internal/CamSpecHM_10_7.dataset')
            if '_TT_' in jobItem.chainRoot:
                like = CamSpec_likelihood(fname, {'use_cl': '100x100 143x143 217x217 143x217'})
            else:
                like = CamSpec_likelihood(fname)

            if want_unbinned:
                return like, like
            else:
                return like


def getPlotter(plot_data=None, chain_dir=None, **kwargs):
    global plotter, rootdir
    if not kwargs.get('settings'):
        kwargs['settings'] = s
    # make sure rc changed as desired  (e.g. module is reloaded)
    setRc()
    if plot_data is not None or use_plot_data:
        plotter = planckPlotter(plot_data or os.path.join(rootdir, 'plot_data'), **kwargs)
    if chain_dir or not use_plot_data:
        plotter = planckPlotter(chain_dir=chain_dir or rootdir, **kwargs)
    return plotter


def getSubplotPlotter(plot_data=None, chain_dir=None, subplot_size=2, **kwargs):
    s.setWithSubplotSize(subplot_size)
    s.axes_fontsize += 2
    s.colorbar_axes_fontsize += 2
    s.legend_fontsize = s.lab_fontsize + 1
    return getPlotter(plot_data, chain_dir, **kwargs)


def getPlotterWidth(size=1, **kwargs):  # size in mm
    inch_mm = 0.0393700787
    if size == 1:
        width = 88 * inch_mm
    elif size == 2:
        width = 120 * inch_mm
    elif size == 3:
        width = 180 * inch_mm
    else:
        width = size * inch_mm
    s.fig_width_inch = width
    s.setWithSubplotSize(kwargs.get('subplot_size', 2))
    s.rcSizes(**kwargs)
    return getPlotter()


def getSinglePlotter(ratio=3 / 4., plot_data=None, chain_dir=None, width_inch=3.464, **kwargs):
    s.setWithSubplotSize(width_inch)
    s.fig_width_inch = width_inch
    s.rcSizes()
    plotter = getPlotter(plot_data, chain_dir, **kwargs)
    plotter.make_figure(1, xstretch=1 / ratio)
    return plotter


class planckStyleTableFormatter(types.NoLineTableFormatter):
    """Planck style guide compliant formatter

    Andrea Zonca (edits by AL for consistent class structure)"""

    tableOpen = r"""
\begingroup
\openup 5pt
\newdimen\tblskip \tblskip=5pt
\nointerlineskip
\vskip -3mm
\scriptsize
\setbox\tablebox=\vbox{
    \newdimen\digitwidth
    \setbox0=\hbox{\rm 0}
    \digitwidth=\wd0
    \catcode`"=\active
    \def"{\kern\digitwidth}
%
    \newdimen\signwidth
    \setbox0=\hbox{+}
    \signwidth=\wd0
    \catcode`!=\active
    \def!{\kern\signwidth}
%
\halign{"""

    tableClose = r"""} % close halign
} % close vbox
\endPlancktable
\endgroup
"""

    def __init__(self):
        super(planckStyleTableFormatter, self).__init__()
        self.aboveHeader = None
        self.belowHeader = r'\noalign{\vskip 3pt\hrule\vskip 5pt}'
        self.aboveTitles = r'\noalign{\doubleline}'
        self.belowTitles = ''
        self.minorDividor = ''
        self.majorDividor = ''
        self.endofrow = r'\cr'
        self.hline = r'\noalign{\vskip 5pt\hrule\vskip 3pt}'
        self.belowFinalRow = self.hline
        self.belowBlockRow = self.hline
        self.belowRow = None
        self.colDividor = '|'
        self.headerWrapper = "\\omit\\hfil %s\\hfil"
        self.noConstraint = r'\dots'
        self.colSeparator = '&'
        self.spacer = ''

    def formatTitle(self, title):
        return types.texEscapeText(title)

    def belowTitleLine(self, colsPerParam, numResults):
        out = r'\noalign{\vskip -3pt}'
        if colsPerParam > 1:
            out += "\n"
            out += r"\omit"
            out += (r"&\multispan" + str(colsPerParam) + r"\hrulefill") * numResults
            out += r"\cr"
        out += self.getLine("belowTitles")
        return out

    def startTable(self, ncol, colsPerResult, numResults):
        tableOpen = self.tableOpen + "\n"
        tableOpen += r"""\hbox to 0.9in{$#$\leaderfil}\tabskip=1.5em&"""
        if numResults > 3 and colsPerResult == 2:
            for res in range(numResults):
                tableOpen += r"\hfil$#$\hfil\tabskip=0.5em&" + "\n"
                if res < numResults - 1:
                    tableOpen += r"\hfil$#$\hfil\tabskip=1.7em&" + "\n"
        else:
            tableOpen += r"$#$\hfil&" * (colsPerResult * numResults - 1)
        tableOpen += r"\hfil$#$\hfil\tabskip=0pt\cr"
        return tableOpen

    def endTable(self):
        return self.tableClose

    def titleSubColumn(self, colsPerResult, title):
        return '\\multispan' + str(colsPerResult) + '\hfil ' + self.formatTitle(title) + '\hfil'

    def textAsColumn(self, txt, latex=False, separator=False, bold=False):
        bold = False
        if latex:
            res = txt  # there should be NO SPACE after a number in latex AZ
        else:
            wid = len(txt)
            res = txt + self.spacer * max(0, 28 - wid)
        if latex:
            if bold:
                res = '{\\boldmath$' + res + '$}'
            else:
                res = res
        if separator:
            if latex:
                res += self.colSeparator  # there should be NO SPACE after a number in latex AZ
            else:
                res += self.colSeparator
        return res
