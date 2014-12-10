import os, ResultObjs, GetDistPlots, sys, batchJob
from matplotlib import rcParams, rc, pylab
import iniFile

# common setup for matplotlib
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'font.size': 10,
          'legend.fontsize': 9,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'text.usetex': True,
          'font.family':'sans-serif',
          # free font similar to Helvetica
          'font.sans-serif':'FreeSans'}

sfmath = os.path.dirname(os.path.abspath(__file__)) + os.sep + 'sfmath'
# use of Sans Serif also in math mode
rc('text.latex', preamble=r'\usepackage{' + sfmath.replace(os.sep, '/') + '}')

rcParams.update(params)

if False:
    non_final = True
    version = 'CamSpec v910HM'
    defdata_root = 'CamSpecHM'
else:
    non_final = False
    version = 'clik10.2'
    defdata_root = 'plikHM'

datalabel = dict()
defdata_TT = defdata_root + '_TT_lowTEB'
datalabel[defdata_TT] = r'\textit{Planck} TT$+$lowP'
defdata_TE = defdata_root + '_TE_lowEB'
datalabel[defdata_TE] = r'\textit{Planck} TE$+$lowP'
defdata_EE = defdata_root + '_EE_lowEB'
datalabel[defdata_EE] = r'\textit{Planck} EE$+$lowP'
defdata_TE_TEB = defdata_root + '_TE_lowTEB'
datalabel[defdata_TE_TEB] = r'\textit{Planck} TE$+$lowT,P'
defdata_EE_TEB = defdata_root + '_EE_lowTEB'
datalabel[defdata_EE_TEB] = r'\textit{Planck} EE$+$lowT,P'


defdata_all = defdata_root + '_TTTEEE_lowTEB'
datalabel[defdata_all] = r'\textit{Planck} TT,TE,EE$+$lowP'
defdata_TTTEEE = defdata_all
defdata_TTonly = defdata_root + '_TT_lowl'
datalabel[defdata_TTonly] = r'\textit{Planck} TT'
defdata_allNoLowE = defdata_root + '_TTTEEE_lowl'
datalabel[defdata_allNoLowE] = r'\textit{Planck} TT,TE,EE'

defdata = defdata_TT
deflabel = datalabel[defdata_TT]

defdata_lensing = defdata_TT + '_lensing'
datalabel[defdata_lensing] = datalabel[defdata_TT] + '$+$lensing'
defdata_all_lensing = defdata_all + '_lensing'
datalabel[defdata_all_lensing] = datalabel[defdata_all] + '$+$lensing'

planck = r'\textit{Planck}'

planckTT = datalabel[defdata_TTonly]
planckTTlowTEB = datalabel[defdata_TT]
planckall = datalabel[defdata_all]
NoLowLE = datalabel[defdata_allNoLowE]
lensing = datalabel[defdata_lensing]
lensingall = datalabel[defdata_all_lensing]
defplanck = datalabel[defdata]

shortlabel = {}
for key, value in datalabel.items():
    shortlabel[key] = value.replace(planck + ' ', '')

NoLowLhighLtau = r'\textit{Planck}$-$lowL+highL+$\tau$prior'
NoLowLhighL = r'\textit{Planck}$-$lowL+highL'
WPhighLlensing = r'\textit{Planck}+lensing+WP+highL'
WP = r'\textit{Planck}+WP'
WPhighL = r'\textit{Planck}+WP+highL'
NoLowL = r'\textit{Planck}$-$lowL'
NoLowLtau = r'\textit{Planck}$-$lowL+$\tau$prior'
lensonly = 'lensing'
HST = r'$H_0$'
BAO = 'BAO'


LCDM = r'$\Lambda$CDM'

s = GetDistPlots.defaultSettings
s.legend_frame = False
s.figure_legend_frame = False
s.prob_label = r'$P/P_{\rm max}$'
s.norm_prob_label = 'Probability density'
s.prob_y_ticks = True
s.param_names_for_labels = 'clik_units.paramnames'
s.alpha_filled_add = 0.85
# s.solid_colors = ['#006FED', '#E03424', 'gray', '#009966' ]
s.solid_contour_palefactor = 0.6

s.solid_colors = [('#8CD3F5', '#006FED'), ('#F7BAA6', '#E03424'), ('#D1D1D1', '#A1A1A1'), 'g', 'cadetblue', 'indianred']
s.axis_marker_lw = 0.6
s.lw_contour = 1

use_plot_data = True

def getRoot():
    global use_plot_data
    codeRoot = batchJob.getCodeRootPath()
    configf = os.path.join(codeRoot, 'python', 'config.ini')
    root = None
    output_base_dir = None
    if os.path.exists(configf):
        ini = iniFile.iniFile(configf)
        root = ini.string('default_grid_root', '')
        output_base_dir = ini.string('output_base_dir', '')
        use_plot_data = ini.bool('use_plot_data', use_plot_data)
    return root or os.path.join(codeRoot, 'main'), output_base_dir or codeRoot

rootdir, output_base_dir = getRoot()

H0_high = [73.9, 2.7]
H0_Freeman12 = [74.3, 2.1]
H0_gpe = [70.6, 3.3]

# various Omegam sigma8 constraints for plots
def PLSZ(omm, sigma):
    # from Anna 18/7/2014 for fixed b=0.8
    return  (0.757 + 0.013 * sigma) * (omm / 0.32) ** (-0.3)

def CFTHlens_Kilbinger(omm, sigma):
    return  (0.79 + 0.03 * sigma) * (omm / 0.27) ** (-0.6)


def CFTHlens_LCDM(omm, sigma):  # 1408.4742
    return  (0.74 + 0.03 * sigma) * (omm / 0.27) ** (-0.47)

def CFTHlens_mnu(omm, sigma):  # 1408.4742, same for sterile case
    return  (0.72 + 0.03 * sigma) * (omm / 0.27) ** (-0.48)


def galaxygalaxy(omm, sigma):  # mandelbaum
    return  (0.8 + 0.05 * sigma) * (omm / 0.25) ** (-0.57)

def planck_lensing(omm, sigma):
    # g60_full
    return  (0.592 + 0.021 * sigma) * omm ** (-0.25)


def plotBounds(omm, data, c='gray'):
    pylab.fill_between(omm, data(omm, -2), data(omm, 2), facecolor=c, alpha=0.15, edgecolor=c, lw=0)
    pylab.fill_between(omm, data(omm, -1), data(omm, 1), facecolor=c, alpha=0.25, edgecolor=c, lw=0)


class planckPlotter(GetDistPlots.GetDistPlotter):

    def getBatch(self):
        if not hasattr(self, 'batch'): self.batch = batchJob.readobject(rootdir)
        return self.batch

    def doExport(self, fname=None, adir=None, watermark=None, tag=None):
        if watermark is None and non_final:
            watermark = version
        if adir:
            if not os.sep in adir: adir = os.path.join(output_base_dir, adir)
        super(planckPlotter, self).export(fname, adir, watermark, tag)

    def export(self, fname=None, tag=None):
        self.doExport(fname, 'outputs', tag=tag)

    def exportExtra(self, fname=None):
        self.doExport(fname, 'plots')

    def getRoot(self, paramtag, datatag, returnJobItem=False):
        return self.getBatch().resolveName(paramtag, datatag, returnJobItem=returnJobItem)

    def getJobItem(self, paramtag, datatag):
        jobItem = self.getRoot(paramtag, datatag, returnJobItem=True)
        jobItem.loadJobItemResults(paramNameFile=self.settings.param_names_for_labels)
        return jobItem

def getPlotter(plot_data=None, grid_dir=None):
    global plotter, rootdir
    if plot_data is not None or use_plot_data:
        plotter = planckPlotter(plot_data or os.path.join(rootdir, 'plot_data'))
    if grid_dir or not use_plot_data:
        plotter = planckPlotter(chain_dir=grid_dir or rootdir)
    return plotter

plotter = getPlotter()


def getSubplotPlotter(plot_data=None, grid_dir=None):
    s.setWithSubplotSize(2)
    s.axes_fontsize += 2
    s.colorbar_axes_fontsize += 2
#    s.lab_fontsize += 2
    s.legend_fontsize = s.lab_fontsize + 1
    return getPlotter(plot_data, grid_dir)

def getPlotterWidth(size=1, **kwargs):  # size in mm
    inch_mm = 0.0393700787
    if size == 1:
        width = 88 * inch_mm
    elif size == 2: width = 120 * inch_mm
    elif size == 3: width = 180 * inch_mm
    else: width = size * inch_mm
    s.fig_width_inch = width
    s.setWithSubplotSize(2)
    s.rcSizes(**kwargs)
    return getPlotter()

def getSinglePlotter(ratio=3 / 4., plot_data=None, grid_dir=None):
    s.setWithSubplotSize(3.5)
    s.rcSizes()
    plotter = getPlotter(plot_data, grid_dir)
    plotter.make_figure(1, xstretch=1 / ratio)
    return plotter


class planckStyleTableFormatter(ResultObjs.noLineTableFormatter):
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
        return ResultObjs.texEscapeText(title)

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
            if bold: res = '{\\boldmath$' + res + '$}'
            else:  res = res
        if separator:
            if latex:
                res += self.colSeparator  # there should be NO SPACE after a number in latex AZ
            else:
                res += self.colSeparator
        return res
