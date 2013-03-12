import os, ResultObjs, GetDistPlots
from matplotlib import rcParams, rc

# common setup for matplotlib
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'font.family':'sans-serif',
          # free font similar to Helvetica
          'font.sans-serif':'FreeSans'}

sfmath = os.path.dirname(os.path.abspath(__file__)) + os.sep + 'sfmath'
# use of Sans Serif also in math mode
rc('text.latex', preamble=r'\usepackage{' + sfmath + '}')

rcParams.update(params)

planck = r'\textit{Planck}'
WP = r'\textit{Planck}+WP'
WPhighL = r'\textit{Planck}+WP+highL'
lensing = r'\textit{Planck}+lensing'
WPhighLlensing = r'\textit{Planck}+lensing+WP+highL'
NoLowL = r'\textit{Planck}$-$lowL'
NoLowLhighL = r'\textit{Planck}$-$lowL+highL'
NoLowLtau = r'\textit{Planck}$-$lowL+$\tau$prior'
NoLowLhighLtau = r'\textit{Planck}$-$lowL+highL+$\tau$prior'
LCDM = r'$\Lambda$CDM'

s = GetDistPlots.defaultSettings
s.legend_frame = False
s.figure_legend_frame = False
s.prob_label = r'$P/P_{\rm max}$'
s.prob_y_ticks = True
s.param_names_for_labels = 'clik_units.paramnames'
s.alpha_filled_add = 0.85
# s.solid_colors = ['#006FED', '#E03424', 'gray', '#009966' ]
s.solid_contour_palefactor = 0.6

s.solid_colors = [('#8CD3F5', '#006FED'), ('#F7BAA6', '#E03424'), ('#D1D1D1', '#A1A1A1')]
s.axis_marker_lw = 0.6
s.lw_contour = 0.2

class planckPlotter(GetDistPlots.GetDistPlotter):
    def export(self, fname):
        if '.' in fname:GetDistPlots.GetDistPlotter.export(self, fname)
        else:
            GetDistPlots.GetDistPlotter.export(self, 'outputs/' + fname + '.pdf')

    def exportExtra(self, fname):
        GetDistPlots.GetDistPlotter.export(self, 'plots/' + fname + '.pdf')


plotter = planckPlotter('main/plot_data')

def getSinglePlotter(ratio=3 / 4.):
    s.setWithSubplotSize(3.5)
    s.legend_font_size = 9
    s.colorbar_axes_fontsize = 8
    s.scatter_size = 2
    s.lab_fontsize = 13
    s.axes_fontsize = 9

    plotter.make_figure(1, ystretch=ratio)
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
            res = txt + ' ' * max(0, 28 - wid)
        if latex:
            if bold: res = '{\\boldmath$' + res + '$}'
            else:  res = res
        if separator:
            if latex:
                res += '& '  # there should be NO SPACE after a number in latex AZ
            else:
                res += ' & '
        return res
