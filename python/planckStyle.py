import ResultObjs, GetDistPlots
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

# use of Sans Serif also in math mode
rc('text.latex')
# , preamble='\usepackage{sfmath}')

rcParams.update(params)

planck = r'\textit{Planck}'
WP = r'\textit{Planck}+WP'
WPhighL = r'\textit{Planck}+WP+highL'
lensing = r'\textit{Planck}+TT lensing'
WPhighLlensing = r'\textit{Planck}+WP+highL+lensing'


class planckPlotter(GetDistPlots.GetDistPlotter):
    def export(self, fname):
        GetDistPlots.GetDistPlotter.export(self, 'outputs/' + fname + '.pdf')
#        GetDistPlots.GetDistPlotter.export(self, 'outputs/' + fname + '.eps')

    def exportExtra(self, fname):
        GetDistPlots.GetDistPlotter.export(self, 'plots/' + fname + '.pdf')


plotter = planckPlotter('main/plot_data')


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
