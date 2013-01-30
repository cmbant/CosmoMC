import paramNames, decimal

class textFile:

    def __init__(self, lines=[]):
        self.lines = lines

    def write(self, outfile):
        textFileHandle = open(outfile, 'w')
        textFileHandle.write("\n".join(self.lines))
        textFileHandle.close()

def texEscapeText(string):
        return string.replace('_', '{\\textunderscore}')


def float_to_decimal(f):
    # http://docs.python.org/library/decimal.html#decimal-faq
    "Convert a floating point number to a Decimal with no loss of information"
    n, d = f.as_integer_ratio()
    numerator, denominator = decimal.Decimal(n), decimal.Decimal(d)
    ctx = decimal.Context(prec=60)
    result = ctx.divide(numerator, denominator)
    while ctx.flags[decimal.Inexact]:
        ctx.flags[decimal.Inexact] = False
        ctx.prec *= 2
        result = ctx.divide(numerator, denominator)
    return result

def numberFigs(number, sigfig):
    # http://stackoverflow.com/questions/2663612/nicely-representing-a-floating-point-number-in-python/2663623#2663623
    assert(sigfig > 0)
    try:
        d = decimal.Decimal(number)
    except TypeError:
        d = float_to_decimal(float(number))
    sign, digits = d.as_tuple()[0:2]
    if len(digits) < sigfig:
        digits = list(digits)
        digits.extend([0] * (sigfig - len(digits)))
    shift = d.adjusted()
    result = int(''.join(map(str, digits[:sigfig])))
    # Round the result
    if len(digits) > sigfig and digits[sigfig] >= 5: result += 1
    result = list(str(result))
    # Rounding can change the length of result
    # If so, adjust shift
    shift += len(result) - sigfig
    # reset len of result to sigfig
    result = result[:sigfig]
    if shift >= sigfig - 1:
        # Tack more zeros on the end
        result += ['0'] * (shift - sigfig + 1)
    elif 0 <= shift:
        # Place the decimal point in between digits
        result.insert(shift + 1, '.')
    else:
        # Tack zeros on the front
        assert(shift < 0)
        result = ['0.'] + ['0'] * (-shift - 1) + result
    if sign:
        result.insert(0, '-')
    return ''.join(result)

class numberFormatter():
    def __init__(self, sig_figs=4):
        self.sig_figs = sig_figs

    def namesigFigs(self, value, limplus, limminus):
        frac = limplus / (abs(value) + limplus)
        err_sf = 2
        if value >= 20 and frac > 0.1 and limplus >= 2: err_sf = 1

        plus_str = self.formatNumber(limplus, err_sf, True)
        minus_str = self.formatNumber(limminus, err_sf, True)
        sf = self.sig_figs
        if frac > 0.1 and value < 100 and value >= 20: sf = 2
        elif frac > 0.01 and value < 1000: sf = 3
#        if abs(value) < 1 and limplus - limminus > abs(value): sf = 2
        res = self.formatNumber(value, sf)
        maxdp = max(self.decimal_places(plus_str), self.decimal_places(minus_str))
        while abs(value) < 1 and maxdp < self.decimal_places(res):
            sf -= 1
            if sf == 0:
                res = ('%.' + str(maxdp) + 'f') % value
                if (float(res) == 0.0): res = ('%.' + str(maxdp) + 'f') % 0
                break
            else: res = self.formatNumber(value, sf)

        while self.decimal_places(plus_str) > self.decimal_places(res):
            sf += 1
            res = self.formatNumber(value, sf)
        return (res, plus_str, minus_str)

    def formatNumber(self, value, sig_figs=None, wantSign=False):
        if sig_figs is None:
            sf = self.sig_figs
        else: sf = sig_figs
        s = numberFigs(value, sf)
        if wantSign:
            if s[0] != '-' and float(s) < 0: s = '-' + s
            if  float(s) > 0: s = '+' + s
        return s

    def decimal_places(self, s):
        i = s.find('.')
        if i > 0: return len(s) - i - 1
        return 0


class tableFormatter():
    def __init__(self):
        self.border = '|'
        self.endofrow = '\\\\'
        self.hline = '\\hline'
        self.paramText = 'Parameter'
        self.aboveTitles = self.hline
        self.majorDividor = '|'
        self.minorDividor = '|'
        self.colDividor = '||'
        self.belowTitles = ''

    def getLine(self, position=None):
        if position is not None and hasattr(self, position): return getattr(self, position)
        return self.hline

    def belowTitleLine(self, colsPerParam):
        return self.getLine("belowTitles")

    def startTable(self, ncol, colsPerResult, numResults):
        part = self.majorDividor + (" c" + self.minorDividor) * (colsPerResult - 1) + ' c'
        return '\\begin{tabular} {' + self.border + " l " + part * numResults + (self.colDividor + " l " + part * numResults) * (ncol - 1) + self.border + '}'

    def endTable(self):
        return '\\end{tabular}'

    def titleSubColumn(self, colsPerResult, title):
        return ' \\multicolumn{' + str(colsPerResult) + '}{' + self.majorDividor + 'c' + self.majorDividor + '}{' + self.formatTitle(title) + '}'

    def formatTitle(self, title):
        return '\\bf ' + texEscapeText(title)

class planckTableFormatter(tableFormatter):

    def __init__(self):
        tableFormatter.__init__(self)
        self.border = ''
        self.aboveTitles = r'\noalign{\vskip 3pt}' + self.hline + r'\noalign{\vskip 1.5pt}' + self.hline + r'\noalign{\vskip 5pt}'
        self.belowTitles = r'\noalign{\vskip 3pt}' + self.hline
        self.aboveHeader = ''
        self.belowHeader = self.hline
        self.minorDividor = ''
        self.belowFinalRow = ''

    def titleSubColumn(self, colsPerResult, title):
        return ' \\multicolumn{' + str(colsPerResult) + '}{' + 'c' + '}{ ' + self.formatTitle(title) + '}'

class planckNoLineTableFormatter(planckTableFormatter):

    def __init__(self):
        planckTableFormatter.__init__(self)
        self.aboveHeader = ''
#        self.belowHeader = r'\noalign{\vskip 5pt}'
        self.minorDividor = ''
        self.majorDividor = ''
        self.belowFinalRow = self.hline
        self.belowBlockRow = self.hline
        self.colDividor = '|'
        self.hline = ''

    def belowTitleLine(self, colsPerParam):
        return r'\noalign{\vskip 3pt}\cline{2-' + str(colsPerParam + 1) + r'}\noalign{\vskip 3pt}'

class resultTable():

    def __init__(self, ncol, results, tableParamNames=None, titles=None, formatter=None, numFormatter=None, blockEndParams=None, paramList=None):
# results is a margeStats or bestFit table
        self.lines = []
        if formatter is None: self.format = planckNoLineTableFormatter()
        else: self.format = formatter
        self.ncol = ncol
        if tableParamNames is None:
            self.tableParamNames = results[0]
        else: self.tableParamNames = tableParamNames
        if paramList is not None: self.tableParamNames = self.tableParamNames.filteredCopy(paramList)
        if numFormatter is None: self.numberFormatter = numberFormatter()
        else: self.numberFormatter = numFormatter

        self.results = results
        self.boldBaseParameters = True
        self.colsPerResult = len(results[0].columns)
        self.colsPerParam = len(results) * self.colsPerResult

        nparams = self.tableParamNames.numParams()
        numrow = nparams / ncol
        if nparams % ncol != 0: numrow += 1
        rows = []
        for par in self.tableParamNames.names[0:numrow]:
            rows.append([par])
        for col in range(1, ncol):
            for i in range(numrow * col, min(numrow * (col + 1), nparams)):
                rows[i - numrow * col].append(self.tableParamNames.names[i]);

        self.lines.append(self.format.startTable(ncol, self.colsPerResult, len(results)))
        if titles is not None: self.addTitlesRow(titles)
        self.addHeaderRow()
        for row in rows[:-1]:
            self.addFullTableRow(row)
            if ncol == 1 and blockEndParams is not None and row[0].name in blockEndParams: self.addLine("belowBlockRow")
            else: self.addLine("belowRow")
        self.addFullTableRow(rows[-1])
        self.addLine("belowFinalRow")
        self.endTable()


    def addFullTableRow(self, row):
        txt = " & ".join(self.paramLabelColumn(param) + self.paramResultsTex(param) for param in row)
        if not self.ncol == len(row):
            txt += ' & ' * ((1 + self.colsPerParam) * (self.ncol - len(row)))
        self.lines.append(txt + self.format.endofrow)

    def addLine(self, position):
        return self.lines.append(self.format.getLine(position))

    def addTitlesRow(self, titles):
        self.addLine("aboveTitles")
        res = self.format.titleSubColumn(1, '') + ' & ' + " & ".join(self.format.titleSubColumn(self.colsPerResult, title) for title in titles)
        self.lines.append((('& ' + res) * self.ncol)[1:] + self.format.endofrow)
        self.lines.append(self.format.belowTitleLine(self.colsPerParam))

    def addHeaderRow(self):
        self.addLine("aboveHeader")
        res = '& ' + self.format.paramText
        for result in self.results:
            res += ' & ' + " & ".join(result.columns)
        self.lines.append((res * self.ncol)[1:] + self.format.endofrow)
        self.addLine("belowHeader")

    def paramResultsTex(self, param):
        return " & ".join(self.paramResultTex(result, param) for result in self.results)

    def paramResultTex(self, result, p):
        values = result.texValues(self.numberFormatter, p)
        if values is not None:
            if len(values) > 1: txt = self.textAsColumn(values[1], True, separator=True)
            else: txt = ''
            txt += self.textAsColumn(values[0], values[0] != '---')
            return txt
        else: return self.textAsColumn('') * len(result.columns)

    def textAsColumn(self, txt, latex=False, separator=False, bold=False):
        wid = len(txt)
        if latex:
            wid += 2
            if bold: wid += 11
        res = txt + ' ' * max(0, 28 - wid)
        if latex:
            if bold: res = '{\\boldmath$' + res + '$}'
            else:  res = '$' + res + '$'
        if separator: res += ' & '
        return res

    def paramLabelColumn(self, param):
        return  self.textAsColumn(param.label, True, separator=True, bold=not param.isDerived)

    def endTable(self):
        self.lines.append(self.format.endTable())

    def tableTex(self):
        return "\n".join(self.lines)

    def writeTable(self, fname):
        textFile(self.lines).write(fname)


class paramResults(paramNames.paramList): pass

class bestFit(paramResults):

    def __init__(self, fileName=None, setParamNameFile=None, want_fixed=False):
        paramResults.__init__(self)
        self.columns = ['Best fit']
        if fileName is not None: self.loadFromFile(fileName, want_fixed=want_fixed)
        if setParamNameFile is not None: self.setLabelsAndDerivedFromParamNames(setParamNameFile)

    def loadFromFile(self, filename, want_fixed=False):
        textFileLines = self.fileList(filename)
        first = textFileLines[0].strip().split('=')
        self.logLike = float(first[1].strip())
        isFixed = False
        isDerived = False
        self.chiSquareds = []
        chunks = 0
        for ix in range(2, len(textFileLines)):
            line = textFileLines[ix]
            if len(line.strip()) == 0:
                chunks += 1
                isFixed = not isFixed
                isDerived = True
                if chunks == 3:
                    if ix + 2 >= len(textFileLines): break
                    for likePart in textFileLines[ix + 2:]:
                        if len(likePart.strip()) != 0:
                            (chisq, name) = [s.strip() for s in likePart.split(None, 2)][1:]
                            name = [s.strip() for s in name.split(':', 1)]
                            if len(name) > 1:
                                (kind, name) = name
                            else: kind = ''
                            self.chiSquareds.append((kind, name, float(chisq)))
                    break
                continue
            if not isFixed or want_fixed:
                param = paramNames.paramInfo()
                param.isDerived = isDerived
                (param.number, param.best_fit, param.name, param.label) = [s.strip() for s in line.split(None, 3)]
                param.number = int(param.number)
                param.best_fit = float(param.best_fit)
                self.names.append(param)

    def sortedChiSquareds(self):
        likes = dict()
        for (kind, name, chisq) in self.chiSquareds:
            if not kind in likes: likes[kind] = []
            likes[kind].append((name, chisq))
        return sorted(likes.iteritems())

    def texValues(self, formatter, p):
        param = self.parWithName(p.name)
        if param is not None: return [formatter.formatNumber(param.best_fit)]
        else: return None


class margeStats(paramResults):

    def loadFromFile(self, filename):
        textFileLines = self.fileList(filename)
        for i in range(1, len(textFileLines)):
            line = textFileLines[i]
            if len(line.strip()) == 0:
                lims = textFileLines[i + 1].split(':')[1]
                self.limits = [float(s.strip()) for s in lims.split(';')]
                for line in textFileLines[i + 2:]:
                    par = self.parWithNumber(int(line.split()[0].strip()))
                    par.name = line.split()[1].strip()
                    tag = line[27:39].strip()
                    par.twotail = tag == 'two tail'
                    par.lim_top = tag == '< one tail' or tag == 'no tail'
                    par.lim_bot = tag == '> one tail' or tag == 'no tail'
                break
            param = paramNames.paramInfo()
            items = [s.strip() for s in line.split(None, 7)]
            param.number = int(items[0])
            param.mean = float(items[1])
            param.err = float(items[2])

            param.limits = [[float(s) for s in items[3:5]], [float(s) for s in items[5:7]]]
            param.label = items[7]
            self.names.append(param)
        self.columns = [str(round(float(self.limits[1]) * 100)) + '\\% limits']


    def addBestFit(self, bf):
        self.columns = ['Best fit'] + self.columns
        self.logLike = bf.logLike
# the next line deletes parameters not in best-fit; this is good e.g. to get rid of yhe from importance sampled result
        self.names = [x for x in self.names if bf.parWithName(x.name) is not None]
        for par in self.names:
            param = bf.parWithName(par.name)
            par.best_fit = param.best_fit
            par.isDerived = param.isDerived

    def texValues(self, formatter, p):
        param = self.parWithName(p.name)
        if not param is None:
            lims = param.limits[1]
            sf = 3
            if param.twotail:
                res, plus_str, minus_str = formatter.namesigFigs(param.mean, lims[1] - param.mean, lims[0] - param.mean)
                res = res + '^{' + plus_str + '}_{' + minus_str + '}'
            elif param.lim_bot and not param.lim_top:
                res = '< ' + formatter.formatNumber(lims[1], sf)
            elif param.lim_top and not param.lim_bot:
                res = '> ' + formatter.formatNumber(lims[0], sf)
            else: res = '---'
            if len(self.columns) > 1:  # add best fit too
                rangew = (lims[1] - lims[0]) / 10
                bestfit = formatter.namesigFigs(param.best_fit, rangew, -rangew)[0]
                return [res, bestfit]
            return [res]
        else: return None


class convergeStats(paramResults):
    def loadFromFile(self, filename):
        textFileLines = self.fileList(filename)
        self.R_eigs = []
        for i in range(len(textFileLines)):
            if textFileLines[i].find('var(mean)') >= 0:
                for line in textFileLines[i + 1:]:
                    if len(line.strip()) == 0:return
                    try: self.R_eigs.append(line.split()[1])
                    except: self.R_eigs.append('1e30')

    def worstR(self):
        return self.R_eigs[len(self.R_eigs) - 1]


