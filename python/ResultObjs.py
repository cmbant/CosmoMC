import paramNames, decimal

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

class resultTable():

    def __init__(self, ncol, caption='', results=None, compareNum=0, numResults=1, border='||', position='!ht'):
        self.lines = []
        self.caption = caption
        self.sig_figs = 4
        self.numResults = numResults
        self.compareNum = compareNum
        self.ncol = ncol
        self.results = results
        self.boldBaseParameters = True
        self.lines.append('\\begin{table}[' + position + ']')
        self.lines.append('\\begin{tabular} {' + (border + " l " + "| c"* numResults) * ncol + border + '}')
        self.endofrow = '\\\\'
        self.paramText = 'Parameter'
        self.bestFitText = 'Best fit'
        self.addLine()

    def addLine(self):
        return self.lines.append('\\hline')

    def margeHeaderRow(self):
        res = '& ' + self.paramText + ' & '
        paramHead = ''
        if self.numResults > 1: paramHead += self.bestFitText + ' & '
        paramHead += str(round(float(self.results.limits[1]) * 100)) + '\\% limits'
        res = res + paramHead * (self.compareNum + 1)
        self.lines.append((res * self.ncol)[1:] + self.endofrow)
        self.addLine()
        self.addLine()

    def addBestFitParamRow(self, row):
        line = " & ".join(self.bestFitParamTex(param) for param in row)
        self.lines.append(line + ' & & ' * (self.ncol - len(row)) + self.endofrow)
        self.addLine()

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

    def labelAndResult(self, param, res, latexResult=False):
        return self.textAsColumn(param.label, True) + ' & ' + self.textAsColumn(res, latexResult)

    def bestFitParamTex(self, param):
        res = self.formatNumber(param.best_fit)
        return self.labelAndResult(param, res, True)

    def paramLabelColumn(self, param):
        return  self.textAsColumn(param.label, True, separator=True, bold=not param.isDerived)

    def margeTexRow(self, param):
        return self.paramLabelColumn(param) + self.margeAndBestFitTexParam(param)

    def margeAndBestFitTexParam(self, param):
        lims = param.limits[1]
        sf = 3
        if param.twotail:
            res, plus_str, minus_str = self.namesigFigs(param.mean, lims[1] - param.mean, lims[0] - param.mean)
            res = res + '^{' + plus_str + '}_{' + minus_str + '}'
        elif param.lim_bot and not param.lim_top:
            res = '< ' + self.formatNumber(lims[1], sf)
        elif param.lim_top and not param.lim_bot:
            res = '> ' + self.formatNumber(lims[0], sf)
        else: res = '---'
        paramResultText = ''
        if self.numResults == 2:  # add best fit too
            rangew = (lims[1] - lims[0]) / 10
            bres = self.namesigFigs(param.best_fit, rangew, -rangew)[0]
            paramResultText += self.textAsColumn(bres, True, separator=True)
        paramResultText += self.textAsColumn(res, res != '---')
        return paramResultText


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

    def addMargeStatParamRow(self, row):
        line = " & ".join(self.margeTexRow(param) for param in row)
        self.lines.append(line + ' & ' * ((1 + self.numResults) * (self.ncol - len(row))) + '\\\\')
        self.lines.append('\\hline')


    def endTable(self):
        self.lines.append('\\end{tabular}')
        if self.caption != '': self.lines.append('\\caption{' + self.caption + '}')
        self.lines.append('\\end{table}')

    def tableTex(self):
        return "\n".join(self.lines)

    def writeTable(self, fname):
        textFileHandle = open(fname, 'w')
        textFileHandle.write(self.tableTex())
        textFileHandle.close()


class paramResults(paramNames.paramList):

    def resultTable(self, ncol, caption='', compareResults=[]):
        numrow = self.numParams() / ncol
        if self.numParams() % ncol != 0: numrow += 1
        rows = []
        for par in self.names[0:numrow]:
            rows.append([par])
        for col in range(1, ncol):
            for i in range(numrow * col, min(numrow * (col + 1), self.numParams())):
                rows[i - numrow * col].append(self.names[i]);
        hasBestFit = hasattr(self, 'logLike')
        r = resultTable(ncol, caption, self, compareNum=len(compareResults) , numResults=(1, 2)[hasBestFit])

        self.addHeader(r)
        for row in rows: self.addParamTableRow(r, row)
        r.endTable()
        return r



class bestFit(paramResults):

    def loadFromFile(self, filename):
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
                continue
            if not isFixed:
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

    def addParamTableRow(self, resultTable, row):
        resultTable.addBestFitParamRow(row)

    def addHeader(self, table):pass



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

    def addParamTableRow(self, resultTable, row):
        resultTable.addMargeStatParamRow(row)

    def addBestFit(self, bf):
        self.logLike = bf.logLike
# the next line deletes parameters not in best-fit; this is good e.g. to get rid of yhe from importance sampled result
        self.names = [x for x in self.names if bf.parWithName(x.name) is not None]
        for par in self.names:
            param = bf.parWithNumber(par.number)
            par.best_fit = param.best_fit
            par.isDerived = param.isDerived

    def addHeader(self, table):
        table.margeHeaderRow()


class convergeStats(paramResults):
    def loadFromFile(self, filename):
        textFileLines = self.fileList(filename)
        self.R_eigs = []
        for i in range(len(textFileLines)):
            if textFileLines[i].find('var(mean)') >= 0:
                for line in textFileLines[i + 1:]:
                    if len(line.strip()) == 0:return
                    self.R_eigs.append(line.split()[1])

    def worstR(self):
        return self.R_eigs[len(self.R_eigs) - 1]

