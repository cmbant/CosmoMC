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

    def __init__(self, ncol, caption=''):
        self.lines = []
        self.caption = caption
        self.sig_figs = 4
        self.lines.append('\\begin{table}')

        self.lines.append('\\begin{tabular} {' + "|| l | c"*ncol + '||}')
        self.lines.append('\\hline')

    def addLine(self):
        self.lines.append('\\hline')


    def addBestFitParamRow(self, row, ncol):
        line = " & ".join(self.bestFitParamTex(param) for param in row)
        self.lines.append(line + ' & & ' * (ncol - len(row)) + '\\\\')
        self.lines.append('\\hline')

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

#        form = "%0." + str(sf - 1) + "e"
#        st = form % value
#        expo = int(st.split('e')[1])
#        if (expo > 0): sf = max(0, sf - expo - 1)
#        return (('%' + ('', '+')[wantSign] + '.' + str(sf) + 'g') % float(st))

    def labelAndLatexResult(self, param, res):
        return '$' + param.label + ' ' * max(0, 22 - len(param.label)) + '$ & $' + res + '$' + ' ' * max(0, 22 - len(res))

    def labelAndResult(self, param, res):
        return '$' + param.label + ' ' * max(0, 22 - len(param.label)) + '$ &  ' + res + ' ' + ' ' * max(0, 22 - len(res))

    def bestFitParamTex(self, param):
        res = self.formatNumber(param.best_fit)
        return self.labelAndLatexResult(param, res)

    def margeTex(self, param):
        lims = param.limits[1]
        if param.twotail:
            res = self.plus_minus(param.mean, lims[1] - param.mean, lims[0] - param.mean)
        elif param.lim_bot and not param.lim_top:
            res = '< ' + self.formatNumber(lims[1], 2)
        elif param.lim_top and not param.lim_bot:
            res = '> ' + self.formatNumber(lims[0], 2)
        else: return self.labelAndResult(param, '---')
        return self.labelAndLatexResult(param, res)


    def plus_minus(self, value, limplus, limminus):
        frac = limplus / value
        err_sf = 2
        if value >= 20 and frac > 0.1 and limplus >= 2: err_sf = 1

        plus_str = self.formatNumber(limplus, err_sf, True)
        minus_str = self.formatNumber(limminus, err_sf, True)
        sf = self.sig_figs
        if frac > 0.1 and value < 100 and value >= 20: sf = 2
        elif frac > 0.01 and value < 1000: sf = 3
#        if abs(value) < 1 and limplus - limminus > abs(value): sf = 2
        res = self.formatNumber(value, sf)
        while abs(value) < 1 and max(self.decimal_places(plus_str), self.decimal_places(minus_str)) \
                < self.decimal_places(res):
            sf -= 1
            res = self.formatNumber(value, sf)

        while self.decimal_places(plus_str) > self.decimal_places(res):
            sf += 1
            res = self.formatNumber(value, sf)
        res = res + '^{' + plus_str + '}_{' + minus_str + '}'
        return res


    def addMargeStatParamRow(self, row, ncol):
        line = " & ".join(self.margeTex(param) for param in row)
        self.lines.append(line + ' & & ' * (ncol - len(row)) + '\\\\')
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


class paramResults():

    def __init__(self, filename=''):

        self.results = []
        self.logLike = None
        if filename != '': self.loadFromFile(filename)

    def numParams(self):
        return len(self.results)

    def numDerived(self):
        return len([1 for info in self.results if info.isDerived])

    def numNonDerived(self):
        return len([1 for info in self.results if not info.isDerived])

    def parWithNumber(self, num):
        for par in self.results:
            if par.num == num:
                return par
        return None

    def resultTable(self, ncol, caption=''):
        numrow = self.numParams() / ncol
        if self.numParams() % ncol != 0: numrow += 1
        rows = []
        for par in self.results[0:numrow]:
            rows.append([par])
        for col in range(1, ncol):
            for i in range(numrow * col, min(numrow * (col + 1), self.numParams())):
                rows[i - numrow * col].append(self.results[i]);
        r = resultTable(ncol, caption)
        for row in rows: self.addParamTableRow(r, row, ncol)
        r.endTable()
        return r



class bestFit(paramResults):

    def loadFromFile(self, filename):
        textFileHandle = open(filename)
        textFileLines = textFileHandle.readlines()
        textFileHandle.close()
        first = textFileLines[0].strip().split('=')
        self.logLike = float(first[1].strip())
        isFixed = False
        isDerived = False
        for line in textFileLines[2:]:
            if len(line.strip()) == 0:
                isFixed = not isFixed
                isDerived = True
                continue
            if not isFixed:
                param = paramNames.paramInfo()
                param.isDerived = isDerived
                (param.number, param.best_fit, param.name, param.label) = [s.strip() for s in line.split(None, 3)]
                param.number = int(param.number)
                param.best_fit = float(param.best_fit)
                self.results.append(param)

    def addParamTableRow(self, resultTable, row, ncol):
        resultTable.addBestFitParamRow(row, ncol)



class margeStats(paramResults):

    def loadFromFile(self, filename):
        textFileHandle = open(filename)
        textFileLines = textFileHandle.readlines()
        textFileHandle.close()
        for i in range(1, len(textFileLines)):
            line = textFileLines[i]
            if len(line.strip()) == 0:
                lims = textFileLines[i + 1].split(':')[1]
                self.limits = [float(s.strip()) for s in lims.split(';')]
                for line in textFileLines[i + 2:]:
                    par = self.parWithNumber(int(line.split()[0].strip()))
                    tag = line[27:39].strip()
                    par.twotail = tag == 'two tail'
                    par.lim_top = tag == '< one tail' or tag == 'no tail'
                    par.lim_bot = tag == '> one tail' or tag == 'no tail'
                break
            param = paramNames.paramInfo()
            items = [s.strip() for s in line.split(None, 7)]
            param.num = int(items[0])
            param.mean = float(items[1])
            param.err = float(items[2])

            param.limits = [[float(s) for s in items[3:5]], [float(s) for s in items[5:7]]]
            param.label = items[7]
            self.results.append(param)

    def addParamTableRow(self, resultTable, row, ncol):
        resultTable.addMargeStatParamRow(row, ncol)

