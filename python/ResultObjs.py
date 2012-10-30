import paramNames

class resultTable():

    def __init__(self, ncol, caption=''):
        self.lines = []
        self.caption = caption
        self.lines.append('\\begin{table}')

        self.lines.append('\\begin{tabular} {' + "|| l | c"*ncol + '||}')
        self.lines.append('\\hline')

    def addLine(self):
        self.lines.append('\\hline')


    def addBestFitParamRow(self, row, ncol):
        line = " & ".join(self.bestFitParamTex(param) for param in row)
        self.lines.append(line + ' & & ' * (ncol - len(row)) + '\\\\')

    def bestFitParamTex(self, param):
        value = float(param.value)
        sf = 4
        form = "%0." + str(sf - 1) + "e"
        st = form % value
        num, expo = st.split('e')
        expo = int(expo)
        if (expo > 0): sf = max(0, sf - expo)
        return '$' + param.label + '$ & ' + (('%0.' + str(sf) + 'f') % float(st))

    def endTable(self):
        self.lines.append('\\hline')
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
        self.logLike = first[1].strip()
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
                (param.number, param.value, param.name, param.label) = [s.strip() for s in line.split(None, 3)]
                self.results.append(param)
#                print param.label, param.value

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
                    par.twotail = line.find('two tail') >= 0
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

