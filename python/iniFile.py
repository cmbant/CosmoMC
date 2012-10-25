# AL Apr 11
import os

class iniFile:


    def __init__(self, filename=''):

        self.params = dict()
        self.readOrder = []
        if filename != '': self.readFile(filename)


    def readFile(self, filename):
        includes = []

        textFileHandle = open(filename)
        # Remove blanck lines and comment lines from the python list of lists.
        for line in textFileHandle:
            s = line.strip()
            if s == 'END':break
            if s.startswith('#'):continue
            elif s.startswith('INCLUDE('):
                includes.append(s[s.find('(') + 1:s.rfind(')')])
            elif s != '':
                eq = s.find('=')
                if eq >= 0:
                    key = s[0:eq].strip()
                    if key in self.params:
                        raise Exception('Error: duplicate key: ' + key)
                    value = s[eq + 1:].strip()
                    self.params[key] = value;
                    self.readOrder.append(key);

        textFileHandle.close()
        for ffile in includes:
            if os.path.isabs(ffile):
                self.readFile(ffile)
            else:
                self.readFile(os.path.join(os.path.dirname(filename), ffile))

        return self.params

    def saveFile(self, filename):
        f = open(filename, 'w')
        f.write("\n".join(self.fileLines()))
        f.close()

    def fileLines(self):

        def asIniText(value):
            if type(value) == type(''): return value
            if type(value) == type(True):
                return str(value)[0]
            return str(value)

        parameterLines = []
        keys = self.params.keys()
        keys.sort()

        for key in self.readOrder:
            if key in keys:
                parameterLines.append(key + '=' + asIniText(self.params[key]));
                keys.remove(key)
        for key in keys:
            parameterLines.append(key + '=' + asIniText(self.params[key]));

        return parameterLines


    def replaceTags(self, placeholder, text):

        for key in self.params:
            self.params[key] = self.params[key].replace(placeholder, text);

            return self.params
