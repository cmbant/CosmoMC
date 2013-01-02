import os, iniFile, paramNames, pickle
import numpy as np


def loadChains(root=None, chain_indices=None, ignore_rows=0, ignore_frac=0, no_cache=False, separate_chains=False):
    c = chains(root, ignore_rows=ignore_rows)
    if root is not None:
        files = c.chainFiles(root, chain_indices)
        c.cachefile = root + '.pysamples'
        if not separate_chains and chain_indices is None and not no_cache and os.path.exists(c.cachefile) and c.lastModified(files) < os.path.getmtime(c.cachefile):
            with open(c.cachefile, 'rb') as inp:
                return pickle.load(inp)
        elif c.loadChains(root, files):
            c.removeBurnFraction(ignore_frac)
            c.deleteFixedParams()
            c.getChainStats()
            if not separate_chains:
                c.makeSingle()
                with open(c.cachefile, 'wb') as output:
                    pickle.dump(c, output, pickle.HIGHEST_PROTOCOL)
            return c
        return None


class chain():
    def __init__(self, filename, ignore_rows=0):
        self.setColData(np.loadtxt(filename, skiprows=ignore_rows))

    def setColData(self, coldata):
        self.coldata = coldata
        self.n = self.coldata.shape[1] - 2
        self.means = []

    def getCov(self, n):
        if len(self.means) == 0: self.getMeans()
        weights = self.coldata[:, 0]
        cov = np.zeros((n, n))
        for i in range(n):
            for j in range(i, n):
                cov[i, j] = np.dot(weights, (self.coldata[:, i + 2] - self.means[i]) * (self.coldata[:, j + 2] - self.means[j])) / self.norm
                cov[j, i] = cov[i, j]
        self.cov = cov
        return cov

    def getMeans(self):
        means = []
        weights = self.coldata[:, 0]
        self.norm = weights.sum()
        for i in range(2, self.n + 2):
            means.append(np.average(self.coldata[:, i], weights=weights))
        self.means = np.asarray(means)
        return means

class chains():

    def __init__(self, root=None, ignore_rows=0):
        self.precision = '%.8e'
        self.ignore_rows = ignore_rows
        self.root = root
        self.chains = []
        self.samples = None
        self.paramNames = paramNames.paramNames(root + '.paramnames')

    def lastModified(self, files):
        return max([os.path.getmtime(fname) for fname in files])

    def chainFiles(self, root, chain_indices=None):
        index = -1
        files = []
        while True:
            index += 1
            fname = root + ('', '_' + str(index))[index > 0] + '.txt'
            if index > 0 and not os.path.exists(fname): break
            if (chain_indices is None or index in chain_indices) and os.path.exists(fname):
                files.append(fname)
        return files

    def loadChains(self, root, files):
        self.chains = []
        for fname in files:
                print fname
                self.chains.append(chain(fname, self.ignore_rows))
        if len(self.chains) == 0: print 'loadChains - no chains found for ' + root
        return len(self.chains) > 0

    def getChainsStats(self):
        nparam = self.paramNames.numNonDerived()
        for chain in self.chains: chain.getCov(nparam)
        self.means = np.zeros(nparam)
        norm = np.sum([chain.norm for chain in self.chains])
        for chain in self.chains:
            self.means = self.means + chain.means[0:nparam] * chain.norm
        self.means /= norm
        meanscov = np.zeros((nparam, nparam))
        for i in range(nparam):
            for j in range(nparam):
                meanscov[i, j] = np.sum([chain.norm * (chain.means[i] - self.means[i]) * (chain.means[j] - self.means[j]) for chain in self.chains])
                meanscov[j, i] = meanscov[i, j]
        meanscov *= len(self.chains) / (len(self.chains) - 1) / norm
        self.meancov = np.zeros((nparam, nparam))
        for chain in self.chains:
            self.meancov += chain.cov * chain.norm
        self.meancov /= norm
        M = self.meancov
        for i in range(nparam):
            norm = np.sqrt(self.meancov[i, i])
            M[i, :] /= norm
            M[:, i] /= norm
            meanscov[:, i] /= norm
            meanscov[i, :] /= norm
        R = np.linalg.inv(np.linalg.cholesky(M))
        D = np.linalg.eigvals(np.dot(R, meanscov).dot(R.T))
        self.GelmanRubin = max(np.real(D))
        print 'R-1 = ', self.GelmanRubin


    def makeSingle(self):
        self.weights = np.hstack((chain.coldata[:, 0] for chain in self.chains))
        self.loglikes = np.hstack((chain.coldata[:, 1] for chain in self.chains))
        self.samples = np.vstack((chain.coldata[:, 2:] for chain in self.chains))
        self.numrows = self.samples.shape[0]
        del(self.chains)
        return self

    def loadWMAPChain(self, chainPath, namesFile, thinfac=1):
        params = paramNames.paramNames(namesFile)
        cols = []
        self.weights = np.loadtxt(chainPath + 'weight', usecols=(1,))
        self.loglikes = np.loadtxt(chainPath + 'neglnlike', usecols=(1,))
        if thinfac <> 1:
            thin_ix = self.thin_indices(thinfac)
            self.weights = np.ones(len(thin_ix))
            self.loglikes = self.loglikes[thin_ix]
        else: thin_ix = None

        usednames = []
        for param in params.names:
            if os.path.exists(chainPath + param.name):
                col = np.loadtxt(chainPath + param.name, usecols=(1,))
                if thin_ix is not None: col = col[thin_ix]
                cols.append(col)
                usednames.append(param)
        params.names = usednames
        self.paramNames = params
        self.samples = np.vstack(cols).transpose()
        self.numrows = self.samples.shape[0]


    def removeBurnFraction(self, ignore_frac):
        for chain in self.chains:
            chain.coldata = chain.coldata[round(chain.coldata.shape[0] * ignore_frac):, :]

    def deleteFixedParams(self):
        fixed = []
        if self.samples is not None:
            for i in range(self.samples.shape[1]):
                if np.all(self.samples[:, i] == self.samples[0, i]): fixed.append(i)
            self.samples = np.delete(self.samples, fixed, 1)
        else:
            chain = self.chains[0]
            for i in range(2, chain.n + 2):
                if np.all(chain.coldata[:, i] == chain.coldata[0, i]): fixed.append(i)
            for chain in self.chains:
                chain.setColData(np.delete(chain.coldata, fixed, 1))
        self.paramNames.deleteIndices([fix - 2 for fix in fixed])
        print 'Parameters: ', [name.name for name in self.paramNames.names[0:self.paramNames.numNonDerived()]]


    def writeSingle(self, root):
        np.savetxt(root + '.txt', np.hstack((self.weights.reshape(-1, 1), self.loglikes.reshape(-1, 1), self.samples)), fmt=self.precision)
        self.paramNames.saveAsText(root + '.paramnames')

    def thin_indices(self, factor):
        thin_ix = []
        tot = 0
        i = 0
        self.numrows = len(self.weights)
        mult = self.weights[i]
        if abs(round(mult) - mult) > 1e-4:
                raise Exception('Can only thin with integer weights')

        while i < self.numrows:
            if (mult + tot < factor):
                tot += mult
                i += 1
                if i < self.numrows: mult = self.weights[i]
            else:
                thin_ix.append(i)
                if mult == factor - tot:
                    i += 1
                    if i < self.numrows: mult = self.weights[i]
                else:
                    mult -= (factor - tot)
                tot = 0
        return thin_ix

    def thin(self, factor):
        thin_ix = self.thin_indices(factor)
        self.samples = self.samples[thin_ix, :]
        self.weights = np.ones(len(thin_ix))
        self.loglikes = self.loglikes[thin_ix]


# c = loadChains('C:\\tmp\\Planck\\chains\\base_nrun_r_planck_CAMspec_lowl_lowLike', ignore_frac=0.3, separate_chains=False)
# c.getChainsStats()

# c = chains()
# c.loadWMAPChain('C:\\tmp\\Planck\\WMAP9\\test\\', 'C:\\Work\\Dist\\git\\cosmomcplanck\\LAMBDA.paramnames', thinfac=50)
# c.thin(30)
# c.writeSingle('z:\\test_thin')
