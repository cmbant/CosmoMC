import os, pickle, random
import numpy as np
from getdist import paramNames

def lastModified(files):
    return max([os.path.getmtime(fname) for fname in files if os.path.exists(fname)])

def chainFiles(root, chain_indices=None, ext='.txt'):
    index = -1
    files = []
    while True:
        index += 1
        fname = root + ('', '_' + str(index))[index > 0] + ext
        if index > 0 and not os.path.exists(fname): break
        if (chain_indices is None or index in chain_indices) and os.path.exists(fname):
            files.append(fname)
    return files


def loadChains(root=None, chain_indices=None, ignore_rows=0, ignore_frac=0, no_cache=False, separate_chains=False, no_stat=False, ext='.txt'):
    c = chains(root, ignore_rows=ignore_rows)
    if root is not None:
        files = chainFiles(root, chain_indices, ext=ext)
        c.cachefile = root + ext + '.pysamples'
        if not separate_chains and chain_indices is None and not no_cache and os.path.exists(c.cachefile) and lastModified(files) < os.path.getmtime(c.cachefile):
            with open(c.cachefile, 'rb') as inp:
                return pickle.load(inp)
        elif c.loadChains(root, files):
            c.removeBurnFraction(ignore_frac)
            c.deleteFixedParams()
            if not no_stat: c.getChainsStats()
            if not separate_chains:
                c.makeSingle()
                c.updateChainBaseStatistics()
                with open(c.cachefile, 'wb') as output:
                    pickle.dump(c, output, pickle.HIGHEST_PROTOCOL)
            return c
        return None

def loadGridChain(rootdir, model, data, importance=None, ignore_frac=0):
    tags = [model, data]
    if importance is not None: tags += ['post', importance]
    path = model + os.sep + data + os.sep + "_".join(tags)
    print path
    return loadChains(rootdir + os.sep + path, ignore_frac=ignore_frac, separate_chains=False, no_stat=True)


def getSignalToNoise(C, noise=None, R=None, eigs_only=False):
    if R is None:
        if noise is None: raise Exception('Must give noise or rotation R')
#        w, U = np.linalg.eigh(noise)  # noise = (U * w).dot(U.T)
#        w = w ** (-0.5)
#        R = (U * w).dot(U.T)
        R = np.linalg.inv(np.linalg.cholesky(noise))

    M = np.dot(R, C).dot(R.T)
    if eigs_only:
        return np.linalg.eigvalsh(M)
    else:
        w, U = np.linalg.eigh(M)
        U = np.dot(U.T, R)
        return w, U


class chain(object):
    def __init__(self, filename=None, ignore_rows=0, samples=None, weights=None, loglikes=None):
        if filename:
            self.setColData(np.loadtxt(filename, skiprows=ignore_rows))
        else:
            self.setSamples(samples)
            self.weights = weights
            self.loglikes = loglikes


    def setColData(self, coldata):
        self.setSamples(coldata[:, 2:])
        self.weights = coldata[:, 0]
        self.loglikes = coldata[:, 1]

    def setSamples(self, samples):
        self.samples = samples
        self.n = self.samples.shape[1]
        self.means = []

    def getCov(self, n):
        if len(self.means) == 0: self.getMeans()
        cov = np.zeros((n, n))
        for i in range(n):
            weightdiff = (self.samples[:, i] - self.means[i]) * self.weights
            for j in range(i, n):
                cov[i, j] = np.dot(weightdiff , self.samples[:, j] - self.means[j]) / self.norm
                cov[j, i] = cov[i, j]
        self.cov = cov
        return cov

    def getMeans(self):
        self.norm = np.sum(self.weights)
        self.means = np.empty(self.n)
        for i in range(self.n):
            self.means[i] = np.dot(self.samples[:, i], self.weights) / self.norm
        return self.means

class paramConfidenceData(object):
    pass

class parSamples(object): pass

class chains(object):

    def __init__(self, root=None, ignore_rows=0, jobItem=None):
        self.jobItem = jobItem
        self.precision = '%.8e'
        self.ignore_lines = ignore_rows
        self.root = root
        self.samples = None
        self.hasNames = os.path.exists(root + '.paramnames')
        self.needs_update = True
        if self.hasNames:
            self.paramNames = paramNames.paramNames(root + '.paramnames')
            self.getParamIndices()

    def getParamIndices(self):
        index = dict()
        for i, name in enumerate(self.paramNames.names):
            index[name.name] = i
        self.index = index

    def setParams(self, obj):
        for i, name in enumerate(self.paramNames.names):
            setattr(obj, name.name, self.samples[:, i])

    def getParams(self):
        pars = parSamples()
        self.setParams(pars)
        return pars

    def valuesForParam(self, param, force_array=False):
        if isinstance(param, np.ndarray): return param
        if isinstance(param, basestring): param = [param]
        results = []
        for par in param:
            if isinstance(par, basestring):
                ix = self.index[par]
                results.append(self.samples[:, ix])
            else: results.append(par)
        if len(results) == 1 and not force_array: return results[0]
        return results

    def weighted_sum(self, paramVec, where=None):
        if where is None: return np.dot(paramVec, self.weights)
        return np.dot(paramVec[where], self.weights[where])

    def get_norm(self, where=None):
        if where is None:
            if not hasattr(self, 'norm'): self.norm = np.sum(self.weights)
            return self.norm
        else:
            return np.sum(self.weights[where])

    def mean(self, paramVec, where=None):
        return self.weighted_sum(paramVec, where) / self.get_norm(where)

    def var(self, paramVec, where=None):
        return self.weighted_sum((paramVec - self.mean(paramVec, where)) ** 2 , where) / self.get_norm(where)

    def std(self, paramVec, where=None):
        return np.sqrt(self.var(paramVec, where))

    def twoTailLimits(self, paramVec, confidence):
        limits = np.array([(1 - confidence) / 2, 1 - (1 - confidence) / 2])
        return self.confidence(paramVec, limits)


    def initParamConfidenceData(self, paramVec, start=0, end=None, weights=None):
        if weights is None: weights = self.weights
        d = paramConfidenceData()
        d.paramVec = self.valuesForParam(paramVec)[start:end]
        d.norm = np.sum(weights[start:end])
        d.indexes = d.paramVec.argsort()
        weightsort = weights[start + d.indexes]
        d.cumsum = np.cumsum(weightsort)
        return d

    def confidence(self, paramVec, limfrac, upper=False, start=0, end=None, weights=None):
        if isinstance(paramVec, paramConfidenceData):
            d = paramVec
        else:
            d = self.initParamConfidenceData(paramVec, start, end, weights)

        if not upper: target = d.norm * limfrac
        else: target = d.norm * (1 - limfrac)
        ix = np.searchsorted(d.cumsum, target)
        return d.paramVec[d.indexes[np.minimum(ix, d.indexes.shape[0] - 1)]]

    def cov(self, paramVecs):
        paramVecs = self.valuesForParam(paramVecs, force_array=True)
        diffs = [vec - self.mean(vec) for vec in paramVecs]
        n = len(paramVecs)
        cov = np.zeros((n, n))
        for i, diff1 in enumerate(diffs):
            weightdiff = diff1 * self.weights
            for j, diff2 in enumerate(diffs[:i + 1]):
                cov[i, j] = np.dot(weightdiff , diff2) / self.norm
                cov[j, i] = cov[i, j]
        return cov

    def corr(self, paramVecs, cov=None):
        if cov is not None:
            corr = np.copy(cov)
        else:
            corr = self.cov(paramVecs)
        diag = [np.sqrt(corr[i, i]) for i in range(len(paramVecs))]
        for i, di in enumerate(diag):
            if di:
                corr[i, :] /= di
                corr[:, i] /= di
        return corr

    def getSignalToNoise(self, params, noise=None, R=None, eigs_only=False):
        """
        Returns w, M, where w is the eigenvalues of the signal to noise (small means better constrained)
        """
        C = self.cov(params)
        return getSignalToNoise(C, noise, R, eigs_only)


    def updateChainBaseStatistics(self):
        self.norm = np.sum(self.weights)
        self.numrows = self.samples.shape[0]
        self.num_vars = self.samples.shape[1]
        self.getParamIndices()
        self.needs_update = False

    def addDerived(self, paramVec, **kwargs):
        self.samples = np.c_[self.samples, paramVec]
        self.needs_update = True
        return self.paramNames.addDerived(**kwargs)
    #    self.updateChainBaseStatistics()

    def loadChains(self, root, files):
        self.chains = []
        for fname in files:
                print fname
                self.chains.append(chain(fname, self.ignore_lines))
        if len(self.chains) == 0: print 'loadChains - no chains found for ' + root
        return len(self.chains) > 0

    def getChainsStats(self, nparam=None):
        nparam = nparam or self.paramNames.numNonDerived()
        for chain in self.chains: chain.getCov(nparam)
        norm = np.sum([chain.norm for chain in self.chains])
        self.means = np.zeros(nparam)
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
        D = np.linalg.eigvalsh(np.dot(R, meanscov).dot(R.T))
        self.GelmanRubin = max(np.real(D))
        print 'R-1 = ', self.GelmanRubin


    def makeSingle(self):
        self.chain_offsets = np.cumsum(np.array([0] + [chain.samples.shape[0] for chain in self.chains]))
        self.weights = np.hstack((chain.weights for chain in self.chains))
        self.loglikes = np.hstack((chain.loglikes for chain in self.chains))
        self.samples = np.vstack((chain.samples for chain in self.chains))
        del(self.chains)
        return self

    def getSeparateChains(self):
        if hasattr(self, 'chains'): return self.chains
        chainlist = []
        for off1, off2 in zip(self.chain_offsets[:-1], self.chain_offsets[1:]):
            chainlist.append(chain(samples=self.samples[off1:off2], weights=self.weights[off1:off2], loglikes=self.loglikes[off1:off2]))
        return chainlist

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
            ix = int(round(chain.samples.shape[0] * ignore_frac))
            chain.samples = chain.samples[ix:, :]
            if chain.weights is not None:
                chain.weights = chain.weights[ix:]
            if chain.loglikes is not None:
                chain.loglikes = chain.loglikes[ix:]

    def deleteFixedParams(self):
        fixed = []
        if self.samples is not None:
            for i in range(self.samples.shape[1]):
                if np.all(self.samples[:, i] == self.samples[0, i]): fixed.append(i)
            self.samples = np.delete(self.samples, fixed, 1)
        else:
            chain = self.chains[0]
            for i in range(chain.n):
                if np.all(chain.samples[:, i] == chain.samples[0, i]): fixed.append(i)
            for chain in self.chains:
                chain.setSamples(np.delete(chain.samples, fixed, 1))
        if self.hasNames:
            self.paramNames.deleteIndices(fixed)
            self.getParamIndices()
#           print 'Non-derived parameters: ', [name.name for name in self.paramNames.names[0:self.paramNames.numNonDerived()]]


    def writeSingle(self, root):
        np.savetxt(root + '.txt', np.hstack((self.weights.reshape(-1, 1), self.loglikes.reshape(-1, 1), self.samples)), fmt=self.precision)
        self.paramNames.saveAsText(root + '.paramnames')

    def thin_indices(self, factor, weights=None):
        tot = 0
        i = 0
        if weights is None:  weights = self.weights
        numrows = len(weights)
        norm1 = np.sum(weights)
        weights = weights.astype(np.int)
        norm = np.sum(weights)

        if abs(norm - norm1) > 1e-4:
                raise Exception('Can only thin with integer weights')
        if factor <> int(factor):
                raise Exception('Thin factor must be integer')

        if factor >= np.max(weights):
            cumsum = np.cumsum(weights) / int(factor)
            _, thin_ix = np.unique(cumsum, return_index=True)
        else:
            thin_ix = np.empty(norm / factor, dtype=np.int)
            ix = 0
            mult = weights[i]
            while i < numrows:
                if mult + tot < factor:
                    tot += mult
                    i += 1
                    if i < numrows: mult = weights[i]
                else:
                    thin_ix[ix] = i
                    ix += 1
                    if mult == factor - tot:
                        i += 1
                        if i < numrows: mult = weights[i]
                    else:
                        mult -= (factor - tot)
                    tot = 0

        return thin_ix

    def singleSamples_indices(self):
        max_weight = np.max(self.weights)
        thin_ix = []
        for i in range(self.numrows):
            P = self.weights[i] / max_weight
            if random.random() < P:
                thin_ix.append(i)
        return np.array(thin_ix, dtype=np.int)

    def makeSingleSamples(self):
        thin_ix = self.singleSamples_indices()
        self.samples = self.samples[thin_ix, :]
        self.weights = np.ones(len(thin_ix))
        self.loglikes = self.loglikes[thin_ix]
        self.norm = np.sum(self.weights)

    def thin(self, factor):
        thin_ix = self.thin_indices(factor)
        self.samples = self.samples[thin_ix, :]
        self.weights = np.ones(len(thin_ix))
        self.loglikes = self.loglikes[thin_ix]
        self.norm = len(thin_ix)

    def filter(self, where):
        self.samples = self.samples[where, :]
        self.weights = self.weights[where]
        self.loglikes = self.loglikes[where]
        self.norm = np.sum(self.weights)

    def reweight_plus_logLike(self, logLikes):
        scale = np.min(logLikes)
        self.loglikes += logLikes
        self.weights *= np.exp(-(logLikes - scale))
        self.norm = np.sum(self.weights)

