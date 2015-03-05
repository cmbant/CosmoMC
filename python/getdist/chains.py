import os, pickle, random
import numpy as np
from getdist.paramNames import paramNames
from scipy.signal import fftconvolve

class WeightedSampleError(Exception):
    pass

def lastModified(files):
    return max([os.path.getmtime(fname) for fname in files if os.path.exists(fname)])

def convolve1D(x, y, mode):
    if min(x.shape[0], y.shape[0]) > 1000:
        return fftconvolve(x, y, mode)
    else:
        return np.convolve(x, y, mode)

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
                c.updateBaseStatistics()
                with open(c.cachefile, 'wb') as output:
                    pickle.dump(c, output, pickle.HIGHEST_PROTOCOL)
            return c
        return None

def getSignalToNoise(C, noise=None, R=None, eigs_only=False):
    if R is None:
        if noise is None: raise WeightedSampleError('Must give noise or rotation R')
        R = np.linalg.inv(np.linalg.cholesky(noise))

    M = np.dot(R, C).dot(R.T)
    if eigs_only:
        return np.linalg.eigvalsh(M)
    else:
        w, U = np.linalg.eigh(M)
        U = np.dot(U.T, R)
        return w, U

def covToCorr(cov, copy=True):
    if copy: cov = cov.copy()
    for i, di in enumerate(np.sqrt(cov.diagonal())):
        if di:
            cov[i, :] /= di
            cov[:, i] /= di
    return cov

class paramConfidenceData(object): pass

class parSamples(object): pass


class WeightedSamples(object):
    def __init__(self, filename=None, ignore_rows=0, samples=None, weights=None, loglikes=None, name_tag=None):
        if filename:
            self.setColData(np.loadtxt(filename, skiprows=ignore_rows))
            self.name_tag = name_tag or os.path.basename(filename)
        else:
            self.setSamples(samples, weights, loglikes)
            self.name_tag = name_tag

    def setColData(self, coldata):
        self.setSamples(coldata[:, 2:], coldata[:, 0], coldata[:, 1])

    def getName(self):
        return self.name_tag

    def setSamples(self, samples, weights=None, loglikes=None):
        self.weights = weights
        self.loglikes = loglikes
        self.samples = samples
        if samples is not None:
            if isinstance(samples, (list, tuple)):
                samples = np.hstack([x.reshape(-1, 1) for x in samples])
            elif len(samples.shape) == 1:
                samples = np.atleast_2d(samples).transpose()
            self.samples = samples
            self.n = self.samples.shape[1]
            self.numrows = self.samples.shape[0]
        self._weightsChanged()

    def changeSamples(self, samples):
        self.setSamples(samples, self.weights, self.loglikes)

    def _weightsChanged(self):
        if self.weights is not None:
            self.norm = np.sum(self.weights)
        elif self.samples is not None:
            self.weights = np.ones(self.numrows)
            self.norm = np.float64(self.numrows)
        self.means = None
        self.mean_loglike = None
        self.diffs = None
        self.fullcov = None
        self.correlationMatrix = None
        self.vars = None
        self.sddev = None


    def _makeParamvec(self, par):
        if isinstance(par, (int, long)):
            if par >= 0 and par < self.n:
                return self.samples[:, par]
            elif par == -1:
                if self.loglikes is None:
                    raise WeightedSampleError('Samples do not have logLikes (par=-1)' % (par))
                return self.loglikes
            elif par == -2:
                return self.weights
            else: raise WeightedSampleError('Parameter %i does not exist' % (par))
        return par

    def getCov(self, nparam=None):
        if self.fullcov is None:
            self.setCov()
        return self.fullcov[:nparam, :nparam]

    def setCov(self):
        self.fullcov = self.cov()
        return self.fullcov

    def getCorrelationMatrix(self):
        if self.correlationMatrix is None:
            self.correlationMatrix = covToCorr(self.getCov())
        return self.correlationMatrix

    def setMeans(self):
        self.means = self.weights.dot(self.samples) / self.norm
        if self.loglikes is not None: self.mean_loglike = self.weights.dot(self.loglikes) / self.norm
        return self.means

    def getMeans(self):
        if self.means is None:
            return self.setMeans()
        return self.means

    def getVars(self):
        if self.means is None: self.setMeans()
        self.vars = np.empty(self.n)
        for i in range(self.n):
            self.vars[i] = self.weights.dot((self.samples[:, i] - self.means[i]) ** 2) / self.norm
        self.sddev = np.sqrt(self.vars)
        return self.vars

    def setDiffs(self):
        self.diffs = self.mean_diffs()
        return self.diffs

    def getWeightedAutocorrelation(self, paramVec, maxOff=None):
        """ get auto covariance in weight units; 
            divide by var to normalize
            multiply by n/norm to get in sample point (row) units
       """
        if maxOff is None: maxOff = self.n - 1
        d = self.mean_diff(paramVec) * self.weights
        corr = convolve1D(d, d[::-1], 'full')[-d.size:]
        return corr[0:maxOff + 1] * d.size / (self.norm * np.arange(d.size, d.size - maxOff - 1, -1))

    def getEffectiveSamples(self, j=0, min_corr=0.05):
        corr = self.getWeightedAutocorrelation(j, self.numrows / 10)
        corr /= self.var(j)
        ix = np.argmin(corr > min_corr * corr[0])
        N = corr[0] + 2 * np.sum(corr[1:ix])
        return self.get_norm() / N

    def weighted_sum(self, paramVec, where=None):
        paramVec = self._makeParamvec(paramVec)
        if where is None: return self.weights.dot(paramVec)
        return np.dot(paramVec[where], self.weights[where])

    def get_norm(self, where=None):
        if where is None:
            if self.norm is None: self.norm = np.sum(self.weights)
            return self.norm
        else:
            return np.sum(self.weights[where])

    def mean(self, paramVec, where=None):
        return self.weighted_sum(paramVec, where) / self.get_norm(where)

    def var(self, paramVec, where=None):
        if where is not None:
            return np.dot(self.mean_diff(paramVec, where) ** 2, self.weights[where]) / self.get_norm(where)
        else:
            return np.dot(self.mean_diff(paramVec) ** 2, self.weights) / self.get_norm()

    def std(self, paramVec, where=None):
        return np.sqrt(self.var(paramVec, where))

    def cov(self, pars=None, where=None):
        diffs = self.mean_diffs(pars, where)
        if pars is None:
            pars = range(self.n)
        n = len(pars)
        cov = np.empty((n, n))
        if where is not None:
            weights = self.weights[where]
        else:
            weights = self.weights
        for i, diff in enumerate(diffs):
            weightdiff = diff * weights
            for j in range(i, n):
                cov[i, j] = weightdiff.dot(diffs[j])
                cov[j, i] = cov[i, j]
        cov /= self.get_norm(where)
        return cov

    def corr(self, pars=None):
        return self.covToCorr(self.cov(pars))

    def mean_diff(self, paramVec, where=None):
        if isinstance(paramVec, (int, long)) and paramVec >= 0 and where is None:
            if self.diffs is not None:
                return self.diffs[paramVec]
            return self.samples[:, paramVec] - self.getMeans()[paramVec]
        paramVec = self._makeParamvec(paramVec)
        if where is None:
            return paramVec - self.mean(paramVec)
        else:
            return paramVec[where] - self.mean(paramVec, where)

    def mean_diffs(self, pars=None, where=None):
        if pars is None: pars = self.n
        if isinstance(pars, (int, long)) and pars >= 0 and where is None:
            means = self.getMeans()
            return [self.samples[:, i] - means[i] for i in range(pars)]
        return [self.mean_diff(i, where) for i in pars]

    def twoTailLimits(self, paramVec, confidence):
        limits = np.array([(1 - confidence) / 2, 1 - (1 - confidence) / 2])
        return self.confidence(paramVec, limits)

    def initParamConfidenceData(self, paramVec, start=0, end=None, weights=None):
        if weights is None: weights = self.weights
        d = paramConfidenceData()
        d.paramVec = self._makeParamvec(paramVec)[start:end]
        d.norm = np.sum(weights[start:end])
        d.indexes = d.paramVec.argsort()
        weightsort = weights[start + d.indexes]
        d.cumsum = np.cumsum(weightsort)
        return d

    def confidence(self, paramVec, limfrac, upper=False, start=0, end=None, weights=None):
        """ 
        Raw sample confidence limits, not using kernel densities
        """
        if isinstance(paramVec, paramConfidenceData):
            d = paramVec
        else:
            d = self.initParamConfidenceData(paramVec, start, end, weights)

        if not upper: target = d.norm * limfrac
        else: target = d.norm * (1 - limfrac)
        ix = np.searchsorted(d.cumsum, target)
        return d.paramVec[d.indexes[np.minimum(ix, d.indexes.shape[0] - 1)]]

    def getSignalToNoise(self, params, noise=None, R=None, eigs_only=False):
        """
        Returns w, M, where w is the eigenvalues of the signal to noise (small means better constrained)
        """
        C = self.cov(params)
        return getSignalToNoise(C, noise, R, eigs_only)

    def thin_indices(self, factor, weights=None):
        """
        Indices to make single weight 1 samples. Assumes intefer weights
        """
        if weights is None:  weights = self.weights
        numrows = len(weights)
        norm1 = np.sum(weights)
        weights = weights.astype(np.int)
        norm = np.sum(weights)

        if abs(norm - norm1) > 1e-4:
                raise WeightedSampleError('Can only thin with integer weights')
        if factor <> int(factor):
                raise WeightedSampleError('Thin factor must be integer')

        if factor >= np.max(weights):
            cumsum = np.cumsum(weights) / int(factor)
            _, thin_ix = np.unique(cumsum, return_index=True)
        else:
            tot = 0
            i = 0
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

    def randomSingleSamples_indices(self):
        max_weight = np.max(self.weights)
        thin_ix = []
        for i in range(self.numrows):
            P = self.weights[i] / max_weight
            if random.random() < P:
                thin_ix.append(i)
        return np.array(thin_ix, dtype=np.int)

    def thin(self, factor):
        thin_ix = self.thin_indices(factor)
        self.setSamples(self.samples[thin_ix, :], loglikes=self.loglikes[thin_ix])

    def filter(self, where):
        self.setSamples(self.samples[where, :], self.weights[where], self.loglikes[where])

    def reweightAddingLogLikes(self, logLikes):
        scale = np.min(logLikes)
        self.loglikes += logLikes
        self.weights *= np.exp(-(logLikes - scale))
        self._weightsChanged()

    def cool(self, cool):
        MaxL = np.max(self.loglikes)
        newL = self.loglikes * cool
        self.weights = self.weights * np.exp(-(newL - self.loglikes) - (MaxL * (1 - cool)))
        self.loglikes = newL
        self._weightsChanged()

    def deleteZeros(self):
        self.filter(self.weights == 0)

    def deleteFixedParams(self):
        fixed = []
        for i in range(self.samples.shape[1]):
            if np.all(self.samples[:, i] == self.samples[0, i]): fixed.append(i)
        self.changeSamples(np.delete(self.samples, fixed, 1))

    def removeBurn(self, remove=0.3):
        if remove >= 1:
            ix = int(remove)
        else:
            ix = int(round(self.numrows * remove))
        if self.weights is not None:
            self.weights = self.weights[ix:]
        if self.loglikes is not None:
            self.loglikes = self.loglikes[ix:]
        self.changeSamples(self.samples[ix:, :])


class chains(WeightedSamples):

    def __init__(self, root=None, jobItem=None, paramNamesFile=None, names=None, **kwargs):
        WeightedSamples.__init__(self, **kwargs)
        self.jobItem = jobItem
        self.precision = '%.8e'
        self.ignore_lines = float(kwargs.get('ignore_rows', 0))
        self.root = root
        paramNamesFile = paramNamesFile or str(root) + '.paramnames'
        self.needs_update = True
        self.chains = None
        self.paramNames = None
        if isinstance(paramNamesFile, paramNames):
            self.paramNames = paramNamesFile
        elif os.path.exists(paramNamesFile):
            self.paramNames = paramNames(paramNamesFile)
        elif names is not None:
            self.paramNames = paramNames(names=names)
        elif self.samples is not None:
            self.paramNames = paramNames(default=self.n)
        if self.paramNames: self.getParamIndices()


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

    def _makeParamvec(self, par):
        if isinstance(par, basestring):
                return self.samples[:, self.index[par]]
        return WeightedSamples._makeParamvec(self, par)

    def updateBaseStatistics(self):
        self.getVars()
        self.mean_mult = self.norm / self.numrows
        self.max_mult = np.max(self.weights)
        self.getParamIndices()
        self.needs_update = False

    def addDerived(self, paramVec, **kwargs):
        self.samples = np.c_[self.samples, paramVec]
        self.needs_update = True
        return self.paramNames.addDerived(**kwargs)

    def loadChains(self, root, files):
        self.chains = []
        self.name_tag = self.name_tag or root
        for fname in files:
                print fname
                self.chains.append(WeightedSamples(fname, self.ignore_lines))
        if len(self.chains) == 0:
            raise WeightedSampleError('loadChains - no chains found for ' + root)
        return len(self.chains) > 0


    def getGelmanRubinEigenvalues(self, nparam=None, chainlist=None):
        # Assess convergence in the var(mean)/mean(var) in the worst eigenvalue
        # c.f. Brooks and Gelman 1997
        if chainlist is None:
            chainlist = self.getSeparateChains()
        nparam = nparam or self.paramNames.numNonDerived()
        meanscov = np.zeros((nparam, nparam))
        means = self.getMeans()[:nparam]
        meancov = np.zeros(meanscov.shape)
        for chain in chainlist:
            diff = chain.getMeans()[:nparam] - means
            meanscov += np.outer(diff, diff)
            meancov += chain.getCov(nparam)
        meanscov /= (len(chainlist) - 1)
        meancov /= len(chainlist)
        w, U = np.linalg.eigh(meancov)
        if np.min(w) > 0:
            U /= np.sqrt(w)
            D = np.linalg.eigvalsh(np.dot(U.T, meanscov).dot(U))
            return D
        else:
            return None

    def getGelmanRubin(self, nparam=None, chainlist=None):
        return np.max(self.getGelmanRubinEigenvalues(nparam, chainlist))

    def makeSingle(self):
        self.chain_offsets = np.cumsum(np.array([0] + [chain.samples.shape[0] for chain in self.chains]))
        weights = np.hstack((chain.weights for chain in self.chains))
        loglikes = np.hstack((chain.loglikes for chain in self.chains))
        self.setSamples(np.vstack((chain.samples for chain in self.chains)), weights, loglikes)
        self.chains = None
        self.needs_update = True
        return self

    def getSeparateChains(self):
        if self.chains is not None:
            return self.chains
        chainlist = []
        for off1, off2 in zip(self.chain_offsets[:-1], self.chain_offsets[1:]):
            chainlist.append(WeightedSamples(samples=self.samples[off1:off2], weights=self.weights[off1:off2], loglikes=self.loglikes[off1:off2]))
        return chainlist

    def removeBurnFraction(self, ignore_frac):
        if self.samples is not None:
            self.removeBurn(ignore_frac)
            self.chains = None
            self.needs_update = True
        else:
            for chain in self.chains:
                chain.removeBurn(ignore_frac)

    def deleteFixedParams(self):
        if self.samples is not None:
            WeightedSamples.deleteFixedParams(self)
            self.chains = None
        else:
            fixed = []
            chain = self.chains[0]
            for i in range(chain.n):
                if np.all(chain.samples[:, i] == chain.samples[0, i]): fixed.append(i)
            for chain in self.chains:
                chain.changeSamples(np.delete(chain.samples, fixed, 1))
        self.paramNames.deleteIndices(fixed)
        self.getParamIndices()

    def writeSingle(self, root):
        np.savetxt(root + '.txt', np.hstack((self.weights.reshape(-1, 1), self.loglikes.reshape(-1, 1), self.samples)), fmt=self.precision)
        self.paramNames.saveAsText(root + '.paramnames')


