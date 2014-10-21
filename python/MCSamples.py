# MCSamples.py

import os
import sys
import glob
import math
import logging
import numpy as np
from scipy.interpolate import interp1d, splrep, splev
from scipy.stats import norm

from chains import chains

# =============================================================================

class Ranges():

    def __init__(self, fileName=None, setParamNameFile=None):
        self.names = []
        self.mins = {}
        self.maxs = {}
        if fileName is not None: self.loadFromFile(fileName)

    def loadFromFile(self, fileName):
        self.filenameLoadedFrom = os.path.split(fileName)[1]
        f = open(fileName)
        for line in f:
            name, mini, maxi = self.readValues(line.strip())
            if name:
                self.names.append(name)
                self.mins[name] = mini
                self.maxs[name] = maxi

    def readValues(self, line):
        name, mini, maxi = None, None, None
        strings = [ s for s in line.split(" ")  if s<>'' ]
        if len(strings)==3:
            name = strings[0]
            mini = float(strings[1])
            maxi = float(strings[2])
        return name, mini, maxi

    def listNames(self):
        return self.names

    def min(self, name, error=False):
        if self.mins.has_key(name):
            return self.mins[name]
        if error: raise Exception("Name not found:" + name)
        return None

    def max(self, name, error=False):
        if self.maxs.has_key(name):
            return self.maxs[name]
        if error: raise Exception("Name not found:" + name)
        return None

# =============================================================================

class Density1D():

    SPLINE_DANGLE = 1.e30

    def __init__(self, n, spacing):
        self.n   = n
        self.X   = np.zeros(n)
        self.P   = np.zeros(n)
        self.spacing = float(spacing)

    def InitSpline(self):
        self.spl = splrep(self.X, self.P, s=0)

    def Prob(self, x):
        # Assuming x is a single value
        return splev([x], self.spl)

    def Limits(self, p):

        # Return values
        mn, mx = 0., 0.
        lim_bot, lim_top = False, False

        factor = 100
        bign = (self.n-1)*factor + 1
        vecx = self.X[0] + np.arange(bign)*self.spacing / factor
        grid = splev(vecx, self.spl)

        norm  = np.sum(grid)
        norm  = norm - (0.5*self.P[-1]) - (0.5*self.P[0])

        try_t = max(grid)
        try_b = 0
        try_last = -1

        while True:
            trial = (try_b + try_t) / 2
            trial_sum = np.sum(grid[grid>=trial])
            if (trial_sum < p*norm):
                try_t = (try_b + try_t) / 2
            else:
                try_b = (try_b + try_t) / 2
            if (math.fabs(trial_sum/try_last - 1) < 1e-4): break
            try_last = trial_sum

        trial = (try_b + try_t) / 2

        lim_bot = (grid[0] >= trial)
        if (lim_bot):
            mn = self.P[0]
        else:
            for i in range(bign):
                if (grid[i] > trial):
                    mn = self.X[0] + (i-1) * self.spacing / factor
                    break

        lim_top = (grid[-1] >= trial)
        if (lim_top):
            mx = self.P[-1]
        else:
            indexes = range(bign)
            indexes.reverse()
            for i in indexes:
                if (grid[i] > trial):
                    mx = self.X[0] + (i-1) * self.spacing / factor
                    break

        return mn, mx, lim_bot, lim_top

# =============================================================================

class MCSamples(chains):

    def __init__(self, root=None, ignore_rows=0):
        chains.__init__(self, root, ignore_rows)

        self.colix = []
        self.limmin = []
        self.limmax = []

        self.has_limits = []
        self.has_limits_bot = []
        self.has_limits_top = []

        self.has_markers = []
        self.markers = []
        self.isused = []

        self.ranges = None
        self.ReadRanges()

        if hasattr(self, 'index'):
            # index is a dict such as { name : index }
            # index2name is a dict such as { index: name }
            self.index2name = dict((v,k) for k, v in self.index.iteritems())
            self.index2name.keys().sort()

        self.density1D = None

        # Other variables
        self.no_plots = False
        self.num_vars = 0
        self.num_bins = 0
        self.num_bins_2D = 0
        self.smooth_scale_1D = -1.0
        self.smooth_scale_2D = -1.0
        self.max_mult = 0
        self.numsamp = 0
        self.plot_data_dir = ""
        self.rootname = ""
        self.rootdirname = ""
        self.num_contours = 0
        self.meanlike = 0.
        self.plot_meanlikes = False
        self.mean_loglikes = False
        self.shade_meanlikes = False
        self.indep_thin = 0
        self.contours = []

        self.credible_interval_threshold = 0.05

        self.subplot_size_inch  = 3.0
        self.subplot_size_inch2 = self.subplot_size_inch
        self.subplot_size_inch3 = self.subplot_size_inch
        self.out_dir = ""

        self.LowerUpperLimits = None
        self.covmat_dimension = 0
        self.max_split_tests = 4
        self.force_twotail = False
        self.mean_mult = 0

        self.corr_length_thin = 0
        self.corr_length_steps = 15


    def ComputeContours(self, ini=None):
        self.num_contours = 2
        if ini:
            self.num_contours = ini.int('num_contours', 2)
        self.contours = []
        # how small the end bin must be relative to max to use two tail
        self.max_frac_twotail = []
        if ini is None: return
        for i in range(1, self.num_contours+1):
            self.contours.append(ini.float('contour'+str(i)))
            max_frac = ini.float('max_frac_twotail'+str(i),
                                 math.exp(-1.0*math.pow(norm.ppf((1-self.contours[i-1])/2), 2)/2))
            self.max_frac_twotail.append(max_frac)



    def ComputeLimits(self, ini=None):

        bin_limits = ""
        if ini: bin_limits = ini.string('all_limits')

        indexes = self.index2name.keys()
        indexes.sort()
        nvars = len(indexes)

        self.limmin = nvars * [0.]
        self.limmax = nvars * [0.]

        self.has_limits = nvars * [False]
        self.has_limits_bot = nvars * [False]
        self.has_limits_top = nvars * [False]

        self.has_markers = nvars * [False]
        self.markers = nvars * [0.]

        for ix in indexes:
            name = self.index2name[ix]
            mini = self.ranges.min(name)
            maxi = self.ranges.max(name)
            if (mini is not None and maxi is not None and mini<>maxi):
                self.limmin[ix] = mini
                self.limmax[ix] = maxi
                self.has_limits_top[ix] = True
                self.has_limits_bot[ix] = True
            if (bin_limits<>''):
                line = bin_limits
            else:
                line = ''
                if ini and ini.params.has_key('limits[%s]'%name.strip()):
                    line = ini.string('limits[%s]'%name.strip())
            if (line<>''):
                limits = [ s for s in line.split(' ') if s<>'' ]
                if len(limits)==2:
                    if limits[0]<>'N':
                        self.limmin[ix] = float(limits[0])
                        self.has_limits_bot[ix] = True
                    if limits[1]<>'N':
                        self.limmax[ix] = float(limits[1])
                        self.has_limits_top[ix] = True
            if ini and ini.params.has_key('marker[%s]'%name.strip()):
                line = ini.string('marker[%s]'%name.strip())
                if (line<>''):
                    self.has_markers = True
                    self.markers[ix] = float(line)


    def AdjustPriors(self):
        sys.exit('You need to write the AdjustPriors function in MCSamples.py first!')
        #print "Adjusting priors"
        #ombh2 = self.samples[0]  # ombh2 prior
        #chisq = (ombh2 - 0.0213)**2/0.001**2
        #self.weights *= np.exp(-chisq/2)
        #self.loglikes += chisq/2

    def MapParameters(self, invars):
        sys.exit('Need to write MapParameters routine first')

    def CoolChain(self, cool):
        print 'Cooling chains by ', cool
        MaxL = np.max(self.loglikes)
        newL = self.loglikes * cool
        newW = self.weights * np.exp(-(newL-self.loglikes) - (MaxL*(1-cool)))
        self.weights = np.hstack(newW)
        self.loglikes = np.hstack(newL)

    def DeleteZeros(self):
        indexes = np.where(self.weights==0)
        self.weights = np.delete(self.weights, indexes)
        self.loglikes = np.delete(self.loglikes, indexes)
        self.samples = np.delete(self.samples, indexes, axis=0)

    def SortColData(self, icol=None):
        """
        Sort coldata in order of likelihood
        """
        if icol is None: return
        if icol==0:
            indexes = self.weights.argsort()
        elif icol==1:
            indexes = self.loglikes.argsort()
        self.weights = self.weights[indexes]
        self.loglikes = self.loglikes[indexes]
        self.samples = self.samples[indexes]

    def ComputeNumSamp(self):
        self.numsamp = np.sum(self.weights)


    def MakeSingleSamples(self, filename="", single_thin=None, writeDataToFile=True):
        """
        Make file of weight-1 samples by choosing samples
        with probability given by their weight.
        """
        if writeDataToFile:
            textFileHandle = open(filename, 'w')
            maxmult = np.max(self.weights)
            for i in range(self.numrows):
                rand = np.random.random_sample()
                if (rand <= self.weights[i]/maxmult/single_thin):
                    textFileHandle.write("%16.7E"%(1.0))
                    textFileHandle.write("%16.7E"%(self.loglikes[i]))
                    for j in range(self.num_vars):
                        textFileHandle.write("%16.7E"%(self.samples[i][j]))
                    textFileHandle.write("\n")
            textFileHandle.close()
        else:
            # return data
            return self.loglikes, self.samples

    def WriteThinData(self, fname, thin_ix, cool):
        nparams = self.samples.shape[1]
        if(cool<>1): print 'Cooled thinned output with temp: ', cool
        MaxL = np.max(self.loglikes)
        textFileHandle = open(fname, 'w')
        i = 0
        for thin in thin_ix:
            if (cool<>1):
                newL = self.loglikes[thin] * cool
                textFileHandle.write("%16.7E"%(
                        math.exp(-(newL - self.loglikes[thin])-MaxL*(1-cool))))
                textFileHandle.write("%16.7E"%(newL))
                for j in nparams:
                    textFileHandle.write("%16.7E"%(self.samples[i][j]))
            else:
                textFileHandle.write("%f"%(1.))
                textFileHandle.write("%f"%(self.loglikes[thin]))
                for j in nparams:
                    textFileHandle.write("%16.7E"%(self.samples[i][j]))
            i += 1

        textFileHandle.close()
        print  'Wrote ', len(thin_ix), ' thinned samples'

    # ThinData(self, fac) => chains.thin_indices(factor)

    def GetCovMatrix(self, writeDataToFile=True):

        nparamNonDerived = self.paramNames.numNonDerived()
        paramVecs = [ self.samples[:, i] for i in range(nparamNonDerived) ]
        self.covmatrix = self.cov(paramVecs)

        if writeDataToFile:
            fname = self.rootdirname + ".covmat"
            textFileHandle = open(fname, "w")
            paramNames = []
            for i in range(nparamNonDerived):
                paramNames.append(self.paramNames.parWithName(self.index2name[i]).name)
            textFileHandle.write("# %s\n"%(" ".join(paramNames)))
            for i in range(nparamNonDerived):
                for j in range(nparamNonDerived):
                    textFileHandle.write("%17.7E"%self.covmatrix[i][j])
                textFileHandle.write("\n")
            textFileHandle.close()

        nparam = self.paramNames.numParams()
        paramVecs = [ self.samples[:, i] for i in range(nparam) ]
        self.corrmatrix = self.corr(paramVecs)
        if writeDataToFile:
            np.savetxt(self.rootdirname + ".corr", self.corrmatrix, fmt="%17.7E")

    # MostCorrelated2D ?

    def GetFractionIndices(self, weights, n):
        nrows = weights.shape[0]
        numsamp = np.sum(weights)

        tot = 0
        aim = numsamp / n
        num = 0

        fraction_indices = []
        fraction_indices.append(0)

        for i in range(nrows):
            tot = tot + weights[i]
            if (tot>aim):
                num = num + 1
                fraction_indices.append(i)
                if (num==n): break
                aim = aim + numsamp / n

        fraction_indices.append(nrows)
        return fraction_indices

    # ConfidVal(self, ix, limfrac, upper) => chains.confidence()

    def PCA(self, pars, param_map, normparam_num):
        """
        Perform principle component analysis. In other words,
        get eigenvectors and eigenvalues for normalized variables
        with optional (log) mapping.
        """

        print 'Doing PCA for ', len(pars),' parameters'
        filename = self.rootdirname + ".PCA"
        textFileHandle = open(filename, "w")
        textFileHandle.write('PCA for parameters:\n')

        if (normparam_num<>0):
            if (normparam_num in pars):
                normparam = pars.index(normparam_num)
        else:
            normparam = -1

        n = len(pars)
        corrmatrix = np.zeros((n, n))
        PCdata = self.samples[:, pars]
        PClabs = []

        PCmean = np.zeros(n)
        sd = np.zeros(n)
        newmean = np.zeros(n)
        newsd = np.zeros(n)

        doexp = False
        for i in range(n):
            name = self.index2name[pars[i]]
            label = self.paramNames.parWithName(name).label
            if (param_map[i]=='L'):
                doexp = True
                PCdata[:, i] = np.log(PCdata[:, i])
                PClabs.append("ln("+label+")")
            elif (param_map[i]=='M'):
                doexp = True
                PCdata[:, i] = np.log(-1.0*PCdata[:, i])
                PClabs.append("ln(-"+label+")")
            else:
                PClabs.append(label)
            textFileHandle.write("%10s :%s\n"%(str(pars[i]+1), str(PClabs[i])))

            PCmean[i] = np.sum(self.weights*PCdata[:, i]) / self.norm
            PCdata[:, i] = PCdata[:, i] - PCmean[i]
            sd[i] = np.sqrt(np.sum(self.weights*np.power(PCdata[:, i], 2)) / self.norm)
            if (sd[i]<>0): PCdata[:, i] = PCdata[:, i] / sd[i]
            corrmatrix[i][i] = 1

        textFileHandle.write('\n')
        textFileHandle.write('Correlation matrix for reduced parameters\n')
        for i in range(n):
            for j in range(n):
                corrmatrix[j][i] = np.sum(self.weights*PCdata[:, i]*PCdata[:, j]) / self.norm
                corrmatrix[i][j] = corrmatrix[j][i]
            textFileHandle.write('%4i :'%(pars[i]+1))
            for j in range(n):
                textFileHandle.write('%8.4f'%corrmatrix[j][i])
            textFileHandle.write('\n')

        u = corrmatrix
        evals, evects = np.linalg.eig(u)
        isorted = evals.argsort()
        u = np.transpose(evects[:, isorted]) # redefining u

        textFileHandle.write('\n')
        textFileHandle.write('e-values of correlation matrix\n')
        for i in range(n):
            isort = isorted[i]
            textFileHandle.write('PC%2i: %8.4f\n'%(i+1, evals[isort]))

        textFileHandle.write('\n')
        textFileHandle.write('e-vectors\n')
        for j in range(n):
            textFileHandle.write('%3i:'%(pars[j]+1))
            for i in range(n):
                isort = isorted[i]
                textFileHandle.write('%8.4f'%(evects[j][isort]))
            textFileHandle.write('\n')

        if (normparam<>-1):
            # Set so parameter normparam has exponent 1
            for i in range(n):
                u[i, :] = u[i, :] / u[i, normparam] * sd[normparam]
        else:
            # Normalize so main component has exponent 1
            for i in range(n):
                maxi = np.abs(u[i, :]).argmax()
                u[i, :] = u[i, :] / u[i, maxi] * sd[maxi]

        nrows = PCdata.shape[0]
        for i in range(nrows):
            PCdata[i, :] = np.dot(u, PCdata[i, :])
            if (doexp): PCdata[i, :] = np.exp(PCdata[i, :])

        textFileHandle.write('\n')
        textFileHandle.write('Principle components\n')

        for i in range(n):
            isort = isorted[i]
            textFileHandle.write('PC%i (e-value: %f)\n'%(i+1, evals[isort]))
            for j in range(n):
                name = self.index2name[pars[j]]
                label = self.paramNames.parWithName(name).label
                if (param_map[j] in ['L', 'M']):
                    expo = "%f"%(1.0/sd[j]*u[i][j])
                    if (param_map[j]=="M"):
                        div = "%f"%(-np.exp(PCmean[j]))
                    else:
                        div = "%f"%(np.exp(PCmean[j]))
                    textFileHandle.write('[%f]  (%s/%s)^{%s}\n'%(u[i][j], label, div, expo))
                else:
                    expo = "%f"%(sd[j]/u[i][j])
                    if (doexp):
                        textFileHandle.write('[%f]   exp((%s-%f)/%s)\n'%(u[i][j], label, PCmean[j], expo))
                    else:
                        textFileHandle.write('[%f]   (%s-%f)/%s)\n'%(u[i][j], label, PCmean[j], expo))

            newmean[i] = np.sum(self.weights*PCdata[:, i]) / self.norm
            newsd[i] = np.sqrt(np.sum(self.weights*np.power(PCdata[:, i]-newmean[i], 2)) / self.norm)
            textFileHandle.write('          = %f +- %f\n'%(newmean[i], newsd[i]))
            textFileHandle.write('ND limits: %9.3f%9.3f%9.3f%9.3f\n'%(
                    np.min(PCdata[0:self.ND_cont1, i]), np.max(PCdata[0:self.ND_cont1, i]),
                    np.min(PCdata[0:self.ND_cont2, i]), np.max(PCdata[0:self.ND_cont2, i])))
            textFileHandle.write('\n')

        # Find out how correlated these components are with other parameters
        textFileHandle.write('Correlations of principle components\n')
        l = [ "%8i"%i for i in range(1,n+1) ]
        textFileHandle.write('%s\n'%("".join(l)))

        for i in range(n):
            PCdata[:, i] = (PCdata[:, i] - newmean[i]) / newsd[i]

        for j in range(n):
            textFileHandle.write('PC%2i'%(j+1))
            for i in range(n):
                textFileHandle.write('%8.3f'%(
                        np.sum(self.weights*PCdata[:, i]*PCdata[:, j]) / self.norm))
            textFileHandle.write('\n')

        for j in range(self.num_vars):
            if not self.isused[j]: continue
            textFileHandle.write('%4i'%(j+1))
            for i in range(n):
                textFileHandle.write('%8.3f'%(
                        np.sum(self.weights*PCdata[:, i]
                               *(self.samples[:, j]-self.means[j])/self.sddev[j]) / self.norm))

            name = self.index2name[j]
            label = self.paramNames.parWithName(name).label
            textFileHandle.write('   (%s)\n'%(label))

        textFileHandle.close()


    def GetUsedCols(self):
        nparams = self.samples.shape[1]
        self.isused = np.ndarray(nparams, dtype=bool)
        self.isused[:] = False
        for ix in range(nparams):
            self.isused[ix] = not np.all(self.samples[:, ix]==self.samples[0][ix])

    def ComputeMultiplicators(self):
        self.mean_mult = self.norm / self.numrows
        self.max_mult = (self.mean_mult*self.numrows) / min(self.numrows/2, 500)
        outliers = len(self.weights[np.where(self.weights>self.max_mult)])
        if (outliers<>0):
            print 'outlier fraction ', float(outliers)/self.numrows
        self.max_mult = np.max(self.weights)
        self.numsamp = np.sum(self.weights)


    def GetUsedColsChains(self):
        if not hasattr(self, 'chains') or len(self.chains)==0: return
        nparams = self.chains[0].coldata.shape[1] - 2
        self.isused = np.ndarray(nparams, dtype=bool)
        self.isused[:] = False
        for ix in range(nparams):
            self.isused[ix] = not np.all(self.chains[0].coldata[:, ix+2]==self.chains[0].coldata[0][ix+2])

    def DoConvergeTests(self, limfrac):
        """
        Do convergence tests.
        """
        if not hasattr(self, 'chains'): return

        self.GetUsedColsChains()

        # Weights
        weights = np.hstack((chain.coldata[:, 0] for chain in self.chains))

        # Get statistics for individual chains, and do split tests on the samples

        fname = self.rootdirname + '.converge'
        textFileHandle = open(fname, 'w')

        num_chains_used = len(self.chains)
        if (num_chains_used>1):
            print 'Number of chains used = ', num_chains_used

        nparam = self.paramNames.numParams()
        for chain in self.chains: chain.getCov(nparam)
        mean = np.zeros(nparam)
        norm = np.sum([chain.norm for chain in self.chains])
        nrows = np.sum([ chain.coldata.shape[0] for chain in self.chains ])
        for chain in self.chains:
            mean = mean + chain.means[0:nparam] * chain.norm
        mean /= norm

        fullmean = np.zeros(nparam)
        for j in range(nparam):
            for chain in self.chains:
                fullmean[j] += np.sum(chain.coldata[:, 0] * chain.coldata[:, j+2]) / norm

        fullvar = np.zeros(nparam)
        for j in range(nparam):
            for chain in self.chains:
                fullvar[j] += np.sum(chain.coldata[:, 0] * np.power(chain.coldata[:, j+2]-fullmean[j], 2)) / norm

        if (num_chains_used>1):
            textFileHandle.write(" \n")
            textFileHandle.write(" Variance test convergence stats using remaining chains\n")
            textFileHandle.write(" param var(chain mean)/mean(chain var)\n")
            textFileHandle.write(" \n")

            between_chain_var = np.zeros(nparam)
            in_chain_var = np.zeros(nparam)
            chain_means = np.zeros((num_chains_used, nparam))

            #norm = np.sum([chain.norm for chain in self.chains])
            for j in range(nparam):
                if not self.isused[j]: continue

                # Get stats for individual chains - the variance of the means over the mean of the variances
                for i in range(num_chains_used):
                    chain = self.chains[i]
                    chain_means[i][j] = np.sum(chain.coldata[:, 0]*chain.coldata[:, j+2]) / chain.norm
                    between_chain_var[j] += np.power(chain_means[i][j] - mean[j], 2)
                    in_chain_var[j] += np.sum(chain.coldata[:, 0] * np.power(chain.coldata[:, j+2]-chain_means[i][j], 2))

                between_chain_var[j] /= (num_chains_used-1)
                in_chain_var[j] /= norm
                label = self.paramNames.names[j].label
                textFileHandle.write("%3i%13.5f  %s\n"%(j+1, between_chain_var[j]/in_chain_var[j], label))

        nparam = self.paramNames.numNonDerived()
        self.covmat_dimension = nparam
        if (num_chains_used>1) and (self.covmat_dimension>0):
            # Assess convergence in the var(mean)/mean(var) in the worst eigenvalue
            # c.f. Brooks and Gelman 1997

            meanscov = np.zeros((nparam, nparam))
            cov = np.zeros((nparam, nparam))

            for j in range(nparam):
                for k in range(nparam):
                    meanscov[k, j] = 0.
                    cov[k, j] = 0.
                    for i in range(num_chains_used):
                        chain = self.chains[i]
                        cov[k, j] += np.sum(
                            chain.coldata[:, 0] *
                            (chain.coldata[:, j+2] - chain_means[i][j]) *
                            (chain.coldata[:, k+2] - chain_means[i][k]) )
                        meanscov[k, j] += (chain_means[i][j] - mean[j]) * (chain_means[i][k] - mean[k])
                    meanscov[j, k] = meanscov[k, j]
                    cov[j, k] = cov[k, j]

            meanscov[:, :] /= (num_chains_used-1)
            cov[:, :] /= norm

            invertible = np.isfinite(np.linalg.cond(cov))
            if (invertible):
                R = np.linalg.inv(np.linalg.cholesky(cov))
                D = np.linalg.eigvalsh(np.dot(R, meanscov).dot(R.T))
                GelmanRubin = max(np.real(D))

                textFileHandle.write("\n")
                textFileHandle.write("var(mean)/mean(var) for eigenvalues of covariance of means of orthonormalized parameters\n")
                for jj in range(nparam):
                    textFileHandle.write("%3i%13.5f\n"%(jj+1, D[jj]))
                print " var(mean)/mean(var), remaining chains, worst e-value: R-1 = %13.5F\n"%GelmanRubin
            else:
                print 'WARNING: Gelman-Rubin covariance not invertible'

        # Do tests for robustness under using splits of the samples
        # Return the rms ([change in upper/lower quantile]/[standard deviation])
        # when data split into 2, 3,.. sets
        textFileHandle.write("\n")
        textFileHandle.write("Split tests: rms_n([delta(upper/lower quantile)]/sd) n={2,3,4}:\n")
        textFileHandle.write("i.e. mean sample splitting change in the quantiles in units of the st. dev.\n")
        textFileHandle.write("\n")

        # Need these values here
        self.weights = np.hstack((chain.coldata[:, 0] for chain in self.chains))
        self.norm = np.sum(self.weights)

        split_tests = {}
        nparam = self.paramNames.numParams()
        for j in range(nparam):
            if not self.isused[j]: continue
            coldata = np.hstack((chain.coldata[:, j+2] for chain in self.chains))
            for endb in [0, 1]:
                for split_n in range(2, self.max_split_tests+1):
                    frac = self.GetFractionIndices(self.weights, split_n)
                    split_tests[split_n] = 0.
                    confid = self.confidence(coldata, (1-limfrac)/2., endb==0)
                    for i in range(split_n):
                        split_tests[split_n] = split_tests[split_n] + math.pow(self.confidence(coldata[frac[i]:frac[i+1]], (1-limfrac)/2., endb==0, start=frac[i], end=frac[i+1])-confid, 2)

                    split_tests[split_n] = math.sqrt(split_tests[split_n]/split_n/fullvar[j])
                if (endb==0):
                    typestr = 'upper'
                else:
                    typestr = 'lower'

                textFileHandle.write("%3i"%(j+1))
                for split_n in range(2, self.max_split_tests+1):
                    textFileHandle.write("%9.4f"%(split_tests[split_n]))
                label = self.paramNames.names[j].label
                textFileHandle.write("  %s %s\n"%(label, typestr))

        # Now do Raftery and Lewis method
        # See http://www.stat.washington.edu/tech.reports/raftery-lewis2.ps
        # Raw non-importance sampled chains only
        thin_fac = np.zeros(num_chains_used)

        tran = np.zeros( (2,2,2), dtype=np.int )
        tran2 = np.zeros( (2,2), dtype=np.int )

        epsilon = 0.001

        if (np.all(np.abs(weights-np.round(np.where(weights>0.6, weights, 0.6)))<1e-4)):

            nburn = np.zeros(num_chains_used, dtype=np.int)
            markov_thin = np.zeros(num_chains_used, dtype=np.int)
            hardest = -1
            hardestend = 0
            for ix in range(num_chains_used):
                chain = self.chains[ix]
                thin_fac[ix] = int(round(np.max(chain.coldata[:, 0])))

                for j in range(self.covmat_dimension):
                    coldata = np.hstack((chain.coldata[:, j] for chain in self.chains))
                    if (self.force_twotail or not self.has_limits[j]):
                        for endb in [0, 1]:
                            # Get binary chain depending on whether above or below confidence value
                            u = self.confidence(chain.coldata[:, j], (1-limfrac)/2, endb==0)
                            while(True):
                                thin_ix = self.thin_indices(thin_fac[ix])
                                thin_rows = len(thin_ix)
                                if (thin_rows<2): break
                                binchain = np.zeros(thin_rows)
                                binchain[:] = 1
                                indexes = np.where(coldata[thin_ix] >= u)
                                binchain[indexes] = 0

                                tran[:,:,:] = 0
                                # Estimate transitions probabilities for 2nd order process
                                for i in range(2, thin_rows):
                                    tran[binchain[i-2]][binchain[i-1]][binchain[i]] += 1

                                # Test whether 2nd order is better than Markov using BIC statistic
                                g2 = 0
                                for i1 in [0, 1]:
                                    for i2 in [0, 1]:
                                        for i3 in [0, 1]:
                                            if (tran[i1][i2][i3]<>0):
                                                fitted = float(
                                                    (tran[i1][i2][0]+tran[i1][i2][1]) *
                                                    (tran[0][i2][i3]+tran[1][i2][i3]) ) \
                                                / float(tran[0][i2][0]+tran[0][i2][1]+
                                                        tran[1][i2][0]+tran[1][i2][1])
                                                focus = float(tran[i1][i2][i3])
                                                g2 += math.log(focus/fitted) * focus
                                g2 = g2 * 2

                                if (g2 - math.log(float(thin_rows-2)) * 2 < 0): break
                                thin_fac[ix] += 1

                            # Get Markov transition probabilities for binary processes
                            if (np.sum(tran[:, 0, 1])==0 or np.sum(tran[:, 1, 0])==0):
                                thin_fac[ix] = 0
                                #goto 203

                            alpha = np.sum(tran[:, 0, 1]) / float(np.sum(tran[:, 0, 0])+np.sum(tran[:, 0, 1]))
                            beta = np.sum(tran[:, 1, 0]) / float(np.sum(tran[:, 1, 0])+np.sum(tran[:, 1, 1]))
                            probsum = alpha + beta
                            tmp1 = math.log(probsum*epsilon/max(alpha,beta)) / math.log(abs(1.0-probsum))
                            if (int(tmp1+1)*thin_fac[ix]>nburn[ix]):
                                nburn[ix] = int(tmp1+1)*thin_fac[ix]
                                hardest = j
                                hardestend = endb

                markov_thin[ix] = thin_fac[ix]

                # Get thin factor to have independent samples rather than Markov
                hardest = max(hardest, 0)
                u = self.confidence(chain.coldata[:, hardest], (1-limfrac)/2, hardestend==0)
                thin_fac[ix] += 1

                while(True):
                    thin_ix = self.thin_indices_chain(chain.coldata[:, 0], thin_fac[ix])
                    thin_rows = len(thin_ix)
                    if (thin_rows<2): break
                    binchain = np.zeros(thin_rows)
                    binchain[:] = 1

                    coldata = np.hstack((chain.coldata[:, hardest] for chain in self.chains))
                    coldata = coldata[thin_ix]
                    indexes = np.where(coldata >= u)
                    binchain[indexes] = 0

                    tran2[:,:] = 0
                    # Estimate transitions probabilities for 2nd order process
                    for i in range(1, thin_rows):
                        tran2[binchain[i-1]][binchain[i]] += 1

                    # Test whether independence is better than Markov using BIC statistic
                    g2 = 0
                    for i1 in [0, 1]:
                        for i2 in [0, 1]:
                            if (tran2[i1][i2]<>0):
                                fitted = float(
                                    (tran2[i1][0]+tran2[i1][1]) *
                                    (tran2[0][i2]+tran2[1][i2]) ) / float(thin_rows-1)
                                focus = float(tran2[i1][i2])
                                if (fitted<=0 or focus<=0):
                                    print 'Raftery and Lewis estimator had problems'
                                    return
                                g2 += math.log(focus/fitted) * focus
                    g2 = g2 * 2

                    if (g2 - math.log(float(thin_rows-1))< 0): break

                    thin_fac[ix] += 1
#goto 203
                if (thin_rows<2): thin_fac[ix] = 0

            textFileHandle.write("\n")
            textFileHandle.write("Raftery&Lewis statistics\n")
            textFileHandle.write("\n")
            textFileHandle.write("chain  markov_thin  indep_thin    nburn\n")

            # Computation of mean_mult
            self.mean_mult = norm / nrows

            for ix in range(num_chains_used):
                if (thin_fac[ix]==0):
                    textFileHandle.write("%4i      Not enough samples\n"%ix)
                else:
                    textFileHandle.write("%4i%12i%12i%12i"%(
                            ix, markov_thin[ix], thin_fac[ix], nburn[ix]))

            if (not np.all(thin_fac!=0)):
                print 'RL: Not enough samples to estimate convergence stats'
            else:
                print 'RL: Thin for Markov: ', np.max(markov_thin)
                self.indep_thin = np.max(thin_fac)
                print 'RL: Thin for indep samples:  ', str(self.indep_thin)
                print 'RL: Estimated burn in steps: ', np.max(nburn), ' (', int(round(np.max(nburn)/self.mean_mult)), ' rows)'

            # Get correlation lengths
            textFileHandle.write("\n")
            textFileHandle.write("Parameter auto-correlations as function of step separation\n")
            textFileHandle.write("\n")
            if (self.corr_length_thin<>0):
                autocorr_thin = self.corr_length_thin
            else:
                if (self.indep_thin==0):
                    autocorr_thin = 20
                elif (self.indep_thin<=30):
                    autocorr_thin = 5
                else:
                    autocorr_thin = 5 * (self.indep_thin/30)

            thin_ix = self.thin_indices(autocorr_thin)
            thin_rows = len(thin_ix)
            maxoff = int(min(self.corr_length_steps, thin_rows/(autocorr_thin*num_chains_used)))

            corrs = np.zeros([maxoff, nparam])
            for off in range(maxoff):
                for i in range(off, thin_rows):
                    for j in range(nparam):
                        if not self.isused[j]: continue
                        coldata = np.hstack((chain.coldata[:, j] for chain in self.chains))
                        corrs[off][j] += (coldata[thin_ix[i]]-fullmean[j]) * (coldata[thin_ix[i-off]]-fullmean[j])
                for j in range(nparam):
                    if not self.isused[j]: continue
                    corrs[off][j] /= (thin_rows-off) / fullvar[j]

            if (maxoff>0):
                for i in range(maxoff):
                    textFileHandle.write("%8i"%((i+1)*autocorr_thin))
                textFileHandle.write("\n")
                for j in range(nparam):
                    if (self.isused[j]):
                        name = self.index2name.get(j, "NOTFOUND")
                        label = self.paramNames.parWithName(name).label
                        textFileHandle.write("%3i"%j+1)
                        for i in range(maxoff):
                            textFileHandle.write("%8.3f"%corrs[i][j])
                        textFileHandle.write("%s\n"%label)
        textFileHandle.close()


    def PreComputeDensity(self, j):

        # Return values
        width, smooth_1D, end_edge = 0, 0, 0

        ix = j # ix = self.colix[j]
        paramVec = self.samples[:, ix]

        self.param_min[j] = np.min(paramVec)
        self.param_max[j] = np.max(paramVec)
        self.range_min[j] = self.confidence(paramVec, 0.0005, upper=False)
        self.range_max[j] = self.confidence(paramVec, 0.0005, upper=True)
        width = (self.range_max[j]-self.range_min[j])/(self.num_bins+1)
        if (width==0):
            return width, smooth_1D, end_edge

        logging.debug("Smooth scale ... ")
        if (self.smooth_scale_1D<=0):
            # Automatically set smoothing scale from rule of thumb for Gaussian
            opt_width = 1.06 / math.pow(max(1.0, self.numsamp/self.max_mult), 0.2) * self.sddev[j]
            smooth_1D = opt_width / width * abs(self.smooth_scale_1D)
            if (smooth_1D<0.5):
                print 'Warning: num_bins not large enough for optimal density'
            smooth_1D = max(1.0, smooth_1D)
        elif (self.smooth_scale_1D<1.0):
            smooth_1D = self.smooth_scale_1D * self.sddev[j]/width
            if (smooth_1D<1):
                print 'Warning: num_bins not large enough to well sample smoothed density'
        else:
            smooth_1D = self.smooth_scale_1D

        end_edge = int(round(smooth_1D * 2))

        logging.debug("Limits ... ")
        if (self.has_limits_bot[ix]):
            if ( (self.range_min[j]-self.limmin[ix]>(width*end_edge)) and
                 (self.param_min[j]-self.limmin[ix]>(width*smooth_1D)) ):
                # long way from limit
                self.has_limits_bot[ix] = False
            else:
                self.range_min[j] = self.limmin[ix]

        if (self.has_limits_top[ix]):
            if ((self.limmax[ix]-self.range_max[j]>(width*end_edge)) and
                (self.limmax[ix]-self.param_max[j]>(width*smooth_1D)) ):
                self.has_limits_top[ix] = False
            else:
                self.range_max[j] = self.limmax[ix]
        self.has_limits[ix] = self.has_limits_top[ix] or self.has_limits_bot[ix]

        if (self.has_limits_top[ix]):
            self.center[j] = self.range_max[j]
        else:
            self.center[j] = self.range_min[j]

        self.ix_min[j] = int(round((self.range_min[j] - self.center[j])/width))
        self.ix_max[j] = int(round((self.range_max[j] - self.center[j])/width))

        if (not self.has_limits_bot[ix]): self.ix_min[j] -= end_edge
        if (not self.has_limits_top[ix]): self.ix_max[j] += end_edge

        return width, smooth_1D, end_edge


    def Get1DDensity(self, j, writeDataToFile=True):

        ix = j # ix = self.colix[j]
        logging.debug("index is %i"%j)

        paramVec = self.samples[:, ix]

        width, smooth_1D, end_edge = self.PreComputeDensity(j)

        if (width==0):
            print "Warning width is 0"
            return

        # Using index mapping for f90 arrays with non standard indexes.

        # In f90, binsraw(ix_min(j):ix_max(j))
        binsraw = np.zeros( self.ix_max[j] - self.ix_min[j] + 1 )

        fine_fac = 10
        winw = int(round(2.5 * fine_fac * smooth_1D))
        fine_edge = winw + fine_fac * end_edge
        fine_width = width / fine_fac

        imin = int(round((self.param_min[j] - self.center[j]) / fine_width))
        imax = int(round((self.param_max[j] - self.center[j]) / fine_width))
        # In f90, finebins(imin-fine_edge:imax+fine_edge)
        finebins = np.zeros( (imax+fine_edge) - (imin-fine_edge) + 1 )

        if (self.plot_meanlikes):
            # In f90, finebinlikes(imin-fine_edge:imax+fine_edge)
            finebinlikes = np.zeros( (imax+fine_edge) - (imin-fine_edge) + 1 )

        for i in range(len(paramVec)):
            ix2 = int(round( (paramVec[i]-self.center[j]) / width ))
            if ( (ix2<=self.ix_max[j]) and (ix2>=self.ix_min[j]) ):
                binsraw[ix2 - self.ix_min[j]] += paramVec[i]
            ix2 = int(round( (paramVec[i]-self.center[j]) / fine_width ))
            finebins[ix2 - (imin-fine_edge)] += self.weights[i]
            if (self.plot_meanlikes):
                if (self.mean_loglikes):
                    finebinlikes[ix2 - (imin-fine_edge)] += self.weights[i] * self.loglikes[i]
                else:
                    finebinlikes[ix2 - (imin-fine_edge)] += self.weights[i] * np.exp(self.meanlike-self.loglikes[i])

        if (self.ix_min[j]<>self.ix_max[j]):
            # account for underweighting near edges
            if ( (not self.has_limits_bot[ix]) and (binsraw[end_edge-1]==0) and
                 (binsraw[end_edge]>np.max(binsraw)/15) ):
                self.EdgeWarning(ix)
            if ( (not self.has_limits_top[ix]) and (binsraw[self.ix_max[j]-end_edge+1-self.ix_min[j]]==0) and
                 (binsraw[self.ix_max[j]-end_edge+1-self.ix_min[j]]>np.max(binsraw)/15) ):
                self.EdgeWarning(ix)

        # In f90, Win(-winw:winw)
        Win = np.zeros( (2*winw) + 1 )
        for i in range(-winw, winw+1):
            Win[i - (-winw)] = math.exp( - math.pow(i, 2) / math.pow(fine_fac*smooth_1D, 2) / 2 )
        Win = Win / np.sum(Win)

        has_prior = self.has_limits_bot[ix] or self.has_limits_top[ix]
        if (has_prior):
            # In f90, prior_mask(imin-fine_edge:imax+fine_edge)
            prior_mask = np.ones( (imax+fine_edge) - (imin-fine_edge) + 1 )

            if (self.has_limits_bot[ix]):
                index = (self.ix_min[j]*fine_fac) - (imin-fine_edge)
                prior_mask[ index ] = 0.5
                prior_mask[ : index ] = 0
            if (self.has_limits_top[ix]):
                index = (self.ix_max[j]*fine_fac) - (imin-fine_edge)
                prior_mask[ index ] = 0.5
                prior_mask[ index+1 : ] = 0

        # High resolution density (sampled many times per smoothing scale)
        if (self.has_limits_bot[ix]): imin = self.ix_min[j] * fine_fac
        if (self.has_limits_top[ix]): imax = self.ix_max[j] * fine_fac

        logging.debug("Density1D ... ")
        self.density1D = Density1D(imax-imin+1, fine_width)
        for i in range(imin, imax+1):
            istart, iend = (i-winw)-(imin-fine_edge), (i+winw+1)-(imin-fine_edge)
            self.density1D.P[i-imin] = np.dot(Win, finebins[istart:iend])
            self.density1D.X[i-imin] = self.center[j] + (i*fine_width)
            if (has_prior and self.density1D.P[i-imin]>0):
                # correct for normalization of window where it is cut by prior boundaries
                edge_fac = 1 / np.dot(Win, prior_mask[istart:iend])
                self.density1D.P[i-imin] *= edge_fac

        maxbin = np.max(self.density1D.P)
        if (maxbin==0):
            print 'no samples in bin, param: ', self.index2name[ix]
            sys.exit()
        self.density1D.P /= maxbin

        logging.debug("InitSpline ...")
        self.density1D.InitSpline()

        logZero = 1e30
        if (not self.no_plots):
            # In f90, binCounts(ix_min(j):ix_max(j))
            bincounts = np.zeros( self.ix_max[j] - self.ix_min[j] + 1 )
            if (self.plot_meanlikes):
                # In f90, binlikes(ix_min(j):ix_max(j))
                binlikes = np.zeros( self.ix_max[j] - self.ix_min[j] + 1 )
                if (self.mean_loglikes): binlikes[:] = logZero

            # Output values for plots
            for ix2 in range(self.ix_min[j], self.ix_max[j]+1):
                istart, iend = (ix2*fine_fac-winw)-(imin-fine_edge), (ix2*fine_fac+winw+1)-(imin-fine_edge)
                bincounts[ix2 - self.ix_min[j]] = np.dot(Win, finebins[istart:iend])
                if (self.plot_meanlikes and (bincounts[ix2 - self.ix_min[j]]>0)):
                    binlikes[ix2 - self.ix_min[j]] = np.dot(Win, finebinlikes[istart:iend]) / bincounts[ix2 - self.ix_min[j]]
                if (has_prior):
                    # correct for normalization of window where it is cut by prior boundaries
                    edge_fac = 1 / np.dot(Win, prior_mask[istart:iend])
                    bincounts[ix2 - self.ix_min[j]] *= edge_fac

            bincounts = bincounts / maxbin
            if (self.plot_meanlikes and self.mean_loglikes):
                maxbin = min(binlikes)
                binlikes = np.where( (binlikes-maxbin)<30, np.exp(-(binlikes-maxbin)), 0 )



            if writeDataToFile:
                logging.debug("Write data to file ...")

                fname = self.rootname + "_p_" + str(self.index2name[j])

                filename = os.path.join(self.plot_data_dir, fname + ".dat")
                textFileHandle = open(filename, 'w')
                for i in range(self.ix_min[j], self.ix_max[j]+1):
                    textFileHandle.write("%16.7E%16.7E\n"%(self.center[j] + i*width, bincounts[i - self.ix_min[j]]))
                if (self.ix_min[j]==self.ix_max[j]):
                    for i in range(self.ix_min.shape[0]):
                        textFileHandle.write("%16.7E"%(self.center[j] + self.ix_min[i]*width))
                    textFileHandle.write("\n")
                textFileHandle.close()

                if (self.plot_meanlikes):
                    maxbin = max(binlikes)
                    filename_like = os.path.join(self.plot_data_dir, fname + ".likes")
                    textFileHandle = open(filename_like, 'w')
                    for i in range(self.ix_min[j], self.ix_max[j]+1):
                        textFileHandle.write("%16.7E%16.7E\n"%(self.center[j] + i*width, binlikes[i - self.ix_min[j]]/maxbin))
                    textFileHandle.close()

            else:

                logging.debug("Return data ...")
                dat, likes = None, None

                ncols = 2
                nrows = self.ix_max[j]+1 - self.ix_min[j]
                if (self.ix_min[j]==self.ix_max[j]): nrows += 1

                dat = np.ndarray((nrows, ncols))
                index = 0
                for i in range(self.ix_min[j], self.ix_max[j]+1):
                    dat[index] = self.center[j] + i*width, bincounts[i - self.ix_min[j]]
                    index += 1
                if (self.ix_min[j]==self.ix_max[j]):
                    dat[index] = self.center[j] + self.ix_min[0]*width, self.center[j] + self.ix_min[1]*width
                logging.debug("dat.shape = %s"%str(dat.shape))

                if (self.plot_meanlikes):
                    nrows = self.ix_max[j]+1 - self.ix_min[j]
                    likes = np.ndarray((nrows, ncols))
                    index = 0
                    for i in range(self.ix_min[j], self.ix_max[j]+1):
                        likes[index] = self.center[j] + i*width, binlikes[i - self.ix_min[j]]/maxbin
                        index += 1
                    logging.debug("likes.shape = %s"%str(likes.shape))

                return dat, likes

        return None, None



    def Get2DPlotData(self, j, j2, writeDataToFile=True):
        """
        Get 2D plot data.
        """
        fine_fac_base = 5

        has_prior = self.has_limits[j] or self.has_limits[j2]

        corr = self.corrmatrix[j][j2]
        # keep things simple unless obvious degeneracy
        if (abs(corr)<0.3): corr = 0.
        corr = max(-0.95, corr)
        corr = min(0.95, corr)

        # for tight degeneracies increase bin density
        nbin2D = min(4*self.num_bins_2D, int(round(self.num_bins_2D/(1-abs(corr)))))

        widthx = (self.range_max[j]-self.range_min[j]) / (nbin2D+1)
        widthy = (self.range_max[j2]-self.range_min[j2]) / (nbin2D+1)
        smooth_scale = (self.smooth_scale_2D*nbin2D) / self.num_bins_2D
        fine_fac = max(2, int(round(fine_fac_base/smooth_scale)))

        ixmin = int(round((self.range_min[j] - self.center[j]) / widthx))
        ixmax = int(round((self.range_max[j] - self.center[j]) / widthx))

        iymin = int(round((self.range_min[j2] - self.center[j2]) / widthy))
        iymax = int(round((self.range_max[j2] - self.center[j2]) / widthy))

        if (not self.has_limits_bot[j]): ixmin -= 1
        if (not self.has_limits_bot[j2]): iymin -= 1
        if (not self.has_limits_top[j]): ixmax += 1
        if (not self.has_limits_top[j2]): iymax += 1

        # Using index mapping for f90 2D arrays with non standard indexes.

        # In f90, bins2D(ixmin:ixmax,iymin:iymax) and bin2Dlikes(ixmin:ixmax,iymin:iymax)
        bins2D = np.zeros( (ixmax-ixmin+1, iymax-iymin+1) )
        bin2Dlikes = np.zeros( (ixmax-ixmin+1, iymax-iymin+1) )

        winw = int(round(fine_fac*smooth_scale))
        imin = (ixmin-3)*winw + 1
        imax = (ixmax+3)*winw - 1
        jmin = (iymin-3)*winw + 1
        jmax = (iymax+3)*winw - 1

        # In f90, finebins(imin:imax,jmin:jmax)
        finebins = np.zeros( (imax-imin+1, jmax-jmin+1) )
        if (self.shade_meanlikes):
            # In f90, finebinlikes(imin:imax,jmin:jmax)
            finebinlikes = np.zeros( (imax-imin+1, jmax-jmin+1) )

        widthj  = widthx / fine_fac
        widthj2 = widthy / fine_fac
        col1 = self.colix[j]
        col2 = self.colix[j2]
        for i in range(self.numrows):
            ix1 = int(round(((self.samples[i, j]-self.center[j]) / widthj)))
            ix2 = int(round(((self.samples[i, j2]-self.center[j2]) / widthj2)))
            if ( (ix1>=imin) and (ix1<=imax) and (ix2>=jmin) and (ix2<=jmax)):
                finebins[ix1-imin][ix2-jmin] += self.weights[i]
                if (self.shade_meanlikes):
                    finebinlikes[ix1-imin][ix2-jmin] += self.weights[i] * (math.exp(self.meanlike-self.loglikes[i]))

        winw = int(round(2*fine_fac*smooth_scale))
        # In f90, Win(-winw:winw,-winw:winw)
        Win = np.ones( ((2*winw)+1, (2*winw)+1) )
        indexes = range(-winw, winw+1)
        for ix1 in indexes:
            for ix2 in indexes:
                Win[ix1-(-winw)][ix2-(-winw)] = math.exp(
                    - ( ((ix1*ix1) + (ix2*ix2) - 2*corr*ix1*ix2)) /
                      (2 * (fine_fac*fine_fac) * (smooth_scale*smooth_scale) * (1-corr*corr)) )

        if (has_prior):
            norm = np.sum(Win)
            # In f90, prior_mask(imin:imax,jmin:jmax)
            prior_mask = np.ones( (imax-imin+1, jmax-jmin+1) )
            if (self.has_limits_bot[j]):
                prior_mask[(ixmin*fine_fac)-imin, :] /= 2
                prior_mask[:(ixmin*fine_fac)-imin, :] = 0
            if (self.has_limits_top[j]):
                prior_mask[(ixmax*fine_fac)-imin, :] /= 2
                prior_mask[(ixmax*fine_fac)-imin:, :] = 0

            if (self.has_limits_bot[j2]):
                prior_mask[:, (iymin*fine_fac)-jmin] /= 2
                prior_mask[:, :(iymin*fine_fac)-jmin] = 0
            if (self.has_limits_top[j2]):
                prior_mask[:, (iymax*fine_fac)-jmin] /= 2
                prior_mask[:, :(iymax*fine_fac)-jmin] = 0

        for ix1 in range(ixmin, ixmax+1):
            for ix2 in range(iymin, iymax+1):
                ix1start, ix1end = ix1*fine_fac-winw-imin, ix1*fine_fac+winw+1-imin
                ix2start, ix2end = ix2*fine_fac-winw-jmin, ix2*fine_fac+winw+1-jmin
                bins2D[ix1-ixmin][ix2-iymin] = np.sum(np.multiply(Win, finebins[ix1start:ix1end, ix2start:ix2end]))
                if (self.shade_meanlikes):
                    bin2Dlikes[ix1-ixmin][ix2-iymin] = np.sum(np.multiply(Win, finebinlikes[ix1start:ix1end, ix2start:ix2end]))

                if (has_prior):
                    # correct for normalization of window where it is cut by prior boundaries
                    denom = np.sum(np.multiply(Win, prior_mask[ix1start:ix1end, ix2start:ix2end]))
                    if denom!=0.:
                        edge_fac = norm / denom
                    else:
                        edge_fac = 0.
                    bins2D[ix1-ixmin][ix2-iymin] *= edge_fac
                    if (self.shade_meanlikes):
                        bin2Dlikes[ix1-ixmin][ix2-iymin] *= edge_fac

        if (self.shade_meanlikes):
            for ix1 in range(ixmin, ixmax+1):
                for ix2 in range(iymin, iymax+1):
                    if (bins2D[ix1-ixmin][ix2-iymin]>0):
                        bin2Dlikes[ix1-ixmin][ix2-iymin] = bin2Dlikes[ix1-ixmin][ix2-iymin] / bins2D[ix1-ixmin][ix2-iymin]

        bins2D = bins2D / np.max(bins2D)

        # Get contour containing contours(:) of the probability
        norm = np.sum(bins2D)

        contour_levels = np.zeros(self.num_contours)
        for ix1 in range(self.num_contours):
            try_t = np.max(bins2D)
            try_b = 0
            lasttry = -1
            while True:
                try_sum = np.sum(bins2D[np.where(bins2D < (try_b + try_t) / 2)])
                if (try_sum > (1-self.contours[ix1])*norm):
                    try_t = (try_b + try_t) / 2
                else:
                    try_b = (try_b + try_t) / 2
                if (try_sum == lasttry): break
                lasttry = try_sum
            contour_levels[ix1] = (try_b + try_t) / 2

        bins2D[np.where(bins2D < 1e-30)] = 0

        if writeDataToFile:

            name = self.index2name.get(j, "NOTFOUND")
            name2 = self.index2name.get(j2, "NOTFOUND")
            plotfile = self.rootname + "_2D_%s_%s"%(name, name2)
            filename = os.path.join(self.plot_data_dir, plotfile)
            textFileHandle = open(filename, 'w')
            for ix1 in range(ixmin, ixmax+1):
                for ix2 in range(iymin, iymax+1):
                    textFileHandle.write("%16.7E"%(bins2D[ix1-ixmin][ix2-iymin]))
                textFileHandle.write("\n")
            textFileHandle.close()

            textFileHandle = open(filename + "_y", 'w')
            for i in range(ixmin, ixmax+1):
                textFileHandle.write("%16.7E\n"%(self.center[j] + i*widthx))
            textFileHandle.close()

            textFileHandle = open(filename + "_x", 'w')
            for i in range(iymin, iymax+1):
                textFileHandle.write("%16.7E\n"%(self.center[j2] + i*widthy))
            textFileHandle.close()

            textFileHandle = open(filename + "_cont", 'w')
            s_levels = [ "%16.7E"%level for level in contour_levels ]
            textFileHandle.write("%s\n"%(" ".join(s_levels)))
            textFileHandle.close()

            if (self.shade_meanlikes):
                textFileHandle = open(filename + "_likes", 'w')
                maxbin = np.max(bin2Dlikes)
                for ix1 in range(ixmin, ixmax+1):
                    for ix2 in range(iymin, iymax+1):
                        textFileHandle.write("%16.7E"%(bin2Dlikes[ix1-ixmin][ix2-iymin]/maxbin))
                    textFileHandle.write("\n")
                textFileHandle.close()
        else:
            dat, likes, cont, x, y = None, None, None, None, None

            ncols = ixmax+1 - ixmin
            nrows = iymax+1 - iymin
            dat = np.ndarray((ncols, nrows))
            irow, icol = 0, 0
            for ix1 in range(ixmin, ixmax+1):
                icol = 0
                for ix2 in range(iymin, iymax+1):
                    dat[irow][icol] = bins2D[ix1-ixmin][ix2-iymin]
                    icol += 1
                irow += 1

            if (self.shade_meanlikes):
                maxbin = np.max(bin2Dlikes)
                likes = np.ndarray((ncols, nrows))
                irow, icol = 0, 0
                for ix1 in range(ixmin, ixmax+1):
                    icol = 0
                    for ix2 in range(iymin, iymax+1):
                        likes[irow][icol] = bin2Dlikes[ix1-ixmin][ix2-iymin]/maxbin
                        icol += 1
                    irow += 1

            nlevels = len(contour_levels)
            cont = np.ndarray(nlevels)
            idx = 0
            for level in contour_levels:
                cont[idx] = level
                idx += 1

            x = np.ndarray(nrows)
            idx = 0
            for i in range(iymin, iymax+1):
                x[idx] = self.center[j2] + i*widthy
                idx += 1

            y = np.ndarray(ncols)
            idx = 0
            for i in range(ixmin, ixmax+1):
                y[idx] = self.center[j] + i*widthx
                idx += 1

            return dat, likes, cont, x, y


    def EdgeWarning(self, i):
        if self.index2name.has_key(i):
            print 'Warning: sharp edge in parameter %i - check limits%i'%(i, i+1)
        else:
            name = self.index2name[i]
            print 'Warning: sharp edge in parameter % - check limits[%s] or limits%i'%(name, name, i+1)

    def GetChainLikeSummary(self, toStdOut=False):
        text = ""
        #maxlike = np.max(self.loglikes)
        maxlike = self.loglikes[0]
        text += "Best fit sample -log(Like) = %f\n"%maxlike
        if ((self.loglikes[self.numrows-1]-maxlike)<30):
            self.meanlike = np.log( np.sum(np.exp(self.loglikes-maxlike)*self.weights) / self.norm ) + maxlike
            text += "Ln(mean 1/like) = %f\n"%(self.meanlike)
        self.meanlike = np.sum(self.loglikes*self.weights) / self.norm
        text += "mean(-Ln(like)) = %f\n"%(self.meanlike)
        self.meanlike = - np.log( np.sum(np.exp(-(self.loglikes-maxlike))*self.weights) / self.norm ) + maxlike
        text += "-Ln(mean like)  = %f\n"%(self.meanlike)
        if toStdOut:
            print text
        else:
            return text

    def ComputeColix(self):
        self.num_vars = self.samples.shape[1]
        colix = [0] * self.num_vars
        for i in range(self.num_vars): colix[i] = i
        self.colix = colix


    def ComputeStats(self):
        """
        Compute mean and std dev.
        """
        nparam = self.samples.shape[1]
        self.means = np.zeros(nparam)
        self.sddev = np.zeros(nparam)
        for i in range(nparam): self.means[i] = self.mean(self.samples[:, i])
        for i in range(nparam): self.sddev[i] = self.std(self.samples[:, i])

    def Init1DDensity(self):
        nparam = self.samples.shape[1]
        self.param_min = np.zeros(nparam)
        self.param_max = np.zeros(nparam)
        self.range_min = np.zeros(nparam)
        self.range_max = np.zeros(nparam)
        self.center = np.zeros(nparam)
        self.ix_min = np.zeros(nparam, dtype=np.int)
        self.ix_max = np.zeros(nparam, dtype=np.int)
        #
        self.marge_limits_bot = np.ndarray([self.num_contours, nparam], dtype=bool)
        self.marge_limits_top = np.ndarray([self.num_contours, nparam], dtype=bool)

    def GetConfidenceRegion(self):
        self.ND_cont1, self.ND_cont2 = -1, -1
        cumsum = np.cumsum(self.weights)
        indexes = np.where(cumsum>self.norm*self.contours[0])
        self.ND_cont1 = indexes[0][0]
        indexes = np.where(cumsum>self.norm*self.contours[1])
        self.ND_cont2 = indexes[0][0]

    def ReadRanges(self):
        ranges_file = self.root + '.ranges'
        if (os.path.isfile(ranges_file)):
            self.ranges = Ranges(ranges_file)
        else:
            print "No file %s"%ranges_file

    def WriteBounds(self, filename=None):
        if not hasattr(self, 'ranges') or self.ranges is None: return

        write_to_file = filename is not None
        if write_to_file:
            textFileHandle = open(filename, 'w')
        else:
            upper = dict()
            lower = dict()
        for i in self.index2name.keys():
            if not self.isused[i]: continue
            if (self.has_limits_bot[i] or self.has_limits_top[i]):
                name = self.index2name[i]
                valMin = self.ranges.min(name)
                if (valMin is not None):
                    lim1 = "%15.7E"%valMin
                else:
                    lim1 = "    N"
                valMax = self.ranges.max(name)
                if (valMax is not None):
                    lim2 = "%15.7E"%valMax
                else:
                    lim2 = "    N"

                if write_to_file:
                    textFileHandle.write("%22s%17s%17s\n"%(name, lim1, lim2))
                else:
                    if lim1.strip() != 'N': lower[name] = float(valMin)
                    if lim2.strip() != 'N': upper[name] = float(valMax)
        if write_to_file:
            textFileHandle.close()
        else:
            return lower, upper



    def OutputMargeStats(self, writeDataToFile=True):
        contours_str = '; '.join([ str(c) for c in self.contours ])


        maxLen = max([ len(name) for name in self.index.keys() ])
        j = max(9, maxLen)

        text  = ""
        text += "Marginalized limits: %s\n"%contours_str
        text += "\n"
        text += "%-15s"%("parameter")
        text += "%-15s "%("mean")
        text += "%-15s "%("sddev")
        for j in range(self.num_contours):
            text += "%-15s "%("lower"+str(j+1))
            text += "%-15s "%("upper"+str(j+1))
            text += "%-7s"%("limit"+str(j+1))
        text += "\n"

        for j in range(self.num_vars):
            if not self.isused[j]: continue
            text += "%-12s"%(self.index2name[j])
            text += "%16.7E%16.7E"%(self.means[j], self.sddev[j])
            for i in range(self.num_contours):
                text += "%16.7E%16.7E"%(self.LowerUpperLimits[j][0][i], self.LowerUpperLimits[j][1][i])
                if (self.marge_limits_bot[i][j] and self.marge_limits_top[i][j]):
                    tag = 'none'
                elif (self.marge_limits_bot[i][j]):
                    tag = '>'
                elif (self.marge_limits_top[i][j]):
                    tag = '<'
                else:
                    tag = 'two'
                text += "  %-5s"%(tag)
            label = self.paramNames.names[j].label
            text += "   %s\n"%(label)

        if writeDataToFile:
            filename = self.rootdirname + '.margestats'
            textFileHandle = open(filename, 'w')
            textFileHandle.write(text)
            textFileHandle.close()
        else:
            return text

    def GetUsedParamNames(self):
        names = []
        for info in self.paramNames.names:
            index = self.index[info.name]
            if self.isused[index]:
                names.append(info.name)
        return names

    def WriteParamNames(self, filename, indices=None, add_derived=None):
        textFileHandle = open(filename, 'w')
        for info in self.paramNames.names:
            index = self.index[info.name]
            if self.isused[index]:
                textFileHandle.write(info.string() + '\n')
        textFileHandle.close()

    # Pass weights as parameter (for chain only)
    def thin_indices_chain(self, weights, factor):
        thin_ix = []
        tot = 0
        i = 0
        numrows = len(weights)
        mult = weights[i]
        if abs(round(mult) - mult) > 1e-4:
                raise Exception('Can only thin with integer weights')

        while i < numrows:
            if (mult + tot < factor):
                tot += mult
                i += 1
                if i < numrows: mult = weights[i]
            else:
                thin_ix.append(i)
                if mult == factor - tot:
                    i += 1
                    if i < numrows: mult = weights[i]
                else:
                    mult -= (factor - tot)
                tot = 0
        return thin_ix


    # Bins functions

    def createLowerUpperLimits(self):
        self.LowerUpperLimits = np.zeros([self.num_vars, 2, self.num_contours])


    def Do1DBins(self, max_frac_twotail=None, writeDataToFile=True):
        if max_frac_twotail is None:
            max_frac_twotail = self.max_frac_twotail
        self.createLowerUpperLimits()

        for j in range(self.num_vars):
            if not self.isused[j]: continue
            self.Get1DDensity(j, writeDataToFile)
            self.setLowerUpperLimits(j, max_frac_twotail)



    def setLowerUpperLimits(self, ix, max_frac_twotail):

        # Get limits, one or two tail depending on whether posterior
        # goes to zero at the limits or not
        for ix1 in range(self.num_contours):

            self.marge_limits_bot[ix1][ix] = self.has_limits_bot[ix] and \
            (not self.force_twotail) and (self.density1D.P[0] > max_frac_twotail[ix1])
            self.marge_limits_top[ix1][ix] = self.has_limits_top[ix] and \
            (not self.force_twotail) and (self.density1D.P[-1] > max_frac_twotail[ix1])

            if ( (not self.marge_limits_bot[ix1][ix]) or \
                 (not self.marge_limits_top[ix1][ix]) ):
                # give limit
                tail_limit_bot, tail_limit_top, marge_bot, marge_top = self.density1D.Limits(self.contours[ix1])
                self.marge_limits_bot[ix1][ix] = marge_bot
                self.marge_limits_top[ix1][ix] = marge_top

                limfrac = 1 - self.contours[ix1]

                if (self.marge_limits_bot[ix1][ix]):
                    # fix to end of prior range
                    tail_limit_bot = self.range_min[ix]
                elif (self.marge_limits_top[ix1][ix]):
                    # 1 tail limit
                    tail_limit_bot = self.confidence(self.samples[:, ix],
                                                     limfrac, upper=False)
                else:
                    # 2 tail limit
                    tail_confid_bot = self.confidence(self.samples[:, ix],
                                                      limfrac/2, upper=False)

                if (self.marge_limits_top[ix1][ix]):
                    tail_limit_top = self.range_max[ix]
                elif (self.marge_limits_bot[ix1][ix]):
                    tail_limit_top = self.confidence(self.samples[:, ix],
                                                     limfrac, upper=True)
                else:
                    tail_confid_top = self.confidence(self.samples[:, ix],
                                                      limfrac/2, upper=True)

                if ( (not self.marge_limits_bot[ix1][ix]) and \
                     (not self.marge_limits_top[ix1][ix]) ):
                    # Two tail, check if limits are at very differen density
                    if (math.fabs( self.density1D.Prob(tail_confid_top) -
                                   self.density1D.Prob(tail_confid_bot))
                        < self.credible_interval_threshold):
                        tail_limit_top = tail_confid_top
                        tail_limit_bot = tail_confid_bot

                self.LowerUpperLimits[ix][1][ix1] = tail_limit_top
                self.LowerUpperLimits[ix][0][ix1] = tail_limit_bot
            else:
                # no limit
                self.LowerUpperLimits[ix][1][ix1] = self.range_max[ix]
                self.LowerUpperLimits[ix][0][ix1] = self.range_min[ix]



    def GetCust2DPlots(self, num_cust2D_plots):

        try_t = 1e5
        x, y = 0, 0
        cust2DPlots = []
        for j in range(num_cust2D_plots):
            try_b = -1e5
            for ix1 in range(self.num_vars):
                for ix2 in range(ix1+1, self.num_vars):
                    if ( not self.isused[ix1] or \
                             not self.isused[ix2]): continue

                    if ( abs(self.corrmatrix[ix1][ix2]) < try_t ) and \
                            ( abs(self.corrmatrix[ix1][ix2]) > try_b) :
                        try_b = abs(self.corrmatrix[ix1][ix2])
                        x, y = ix1, ix2
            if (try_b==-1e5):
                num_cust2D_plots = j-1
                break
            try_t = try_b
            cust2DPlots.append(x + y*1000)

        return cust2DPlots, num_cust2D_plots


    # Write functions

    def WriteScriptPlots1D(self, filename):
        textFileHandle = open(filename, 'w')
        textInit = WritePlotFileInit()
        textFileHandle.write(textInit%(
                self.plot_data_dir, self.subplot_size_inch,
                self.out_dir, self.rootname))
        text = 'g.plots_1d(roots)\n'
        textFileHandle.write(text)
        textExport = WritePlotFileExport()
        fname = self.rootname + '.' + 'pdf'
        textFileHandle.write(textExport%(fname))
        textFileHandle.close()


    def WriteScriptPlots2D(self, filename, plot_2D_param, num_cust2D_plots, cust2DPlots, plots_only):
        self.done2D = np.ndarray([self.num_vars, self.num_vars], dtype=bool)
        self.done2D[:, :] = False

        textFileHandle = open(filename, 'w')
        textInit = WritePlotFileInit()
        textFileHandle.write(textInit%(
                self.plot_data_dir, self.subplot_size_inch2,
                self.out_dir, self.rootname))
        textFileHandle.write('pairs=[]\n')
        plot_num = 0
        for j in range(self.num_vars):
            if (self.ix_min[j]<>self.ix_max[j]):
                if ( (plot_2D_param<>0) or (num_cust2D_plots<>0) ):
                    if (j==plot_2D_param): continue
                    j2min = 0
                else:
                    j2min = j + 1

                for j2 in range(j2min, self.num_vars):
                    if (self.ix_min[j2]<>self.ix_max[j2]):
                        if ( (plot_2D_param<>0) and (j2<>plot_2D_param) ): continue
                        if ( (num_cust2D_plots<>0) and (cust2DPlots.count(j*1000+j2)==0) ): continue
                        plot_num += 1
                        self.done2D[j][j2] = True
                        if (not plots_only): self.Get2DPlotData(j, j2)
                        name1 = self.index2name[j]
                        name2 = self.index2name[j2]
                        textFileHandle.write("pairs.append(['%s','%s'])\n"%(name1, name2))
        textFileHandle.write('g.plots_2d(roots,param_pairs=pairs)\n')
        textExport = WritePlotFileExport()
        fname = self.rootname + '_2D.' + 'pdf'
        textFileHandle.write(textExport%(fname))
        textFileHandle.close()


    def WriteScriptPlotsTri(self, filename, triangle_params):
        textFileHandle = open(filename, 'w')
        textInit = WritePlotFileInit()
        textFileHandle.write(textInit%(
                self.plot_data_dir, self.subplot_size_inch,
                self.out_dir, self.rootname))
        names = [ self.index2name[i] for i in triangle_params if self.isused[i] ]
        text = 'g.triangle_plot(roots, %s)\n'%str(names)
        textFileHandle.write(text)
        textExport = WritePlotFileExport()
        fname = self.rootname + '_tri.' + 'pdf'
        textFileHandle.write(textExport%(fname))
        textFileHandle.close()


    def WriteScriptPlots3D(self, filename, num_3D_plots, plot_3D):
        textFileHandle = open(filename, 'w')
        textInit = WritePlotFileInit()
        textFileHandle.write(textInit%(
                self.plot_data_dir, self.subplot_size_inch3,
                self.out_dir, self.rootname))
        textFileHandle.write('sets=[]\n')
        text = ""
        for j in range(num_3D_plots):
            v1, v2, v3 = plot_3D[j]
            text += "sets.append(['%s','%s','%s'])\n"%(v1, v2, v3)
        text += 'g.plots_3d(roots,sets)\n'
        textFileHandle.write(text)
        fname = self.rootname + '_3D.' + 'pdf'
        textExport = WritePlotFileExport()
        textFileHandle.write(textExport%(fname))
        textFileHandle.close()


    def WriteGlobalLikelihood(self, filename):
        bestfit_ix = 0 # Since we have sorted the lines
        textFileHandle = open(filename, 'w')
        textInit = self.GetChainLikeSummary(toStdOut=False)
        textFileHandle.write(textInit)
        textFileHandle.write("\n")
        textFileHandle.write('param  bestfit        lower1         upper1         lower2         upper2\n')
        for j in range(self.num_vars):
            if not self.isused[j]: continue
            best = self.samples[bestfit_ix][j]
            min1 = min(self.samples[0:self.ND_cont1, j])
            max1 = max(self.samples[0:self.ND_cont1, j])
            min2 = min(self.samples[0:self.ND_cont2, j])
            max2 = max(self.samples[0:self.ND_cont2, j])
            label = self.paramNames.names[j].label
            textFileHandle.write('%5i%15.7E%15.7E%15.7E%15.7E%15.7E   %s\n'%(j+1, best, min1, max1, min2, max2, label))
        textFileHandle.close()


def WritePlotFileInit():
    text = """import GetDistPlots, os
g=GetDistPlots.GetDistPlotter('%s')
g.settings.setWithSubplotSize(%f)
outdir='%s'
roots=['%s']
"""
    return text

def WritePlotFileExport():
    text = "g.export(os.path.join(outdir,'%s'))\n"
    return text



# ==============================================================================

# Usefull functions

def GetParamNamesFiles(rootdir):
    pattern = os.path.join(rootdir, '*.paramnames')
    files = glob.glob(pattern)
    files.sort()
    return files

def GetRootFileName(rootdir):
    rootFileName = ""
    pattern = os.path.join(rootdir, '*_*.txt')
    chain_files = glob.glob(pattern)
    chain_files.sort()
    if chain_files:
        chain_file0 = chain_files[0]
        rindex = chain_file0.rindex('_')
        rootFileName = chain_file0[:rindex]
    return rootFileName

def GetChainFiles(rootdir):
    chain_files = []
    iFile = 1
    while 1:
        fname = rootdir + '_' + str(iFile) + '.txt'
        if os.path.isfile(fname):
            chain_files.append(fname)
            iFile += 1
        else: break
    return chain_files

