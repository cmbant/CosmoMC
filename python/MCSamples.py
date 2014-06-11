# MCSamples.py

import os
import sys
import numpy as np
from chains import chains

class MCSamples(chains):

    def __init__(self, root=None, ignore_rows=0):
        chains.__init__(self, root, ignore_rows)
        self.ReadRanges()


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
        self.samples = np.delete(self.samples, indexes)
        
    def SortColData(self, icol=None):
        if icol is None: return
        if icol==0: 
            indexes = self.weights.argsort()
        elif icol==1:
            indexes = self.loglikes.argsort()
        self.weights = self.weights[indexes]
        self.loglikes = self.loglikes[indexes]
        #todo 
        #self.samples[] = ...


    def WriteSingleSamples(self, filename, single_thin):
        """
        Make file of weight-1 samples by choosing samples 
        with probability given by their weight.
        """
        textFileHandle = open(filename, 'w')
        maxmult = np.max(self.weights)
        for i in range(self.numrows):
            rand = np.random.random_sample()
            if (rand <= self.weights[i]/maxmult/single_thin):
                textFileHandle.write("%16.7E"%(1.0))
                textFileHandle.write("%16.7E"%(self.loglikes[i]))
                #todo write samples
        textFileHandle.close()


    def WriteThinData(self, fname, cool):
        """
        Write thin data.
        """
        if(cool<>1): print 'Cooled thinned output with temp: ', cool
        textFileHandle = open(fname, 'w')
        #Data for thin ...
        thin_rows = 0
        textFileHandle.close()
        print  'Wrote ', thin_rows, ' thinned samples'

    # ThinData(self, fac) => chains.thin_indices(factor)


    def GetCovMatrix(self):
        """
        Get covariance matrix.
        """
        # use chains.cov and chains.corr
        #textFileHandle = open(fname)
        #textFileHandle.close()
        pass


    def MostCorrelated2D(self, i1, i2, direc):
        """
        Find which parameter is most correllated with the degeneracy in ix1, ix2
        """
        if not direc in [0, -1]:
            sys.exit('Invalid 3D color parameter')
        # Used ?


    def GetFractionIndices(self, fraction_indices, n):
        """
        Get fraction indices.
        """
        # Used ?
        pass


    # ConfidVal(self, ix, limfrac, upper) => chains.confidence()


    def PCA(self, npars, n, param_map, normparam_num):
        """
        Perform principle component analysis. In other words, 
        get eigenvectors and eigenvalues for normalized variables
        with optional (log) mapping.
        """
        val = 0.

        print 'Doing PCA for ',n,' parameters'
        filename = rootdirname + ".PCA"
        textFileHandle = open(filename)
        textFileHandle.write('PCA for parameters: ')
        
        if (normparam_num<>0):
            if (normparam_num in pars):
                normparam = pars.index(normparam_num)+1
            else:
                sys.exit('Invalid PCA normalization parameter')
        else:
            normparam = 0

        # ...
        # matplotlib.mlab import PCA
        # scipy.linalg.svd

        for i in range(n):
            # ...
            textFileHandle.write("%i:%f"%(pars[i], val))
        
        textFileHandle.write('\n')
        textFileHandle.write('Correlation matrix for reduced parameters\n')
        
        for i in range(n):
            # ...
            textFileHandle.write('\n') # ...

        # ...

        textFileHandle.write('\n')
        textFileHandle.write('e-values of correlation matrix\n')
            
        # ...

        textFileHandle.write('\n')
        textFileHandle.write('e-vectors\n')

        # ...

        textFileHandle.write('\n')
        textFileHandle.write('Principle components\n')

        for i in range(n):
            # ...
            textFileHandle.write('PC%i (e-value: %f)\n'%(i, val))

            for j in range(n):
                # ...
                pass


        textFileHandle.write('          = %f +- %f\n'%(val, val))
        textFileHandle.write('\n')
                
        textFileHandle.write('Correlations of principle components\n')


        # ...

        textFileHandle.close()


    def GetUsedCols(self):
        """
        Get used columns.
        """
        #Still used ?
        pass
        
    
        
    def DoConvergeTests(self, limfrac):
        """
        Do convergence tests.
        """
        return

        # Get statistics for individual chains, and do split tests on the samples


        rootdirname = ""
        num_chains_used = 0
        nrows = 0
        chain_indices = []


        fname = rootdirname.strip() + '.converge'
        textFileHandle = open(fname, 'w')
        
        if (num_chains_used>1):
            print 'Number of chains used =  ', num_chains_used

        chain_indices[num_chains_used+1] = nrows
        for i in range(num_chains_used):
            #chain_start(i) =  chain_indices(i) + nint((chain_indices(i+1)-chain_indices(i))*cutfrac)
            #chain_samp(i) = sum(coldata(1,chain_start(i):chain_indices(i+1)-1))
            pass

        #usedsamps = sum(chain_samp(1:num_chains_used))
        #maxsamp = maxval(chain_samp(1:num_chains_used))
        num = 0
        ncols = 0

        for j in range(3, ncols+1):
            #
            pass
        

        for j in range(3, ncols+1):
            #if (isused(j)) &
            #fullvar(j)=  sum(coldata(1,0:nrows-1)*(coldata(j,0:nrows-1)-fullmean(j))**2)/numsamp
            pass

        if (num_chains_used>1):
            textFileHandle.write("\n")
            textFileHandle.write("Variance test convergence stats using remaining chains\n")
            textFileHandle.write("param var(chain mean)/mean(chain var)\n")
            textFileHandle.write("\n")

            for j in range(3, ncols+1):
                # 
                pass


        if (num_chains_used>1) and (covmat_dimension>0):
            # Assess convergence in the var(mean)/mean(var) in the worst eigenvalue
            # c.f. Brooks and Gelman 1997
        
            while(usedvars[num]>(covmat_dimension+2)):
                num -= 1
            
            #allocate(meanscov(num,num))
            #allocate(cov(num,num))

            for jj in range(0, num):
                j = usedvars[jj]
                for kk in range(jj, num):
                    k = usedvars[kk]
                    # ...

            #meanscov = meanscov/(num_chains_used-1) !(usedsamps/maxsamp -1)
            #cov = cov / usedsamps

            #invertible = GelmanRubinEvalues(cov, meanscov, evals, num)
            if (invertible):
                textFileHandle.write("\n")
                textFileHandle.write("var(mean)/mean(var) for eigenvalues of covariance of means of orthonormalized parameters\n")
                R = 0
                # ...
                for jj in range(0, num):
                    #write (F%unit,'(1I3,f13.5)') jj,evals(jj)
                    #R = max(R,evals(jj))
                    textFileHandle.write("\n")
                textFileHandle.write(" var(mean)/mean(var), remaining chains, worst e-value: R-1 = %13.5F\n"%R)
                #deallocate(cov,meanscov)
            else:
                print 'WARNING: Gelman-Rubin covariance not invertible'


        # Do tests for robustness under using splits of the samples
        # Return the rms ([change in upper/lower quantile]/[standard deviation])
        # when data split into 2, 3,.. sets
        textFileHandle.write("\n")
        textFileHandle.write("Split tests: rms_n([delta(upper/lower quantile)]/sd) n={2,3,4}:\n")
        textFileHandle.write("i.e. mean sample splitting change in the quantiles in units of the st. dev.\n")
        textFileHandle.write("\n")
        for j in range(3, self.ncols+1):
            # ...
            pass


        # Now do Raftery and Lewis method
        # See http://www.stat.washington.edu/tech.reports/raftery-lewis2.ps
        # Raw non-importance sampled chains only

        if (1):

            nburn = 0
            hardest=-1
            hardestend=0

            for ix in range(num_chains_used):
                #thin_fac(ix) = nint(maxval(coldata(1,chain_indices(ix):chain_indices(ix+1)-1)))

                for j in range(2, covmat_dimension+2):
                    #...
                    pass

                # Get thin factor to have independent samples rather than Markov
                #hardest = max(hardest,1)
                #u = ConfidVal(hardest,(1-limfrac)/2,hardestend==0)
                #thin_fac(ix) = thin_fac(ix) + 1

                while(True):
                    
                    mc.ThinData()
                    if (thin_rows<2): break
                    #binchain 

                    tran2 = 0
                    # Estimate transitions probabilities for 2nd order process
                    
                    for i in range(thin_rows-1):
                        #tran2(binchain(i-1),binchain(i)) = tran2(binchain(i-1),binchain(i)) +1 
                        pass
                    
                
                    # Test whether independence is better than Markov using BIC statistic
                    g2 = 0
                    # ...
                    g2 = g2 * 2


                    #if (g2 - log( dble(thin_rows-1) ) < 0) exit
                    thin_fac[ix] = thin_fac[ix] + 1

                if (thin_rows<2): 
                    thin_fac[ix] = 0

            textFileHandle.write("\n")
            textFileHandle.write("Raftery&Lewis statistics\n")
            textFileHandle.write("\n")
            textFileHandle.write("chain  markov_thin  indep_thin    nburn\n")
            for ix in range(num_chains_used):
                if (1): # thin_fac(ix)==0
                    textFileHandle.write("%4i      Not enough samples\n"%chain_numbers[ix])
                else:
                    textFileHandle.write("%4i%12i%12i%12i"%(
                            chain_numbers[ix], markov_thin[ix], thin_fac[ix], nburn[ix]))

            if (1):
                print 'RL: Not enough samples to estimate convergence stats'
            else:
                print 'RL: Thin for Markov: '
                print 'RL: Thin for indep samples:  '
                print 'RL: Estimated burn in steps: '

            # Get correlation lengths
            textFileHandle.write("\n")
            textFileHandle.write("Parameter auto-correlations as function of step separation\n")
            textFileHandle.write("\n")



            mc.ThinData(autocorr_thin)
            #maxoff = min(corr_length_steps,thin_rows/(autocorr_thin*num_chains_used))

            #corrs = 0
            for j in range(maxoff):
                pass

            if (maxoff>0):
                pass



        textFileHandle.close()




    
    def Get1DDensity(self, paramVec):
        """
        Get 1D density.
        """

        # input parameters
        num_bins = 1
        smooth_scale_1D = 1
        filename = "" # plot_data_dir, rootname, j, .dat
        filename_like = "" # plot_data_dir, rootname, j, .likes
        #

        fine_fac = 10
        logZero = math.pow(1, 30) # ?

        param_min = np.min(paramVec)
        param_max = np.max(paramVec)
        range_min = self.confidence(paramVec, 0.0005, upper=False) 
        range_max = self.confidence(paramVec, 0.0005, upper=True) 
        width = (range_max-range_min)/(num_bins+1)
        if (width==0):
            print "Warning width is 0"
            return

        if (smooth_scale_1D<=0):
            # Automatically set smoothing scale from rule of thumb for Gaussian
            pass
        else:
            pass


        end_edge = round(smooth_1D*2)
        
        #if (has_limits_bot(ix)) then
        #...

        #if (has_limits_top(ix)) then
        #...

        if (has_limits_top(ix)):
            center = range_max
        else:
            center = range_min

        ix_min = round((range_min-center)/width)
        ix_max = round((range_max-center)/width)

        if (not has_limits_bot(ix)): ix_min -= end_edge
        if (not has_limits_top(ix)): ix_max -= end_edge

        binsraw = np.zeros(ix_max+1-ix_min)
        
        winw = round(2.5*fine_fac*smooth_1D)
        fine_edge = winw + fine_fac*end_edge
        fine_width = width/fine_fac
        
        imin = round((param_min-center)/fine_width)
        imax = round((param_max-center)/fine_width)

        #allocate(finebins(imin-fine_edge:imax+fine_edge))
        #finebins=0
        #if (plot_meanlikes) allocate(finebinlikes(imin-fine_edge:imax+fine_edge))
        #if (plot_meanlikes) finebinlikes=0
        
        for i in range(len(paramVec)):
            ix2 = round((paramVec[i]-center)/width)
            if (ix2<=ix_max and ix2>=ix_min): 
                binsraw[ix2-1] -= paramVec[i]
            ix2 = round((paramVec[i]-center)/fine_width)
            finebins[ix2-1] += paramVec[i]
            if (plot_meanlikes):
                finebins[ix2-1] += self.weights[i]*self.loglikes[i]
            else:
                 finebins[ix2-1] += self.weights[i]*np.exp(meanlike-self.loglikes[i])
               
        if (ix_min<>ix_max):
            # account for underweighting near edges
            if (not has_limits_bot(ix) and binsraw[ix_min+end_edge-1]==0 and  
                binsraw[ix_min+end_edge]>np.max(binsraw)/15):
                # call EdgeWarning(ix-2)
                pass
            if (not has_limits_top(ix) and binsraw[ix_max+end_edge+1]==0 and  
                binsraw[ix_max-end_edge]>np.max(binsraw)/15):
                # call EdgeWarning(ix-2)
                pass
            
        #allocate(Win(-winw:winw))
        # ...

        has_prior = has_limits_bot(ix) or has_limits_top(ix)
        if (has_prior):
            # ...
            pass

        

        # High resolution density (sampled many times per smoothing scale)
        if (has_limits_bot(ix)): imin = ix_min*fine_fac
        if (has_limits_top(ix)): imax = ix_max*fine_fac

        #call Density1D%Init(imax-imin+1,fine_width)
        for i in range(imin, imax+1):
            # ...
            pass



        #maxbin = maxval(Density1D%P)
        if (maxbin==0):
            print 'no samples in bin, param: '
            sys.exit()

        #Density1D%P=Density1D%P/maxbin
        #call Density1D%InitSpline()

        if (not no_plots):
            binCounts = np.zeros(ix_max+1-ix_min)
            if (plot_meanlikes):
                binlikes = np.zeros(ix_max+1-ix_min)
                if (mean_loglikes): pass # binlikes=logZero

            # Output values for plots
            for ix2 in range(ix_min, ix_max+1):
                # ...
                pass
            

            bincounts = bincounts/maxbin
            if (plot_meanlikes and mean_loglikes):
                maxbin = np.min(binlikes)
                binlikes = np.where(binlikes-maxbin<30, np.exp(-(binlikes-maxbin)), 0)

            textFileHandle = open(filename, 'w')
            for i in range(ix_min, ix_max+1):
                textFileHandle.write("%f%16.7E\n"%(center+i*width, bincounts[i]))
                if (ix_min==ix_max): 
                    textFileHandle.write("%16.7E\n"%(center+ix_min*width))
            textFileHandle.close()
        
        if (plot_meanlikes):
            maxbin = max(binlikes)
            textFileHandle = open(filename_like, 'w')
            for i in range(ix_min, ix_max+1):
                textFileHandle.write("%f%16.7E\n"%(center+i*width, binlikes[i]/maxbin))
            textFileHandle.close()

        

    def Get2DPlotData(self, j, j2):
        """
        Get 2D plot data.
        """

        # ...


        plotfile = self.dat_file_2D( self.rootname, j, j2)
        filename = os.path.join(self.plot_data_dir, plotfile)
        textFileHandle = open(filename, 'w')
        # ...
        textFileHandle.close()


        textFileHandle = open(filename + "_y", 'w')
        # ...
        textFileHandle.close()


        textFileHandle = open(filename + "_x", 'w')
        # ...
        textFileHandle.close()


        textFileHandle = open(filename + "_cont", 'w')
        # ...
        textFileHandle.close()


        textFileHandle = open(filename + "_likes", 'w')
        # ...
        textFileHandle.close()



    def EdgeWarning(self, param):
        # FIXME not many calls. can  be replaced easily 
        pass


    # 
    def dinvnorm(self, p):
        """
        Normal inverse translate from
        http://home.online.no/~pjacklam/notes/invnorm
        a routine written by john herrero.

        :param p: value for 
        :type p: real type
        """
        a1=-39.6968302866538
        a2=220.946098424521
        a3=-275.928510446969
        a4=138.357751867269
        a5=-30.6647980661472
        a6=2.50662827745924
        b1=-54.4760987982241
        b2=161.585836858041
        b3=-155.698979859887
        b4=66.8013118877197
        b5=-13.2806815528857
        c1=-0.00778489400243029
        c2=-0.322396458041136
        c3=-2.40075827716184
        c4=-2.54973253934373
        c5=4.37466414146497
        c6=2.93816398269878
        d1=0.00778469570904146
        d2=0.32246712907004
        d3=2.445134137143
        d4=3.75440866190742
        p_low=0.02425
        p_high=1-p_low
        dinvnorm = 0
        if(p<p_low):
            q=dsqrt(-2*dlog(p))
            dinvnorm=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
        elif (p>=p_high):
            q=p-0.5
            r=q*q
            dinvnorm=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
        else:
            q=dsqrt(-2*dlog(1-p))
            dinvnorm=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/ ((((d1*q+d2)*q+d3)*q+d4)*q+1)
        return dinvnorm


    def GetChainLikeSummary(self, toStdOut=False):
        """
        """
        text = ""

        maxlike = np.max(self.loglikes)
        text += "Best fit sample -log(Like) = %f\n"%maxlike
        
        if ((self.loglikes[self.numrows-1]-maxlike)<30):
            meanlike = np.log(np.sum(np.exp(self.loglikes-maxlike)*self.weights)/self.norm)+maxlike
            text += "Ln(mean 1/like) = %f\n"%(meanlike)

        meanlike = np.sum(self.loglikes*self.weights)/self.norm
        text += "mean(-Ln(like)) = %f\n"%(meanlike)

        meanlike = -np.log(np.sum(np.exp(self.loglikes-maxlike)*self.weights)/self.norm)+maxlike
        text = "-Ln(mean like)  = %f\n"%(meanlike)

        if toStdOut:
            print text
        else:
            return text


    def SetContours(self, contours=[]):
        self.contours = contours

    def GetConfidenceRegion(self):
        ND_cont1, ND_cont2 = -1, -1
        #indexes = np.where(paramVec[:]>self.get_norm(paramVec)*self.contours[0])
        #ND_cont1 = indexes[0][0]
        #indexes = np.where(paramVec[:]>self.get_norm(paramVec)*self.contours[1])
        #ND_cont2 = indexes[0][0]
        return ND_cont1, ND_cont2


    def ReadRanges(self):
        ranges_file = self.root + '.ranges'
        if (os.path.isfile(ranges_file)):
            self.ranges = Ranges(ranges_file)
        else:
            self.ranges = None


    def WriteBounds(self, filename):
        if not hasattr(self, 'ranges') or self.ranges is None: return 

        textFileHandle = open(filename, 'w')
        names = dict((v,k) for k, v in self.index.iteritems()) #{index:name}
        names.keys().sort()
        indexes = names.keys()
        for i in indexes:
            name = names[i]
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
            textFileHandle.write("%22s%17s%17s\n"%(name, lim1, lim2))
        textFileHandle.close()


    def OutputMargeStats(self):
        pass
            

    def WriteParamNames(self, filename):
        pass



#FIXME?: pass values as parameters 
def WritePlotFileInit():
    text = """
import GetDistPlots, os
g=GetDistPlots.GetDistPlotter('%s')
g.settings.setWithSubplotSize(%f)
outdir='%s'
roots=['%s']
"""
    return text


def WritePlotFileExport():
    text = "g.export(os.path.join(outdir,'%s'))"
    return text


# 

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

