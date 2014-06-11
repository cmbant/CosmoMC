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
        self.samples = np.delete(self.samples, indexes, axis=0)
        
    def SortColData(self, icol=None):
        if icol is None: return
        if icol==0: 
            indexes = self.weights.argsort()
        elif icol==1:
            indexes = self.loglikes.argsort()
        self.weights = self.weights[indexes]
        self.loglikes = self.loglikes[indexes]
        self.samples = self.samples[indexes]

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
        #fixme: see chain.writeSingle
        if(cool<>1): print 'Cooled thinned output with temp: ', cool
        MaxL = np.max(self.loglikes)
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

    # MostCorrelated2D ? 

    # GetFractionIndices ?

    # ConfidVal(self, ix, limfrac, upper) => chains.confidence()

    def PCA(self, pars, param_map, normparam_num):
        """
        Perform principle component analysis. In other words, 
        get eigenvectors and eigenvalues for normalized variables
        with optional (log) mapping.
        """
        print 'Doing PCA for ', len(pars),' parameters'
        filename = self.root + ".PCA"
        textFileHandle = open(filename)
        textFileHandle.write('PCA for parameters:\n')
        
        if (normparam_num<>0):
            if (normparam_num in pars):
                normparam = pars.index(normparam_num)+1
            else:
                sys.exit('Invalid PCA normalization parameter')
        else:
            normparam = 0

        PCdata = self.samples[pars]
        PClabs = []
        doexp = False
        for i in range(len(pars)):
            if (param_map[i]=='L'):
                doexp = True
                PCdata[:, i] = np.log(PCdata[:, i])
                PClabs.append("ln("+self.index2name[i]+")")
            elif (param_map[i]=='M'):
                doexp = True
                PCdata[:, i] = np.log(-1.0*PCdata[:, i])
                PClabs.append("ln(-"+self.index2name[i]+")")
            else:
                PClabs.append(self.index2name[i])
            textFileHandle.write("%s:%s\n"%(str(pars[i]), str(PClabs[i])))

        #PCmean = np.sum(self.weights
        
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


    # GetUsedCols ?
        
    
        
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

        has_prior = has_limits[colix(j)] or has_limits[colix(j2])

        #corr = corrmatrix(colix(j)-2,colix(j2)-2)
        # keep things simple unless obvious degeneracy
        if (abs(corr)<0.3): corr=0. 
        corr = max(-0.95, corr)
        corr = min(0.95, corr)
        
        # for tight degeneracies increase bin density
        nbin2D = min(4*num_bins_2D, round(num_bins_2D/(1-abs(corr))))
        
        widthx = (range_max[j]-range_min[j])/(nbin2D+1)
        widthy = (range_max[j2]-range_min[j2])/(nbin2D+1)
        smooth_scale = (smooth_scale_2D*nbin2D)/num_bins_2D
        fine_fac = max(2, round(fine_fac_base/smooth_scale))

        ixmin = round((range_min[j] - center[j])/widthx)
        ixmax = round((range_max[j] - center[j])/widthx)

        iymin = round((range_min[j2] - center[j2])/widthy)
        iymax = round((range_max[j2] - center[j2])/widthy)

        if (not has_limits_bot(colix[j])): ixmin -= 1
        if (not has_limits_bot(colix[j2])): iymin -= 1
        if (not has_limits_top(colix[j])): ixmax += 1
        if (not has_limits_top(colix[j2])): iymax += 1
        
        #allocate(bins2D(ixmin:ixmax,iymin:iymax))
        #allocate(bin2Dlikes(ixmin:ixmax,iymin:iymax))
        #bins2D = 0
        #bin2Dlikes = 0

        winw = round(fine_fac*smooth_scale)
        imin = (ixmin-3)*winw+1
        imax = (ixmax+3)*winw-1
        jmin = (iymin-3)*winw+1
        jmax = (iymax+3)*winw-1
        #allocate(finebins(imin:imax,jmin:jmax))
        #if (shade_meanlikes) allocate(finebinlikes(imin:imax,jmin:jmax))
        #finebins = 0
        #if (shade_meanlikes) finebinlikes=0
        
        widthj = widthx/fine_fac
        widthj2 = widthy/fine_fac
        #col1 = colix(j)
        #col2 = colix(j2)
#    do i = 0, nrows-1
#        ix1=nint((coldata(col1,i)-center(j))/widthj)
#        ix2=nint((coldata(col2,i)-center(j2))/widthj2)
#        if (ix1>=imin .and. ix1<=imax .and. ix2>=jmin .and. ix2 <=jmax) then
#            finebins(ix1,ix2) = finebins(ix1,ix2) + coldata(1,i)
#            if (shade_meanlikes) finebinlikes(ix1,ix2) = finebinlikes(ix1,ix2) + coldata(1,i)*exp(meanlike - coldata(2,i))
#        end if
#    end do

        winw = round(2*fine_fac*smooth_scale)
        #allocate(Win(-winw:winw,-winw:winw))

#    do ix1=-winw,winw
#        do ix2=-winw,winw
#            !                Win(ix1,ix2) = exp(-(ix1**2+ix2**2)/real(fine_fac**2,mcp)/2)
#            Win(ix1,ix2) = exp(-(ix1**2+ix2**2 - 2*corr*ix1*ix2)/(2*fine_fac**2*smooth_scale**2*(1-corr**2)))
#        end do
#    end do

        if (has_prior):
            norm = sum(win)
            #allocate(prior_mask(imin:imax,jmin:jmax))
            #prior_mask =1
            if (has_limits_bot(colix(j))):
                prior_mask(ixmin*fine_fac,:) = prior_mask(ixmin*fine_fac,:)/2
                prior_mask(imin:ixmin*fine_fac-1,:) = 0
            if (has_limits_top(colix(j))):
                prior_mask(ixmax*fine_fac,:) = prior_mask(ixmax*fine_fac,:)/2
                prior_mask(ixmax*fine_fac+1:imax,:) = 0
            if (has_limits_bot(colix(j2))):
                prior_mask(:,iymin*fine_fac) = prior_mask(:,iymin*fine_fac)/2
                prior_mask(:,jmin:iymin*fine_fac-1) = 0
            if (has_limits_top(colix(j2))):
                prior_mask(:,iymax*fine_fac) = prior_mask(:,iymax*fine_fac)/2
                prior_mask(:,iymax*fine_fac+1:jmax) = 0


#    do ix1=ixmin, ixmax
#        do ix2=iymin,iymax
# ... 


#    if (shade_meanlikes) then
#        deallocate(finebinlikes)
#        do ix1=ixmin,ixmax
#            do ix2 =iymin,iymax
#                if (bins2D(ix1,ix2) >0) bin2Dlikes(ix1,ix2) = bin2Dlikes(ix1,ix2)/bins2D(ix1,ix2)
#            end do
#        end do
#    end if


        bins2D = bins2D/np.max(bins2D)
        # Get contour containing contours(:) of the probability
        norm = np.sum(bins2D)

        for ix1 in range(num_contours):

            try_t = np.max(bins2D)
            try_b = 0

            lasttry = -1
            while True:
                try_sum = np.sum(bin2D[np.where(bins2D < (try_b + try_t) / 2)])
                if (try_sum > (1-contours[ix1])*norm):
                    try_t = (try_b + try_t) / 2
                else:
                    try_b = (try_b + try_t) / 2
                if (try_sum == try_last): break
                try_last = try_sum
            contour_levels[ix1] = (try_b + try_t) / 2
            
        bind2D[np.where(bins2D < 1e-30)] = 0

        #fixme: use np.savetxt
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
        pass

    def GetChainLikeSummary(self, toStdOut=False):
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
        self.index2name = dict((v,k) for k, v in self.index.iteritems()) 
        self.index2name.keys().sort()
        for i in self.index2name.keys():
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
