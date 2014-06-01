# MCSamples.py

import os
import sys
import numpy

class MCSamples:
    """
    
    """

    def __init__(self):
        
        # Attributes from MCSamples module
        self.nrows = 0
        self.ncols = 0
        self.thin_rows = 0
        self.num_chains_used = 0
        self.covmat_dimension = 0

        self.plot_ext = "py"

        self.out_dir = ""


        self.thin_ix = [] # TODO numpy int array
        self.colix = [] # TODO numpy int array
        self.isued = [] # TODO numpy bool array
        self.chain_indices = [] # TODO numpy int array
        self.usedvars = [] # TODO numpy int array

        
        # New attributes needed (temporary ?)
        self.rootname = ""
        self.plot_data_dir = ""
        self.rootdirname = ""


    def AdjustPriors(self):
        """
        """
        sys.exit('You need to write the AdjustPriors function in MCSamples.py first!')
    

    def MapParameters(self, invars):
        """
        """
        sys.exit('Need to write MapParameters routine first')


    def CoolChain(self, cool):
        """
        Cool chains.

        :param cool: value for 
        :type cool: real type
        """
        
        print 'Cooling chains by ', cool
        #MaxL = minval(coldata(2,0:nrows-1))
        for i in range(self.nrows):
            #newL = coldata(2,i)*cool
            #coldata(1,i) = coldata(1,i)*exp(-(newL - coldata(2,i)) - MaxL*(1-cool) )
            #coldata(2,i) = newL
            pass
        
    
    def DeleteZeros(self):
        """
        Delete zeros. 
        """
        ii = 0
        for i in range(self.nrows):
            #if (coldata(1,i)/=0) then
            #coldata(:,ii) = coldata(:,i)
            #ii=ii+1
            pass        
        if (ii==0):
            sys.exit('Prior has removed all models!')
        if (ii!=self.nrows):
            print 'Prior removed ', nrows-ii, ' models'
        self.nrows = ii

        
    def SortColData(self, bycol):
        """
        Sort column data.

        :param bycol: value for 
        :type bycol: integer value.
        """
        pass


    def MakeSingleSamples(self, single_thin):
        """
        Make file of weight-1 samples by choosing samples 
        with probability given by their weight.

        :param single_thin: value for 
        :type single_thin: integer type
        """
        #call initRandom()
        fname = os.path.join(self.plot_data_dir, self.rootname.strip()+'_single.txt') 
        textFileHandle = open(fname, 'w')
        #maxmult = maxval(coldata(1,0:nrows-1))
        for i in range(self.nrows):
            #if (ranmar() <= coldata(1,i)/maxmult/single_thin) write (F%unit,float_format) 1.0, coldata(2,i), coldata(colix(1:num_vars),i)
            pass
        textFileHandle.close()


    def WriteThinData(self, fname, cool):
        """
        Write thin data.

        :param fname: file name
        :type fname: string value
        
        :param cool: value for 
        :type cool: real
        """
        if(cool!=1):
            print 'Cooled thinned output with temp: ', cool

        #MaxL = minval(coldata(2,0:nrows-1))

        textFileHandle = open(fname)
        for i in range(self.thin_rows):
            if(cool!=1):
                #newL = coldata(2,thin_ix(i))*cool
                #write (F%unit,float_format) exp(-(newL - coldata(2,thin_ix(i))) - MaxL*(1-cool) ), newL, coldata(3:ncols,thin_ix(i))
                pass
            else:
                #write (F%unit,float_format) 1., coldata(2:ncols,thin_ix(i))
                pass
        textFileHandle.close()
        print  'Wrote ', self.thin_rows, ' thinned samples'


    def ThinData(self, fac, ix1=None, ix2=None):
        """
        Make thinned samples.

        :param fac: value for 
        :type fac: integer type

        :param ix1: optional value for initial value
        :type ix1: integer type

        :param ix2: optional value for ending value
        :type ix2: integer type

        """
        #if (allocated(thin_ix)) deallocate(thin_ix)
        #allocate(thin_ix(0:nint(numsamp)/fac))
        
        tot = 0
        nout = 0
        i = 0
        if (ix1 is not None): i = ix1
        nend = self.nrows
        if (ix2 is not None): nend = ix2
        
        #mult = coldata(1,i)
        while (i<nend):
            #if (abs(nint(coldata(1,i)) - coldata(1,i)) > 1e-4) &
            #stop 'non-integer weights in ThinData'
            
            if (mult+tot<fac):
                tot = tot + mult
                i = i+1
                #if (i<nend): mult = nint(coldata(1,i))
            else:
                #thin_ix(nout) = i
                nout = nout+1
                if (mult==(fac-tot)):
                    i = i+1
                    #if (i<nend) mult = nint(coldata(1,i))
                else:
                    mult = mult - (fac -tot)
                tot = 0

        self.thin_rows = nout
        
        #close(50) ??? 


    def GetCovMatrix(self):
        """
        Get covariance matrix.
        """
        pass


    def MostCorrelated2D(self, i1, i2, direc):
        """
        Most correlated ...

        :param i1: value for 
        :type i1: integer type

        :param i2: value for 
        :type i2: integer type 

        :param direc: value for 
        :type direc: integer type 
        """

        if (direc!=0) and (direc!=-1):
            sys.exit('Invalid 3D color parameter')

        # Use numpy here ?
        pars = [i1, i2]
        
        # Matrix operations ...
        
        #return 


    def GetFractionIndices(self, fraction_indices, n):
        """
        Get fraction indices.

        :param fraction_indices: fraction indices 
        :type fraction_indices: array of integer type

        :param n: value
        :type n: integer type
        """

        tot = 0
        aim = numsamp/n
        num = 1
        fraction_indices[0] = 0
        fraction_indices[n] = self.nrows
        for i in range(self.nrows):
            #tot = tot + coldata(1,i)
            if (tot>aim):
                num = num+1
                fraction_indices[num-1] = i
                if (num==n): break
                aim = aim + (numsamp/n)


    def ConfidVal(self, ix, limfrac, upper, ix1=None, ix2=None):
        """
        
        :param ix: value for 
        :type ix: integer type

        :param limfrac: value for 
        :type limfrac:  type
        
        :param upper: value for 
        :type upper: integer type

        :param ix1: optional value for 
        :type ix1: integer type

        :param ix2: optional value for 
        :type ix2: integer type
        """

        l = 0
        # ...


    def PCA(self, npars, n, param_map, normparam_num):
        """
        Perform principle component analysis. In other words, 
        get eigenvectors and eigenvalues for normalized variables
        with optional (log) mapping.

        ...
        """
        pass 


    def GetUsedCols(self):
        """
        Get used columns.
        """
        #isused(j) = any(coldata(j,0:nrows-1)/=coldata(j,0))
        pass
    
        
    def DoConvergeTests(self, limfrac):
        """
        Do convergence tests.

        :param limfrac: limit to use for split tests and Raftery-Lewis 
        :type : real type
        """

        # Get statistics for individual chains, and do split tests on the samples
        
        fname = self.rootdirname.strip() + '.converge'
        textFileHandle = open(fname, 'w')
        
        if (self.num_chains_used>1):
            print 'Number of chains used =  ', self.num_chains_used

        self.chain_indices[self.num_chains_used+1] = self.nrows
        for i in range(self.self.num_chains_used):
            #chain_start(i) =  chain_indices(i) + nint((chain_indices(i+1)-chain_indices(i))*cutfrac)
            #chain_samp(i) = sum(coldata(1,chain_start(i):chain_indices(i+1)-1))
            pass

        #usedsamps = sum(chain_samp(1:num_chains_used))
        #maxsamp = maxval(chain_samp(1:num_chains_used))
        num = 0


        for j in range(3, self.ncols+1):
            #
            pass
        

        for j in range(3, self.ncols+1):
            #if (isused(j)) &
            #fullvar(j)=  sum(coldata(1,0:nrows-1)*(coldata(j,0:nrows-1)-fullmean(j))**2)/numsamp
            pass

        if (self.num_chains_used>1):
            textFileHandle.write("\n")
            textFileHandle.write("Variance test convergence stats using remaining chains\n")
            textFileHandle.write("param var(chain mean)/mean(chain var)\n")
            textFileHandle.write("\n")

            for j in range(3, self.ncols+1):
                # 
                pass


        if (self.num_chains_used>1) and (self.covmat_dimension>0):
            # Assess convergence in the var(mean)/mean(var) in the worst eigenvalue
            # c.f. Brooks and Gelman 1997
        
            while(self.usedvars[num]>(covmat_dimension+2)):
                num -= 1
            
            #allocate(meanscov(num,num))
            #allocate(cov(num,num))

            for jj in range(0, num):
                j = self.usedvars[jj]
                for kk in range(jj, num):
                    k = self.usedvars[kk]
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

        # ...


        textFileHandle.close()

    
    def Get1DDensity(self, j):
        """
        Get 1D density.

        :param j: value for 
        :type j: integer type
        """

        # ...
        pass
        

    def Get2DPlotData(self, j, j2):
        """
        Get 2D plot data.

        :param j: value for 
        :type j: integer type

        :param j2: value for 
        :type j2: integer type
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

    # FIXME: file unit / file handler as parameter
    def WritePlotFileInit(self, sm, subplot_size):
        """
        Write plot file.

        :param sm: value for 
        :type sm: bool type

        :param subplot_size: size of subplot
        :type subplot_size: integer type
        """
        
        if (self.plot_ext=="py"):
            textFileHandle.write("\n")
            textFileHandle.write("import GetDistPlots, os\n")
            textFileHandle.write("g=GetDistPlots.GetDistPlotter('%s')\n"%self.plot_data_dir)
            textFileHandle.write("g.settings.setWithSubplotSize(%f)\n"%subplot_size)
            textFileHandle.write("outdir='%s'\n"%self.out_dir)
            #write(unit,'(a)', advance='NO') 'roots=['''//trim(rootname)//''''
            #do i = 1, ComparePlots%Count
            #   write(unit,'(a)', advance='NO') ','''// ComparePlots%Item(i)//''''
            #end do
            #write(unit,'(a)') ']'
            textFileHandle.write("\n")
        else:
            # matlab ???
            pass


    # FIXME: file unit / file handler as parameter
    def WritePlotFileExport(self, tag, plot_col, plot_row):
        """
        Write plot file.

        :param tag: tag
        :type tag: string type

        :param plot_col: value for column
        :type plot_col: integer type

        :param plot_row: value for row
        :type plot_row: integer type
        """
        
        # matlab only 
        pass


    def quoted_param_name(self, j):
        """
        Get param name.

        :param j: index
        :type j: integer type
        """
        #res=''''//trim(NameMapping%NameOrNumber(j))//''''
        return ""
        # FIXME: usefull ?

    def python_param_array(self, params,num):
        # FIXME: can be replaced easily 
        pass


    def EdgeWarning(self, param):
        # FIXME not many calls. can  be replaced easily 
        pass


    # FIXME: put in another file 
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


    # FIXME: file unit / file handler as parameter
    def GetChainLikeSummary(self):
        """
        """
        # ...
        pass


if __name__ == "__main__":
    main()

