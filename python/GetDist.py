# GetDist.py


import sys


from MCSamples import MCSamples

import iniFile


# inherits from MCSamples ? 
class GetDist():


    def __init__(self):
        pass


    def run(self, args):

    

        # Declarations of variables
        # -------------------------

        # Needed ?
        # Use file constants.py for constant definition ?


    
        # Initialization
        # --------------

        if (len(args)<2):
            sys.exit('No parameter input file')

        InputFile = args[1]
        ini = iniFile.iniFile()
        ini.readFile(InputFile)
        dataset = ini.params

        # ...


        # TODO: when reading an index variable, use 0 as default value
        #first_chain = Ini%Read_Int('first_chain',1)




        # Chains 
        # ------

        self.nrows = 0
        self.num_chains_used = 0

        for chain_ix in range(first_chain, first_chain + max(1, chain_num)): # ?
            
            #if (any(chain_exclude(1:num_exclude)==chain_ix)) cycle
            
            self.num_chains_used += 1

            if (self.num_chains_used > self.max_chains):
                sys.exit('Increase max_chains in GetDist')

            #chain_indices(num_chains_used) = nrows
            #chain_numbers(num_chains_used) = chain_ix


            if (single_column_chain_files):
                # Use used for WMAP 5-year chains suppled on LAMBDA; code from Mike Nolta
                # Standard CosmoMC case below
                
                first_haschain = 0
                for ip in range(self.ncols):
                    #infile = concat(File%CheckTrailingSlash(concat(in_root,chain_ix)), pname(ip))
                    if not os.path.exists(infile):
                        print 'skipping missing ', infile 
                        #coldata(ip,:) = 0
                        #nrows2(ip) = -1
                    else:
                        print 'reading ', infile 
                        #call ChainFile%Open(infile)
                        if (first_haschain==0): first_haschain = ip
                        # ...
                        idx = 0
                        # ...
                        # TODO: read data from file


                if (first_haschain==0):
                    sys.exit('no chain parameter files read!')

                for ip in range(2, self.ncols+1): # ?
                    # ...
                    pass

                #nrows = nrows2(first_haschain)
                print 'all columns match, nrows = ', self.nrows

            else:
                # Not single column chain files (usual cosmomc format)
                # This increments nrows by number read in
                
                # ...
                pass
                

            if (map_params):
                # ...
                pass
                

            if (ignorerows<1) and (ignorerows!=0):
                # ... 
                pass


        if (self.nrows==0):
            sys.exit('No un-ignored rows! (check number of chains/burn in)')

        if (cool!=1):
            self.CoolChain(cool)

        # Adjust weights if requested
        if (adjust_priors):
            self.AdjustPriors()  

        
        # See which parameters are fixed
        self.GetUsedCols()


        # ...


        if (not no_tests):
            self.DoConvergeTests(converge_test_limit)

        if (adjust_priors):
            self.DeleteZeros()
        
        print 'mean input multiplicity = ', mean_mult

        
        # Output thinned data if requested
        # Must do this with unsorted output
        if (thin_factor!=0):
            self.ThinData(thin_factor)
            self.WriteThinData(self.rootdirname+'_thin.txt', thin_cool) # ?

    
        # Produce file of weight-1 samples if requested

        # ...


        # Only use variables whose labels are not empty (and in list of plotparams if plotparams_num /= 0)
        num_vars = 0

        if (plotparams_num!=0):

            # ...
            pass

        else:
            
            # ...
            pass
           

        for j in range(num_vars):
            #mean(j) = sum(coldata(1,0:nrows-1)*coldata( colix(j),0:nrows-1))/numsamp
            #sddev(j)  = sqrt(sum(coldata(1,0:nrows-1)*(coldata(colix(j),0:nrows-1) -mean(j))**2)/numsamp)
            pass


        if (make_single_samples):
            self.MakeSingleSamples(single_thin)

            
        # IO_WriteBounds

        # Sort data in order of likelihood of points
        #call SortColData(2)
   
        #numsamp = sum(coldata(1,0:nrows-1))

    
        # Get ND confidence region (index into sorted coldata)
        counts = 0
        ND_cont1, ND_cont2 = -1, -1

        for j in range(0, self.nrows-1):
            # ...
            pass


        triangle_plot = triangle_plot and (num_vars>1)
        if (triangle_plot):
            # ...
            pass

        print 'using ',nrows,' rows, processing ',num_vars,' parameters'
        if (indep_thin!=0):
            print 'Approx indep samples: ', int(numsamp/indep_thin) # equiv. to nint ?
        else:
            print  'effective number of samples (assuming indep): ', int(numsamp/max_mult)  # equiv. to nint ?

        # Get covariance matrix and correlation matrix
        self.GetCovMatrix()

        if (PCA_num>0) and not plots_only:
            self.PCA(PCA_params,PCA_num,PCA_func, PCA_NormParam)


        # Find best fit, and mean likelihood
        self.GetChainLikeSummary(stdout)
   
        if (not no_plots):
            # Matlab only ? 
            pass


        LowerUpperLimits = 0


        # 1D 
        # --
        
        

        






        # 2D 
        # --



        # 3D 
        # --



        # Statistics
        # ----------
    



        import pdb
        pdb.set_trace()


        print "done"


if __name__ == "__main__":   
 
    # FIXME: Use argparse here

    gd = GetDist()
    gd.run(sys.argv)
