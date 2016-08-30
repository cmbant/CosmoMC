BICEP2/Keck Array October 2015 Data Products
BICEP2/Keck Array VI: Improved Constraints On Cosmology and Foregrounds When Adding 95 GHz Data From Keck Array
http://bicepkeck.org/

File: BK14_README.txt
Date: 2015-10-29

This README file gives an overview of the BK14 data release products for CosmoMC.
Additional information can be found in the header comments of each file.

Contents of this tarball (file locations as normally organized in CosmoMC):
1. data/BK14/BK14*: These files contain the data (bandpowers, bandpower covariance, and ancillary data) needed to 
   use the BK14 dataset in CosmoMC (including WMAP and Planck polarization data in the BICEP field). 
   BK14_dust.dataset is the main file. Furthur information can be found in comments therein.
2. data/BK14/windows/: This directory contains bandpower window functions.
3. data/BK14/bandpass*: These files contain instrument frequency response.
4. batch2/BK14.ini: This file sets the baseline data selection and foreground nuisance parameters used in BK-VI.
5. batch2/BK14only.ini: For CosmoMC runs where you are using *only* the BK14 data set, you should include it via
   this file, which sets scalar cosmological parameters to nominal values. These parameters are otherwise not
   well constrained by BK14 data. If you are running chains using BK14 alongside CMB TT data or similar, then it 
   is not necessary to fix these parameters.
6. batch2/BK14/BK14_01_baseline.ini, BK14_01_baseline_dist.ini: These files run CosmoMC and getdist to recompute
   the results of the BK14 baseline analysis, as seen in Figure 4 of BK-VI.
7. batch2/BK14/BK14_02_nobetadprior.ini, BK14_03_nobetasprior.ini, ...: These files run CosmoMC and getdist to
   recompute results from the likelihood variations plot, Figure 17 of BK-VI.
8. planck_covmats/BK14.covmat: This file is included by the various ini files found in batch2/BK14/. Using it 
   should speed up the convergence of your chains.
9. source/CMB_BK_Planck.f90: This file is a modified version of the source code for the foreground model used in 
   BK-VI. It is intended to be compiled and run as part of CosmoMC, but may also be a useful reference for 
   technical details of the model. It was originally distributed with the BICEP2/Keck/Planck joint result and 
   data release. With our new publication and data release, this file has been modified to change synchrotron 
   pivot frequency from 150 GHz to 23 GHz. 

