# BICEP2/Keck Array and Planck Joint Analysis January 2015 Data Products
# The BICEP2/Keck and Planck Collaborations, A Joint Analysis of BICEP2/Keck Array and Planck Data
# http://bicepkeck.org/
#
# File: BKPlanck_README.txt
# Date: 2015-01-30 
#
This README gives an overview of the files contained in the BICEP2/Keck Array and Planck Joint Analysis CosmoMC data release.
Additional information can be found in the comments at the top of each data file as well as in the Joint Analysis paper.

The files listed below are available on http://bicepkeck.org/ and in CosmoMC 
(http://cosmologist.info/cosmomc/, as of version Jan 2015).

The files are tested to produce the likelihood distributions given in the Joint Analysis paper using CosmoMC, version Jan 2015.

IMPORTANT NOTE REGARDING smooth_scale_1D used in getdist: 
The smooth_scale_1D parameter can significantly impact the likelihood of r close to zero. 
Make sure to choose a sufficiently small value of smooth_scale_1D. 
See the comments in BKPlanck_01_fiducial_dist.ini for more details.

Brief summary of contents (file locations are as distributed in CosmoMC):
1. data/BKPlanck/BKPlanck_detset*: These files contain the data (bandpowers, noise estimates, and ancillary data) needed to use the 
   BICEP2/Keck Array and Planck joint dataset in CosmoMC. BKPlanck_detset_comb_dust.dataset is the main file. 
   Further information can be found in comments therein.
2. data/BKPlanck/BKPlanck_year*: These files contain an alternate version of the BICEP2/Keck Array and Planck joint dataset with a different 
   set of data splits for Planck single-frequency spectra. The 'detset' split is the fiducial and preferred version. 
   The 'year' split may be useful because it includes lower-frequency bandpowers.
3. data/BKPlanck/windows/: This directory contains bandpower window functions.
4. data/BKPlanck/bandpass*: These files contain instrument frequency response.
5. batch2/BKPlanck/BKPlanck_01_fiducial.ini, BKPlanck_01_fiducial_dist.ini: These files run CosmoMC/getdist to recompute the likelihood results 
   of the fiducial analysis in the Joint Analysis paper.
6. batch2/BKPlanck/BKPlanck_02_y1y2.ini, BKPlanck_03...ini, etc.: These files run CosmoMC/getdist to recompute likelihood results that correspond 
   to analysis alternatives discussed in the Joint Analysis paper, Section III C.
7. batch2/BKPlanck.ini: This file sets the default data selection and foreground parameters as used in the fiducial analysis 
   in the Joint Analysis paper.
8. batch2/BKPlanckonly.ini: This file sets the default data selection, foreground parameters (through BKPlanck.ini) and the cosmological
   parameters as used in the fiducial analysis in the Joint Analysis paper.
9. source/CMB_BK_Planck.f90: The code in this file defines the multi-component CMB+foregrounds model described in Section III 
   of the Joint Analysis paper. This code is not intended to run independently, but as part of CosmoMC. 
   However, it may be a useful reference for technical details of the model.
