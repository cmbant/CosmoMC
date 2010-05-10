#You may need to edit the library paths for MKL for Intel
#Beware of using optmizations that lose accuracy - may give errors when running

#Edit for CFITSIO directories
cfitsio = /usr/local/cfitsio/intel10/64/3.040/lib
#/usr/local/Cluster-Users/cpac/cmb/2.1.0/cfitsio

WMAP = /home/aml1005/WMAP5/likelihood_v3

IFLAG = -I
INCLUDE=

#Intel MPI (assuming mpif77 set to point to ifort)
#change -lmkl_ia32 to -lmkl_p3 for MKL versions earlier than 6 (6 needed for ifc 8+) 
F90C     = mpif90
FFLAGS = -O2 -ip -W0 -WB -openmp -fpp -DMPI -vec_report0
LAPACKL = -L/usr/local/intel/mkl/9.1/lib/em64t -lmkl_lapack -lmkl_em64t -lguide -lpthread

#HPCF settings. Use Inteal 9 or 10.1, not 10.0
#F90C     = mpif90
#FFLAGS = -O2 -Vaxlib -W0 -WB -openmp -fpp -DMPI -vec_report0
#LAPACKL = -L/usr/local/Cluster-Apps/intel/mkl/9.0.018/lib/em64t -lmkl_lapack -lmkl -lguide -lpthread
#INCLUDE=

# Use with Intel 9 on the Altix (eg. COSMOS) [this is old]
# (do "module load icomp90" before compiling, or add to .bashrc)
#F90C = ifort
#FFLAGS = -O3 -w -fpp2 -DMPI 
#LAPACKL = -lscs -lguide -lmpi -ldl

#Intel fortran 8, check you have the latest update from the Intel web pages
#See Makefile_intel for ifc 7.1 or lower (some versions have problems)
#F90C     = ifort
#FFLAGS = -O2 -Vaxlib -ip -W0 -WB -openmp -fpp
#LAPACKL = -L/opt/intel/mkl61/lib/32 -lmkl_lapack -lmkl_ia32 -lguide -lpthread

#Intel fortran 8 on COSMOS
#F90C     = ifort
#FFLAGS = -O2 -Vaxlib -W0 -WB -openmp -fpp2 -lmpi -lscs  -ldl -DMPI
#LAPACKL = -L/opt/intel/mkl80/lib/64 -lmkl_lapack -lmkl_ipf -limf

#SGI, -mp toggles multi-processor. Use -O2 if -Ofast gives problems.
#Not various versions of the compiler are buggy giving erroneous seg faults with -mp.
#Version 7.3 is OK, currently version 7.4 is bugged, as are some earlier versions.
#F90C     = f90
#LAPACKL = -lcomplib.sgimath
#FFLAGS  = -Ofast -mp

#Digital/Compaq fortran, -omp toggles multi-processor
#F90C    = f90
#FFLAGS  = -omp -O4 -arch host -math_library fast -tune host -fpe1
#LAPACKL  = -lcxml

#G95; make sure LAPACK and MPI libs also compiled with g95
#F90C     = mpif90
#FFLAGS = -O2 -DMPI
#LAPACKL = /LAPACK/lapack_LINUX.a /LAPACK/blas_LINUX.a

#GFortran: if pre v4.3 add -D__GFORTRAN__
#F90C     = gfc
#FFLAGS = -O2 -ffree-form 
#LAPACKL = -Wl,-framework -Wl,accelerate

#Absoft ProFortran, single processor, set -cpu:[type] for your local system 
#F90C     = f95
#FFLAGS = -O2 -s -cpu:athlon -lU77 -w -YEXT_NAMES="LCS" -YEXT_SFX="_" 
#LAPACKL =  -llapack -lblas -lg2c 
#IFLAG = -p

#NAGF95, single processor:
#F90C     = f95
#FFLAGS = -DNAGF95 -O3
#LAPACKL = -llapack -lblas -lg2c

#Sun, single processor:
#F90C     = f90
#FFLAGS = -fast -ftrap=%none
#LAPACKL = -lsunperf -lfsu
#LAPACKL = -lsunperf -lfsu -lsocket -lm
#IFLAG = -M

#Sun MPI
#F90C = mpf90
#FFLAGS =  -O4 -openmp -ftrap=%none -dalign -DMPI
#LAPACKL = -lsunperf -lfsu -lmpi_mt
#IFLAG = -M

#Sun parallel enterprise:
#F90C     = f95
#FFLAGS =  -O4 -xarch=native64 -openmp -ftrap=%none
#LAPACKL = -lsunperf -lfsu
#IFLAG = -M


#IBM XL Fortran, multi-processor (run "module load lapack" then run "gmake")
# See also http://cosmocoffee.info/viewtopic.php?t=326
#F90C     = xlf90_r $(LAPACK)
#FFLAGS  = -WF,-DIBMXL -qsmp=omp -qsuffix=f=f90:cpp=F90 -O3 -qstrict -qarch=pwr3 -qtune=pwr3
#INCLUDE = -lessl
#LAPACKL =

PROPOSE = propose.o
CLSFILE = CMB_Cls_simple.o

#Can use params_H if you prefer more generic parameters
PARAMETERIZATION = params_CMB.o

WMAP_inc = $(IFLAG)$(cfitsio)/include $(IFLAG)$(WMAP)
WMAP_lib = -L$(cfitsio)/lib -L$(WMAP) -lcfitsio

CLSLIB  = $(WMAP_lib) -L../camb -lcamb $(LAPACKL) 

F90FLAGS = -DMATRIX_SINGLE $(FFLAGS) $(WMAP_inc) $(IFLAG)../camb $(INCLUDE)

DISTFILES = utils.o Matrix_utils.o settings.o GetDist.o

WMAPOBJS = $(WMAP)/read_archive_map.o \
	$(WMAP)/read_fits.o \
	$(WMAP)/WMAP_5yr_options.o \
	$(WMAP)/WMAP_5yr_util.o \
	$(WMAP)/WMAP_5yr_tt_pixlike.o \
	$(WMAP)/WMAP_5yr_teeebb_pixlike.o \
	$(WMAP)/WMAP_5yr_likelihood.o \
	$(WMAP)/WMAP_5yr_gibbs.o \
	$(WMAP)/WMAP_5yr_tt_beam_ptsrc_chisq.o \
	$(WMAP)/br_mod_dist.o 

OBJFILES= $(WMAPOBJS) utils.o Matrix_utils.o settings.o cmbtypes.o Planck_like.o \
	cmbdata.o WeakLen.o mpk.o supernovae.o SDSSLy-a-v3.o\
	$(CLSFILE) paramdef.o $(PROPOSE) $(PARAMETERIZATION) calclike.o \
	conjgrad_wrapper.o EstCovmat.o postprocess.o MCMC.o driver.o

default: cosmomc

all : cosmomc getdist

settings.o: utils.o
cmbtypes.o: settings.o
Planck_like.o: cmbtypes.o
cmbdata.o: Planck_like.o $(WMAPOBJS)
WeakLen.o: cmbtypes.o
mpk.o: cmbtypes.o
supernovae.o: cmbtypes.o
SDSSLy-a-v3.o: cmbtypes.o
$(CLSFILE): cmbtypes.o
paramdef.o: $(CLSFILE)
$(PROPOSE): paramdef.o
$(PARAMETERIZATION): paramdef.o
calclike.o: $(PARAMETERIZATION)
conjgrad_wrapper.o: calclike.o
EstCovmat.o: conjgrad_wrapper.o
postprocess.o: calclike.o
MCMC.o: calclike.o
driver.o: MCMC.o

.f.o:
	f77 $(F90FLAGS) -c $<

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $*.F90


cosmomc: $(OBJFILES) ../camb/libcamb.a 	
	$(F90C) -o ../cosmomc $(OBJFILES) $(CLSLIB) $(F90FLAGS)


clean:
	rm -f *.o *.mod *.d *.pc *.obj ../core

getdist: $(DISTFILES)
	$(F90C) -o ../getdist $(DISTFILES) $(CLSLIB) $(F90FLAGS) 
