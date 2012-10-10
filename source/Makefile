#use "make RECOMBINATION=cosmorec" to build with CosmoRec rather than RECFAST default

#You may need to edit the library paths for MKL for Intel

##Uncomment the next line to include dr7 LRG
EXTDATA =
#EXTDATA = LRG

RECOMBINATION ?=recfast

#set WMAP empty not to compile with WMAP
WMAP = /home/aml1005/WMAP7/likelihood_v4

#Only needed for WMAP
cfitsio = /usr/local/cfitsio/intel10/64/3.040

#GSL only needed for DR7 LRG
GSLPATH ?= /home/aml1005/libs/gsl

#empty not to compile with CLIK, otherwise set to name of cliklike.f90 file
PLANCKLIKE ?= cliklike_CamSpec
#PLANCKLIKE =cliklike_CamSpec
#PLANCKLIKE =cliklike


IFLAG = -I
INCLUDE =

#clik path and library
ifeq ($(PLANCKLIKE),cliklike)
CLIKPATH = /path/to/clik
CLIKL = -L$(CLIKPATH)/lib -lclik_f90
INCLUDE = $(IFLAG)$(CLIKPATH)/include
endif



#Intel MPI (assuming mpif77 set to point to ifort)
#these settings for ifort 11.1 and higher; may need to add explicit link directory otherwise
#Can add -xHost if your cluster is uniform, or specify specific processor optimizations -x...
#If getdist gives segfaults remove openmp when compiling getdist
F90C     = mpif90
FFLAGS = -O3 -W0 -WB -openmp -fpp -DMPI -vec_report0 -mkl=parallel
LAPACKL = -lmpi

#GFortran: defaults for v4.5; if pre v4.3 add -D__GFORTRAN__
#F90C     = gfortran
#in earlier versions use FFLAGS = -O2 -ffree-form -x f95-cpp-input -D__GFORTRAN__
#FFLAGS =  -O3 -fopenmp -ffree-form -x f95-cpp-input -ffast-math -march=native -funroll-loops
#LAPACKL = -Wl,-framework -Wl,accelerate
#commented above is (I think) for Mac; this is standard linux (sudo apt-get install liblapack-dev)
#LAPACKL = -lblas -llapack

#HPCF settings. Use Intel 9 or 10.1+, not 10.0
#F90C     = mpif90
#FFLAGS = -O2 -Vaxlib -W0 -WB -openmp -fpp -DMPI -vec_report0
#LAPACKL = -L/usr/local/Cluster-Apps/intel/mkl/10.2.2.025/lib/em64t -lmkl_lapack -lmkl -lguide -lpthread
#GSLPATH = /usr/local/Cluster-Apps/gsl/1.9
#cfitsio = /usr/local/Cluster-Users/cpac/cmb/2.1.0/cfitsio
#INCLUDE=

#COSMOS: use "module load cosmolib latest"
#use "runCosmomc" (globally installed) to run, defining required memory usage
ifeq ($(COSMOHOST),cosmos)
F90C = ifort
FFLAGS = -openmp -fast -w -fpp2 -DMPI
LAPACKL = -mkl=sequential -lmkl_lapack -lmpi
cfitsio = $(CFITSIO)
WMAP = $(COSMOLIB)/WMAP7
GSLPATH = $(GSL_ROOT)
endif

#G95; make sure LAPACK and MPI libs also compiled with g95
#F90C     = mpif90
#FFLAGS = -O2 -cpp -DMPI
#LAPACKL = /LAPACK/lapack_LINUX.a /LAPACK/blas_LINUX.a

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

#Absoft ProFortran, single processor, set -cpu:[type] for your local system
#F90C     = f95
#FFLAGS = -O2 -s -cpu:athlon -lU77 -w -YEXT_NAMES="LCS" -YEXT_SFX="_"
#LAPACKL =  -llapack -lblas -lg2c
#IFLAG = -p

#NAGF95, single processor:
#F90C     = f95
#FFLAGS = -DNAGF95 -O3
#LAPACKL = -llapack -lblas -lg2c

#PGF90
#F90C = pgf90
#FFLAGS = -O2 -DESCAPEBACKSLASH -Mpreprocess
#LAPACKL = -llapack -lblas


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

F90FLAGS = -DMATRIX_SINGLE $(FFLAGS) $(IFLAG)../camb $(INCLUDE)
LINKFLAGS = -L../camb -lcamb_$(RECOMBINATION) $(LAPACKL) $(F90CRLINK) $(CLIKL)

DISTFILES = ParamNames.o Matrix_utils.o settings.o IO.o GetDist.o

ifneq ($(PLANCKLIKE),)
FFLAGS += -DCLIK
endif

PLANCKLIKEFILES=

ifeq ($(PLANCKLIKE),cliklike_CamSpec)
PLANCKLIKEFILES += temp_like.o
endif

LIKEFILES = $(PLANCKLIKEFILES) $(PLANCKLIKE).o calclike.o

OBJFILES = ParamNames.o Matrix_utils.o settings.o IO.o cmbtypes.o Planck_like.o  \
	cmbdata.o WeakLen.o bbn.o lrggettheory.o mpk.o bao.o supernovae.o HST.o SDSSLy-a-v3.o \
	$(CLSFILE) paramdef.o $(PROPOSE) $(PARAMETERIZATION) $(LIKEFILES) \
	EstCovmat.o PowellConstrainedMinimize.o minimize.o postprocess.o MCMC.o driver.o


ifeq ($(EXTDATA),LRG)
F90FLAGS += -DDR71RG
OBJFILES += bsplinepk.o
GSLINC = -I$(GSLPATH)/include
LINKFLAGS += -L$(GSLPATH)/lib -lgsl -lgslcblas
endif

F90CRLINK =

ifeq ($(RECOMBINATION),cosmorec)
## This is flag is passed to the Fortran compiler allowing it to link C++ (uncomment the right one).
# GCC (gfortran/g++)
COSMOREC_PATH ?= ../CosmoRec/
F90CRLINK = -lstdc++ -L$(COSMOREC_PATH) -lCosmoRec -L$(GSLPATH)/lib -lgsl -lgslcblas
# Intel Compilers (ifort/icpc)
#F90CRLINK = -cxxlib -L$(COSMOREC_PATH) -lCosmoRec -L$(GSLPATH)/lib -lgsl -lgslcblas
FFLAGS +=  -DCOSMOREC
endif

ifeq ($(RECOMBINATION),hyrec)
HYREC_PATH ?= ../HyRec/
F90CRLINK += -L$(HYREC_PATH) -lhyrec
endif

ifneq ($(WMAP),)
F90FLAGS += $(IFLAG)$(cfitsio)/include $(IFLAG)$(WMAP)
LINKFLAGS +=  -L$(cfitsio)/lib -L$(WMAP) -lcfitsio

OBJFILES += $(WMAP)/read_archive_map.o \
	$(WMAP)/read_fits.o \
	$(WMAP)/healpix_types.o \
	$(WMAP)/WMAP_7yr_options.o \
	$(WMAP)/WMAP_7yr_util.o \
	$(WMAP)/WMAP_7yr_tt_pixlike.o \
	$(WMAP)/WMAP_7yr_teeebb_pixlike.o \
	$(WMAP)/WMAP_7yr_likelihood.o \
	$(WMAP)/WMAP_7yr_gibbs.o \
	$(WMAP)/WMAP_7yr_tt_beam_ptsrc_chisq.o \
	$(WMAP)/br_mod_dist.o
else
F90FLAGS += -DNOWMAP
endif

default: cosmomc

all : cosmomc getdist

cliklike_CamSpec.o: temp_like.o
settings.o: ../camb/libcamb_$(RECOMBINATION).a
cmbtypes.o: settings.o
Planck_like.o: cmbtypes.o
cmbdata.o: Planck_like.o  $(WMAPOBJS)
WeakLen.o: cmbtypes.o
bbn.o: settings.o
mpk.o: cmbtypes.o lrggettheory.o
bao.o: cmbtypes.o mpk.o
HST.o: cmbtypes.o
supernovae.o: cmbtypes.o
SDSSLy-a-v3.o: cmbtypes.o
$(CLSFILE): cmbtypes.o mpk.o HST.o bao.o IO.o
paramdef.o: $(CLSFILE)
$(PROPOSE): paramdef.o
$(PARAMETERIZATION): paramdef.o
calclike.o: $(PARAMETERIZATION)
postprocess.o: calclike.o
MCMC.o: calclike.o
driver.o: EstCovmat.o MCMC.o minimize.o
minimize.o: PowellConstrainedMinimize.o calclike.o

export FFLAGS
export F90C

.f.o:
	f77 $(F90FLAGS) -c $<

%.o: %.c
	$(CC) $(GSLINC) -c $*.c

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $*.F90


cosmomc: camb $(OBJFILES)
	$(F90C) -o ../cosmomc $(OBJFILES) $(LINKFLAGS) $(F90FLAGS)


clean: cleancosmomc
	rm -f ../camb/*.o ../camb/*.obj ../camb/*.mod

cleancosmomc:
	rm -f *.o *.mod *.d *.pc *.obj ../core


getdist: camb $(DISTFILES)
	$(F90C) -o ../getdist $(DISTFILES) $(LINKFLAGS) $(F90FLAGS)

camb:
	cd ../camb && $(MAKE) --file=Makefile_main libcamb_$(RECOMBINATION).a RECOMBINATION=$(RECOMBINATION) EQUATIONS=equations_ppf
