
#Intel MPI
#these settings for ifort 13 and higher
#Can remove -xHost if your cluster is not uniform, or specify specific processor optimizations -x...
F90C     = mpif90
FFLAGS = -mkl -openmp -O3 -xHost -no-prec-div -fpp -DMPI
LAPACKL = 

#gfortran; does not yet work yet, open GCC 4.9 compiler bugs are
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=60232
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=60234
#F90C = gfortran
#FFLAGS = -cpp

#use "make RECOMBINATION=cosmorec" to build with CosmoRec rather than RECFAST default
RECOMBINATION ?=recfast

EXTDATA=

#use PRECISION=SINGLE to use single precision
PRECISION ?=

#set WMAP empty not to compile with WMAP, e.g. WMAP = /scratch/../WMAP9/likelihood_v5
WMAP ?=

#Needed for WMAP; set from ~/.bashrc definition or edit here
cfitsio ?=

#empty not to compile with CLIK, set from ~/.bashrc
PLANCKLIKE ?=

#These is not used in public cosmomc
highL ?=
#highL = ../highL

IFLAG = -I
INCLUDE =

#clik path and library
ifeq ($(PLANCKLIKE),cliklike)
CLIKL = -L$(CLIKPATH)/lib -lclik_f90
INCLUDE = $(IFLAG)$(CLIKPATH)/include
endif


#COSMOS: use "module load cosmolib latest"
#use "runCosmomc" (globally installed) to run, defining required memory usage
ifeq ($(COSMOHOST),cosmos)
F90C = ifort
FFLAGS = -openmp -fast -w -fpp2 -DMPI
LAPACKL = -mkl=sequential -lmkl_lapack -lmpi
cfitsio = $(CFITSIO)
WMAP = $(COSMOLIB)/WMAP9
GSLPATH = $(GSL_ROOT)
endif

ifeq ($(COSMOHOST),darwin)
WMAP = /scratch/aml1005/wmap_likelihood_v5
cfitsio = /usr/local/Cluster-Apps/cfitsio/intel/3.300
FFLAGS = -mkl -openmp -O3 -xHost -no-prec-div -fpp -DMPI
#FFLAGS = -mkl -openmp -g -check all -traceback -no-prec-div -fpp -DMPI -fpe0
endif


#would like to embed this somehow..
#GIT_HASHTAG = "git log --pretty=format:'%H:%cd' -n 1"

PROPOSE = propose.o

CALCFILE = Calculator_CAMB.o

#Can use params_H if you prefer more generic parameters
PARAMETERIZATION = CosmologyParameterizations.o

F90FLAGS = $(FFLAGS) $(IFLAG)../camb $(INCLUDE)
LINKFLAGS = -L../camb -lcamb_$(RECOMBINATION) $(LAPACKL) $(F90CRLINK) $(CLIKL)

DISTFILES = StringUtils.o MpiUtils.o FileUtils.o IniObjects.o ObjectLists.o ObjectParamNames.o \
    Matrix_utils_new.o Random.o settings.o samples.o IO.o GetDist.o



ifneq ($(PRECISION),)
FFLAGS += -D$(PRECISION) -DMATRIX_$(PRECISION)
endif

ifneq ($(PLANCKLIKE),)
FFLAGS += -DCLIK
endif

ifneq ($(highL),)
FFLAGS += -DhighL
LINKFLAGS += -L$(highL) -lhigh
INCLUDE += $(IFLAG)$(highL)
endif

PLANCKLIKEFILES=


ifneq ($(PLANCKLIKE),)
PLANCKLIKEFILES += $(PLANCKLIKE).o
endif


SUPERNOVAE = supernovae_Union2.o supernovae_SNLS.o

DATAMODULES = $(PLANCKLIKEFILES) power_spec.o wigglez.o mpk.o bao.o $(SUPERNOVAE) supernovae.o HST.o

LIKEFILES =  DataLikelihoods.o calclike.o CalcLike_Cosmology.o

OBJFILES = ArrayUtils.o StringUtils.o MpiUtils.o FileUtils.o IniObjects.o ObjectLists.o ObjectParamNames.o \
    Matrix_utils_new.o Interpolation.o Random.o settings.o samples.o IO.o likelihood.o GeneralTypes.o cmbtypes.o \
	Planck_like.o bbn.o Calculator_Cosmology.o $(CALCFILE) $(DATAMODULES) cmbdata.o $(PROPOSE) BaseParameters.o \
	ParamSet.o $(PARAMETERIZATION) $(LIKEFILES) EstCovmat.o PowellConstrainedMinimize.o ImportanceSampling.o \
	MCMC.o minimize.o SampleCollector.o GeneralSetup.o CosmologyConfig.o driver.o

F90CRLINK =

ifeq ($(RECOMBINATION),cosmorec)
## This is flag is passed to the Fortran compiler allowing it to link C++ (uncomment the right one).
# GCC (gfortran/g++)
COSMOREC_PATH ?= ../CosmoRec/
F90CRLINK = -L$(COSMOREC_PATH) -lCosmoRec -L$(GSLPATH)/lib -lgsl -lgslcblas -lstdc++
# Intel Compilers (ifort/icpc)
#F90CRLINK = -cxxlib -L$(COSMOREC_PATH) -lCosmoRec -L$(GSLPATH)/lib -lgsl -lgslcblas
FFLAGS +=  -DCOSMOREC
endif

ifeq ($(RECOMBINATION),hyrec)
HYREC_PATH ?= ../HyRec/
F90CRLINK += -L$(HYREC_PATH) -lhyrec
endif

default: cosmomc

all : cosmomc getdist

GetDist.o: IO.o samples.o
supernovae.o: $(SUPERNOVAE)
IO.o: ObjectParamNames.o settings.o
samples.o: ObjectLists.o settings.o
likelihood.o: ObjectLists.o settings.o ObjectParamNames.o
DataLikelihoods.o: likelihood.o ParamSet.o $(DATAMODULES)
cliklike_CamSpec.o: temp_like.o cmbtypes.o
cliklike.o: IniObjects.o cmbtypes.o
settings.o: ObjectParamNames.o 
cmbtypes.o: settings.o likelihood.o GeneralTypes.o
Planck_like.o: cmbtypes.o
cmbdata.o: Planck_like.o
bbn.o: settings.o likelihood.o Interpolation.o
powerspec.o: cmbtypes.o
wigglez.o: cmbtypes.o power_spec.o
mpk.o: cmbtypes.o power_spec.o
bao.o: cmbtypes.o
HST.o: cmbtypes.o
supernovae.o: cmbtypes.o
Calculator_Cosmology.o: cmbtypes.o
$(CALCFILE): Calculator_Cosmology.o
$(PROPOSE): settings.o
BaseParameters: setings.o 
ParamSet.o: $(CALCFILE) $(PROPOSE) samples.o
$(PARAMETERIZATION): ParamSet.o bbn.o
calclike.o: DataLikelihoods.o
postprocess.o: calclike.o
MCMC.o: calclike.o
driver.o: EstCovmat.o MCMC.o minimize.o $(PARAMETERIZATION)
minimize.o: PowellConstrainedMinimize.o calclike.o

ifneq ($(highL),)
cliklike_CamSpec.o: $(highL)/libhigh.a
endif


ifneq ($(WMAP),)
cmbdata.o: $(WMAP)/libwmap9.a
F90FLAGS += $(IFLAG)$(cfitsio)/include $(IFLAG)$(WMAP)
LINKFLAGS +=  -L$(cfitsio)/lib -L$(WMAP) -lcfitsio -lwmap9
else
F90FLAGS += -DNOWMAP
endif


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
	rm -f ../camb/*.o ../camb/*.obj ../camb/*.mod ../camb/libcamb_$(RECOMBINATION).a

cleancosmomc:
	rm -f *.o *.mod *.d *.pc *.obj ../core


getdist: camb $(DISTFILES)
	$(F90C) -o ../getdist $(DISTFILES) $(LINKFLAGS) $(F90FLAGS)

#If you want to use the ppf module, change equations_ppf to equations
#and halofit_ppf to halofit also comment out wa_ppf in CMB_Cls_simple.f90
camb:
	cd ../camb && \
	$(MAKE) --file=Makefile_main libcamb_$(RECOMBINATION).a \
	RECOMBINATION=$(RECOMBINATION) EQUATIONS=equations_ppf NONLINEAR=halofit_ppf

$(highL)/libhigh.a:
	cd $(highL); make libhigh.a;

$(WMAP)/libwmap9.a:
	cd $(WMAP); make libwmap9.a;
