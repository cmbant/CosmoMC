#Makefile --- High ell Likelihood code...

HIGHELLLIB = libhigh.a
OBJS =  Highell_options.o \
	Highell_subroutines.o \
	Foregrounds_loading.o \
	ACT_south_likelihood.o \
	ACT_equa_likelihood.o \
	SPT_keisler_likelihood.o \
	SPT_reichardt_likelihood.o \
	Highell_likelihood.o

RM = rm -f

LIBS ?= $(LAPACKL)

F90C    ?= mpif90intel 
#problems with -openmp in v12 intel compilers
FFLAGS  =  -mkl -O2 -ip -W0 -WB -fpp2 -vec_report0 -openmp
LIBS  ?=

#			Rules.
#
all: $(HIGHELLLIB) test

check: test
	./test

$(HIGHELLLIB): $(OBJS)
	ar r $@ $^

test: $(OBJS) test.o
	$(F90C) $(FFLAGS) -o $@ test.o $(OBJS) $(LIBS)

%.o: %.f90
	$(F90C) $(FFLAGS) $(INCS) -c -o $@ $<

%.o: %.F90
	$(F90C) $(FFLAGS) $(INCS) -c -o $@ $<

clean:
	$(RM) *.o *.mod *.log *~ *.a test

distclean: clean
	$(RM) test

