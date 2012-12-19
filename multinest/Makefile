AR = ar r  
LINKLIB = ld -shared  
LIBDIR = ./

FFLAGS = -O2 -ip -openmp -fpp -DMPI -vec_report0 -mkl=parallel 


 
NSOBJECTS = utils.o utils1.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o nested.o

%.o: %.f90
	$(F90C) $(FFLAGS) -c -o $@ $^ 

%.o: %.F90
	$(F90C) $(FFLAGS) -c -o $@ $^ 

 
all: libnest3.a
 
libnest3.a: $(NSOBJECTS) 
	$(AR) $@ $^ 

clean: 
	-rm $(LIBDIR)/libnest3.*  *.o *.mod
