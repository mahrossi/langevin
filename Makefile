.PHONY: all clean
FC := gfortran -x f95-cpp-input 
LD := gfortran
FFLAGS := -O3 -mtune=native -DUSELIBS
OBJS := $(patsubst %.f90,%.o,$(wildcard *.f90))

all: langevin.x

langevin.x: $(OBJS)
	$(LD) $(FFLAGS) -o langevin.x $(OBJS) -lblas -llapack

md-gle.o: module_asplines.o md-tools.o
langevin.o: module_asplines.o md-tools.o md-gle.o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm *~ *.o *.mod langevin.x
