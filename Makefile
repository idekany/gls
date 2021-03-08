FC = gfortran
FFLAGS = -fopenmp -O3 -fno-range-check
LFLAGS =
OBJECTS = gls_mod.o mc_mod.o gls.o

.PHONY: clean

FNAME ?= gls

$(FNAME) : $(OBJECTS)
	$(FC) $(FFLAGS) $(LFLAGS) $(OBJECTS) -o $(FNAME)

%.o : %.f90
	$(FC) $(FFLAGS) $(LFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(FNAME)
