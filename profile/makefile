#                      Makefile for ham2d

F90 = ifort -pg
#F90 = gfortran
CC  = icc -pg
#CC  = gcc
FFLAGS =  -w -O2 -r8 
#FFLAGS  =  -w -O2 -lgfortran -fdefault-real-8 # uncomment this for gfortran
CFLAGS =  -w -O2 #-fbacktrace -fbounds-check # -g #-traceback  -check uninit  #-warn-all#-g #-O2 #-g
LFLAGS =$(CFLAGS) -lm
OBJECTS = ham2d.o readGrid.o preprocess.o find_faces.o\
	  initflow.o stepSolution.o updateSoln.o computeRHS.o\
  computeRHSk.o DualcomputeRHSk.o computeRHSkv.o computeLinearRHS.o flux_roe2d.o wallFlux.o periodic_bc.o apply_periodic.o apply_periodic_LHS.o\
	  muscld.o muscld_deriv.o roeflx.o computeForce.o outputSolution.o wrest.o\
	  smoothGrid.o jac_roe.o flux_visc.o triSolvers.o mathops.o findDiagonals.o\
	  ADI.o DualADI.o DADI.o jac_visc.o gaussSeidel.o lineGaussSeidel.o \
          weno3.o weno.c weno_deriv.c\

# Link Instruction
ham2d: $(OBJECTS)
	$(CC) $(OBJECTS) $(LFLAGS) -o ham2d

clean:
	@rm -rf *.o *.mod *.*~ ham2d

%.o:%.F90
	$(F90) $(FFLAGS) -c $< -o $*.o

%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.c
	$(CC) $(CFLAGS) -c $< -o $*.o

%.o:%.mod
