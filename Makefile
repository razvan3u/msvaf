#FC=gfortran
#LDFLAGS=-llapack -lgomp
#FCFLAGS=-O2 -c -Wall -fopenmp

FC=ifort
LDFLAGS=-llapack 
#-openmp
FCFLAGS=-O2 -c -Wall
#FCFLAGS=-O2 -c 
#-openmp 


all: program.x

program.x: program.o functii.o
		$(FC) $(LDFLAGS) program.o functii.o jacobi.o -o program.x
		
functii.o: functii.f90 jacobi.o
		$(FC) $(FCFLAGS) functii.f90

jacobi.o: jacobi.f90 kinds.o
		$(FC) $(FCFLAGS) jacobi.f90

kinds.o: kinds.f90
		$(FC) $(FCFLAGS) kinds.f90

program.o: program.f90 functii.o
		$(FC) $(FCFLAGS) program.f90

clean:
		rm -f functii.o program.o kinds.o jacobi.o functii.mod kinds.mod jacobi.mod program.x

dist:
		tar -cvf fit.tar *.f90 Makefile	

