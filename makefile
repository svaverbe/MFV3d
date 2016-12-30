SHELL=/bin/sh 

FCC=mpif90 # Name of the compiler

# FCC=mpiifort

PROG=MFV3d # Name and path of the executable

FFLAGS=-free -O3 -g2 -debug -traceback -check all -C -heap-arrays -mcmodel=large
FFLAGS=-free -mcmodel=large -O3 -g2 -g -debug -C -fcheck=all -Wall
# FFLAGS=-free -mcmodel=large -O3 -g2 -debug -fcheck=all -C 

# FFLAGS=-free -O3 -g2 -C -mcmodel=large

# FFLAGS=-O3 -mcmodel=large

OBJS=main.o pararange.o advance.o gradients.o limiters.o init.o gravity.o hvol.o FVeqns.o \
tstep.o eos.o reconstruction.o linkedlist.o Riemannsolvers.o RMD.o output.o

# OBJS=MFV3d.o

# Object files
# Lines from here on are the rules that make uses to build the executable

all: $(PROG)

%.o : %.f90
	$(FCC) -c  $<

$(PROG): $(OBJS)
	$(FCC) $(FFLAGS) -o  $@ $(OBJS)
	
clean: 
	rm *.o
