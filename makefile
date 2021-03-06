# Start of the makefile
# Defining variables
#.SUFFIXES:
#.SUFFIXES: .f90 .o
objects = convert.o main_BTO.o Sub_B.o  ff.o
ic = ifort
#CFLAGS=-O3
# Makefile
execname: $(objects)
		$(ic) $(objects)
convert.mod: convert.o convert.f90
		$(ic) -c convert.f90
convert.o: convert.f90
		$(ic) -c convert.f90
main_BTO.o: convert.mod main_BTO.f90
		$(ic) -c main_BTO.f90
Sub_B.o: Sub_B.f90
		$(ic) -c Sub_B.f90
ff.o: ff.f90
		$(ic) -c ff.f90

#%.o : %.f90
#		@$(ic) $(CFLAGS) -c –o $@ $<

clean:
		rm $(objects)
		rm a.out convert.mod

