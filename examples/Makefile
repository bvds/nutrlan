#
# Makefile for examples, using Nutrlan
#
TOP = ${PWD}/..
include ../Make.inc
NUTRLAN_UO = user_ortho.o $(NUTRLAN)
#
# a simple example (serial version)
simple: simple.o user_ortho.o
	$(LOADER) $(OPT) simple.o $(NUTRLAN_UO) $(LIB) -o simple

simple.o : simple.c
	$(LOADER) $(OPT) simple.c -c $(INC) -o simple.o
user_ortho.o : user_ortho.c
	$(LOADER) $(OPT) user_ortho.c -c $(INC) -o user_ortho.o

#
# a simple example (parallel version)
psimple : psimple.o
	$(PLOADER) $(OPT) psimple.o $(NUTRLAN_UO) $(PLIB) -o psimple

psimple.o : psimple.c
	$(PLOADER) $(OPT) psimple.c -c $(INC) -o psimple.o
#
# test nutrlan 
test_ctrlan : test_ctrlan.o
	$(LOADER) $(OPT) test_ctrlan.o $(NUTRLAN_UO) $(LIB) -o test_ctrlan $(DEBUG)

test_ctrlan.o : test_ctrlan.c
	$(LOADER) $(OPT) test_ctrlan.c -c $(INC) -o test_ctrlan.o $(DEBUG)
#
# test ivj
test_ijv : test_ijv.o
	$(LOADER) $(OPT) test_ijv.o $(NUTRLAN_UO) $(LIB) -o test_ijv 

  test_ijv.o : test_ijv.c
	$(LOADER) $(OPT) test_ijv.c -c $(INC) -o test_ijv.o 
#
# test binary
test_pbinary : test_pbinary.o read_aij.o
	$(LOADER) $(OPT) test_pbinary.o read_aij.o $(NUTRLAN_UO) $(LIB) -o test_pbinary

test_pbinary.o : test_pbinary.c
	$(LOADER) $(OPT) test_pbinary.c -c $(INC) -o test_pbinary.o

read_aij.o : Read_Aij.C
	$(LOADER) $(OPT) Read_Aij.C -c $(INC) -o read_aij.o
#
# test ztrlan
test_ztrlan : test_ztrlan.o $(NUTRLAN_UO)
	$(LOADER) $(OPT) test_ztrlan.o $(NUTRLAN_UO) $(LIB) -o test_ztrlan $(DEBUG)

test_ztrlan.o : test_ztrlan.c
	$(LOADER) $(OPT) test_ztrlan.c -c $(INC) -o test_ztrlan.o $(DEBUG)
#
# test ctrlan (parallel version)
test_ptrlan : test_ptrlan.o 
	$(PLOADER) $(OPT) test_ptrlan.o $(NUTRLAN_UO) $(PLIB) -o test_ptrlan

test_ptrlan.o : test_ptrlan.c
	$(PLOADER) $(OPT) test_ptrlan.c -c $(INC) -o test_ptrlan.o
#
# test ztrlan (parallel version)
test_pztrlan : test_pztrlan.o
	$(PLOADER) $(OPT) test_pztrlan.o $(NUTRLAN_UO) $(PLIB) -o test_pztrlan

  test_pztrlan.o : test_pztrlan.c
	$(PLOADER) $(OPT) test_pztrlan.c -c $(INC) -o test_pztrlan.o
#
# simple.f
fsimple : fsimple.o
	$(FC) $(OPF) fsimple.o $(NUTRLAN_UO) $(LIB) -o fsimple

fsimple.o : simple.f
	$(FC) $(OPT) simple.f -c $(FINC) -o fsimple.o
#
# zsimple.f
fzsimple : fzsimple.o
	$(FC) $(OPF) fzsimple.o $(NUTRLAN_UO) $(LIB) -o fzsimple

fzsimple.o : zsimple.f
	$(FC) $(OPT) zsimple.f -c $(FINC) -o fzsimple.o

#(NUTRLAN)::
#	 cd $(HOME)/SRC && $(MAKE) lib_serial
#
# clean
  clean :
	(rm -f *.o *~ LOG* simple psimple test_ctrlan test_ztrlan test_ptrlan test_pztrlan *.exe)
cleanall: clean
	cd .. && $(MAKE) clean
