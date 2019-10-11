#
# Makefile for examples, using Nutrlan
#
TOP = ${PWD}/..
include ../Make.inc
#
# a simple example (serial version)
  simple : simple.o
	$(LOADER) $(CFLAGS) simple.o $(NUTRLAN) $(LIB) -o simple

  simple.o : simple.c
	$(LOADER) $(CFLAGS) simple.c -c $(INC) -o simple.o
#
# a simple example (parallel version)
  psimple : psimple.o
	$(PLOADER) $(CFLAGS) psimple.o $(PNUTRLAN) $(PLIB) -o psimple

  psimple.o : psimple.c
	$(PLOADER) $(CFLAGS) psimple.c -c $(INC) -o psimple.o
#
# test nutrlan 
  test_ctrlan : test_ctrlan.o
	$(LOADER) $(CFLAGS) test_ctrlan.o $(NUTRLAN) $(LIB) -o test_ctrlan

  test_ctrlan.o : test_ctrlan.c
	$(LOADER) $(CFLAGS) test_ctrlan.c -c $(INC) -o test_ctrlan.o
#
# test ivj
  test_ijv : test_ijv.o
	$(LOADER) $(CFLAGS) test_ijv.o $(NUTRLAN) $(LIB) -o test_ijv 

  test_ijv.o : test_ijv.c
	$(LOADER) $(CFLAGS) test_ijv.c -c $(INC) -o test_ijv.o 
#
# test binary
  test_pbinary : test_pbinary.o read_aij.o
	$(LOADER) $(CFLAGS) test_pbinary.o read_aij.o $(NUTRLAN) $(LIB) -o test_pbinary

  test_pbinary.o : test_pbinary.c
	$(LOADER) $(CFLAGS) test_pbinary.c -c $(INC) -o test_pbinary.o

  read_aij.o : Read_Aij.C
	$(LOADER) $(CFLAGS) Read_Aij.C -c $(INC) -o read_aij.o
#
# test ztrlan
  test_ztrlan : test_ztrlan.o $(NUTRLAN)
	$(LOADER) $(CFLAGS) test_ztrlan.o $(NUTRLAN) $(LIB) -o test_ztrlan

  test_ztrlan.o : test_ztrlan.c
	$(LOADER) $(CFLAGS) test_ztrlan.c -c $(INC) -o test_ztrlan.o
#
# test ctrlan (parallel version)
  test_ptrlan : test_ptrlan.o 
	$(PLOADER) $(CFLAGS) test_ptrlan.o $(PNUTRLAN) $(PLIB) -o test_ptrlan

  test_ptrlan.o : test_ptrlan.c
	$(PLOADER) $(CFLAGS) test_ptrlan.c -c $(INC) -o test_ptrlan.o
#
# test ztrlan (parallel version)
  test_pztrlan : test_pztrlan.o
	$(PLOADER) $(CFLAGS) test_pztrlan.o $(PNUTRLAN) $(PLIB) -o test_pztrlan

  test_pztrlan.o : test_pztrlan.c
	$(PLOADER) $(CFLAGS) test_pztrlan.c -c $(INC) -o test_pztrlan.o
#
# simple.f
  fsimple : fsimple.o
	$(FC) $(OPF) fsimple.o $(NUTRLAN) $(LIB) -o fsimple

  fsimple.o : simple.f
	$(FC) $(CFLAGS) simple.f -c $(FINC) -o fsimple.o
#
# zsimple.f
  fzsimple : fzsimple.o
	$(FC) $(OPF) fzsimple.o $(NUTRLAN) $(LIB) -o fzsimple

  fzsimple.o : zsimple.f
	$(FC) $(CFLAGS) zsimple.f -c $(FINC) -o fzsimple.o

#(NUTRLAN)::
#	 cd $(HOME)/SRC && $(MAKE) lib_serial
#
# clean
  clean :
	(rm -f *.o *~ LOG* simple psimple test_ctrlan test_ztrlan test_ptrlan test_pztrlan *.exe)
cleanall: clean
	cd .. && $(MAKE) clean
