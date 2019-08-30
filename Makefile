TOP    = ${PWD}
include Make.inc

################################
# Makefile for nuTRLan library
################################
lib:
	(cd ./SRC; $(MAKE) TOP="$(TOP)" lib_serial)

plib:
	(cd ./SRC; $(MAKE) TOP="$(TOP)" lib_parallel)

blib:
	(cd ./SRC; $(MAKE) TOP="$(TOP)" lib_serial; cd ../CBLAS; $(MAKE) TOP="$(TOP)" all )

pblib:
	(cd ./SRC; $(MAKE) TOP="$(TOP)" lib_parallel; cd ../CBLAS; $(MAKE) TOP="$(TOP)" all )

ftrlan:
	(cd ./FORTRAN; $(MAKE) TOP="$(TOP)" all; cd ../SRC; $(MAKE) TOP="$(TOP)" fortran)

simple: lib
	cd examples && $(MAKE) TOP="$(TOP)" simple

psimple: plib
	cd examples && $(MAKE) TOP="$(TOP)" psimple

clean:
	rm -f $(NUTRLAN)
	cd SRC && $(MAKE) TOP="$(TOP)" clean
	cd CBLAS && $(MAKE) TOP="$(TOP)" clean
	cd examples && $(MAKE) TOP="$(TOP)" clean
	cd FORTRAN && $(MAKE) TOP="$(TOP)" clean

.PHONY: lib plib blib pblib ftrlan simple psimple clean
