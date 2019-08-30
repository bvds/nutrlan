
include ../Make.inc

all: trlan_interface.o trlan_info.o 

clean:
	rm -f *.mod *.o *~

.f.o:
	$(FC) -c $(FOPT) $< $(DEBUG)

