TOP = ${PWD}/..
include ../Make.inc

lib_serial: $(OBJ) $(OBJ_SERIAL)
	ar -cru $(NUTRLAN) $(OBJ) $(OBJ_SERIAL)

lib_parallel: $(OBJ) $(OBJ_PARALLEL)
	ar -cru $(PNUTRLAN) $(OBJ) $(OBJ_PARALLEL)

fortran: $(FOBJ) 
	ar -cru $(NUTRLAN) $(FOBJ) 

clean:
	rm -f *.o *~

#$(OBJ_PARALLEL):
#	$(MPICC) $(INC) -c $(OPT) $< $(DEBUG)

# on seaborg, replace $(CC) with $(MPICC)
.c.o:
	$(CC) $(INC) -c $(OPT) $< $(DEBUG) $(CFLAGS)
