
include ../Make.inc

all: $(OBJ_CBLAS)
	ar -cru $(CTRLAN) $(OBJ_CBLAS)

clean: 
	(rm -f *.o)

.c.o:
	$(CC) $(INC) -c $(OPT) $< $(DEBUG)

