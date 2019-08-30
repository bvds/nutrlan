############################################################
# nuTRLan make include file
############################################################
#
# C compiler
  FC      = ftn
  CC      = cc
  MPICC   = cc
  LOADER  = cc
  PLOADER = cc
#
# Debuging
  DEBUG = 
  OPT = -fastsse -DBLAS
#  If linking to CLAPACK
#  OPT = -O3
#
# Make command
  MAKE = make
#  MAKEFLAGS = TOP=${PWD}
#
# nuTRLan home directory
#  TOP = ${PWD}
  NUTRLAN = $(TOP)/libnutrlan.a
#
# External directories
#  DIR_CLAPACK = /root/ichi/ricpack-1.0/CLAPACK
#
# External libraries
  LIB_DIR = 
  MATH    = 
  BLAS    =
  LAPACK  = 
  SCILIB  = \
            -L/opt/acml/4.0.1a/pgi64/lib -lacml -lacml_mv /opt/pgi/7.2.4/linux86-64/7.2-4/lib/pgi.ld \
            -L/opt/mpt/3.0.2/xt/mpich2-pgi/lib -L/opt/acml/4.0.1a/pgi64/lib -L/opt/xt-libsci/10.2.1/pgi/snos64/lib \
            -L/opt/mpt/3.0.2/xt/sma/lib -L/opt/mpt/3.0.2/xt/util/lib -L/opt/mpt/3.0.2/xt/pmi/lib \
            -L/opt/xt-catamount/2.0.53b/lib/cnos64 -L/opt/xt-pe/2.0.53b/lib/snos64 -L/opt/xt-pe/2.0.53b/cnos/linux/64/lib \
            -L/opt/xt-service/2.0.53b/lib/snos64 -L/opt/pgi/7.2.4/linux86-64/7.2-4/libso -L/opt/pgi/7.2.4/linux86-64/7.2-4/lib \
            -L/usr/lib64 -L/usr/lib64/gcc-lib/x86_64-suse-linux/3.3.3 \
            -ldl -lsci_quadcore -lacml -lpgftnrtl -lsci -lsma -lmpich -lrt -lpct -lpmi -lalpslli -lalpsutil -lportals -lpthread \
            -Wl,-rpath,/opt/pgi/7.2.4/linux86-64/7.2-4/libso -Wl,-rpath,/opt/pgi/7.2.4/linux86-64/7.2-4/lib -lnspgc -lpgc -lgcc_eh -lmpichf90 \
            -Wl,-rpath,/opt/pgi/7.2.4/linux86-64/7.2-4/libso -Wl,-rpath,/opt/pgi/7.2.4/linux86-64/7.2-4/lib \
            -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgf90rtl -lm \
            -L/opt/mpt/3.0.2/xt/mpich2-pgi/lib -L/opt/acml/4.0.1a/pgi64/lib -L/opt/xt-libsci/10.2.1/pgi/snos64/lib \
            -L/opt/mpt/3.0.2/xt/sma/lib -L/opt/mpt/3.0.2/xt/util/lib -L/opt/mpt/3.0.2/xt/pmi/lib -L/opt/xt-catamount/2.0.53b/lib/cnos64 \
            -L/opt/xt-pe/2.0.53b/lib/snos64 -L/opt/xt-pe/2.0.53b/cnos/linux/64/lib -L/opt/xt-service/2.0.53b/lib/snos64 \
            -L/opt/pgi/7.2.4/linux86-64/7.2-4/libso -L/opt/pgi/7.2.4/linux86-64/7.2-4/lib -L/usr/lib64 -L/usr/lib64/gcc-lib/x86_64-suse-linux/3.3.3 \
            -ldl -lsci_quadcore -lacml -lpgftnrtl -lsci -lsma -lmpich -lrt -lpct -lpmi -lalpslli -lalpsutil -lportals -lpthread \
            -Wl,-rpath,/opt/pgi/7.2.4/linux86-64/7.2-4/libso -Wl,-rpath,/opt/pgi/7.2.4/linux86-64/7.2-4/lib -lnspgc -lpgc -lgcc_eh -ldl
#  CLAPACK = $(DIR_CLAPACK)/lapack_IA64.a
#  F2C     = -lF77 -lI77
  LIB     = $(LIB_DIR) $(LAPACK) $(BLAS) $(CLAPACK) $(F2C) $(MATH) $(SCILIB)
  PLIB    = $(LIB)
#
# Include files
#  INC_CLAPACK = -I$(DIR_CLAPACK)/F2CLIBS \
#                -I$(DIR_CLAPACK)/BLAS/WRAP
  INC_NUTRLAN  = -I$(TOP)/INC
  INC         = $(INC_CLAPACK) $(INC_NUTRLAN)
  FINC        = -I$(TOP)/FORTRAN
  CFLAGS      = -DTRL_FORTRAN_COMPATIBLE
#
# Object files
  OBJ_REAL = dsort2.o dstqrb.o trlan.o trlaux.o trlcore.o restart.o trlmap.o
  OBJ_CPLX = zdgemm.o zdgemv.o zdaxpy.o ztrlan.o ztrlaux.o ztrlcore.o
  OBJ = $(OBJ_REAL) $(OBJ_CPLX)
  OBJ_SERIAL = $(OBJ) trl_comm_none.o ztrl_comm_none.o
  OBJ_PARALLEL = $(OBJ) trl_comm_mpi.o ztrl_comm_mpi.o
  FOBJ = ftrlan.o
