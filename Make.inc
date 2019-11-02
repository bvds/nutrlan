############################################################
# nuTRLan make include file for a plain linux installation with LAPACK,
# BLAS and MPI in /usr/local
############################################################
#
# C compiler
  CC      = cc
  MPICC   = mpicc
  LOADER  = cc
  PLOADER = mpicc
#
CFLAGS = -DBLAS
# Debuging
  CFLAGS += -O -Wall -g # -Wextra
# Fast
#  CFLAGS += -O3
#
#  On AMD Epyc, add flags suggested by
#  http://www.prace-ri.eu/best-practice-guide-amd-epyc
#
ifeq ($(shell lscpu|grep -m 1 -c -i epyc),1)
CFLAGS += -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -fomit-frame-pointer
endif
#
# Make command
  MAKE = make
#  MAKEFLAGS = TOP=${PWD}
#
# nuTRLan home directory
#  TOP    = ${PWD}
  NUTRLAN = $(TOP)/libnutrlan.a
  PNUTRLAN = $(TOP)/libnutrlan_mpi.a
#
# External directories
  DIR_CLAPACK = 
#
# External libraries
  LIB_DIR =
  MATH    = -lm
ifeq (1,0) # Use Blis
  BLAS    = -lblis
  LAPACK  = -lflame
else
ifeq (1,1) # Use OpenBLAS
  BLAS    = -L/opt/OpenBLAS/lib -lopenblas
  LAPACK  = -lflame
else   # Use MKL
  MKLROOT = /opt/intel/mkl
  BLAS    = 
  LAPACK  = 
  CLAPACK = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
endif
endif
  F2C     = 
  LIB     = $(LIB_DIR) $(LAPACK) $(BLAS) $(CLAPACK) $(F2C) $(MATH) -lgomp
  PLIB    = $(LIB) -lmpi 
#
# Include files
  INC_CLAPACK =
#
#   /usr/include/mpich or /usr/include/openmpi
  INC_NUTRLAN  = -I$(TOP)/INC -I/usr/include/openmpi-x86_64/
  INC         = $(INC_CLAPACK) $(INC_NUTRLAN)
#
# Object files
  OBJ_REAL = dsort2.o dstqrb.o trlan.o trlaux.o trlcore.o restart.o trlmap.o
  OBJ_CPLX = zdgemm.o zdgemv.o zdaxpy.o ztrlan.o ztrlaux.o ztrlcore.o
  OBJ = $(OBJ_REAL) $(OBJ_CPLX)
  OBJ_SERIAL = $(OBJ) trl_comm_none.o ztrl_comm_none.o
  OBJ_PARALLEL = $(OBJ) trl_comm_mpi.o ztrl_comm_mpi.o
