
#ifndef __TRLMAP_H
#define __TRLMAP_H

#ifdef BLAS
#define		STD_DATA
#define		BLAS_EXT
#endif

#ifdef F2C
#include "f2c.h"
typedef doublecomplex trl_dcomplex;
#define		BLAS_EXT
#endif

#ifdef ESSL
#include "essl.h"
#define STD_DATA
#endif


#ifdef STD_DATA
typedef int integer;
typedef double doublereal;
typedef int logical;
typedef struct {
    double r, i;
} trl_dcomplex;
#define min(a,b) ( a < b ? a : b )
#define max(a,b) ( a > b ? a : b )
#endif

#ifdef CBLAS
typedef struct {
    double r, i;
} doublecomplex;
#define VOID			void
#endif


#ifdef NOCHANGE

#define lsame_		lsame
#define xerbla_		xerbla

#define dlamch_		dlamch
#define dlanst_		dlanst
#define dlascl_		dlascl
#define dlaev2_		dlaev2
#define dlapy2_		dlapy2
#define dlartg_		dlartg
#define dlasr_		dlasr
#define dorgtr_		dorgtr
#define dportrf_	dportrf
#define dportrs_	dportrs
#define dstein_		dstein
#define dtrtrs_         dtrtrs
#define dsyev_		dsyev
#define dsytrd_		dsytrd

#endif

double trl_ddot(integer, const double *, integer, const double *, integer);
void trl_zdotc(trl_dcomplex * ret_val, int n, trl_dcomplex * zx, int incx,
	       trl_dcomplex * zy, int incy);

void trl_daxpy(int n, double da, double *dx, int incx, double *dy,
	       int incy);
void trl_dcopy(int n, double *dx, int incx, double *dy, int incy);
void trl_dgemm(char *transa, char *transb, int m, int n, int k,
	       double alpha, double *a, int lda, double *b, int ldb,
	       double beta, double *c, int ldc);
void trl_dgemv(char *trans, int m, int n, double alpha, double *a, int lda,
	       double *x, int incx, double beta, double *y, int incy);
void trl_dscal(int n, double da, double *dx, int incx);

void trl_zaxpy(int n, trl_dcomplex za, trl_dcomplex * zx, int incx,
	       trl_dcomplex * zy, int incy);
void trl_zcopy(int n, trl_dcomplex * zx, int incx, trl_dcomplex * zy,
	       int incy);
void trl_zgemv(char *trans, int m, int n, trl_dcomplex alpha,
	       trl_dcomplex * a, int lda, trl_dcomplex * x, int incx,
	       trl_dcomplex beta, trl_dcomplex * y, int incy);
void trl_zscal(int n, trl_dcomplex za, trl_dcomplex * zx, int incx);
void trl_zdscal(int n, double da, trl_dcomplex * zx, int incx);
void trl_zgemm(char *transa, char *transb, int m, int n, int k,
	       trl_dcomplex alpha, trl_dcomplex * a, int lda,
	       trl_dcomplex * b, int ldb, trl_dcomplex beta,
	       trl_dcomplex * c, int ldc);

#endif

