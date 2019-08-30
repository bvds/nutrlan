/*
// ZTRLan routine 
// Lawrence Berkeley National Lab.
//
*/
#include <stdio.h>

#ifdef F2C
#include "f2c.h"
#endif

#ifdef ESSL
#include "essl.h"
#endif


#include "trl_map.h"
#include "trlan.h"
#include "trlan_i.h"

//
double trl_ddot(integer n, const double *dx, integer incx,
		const double *dy, integer incy)
{

#ifdef BLAS_EXT
    extern double ddot_();
#endif

    return ddot_(&n, dx, &incx, dy, &incy);
}

//
//
void trl_zdotc(trl_dcomplex * ret_val, int n, trl_dcomplex * zx, int incx,
	       trl_dcomplex * zy, int incy)
{

#ifdef BLAS_EXT
    //extern trl_dcomplex* zdotc();
    //ret_val = zdotc( &n, zx, &incx, zy, &incy );
    extern void zdotc_();
#endif

    /*
       int i;
       for( i=0; i<n; i++ ) {
       printf( "zx[%d]=%e+%ei; zy[%d]=%e+%ei;\n",i,zx[i].r,zx[i].i,i,zy[i].r,zy[i].i );
       }
     */
    zdotc_(ret_val, &n, zx, &incx, zy, &incy);
    //printf( "%e+%ei\n", ret_val->r, ret_val->i );
}

//
//
void trl_dcopy(int n, double *dx, int incx, double *dy, int incy)
{

#ifdef BLAS_EXT
    extern int dcopy_();
#endif

    dcopy_(&n, dx, &incx, dy, &incy);
}

//
//
void trl_dgemv(char *trans, int m, int n, double alpha, double *a, int lda,
	       double *x, int incx, double beta, double *y, int incy)
{

#ifdef BLAS_EXT
    extern int dgemv_();
#endif

    dgemv_(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

//
//
void trl_daxpy(int n, double da, double *dx, int incx, double *dy,
	       int incy)
{

#ifdef BLAS_EXT
    extern int daxpy_();
#endif

    daxpy_(&n, &da, dx, &incx, dy, &incy);
}

//
//
void trl_dgemm(char *transa, char *transb, int m, int n, int k,
	       double alpha, double *a, int lda, double *b, int ldb,
	       double beta, double *c, int ldc)
{

#ifdef BLAS_EXT
    extern int dgemm_();
#endif

    dgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c,
	   &ldc);
}

//
//
void trl_dscal(int n, double da, double *dx, int incx)
{

#ifdef BLAS_EXT
    extern int dscal_();
#endif

    dscal_(&n, &da, dx, &incx);

}

//
//
void trl_zaxpy(int n, trl_dcomplex za, trl_dcomplex * zx, int incx,
	       trl_dcomplex * zy, int incy)
{

#ifdef BLAS_EXT
    extern int zaxpy_();
#endif

    zaxpy_(&n, &za, zx, &incx, zy, &incy);

}

//
//
void trl_zgemv(char *trans, int m, int n, trl_dcomplex alpha,
	       trl_dcomplex * a, int lda, trl_dcomplex * x, int incx,
	       trl_dcomplex beta, trl_dcomplex * y, int incy)
{
#ifdef BLAS_EXT
    extern int zgemv_();
#endif

    zgemv_(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);

}

//
//
void trl_zgemm(char *transa, char *transb, int m, int n, int k,
	       trl_dcomplex alpha, trl_dcomplex * a, int lda,
	       trl_dcomplex * b, int ldb, trl_dcomplex beta,
	       trl_dcomplex * c, int ldc)
{

#ifdef BLAS_EXT
    extern int zgemm_();
#endif

    zgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c,
	   &ldc);

}

//
//
void trl_zscal(int n, trl_dcomplex za, trl_dcomplex * zx, int incx)
{

#ifdef BLAS_EXT
    extern int zscal_();
#endif

    zscal_(&n, &za, zx, &incx);

}

//
//
void trl_zdscal(int n, double da, trl_dcomplex * zx, int incx)
{

#ifdef BLAS_EXT
    extern int zdscal_();
#endif

    zdscal_(&n, &da, zx, &incx);

}

//
//
void trl_zcopy(int n, trl_dcomplex * zx, int incx, trl_dcomplex * zy,
	       int incy)
{

#ifdef BLAS_EXT
    extern int zcopy_();
#endif

    zcopy_(&n, zx, &incx, zy, &incy);

}
