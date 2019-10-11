/*
// ZTRLan routine 
// Lawrence Berkeley National Lab.
//
*/
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <float.h>
#include <mpi.h>

#include "trl_map.h"
#include "trlan.h"
#include "trlan_i.h"
#include "ztrl_comm_i.h"

////
void ztrl_g_dot_(void *mpicomp, int nrow, trl_dcomplex * v1, int ld1, int m1,
		 trl_dcomplex * v2, int ld2, int m2, trl_dcomplex * rr,
		 trl_dcomplex * wrk)
{
//
// Purpose:
// ========
// Implements a distributed version of BLAS routine dgemv, which is used to compute
// dot-products by TRLAN, i.e., wrk = [V1, V2]'*rr.
//
// Arguments:
// ==========
// mpicomp    (input) pointer to MPI_Comm
//             On entry, specifies MPI communicator.
//
// nrow       (input) integer
//             On entry, specifies, the number of rows on the local processor.
//
// v1         (input) double precision array of dimension (ld1,m1)
//             On entry, contains the first part of the matrix.
//
// ld1        (input) integer
//             On entry, specifies the leading dimension of v1.
//
// m1         (input) integer
//             On entry, specifies the number of columns in v1.
//
// v2         (input) double precision array of dimension (ld2,m2)
//             On entry, contains the second part of the matrix.
//
// ld2        (input) integer
//             On entry, specifies the leading dimension of v2.
//
// m2         (input) integer
//             On entry, specifies the number of columns in v2.
//
// rr         (input) double precision array of length (nrow)
//             On entry, contains the vector to be multiplied.
//
// wrk        (output) double precision array of length (m1+m2)
//             On exit, contains th results of this operation.  !! size not checked !!
//
//
// ..
// .. local parameters ..
    char trans = 'C';
    trl_dcomplex one, zero;
    integer c__1 = 1;
    MPI_Comm mpicom = *((MPI_Comm *) mpicomp);
//
// ..
// .. local scalars ..
    int i, nd, m1p1, npe;
//
// ..
// .. executable statements ..
    one.r = 1.0;
    one.i = 0.0;
    zero.r = 0.0;
    zero.i = 0.0;

    nd = m1 + m2;
    // nothing to do if both m1 and m2 are zero
    if (nd <= 0)
	return;
    // make sure the array sizes are correct
    if (ld1 < nrow || ld2 < nrow) {
	fprintf(stderr, "trl_g_dot: incorrect array sizes\n");
	exit(0);
    }
    m1p1 = m1 + 1;
    if (m1 > 2) {
	trl_zgemv(&trans, nrow, m1, one, v1, ld1, rr, c__1, zero, &wrk[nd],
		  c__1);
    } else if (m1 == 2) {
	wrk[nd].r = 0.0;
	wrk[nd].i = 0.0;
	wrk[nd + 1].r = 0.0;
	wrk[nd + 1].i = 0.0;
	for (i = 0; i < nrow; i++) {
	    wrk[nd].r += (v1[i].r * rr[i].r) + (v1[i].i * rr[i].i);
	    wrk[nd].i += (v1[i].r * rr[i].i) - (v1[i].i * rr[i].r);

	    wrk[nd + 1].r +=
		(v1[ld1 + i].r * rr[i].r) + (v1[ld1 + i].i * rr[i].i);
	    wrk[nd + 1].i +=
		(v1[ld1 + i].r * rr[i].i) - (v1[ld1 + i].r * rr[i].i);
	}
    } else if (m1 == 1) {
	trl_zdotc(&wrk[nd], nrow, v1, c__1, rr, c__1);
    }
    if (m2 > 2) {
	trl_zgemv(&trans, nrow, m2, one, v2, ld2, rr, c__1, zero,
		  &wrk[nd + m1p1 - 1], c__1);
    } else if (m2 == 2) {
	wrk[nd + m1p1 - 1].r = 0.0;
	wrk[nd + m1p1 - 1].i = 0.0;
	wrk[nd + m1p1].r = 0.0;
	wrk[nd + m1p1].i = 0.0;
	for (i = 0; i < nrow; i++) {
	    wrk[nd + m1p1 - 1].r +=
		(v2[i].r * rr[i].r) + (v2[i].i * rr[i].i);
	    wrk[nd + m1p1 - 1].i +=
		(v2[i].r * rr[i].i) - (v2[i].i * rr[i].r);

	    wrk[nd + m1p1].r +=
		(v2[ld2 + i].r * rr[i].r) + (v2[ld2 + i].i * rr[i].i);
	    wrk[nd + m1p1].r +=
		(v2[ld2 + i].r * rr[i].i) - (v2[ld2 + i].i * rr[i].r);
	}
    } else if (m2 == 1) {
	trl_zdotc(&(wrk[nd + m1p1 - 1]), nrow, v2, c__1, rr, c__1);
    }
    if (mpicom == MPI_COMM_SELF) {
	npe = 1;
    } else {
	i = MPI_Comm_size(mpicom, &npe);
    }
    if (npe > 1) {
	i = MPI_Allreduce(&wrk[nd], wrk, 2 * nd, MPI_DOUBLE, MPI_SUM,
			  mpicom);
	if (i != MPI_SUCCESS) {
	    printf("TRL_G_DOT: MPI_ALLREDUCE failed with error code %d.\n",
		   i);
	    MPI_Abort(mpicom, i);
	}
    } else {
	trl_zcopy(nd, &wrk[nd], c__1, wrk, c__1);
    }
//
//
}

////
void ztrl_g_sum(void *mpicomp, int nelm, trl_dcomplex * x, trl_dcomplex * y)
{
//
// Purpose:
// ========
// Performs global sum in the parallel environment, nothing is done here.
//
// Arguments:
// ==========
// mpicomp    (input) pointer to MPI_Comm
//            On entry, specifies the MPI communicator.
//
// nelm      (input) integer
//            On entry, specifies the number of elements in x and y.
//
// x         (input/output) double precision array of dimension nelm
//            On entry, contains the resulting vector on this processor.
//            On exit, contain the resulting vector of global sum.
//
// y         (workspace) double precision array of dimensioni nelm
//
    int c__1 = 1;
    int ierr;
    MPI_Comm mpicom = *((MPI_Comm *) mpicomp);
    extern int zcopy_();
    ierr =
	MPI_Allreduce(x, y, 2 * nelm, MPI_DOUBLE, MPI_SUM,
		      mpicom);
    zcopy_(&nelm, y, &c__1, x, &c__1);
    if (ierr != MPI_SUCCESS) {
	printf("TRL_G_SUM: MPI_ALLREDUCE failed with error code %d.\n",
	       ierr);
	MPI_Abort(mpicom, ierr);
    }
    return;
}
