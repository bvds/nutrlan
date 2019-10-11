/*
// ZTRLan routine 
// Lawrence Berkeley National Lab.
//
//
// All communication routines used by TRLAN are in this file. This file contains routines to be
// used on sequential or shared memory environments.   No actual data exchange is performed.
*/
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>

#include "trl_map.h"
#include "trlan.h"
#include "trlan_i.h"
#include "trlaux_i.h"
#include "trl_comm_i.h"

////
//void trl_init_info(trl_info *info, int nrow, int mxlan, int lohi, int ned, 
//                    int nopts, ... ) {
void trl_init_info(trl_info * info, int nrow, int mxlan, int lohi,
		    int ned, double tol, int restart, int maxmv,
		    void *mpicomp)
{
//
// Purpose:
// ========
// Initializes a TRL_INFO variable. This function must be called before calling
// any other user level routine in TRLAN package.
//
// Arguments:
// ==========
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the information
//          about the eigenvalue problem and the progress of TRLAN.
//          On exit, points to the initialized data structure.
//
// nrow    (input) integer
//          On entry, specifies the local dimension of the problem.
//
// mxlan   (input) integer
//          On entry, specifies the maximum number of basis vectors to be used.
//
// lohi    (input)  integer
//          On entry, specifies, which end of the spectrum to compute:
//           lohi < 0, then lower end, the smallest eigenvalues
//           lohi > 0, then high end, the largest eigenvalues
//           lohi = 0, then either lower and or high end, whoever converges first
//          are computed.
//
// ned      (input) integer
//           On entry, specifies the number of wanted eigenvalues and eigenvectors.
//
// tol      (optional) double precision
//           If provided, specifies the tolerance on residual norm. By default,
//           tol is set to be sqrt(epsilon).
//
// trestart (optional) integer
//           If provided, specifies the thick-restarting scheme, 1-4. Default is 1.
//
// mxmv     (optionial) integer
//           If provided, specifies the maximum number of matrix-vector 
//           multiplication allowed. By default, mxmv is set to be 
//           (info->ntot)*(info->ned).
//
// mpicomp  (ignored) pointer to MPI_Comm
//           If provided, specifites the MPI communicator. By default, it is a 
//           duplicate of MPI_COMM_WORLD. In sequential case, this is set to 0.
//
//
// ..
// .. executable statements ..
//  va_list argptr;
//  va_start( argptr,nopts );
//  if( nopts > 0 ) {
    if (tol > 0) {
	//info->tol = va_arg( argptr, double );
	info->tol = tol;
	if (info->tol <= DBL_MIN) {
	    info->tol = DBL_EPSILON;
	} else if (info->tol > 1.0) {
	    info->tol = min(0.1, 1.0 / (info->tol));
	}
    } else {
	info->tol = sqrt(DBL_EPSILON);
    }
    //if( nopts > 1 ) {
    if (restart > 0) {
	//info->restart = va_arg( argptr,int );
	info->restart = restart;
    } else {
	info->restart = 0;
    }
    //if( nopts > 2 ) {
    if (maxmv > 0) {
	//info->maxmv = va_arg( argptr,int );
	info->maxmv = maxmv;
    } else {
	info->maxmv = min(max(info->ntot, 1000), 1000 * info->ned);
    }
    if (mpicomp != NULL) {
        fprintf(stderr, "TRL_INIT_INFO:  ignoring MPI Comm.\n");
    }
    info->mpicomp = NULL;
    //va_end( argptr );
    // setup the rest of arguments
    info->maxlan = mxlan;
    if (mxlan <= ned) {
	info->maxlan = ned + max(ned, 6);
    }
    info->lohi = lohi;
    info->ned = ned;
    info->nloc = nrow;
    info->ntot = nrow;
    info->guess = 0;
    info->nec = 0;
    info->locked = info->nec;
    info->matvec = 0;
    info->nloop = 0;
    info->north = 0;
    info->nrand = 0;
    info->flop = 0;
    info->rflp = 0;
    info->flop_h = 0;
    info->rflp_h = 0;
    info->flop_r = 0;
    info->rflp_r = 0;
    info->clk_rate = CLOCKS_PER_SEC;
#ifdef __64INT
    info->clk_max = 9223372036854775807LL;
#else
    info->clk_max = (clock_t) (pow(2.0, (8.0 * sizeof(clock_t) - 1.0))-1.0);
    if( info->clk_max < 0 ) {
       if( sizeof(clock_t) == 8 ) {
         info->clk_max = 9223372036854775807LL;
       } else {
         printf( "error initializing clock.\n" );
       }
    }
#endif
    if( (double)(info->clk_max) <= 0 ) printf( "??\n" );
    //info->clk_max = -1;
    info->clk_tot = 0;
    info->clk_op = 0;
    info->clk_orth = 0;
    info->clk_res = 0;
    info->tick_t = 0;
    info->tick_o = 0;
    info->tick_h = 0;
    info->tick_r = 0;
    info->clk_in = 0;
    info->clk_out = 0;
    info->wrds_in = 0;
    info->wrds_out = 0;
    info->verbose = 0;
    info->stat = 0;
    info->anrm = 0;
    info->tmv = -1;
    info->trgt = -DBL_MAX;
    info->tres = -1.0;
    info->crat = 1.0;

    info->predicted_crate = 0.0;
    info->old_target = 0.0;
    info->target_id = 0;
    info->ref = 0.0;
    info->avgm = 0.0;
    info->k1 = 0;
    info->k2 = 0;
    info->k = 0;
    info->mvparam = 0;

    info->my_pe = 0;
    info->npes = 1;
    info->cpflag = 0;
    strcpy(info->oldcpf, "");
    // log file pointer
    info->log_io = 99;
    strcpy(info->log_file, "");
    // checkpoint file pointer
    info->cpio = 98;
    strcpy(info->cpfile, "");
//
// .. end of trl_init_info_ ..
//
}

////
void trl_g_sum(void *mpicomp, int nelm, double *x, double *y)
{
//
// Purpose:
// ========
// Performs global sum in the parallel environment, nothing is done here.
//
// Arguments:
// ==========
// mpicomp   (ignored) pointer to MPI_Comm
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
}

////
int trl_sync_flag(void *mpicomp, int inflag)
{
//
// Purpose:
// ========
// Given an integer value, returns the minimum value of all the PEs
//
// Arguments:
// ==========
// mpicomp   (ignored) pointer to MPI_Comm
//            On entry, specifies the MPI communicator.
//
// inflag    (inpuut) integer
//            On entry, specifies the integer value from this processor.
//
    return inflag;
}

////
void trl_g_dot_(void *mpicomp, int nrow, double *v1, int ld1, int m1,
		double *v2, int ld2, int m2, double *rr, double *wrk)
{
//
// Purpose:
// ========
// Implements a distributed version of BLAS routine dgemv, which is used to compute
// dot-products by TRLAN, i.e., wrk = [V1, V2]'*rr.
//
// Arguments:
// ==========
// mpicomp    (ignored) pointer to MPI_Comm
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
// ..
// .. local parameters ..
    char trans = 'T';
    double one = 1.0, zero = 0.0;
    integer c__1 = 1;
//
// ..
// .. local scalars ..
    int i, nd;
//
// ..
// .. executable statements ..
    nd = m1 + m2;
    // nothing to do if both m1 and m2 are zero
    if (nd <= 0)
	return;
    // make sure the array sizes are correct
    if (ld1 < nrow || ld2 < nrow) {
	fprintf(stderr, "trl_g_dot: incorrect array sizes\n");
	exit(0);
    }
    if (m1 > 2) {
	trl_dgemv(&trans, nrow, m1, one, v1, ld1, rr, c__1, zero, wrk,
		  c__1);
    } else if (m1 == 2) {
	wrk[0] = zero;
	wrk[1] = zero;
	for (i = 0; i < nrow; i++) {
	    wrk[0] += v1[i] * rr[i];
	    wrk[1] += v1[ld1 + i] * rr[i];
	}
    } else if (m1 == 1) {
	wrk[0] = trl_ddot(nrow, v1, c__1, rr, c__1);
    }
    if (m2 > 2) {
	trl_dgemv(&trans, nrow, m2, one, v2, ld2, rr, c__1, zero, &wrk[m1],
		  c__1);
    } else if (m2 == 2) {
	wrk[m1] = zero;
	wrk[nd - 1] = zero;
	for (i = 0; i < nrow; i++) {
	    wrk[m1]     += v2[i]        * rr[i];
	    wrk[nd - 1] += v2[ld2 + i] * rr[i];
	}
    } else if (m2 == 1) {
	wrk[m1] = trl_ddot(nrow, v2, c__1, rr, c__1);
    }
//
// .. end of trl_g_dot_ ..
//
}
