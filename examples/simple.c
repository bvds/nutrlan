///
/// This file contains a simple example to call nuTRLan to compute
/// eigenpairs of a symmetric matrix.
/// @author Ichtaro Yamazaki, Kesheng Wu
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <math.h>
//#include "f2c.h"

#include "trlan.h"
#include "trl_map.h"
#include "trl_comm_i.h"
//
/// A simple matrix-vector multiplications routine.
/// Defines a diagonal matrix with values (1, 4, 9, 16, 25, ...).
///
#ifdef TRL_FOTRAN_COMPATIBLE
void diag_op(int *pnrow, int *pncol, double *xin, int *pldx,
	     double *yout, int *pldy) {
    //
    // ..
    // .. local variables ..
    int i, j, nrow, ncol, ldx, ldy;
    //
    // ..
    // .. executable statements ..
    nrow = *pnrow;
    ncol = *pncol;
    ldx  = *pldx;
    ldy  = *pldy;

    for( j=0; j<ncol; j++ ) {
	for( i=0; i<nrow; i++ ) {
	    yout[j*ldy+i] = (i+1.0)*(i+1.0)*xin[j*ldx+i];
	}
    }
}
#else
/* The extra parameter mvparam is not used in this case. */
void diag_op(const int nrow, const int ncol, const double *xin, const int ldx,
	     double *yout, const int ldy, void* mvparam) {
    //
    // ..
    // .. local variables ..
    int i, j;
    //
    // ..
    // .. executable statements ..
    for( j=0; j<ncol; j++ ) {
	for( i=0; i<nrow; i++ ) {
	    yout[j*ldy+i] = (i+1.0)*(i+1.0)*xin[j*ldx+i];
	}
    }
}
#endif

/// A really simple example of how to use TRLAN on a single processor.
int main( int argn, char **argv ) {
    static const int nrow=1897, lohi=-1, ned=5, maxlan=200, mev=100;
    int lwrk;
    // local variable declaration
    double eval[mev], evec[mev*nrow], exact[mev];
    double *res, *wrk;
    trl_info info;
    int i, j, k, fp, check;
    char name2[133], name[150];
    int tmp1, tmp2, nlen;
    // initialize info -- tell TRLAN to compute NEV smallest eigenvalues
    // of the diag_op
    lwrk=maxlan*(maxlan+10);
    if( lwrk > 0 ) {
	res = (double*)malloc(lwrk*sizeof(double));
	wrk = (double*)malloc(lwrk*sizeof(double));
    }
    trl_init_info( &info, nrow, maxlan, lohi, ned, 1.4901e-8, 1, 2000000, -1 );
    trl_set_iguess( &info, 0, 1, 0, NULL );
    // the Lanczos recurrence is set to start with [1,1,...,1]^T
    memset(eval, 0, mev*sizeof(double) );
    for( i=0; i<nrow; i++ ) {
	evec[i] = 1.0;
    }
    //info.verbose =  8;
    // call TRLAN to compute the eigenvalues
    trlan(diag_op, &info, nrow, mev, eval, evec, nrow, lwrk, res );
    trl_print_info(&info, 3*nrow);
    for( i=0; i<mev; i++ ) { // need to match with the definition in diag_op
	exact[i] = (i+1)*(i+1);
    }
    if( info.nec > 0 ) {
	i = info.nec;
    } else {
	i = mev - 1;
    }
    if( info.verbose >= 0) {
	trl_check_ritz( diag_op, &info, nrow, i, evec, nrow, eval,
			&check, res, exact, wrk, lwrk);
    } else {
	trl_terse_info( &info, 0 );
    }
    if( info.verbose > 1 ) {
	trl_rayleigh_quotients( diag_op, &info, i, evec, nrow, res, NULL );
	trl_check_ritz( diag_op, &info, nrow, i, evec, nrow, res, &check,
			&res[i], exact, wrk, lwrk );
    }
    if( info.nec == 0 ) info.nec = min(info.ned, mev-1);
    if( lwrk > 0 ) {
	free(res);
	free(wrk);
    }
    return 0;
} // main
