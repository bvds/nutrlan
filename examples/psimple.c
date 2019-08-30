#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"
#include "trlan.h"

//
// This file contains an example to compute parallel nuTRLan to compute
// eigenpairs of a symmetric matrix.
//
////
// a simple matrix-vector multiplications routine with
// a diagonal matrix with values (1, 4, 9, 16, 25, ....)
#ifdef TRL_FOTRAN_COMPATIBLE
void diag_op(int *pnrow, int *pncol, double *xin, int *pldx,
	     double *yout, int *pldy) {
    // local variables
    int i, j, ioff, joff, doff, nrow, ncol, ldx, ldy;

    nrow = *pnrow;
    ncol = *pncol;
    ldx  = *pldx;
    ldy  = *pldy;

    MPI_Comm_rank(MPI_COMM_WORLD, &doff);
    doff *= nrow;
    for( j=0; j<ncol; j++ ) {
	ioff = j*ldx;
	joff = j*ldy;
	for( i=0; i<nrow; i++ ) {
	    yout[joff+i] = (doff+i+1.0)*(doff+i+1.0)*xin[ioff+i];
	}
    }
}
#else
void diag_op(const int nrow, const int ncol, const double *xin, const int ldx,
	     double *yout, const int ldy, void *mvparam) {
    // local variables
    int i, j, ioff, joff, doff;

    MPI_Comm_rank(MPI_COMM_WORLD, &doff);
    doff *= nrow;
    for( j=0; j<ncol; j++ ) {
	ioff = j*ldx;
	joff = j*ldy;
	for( i=0; i<nrow; i++ ) {
	    yout[joff+i] = (doff+i+1.0)*(doff+i+1.0)*xin[ioff+i];
	}
    }
}
#endif
//
/// a really simple example of how to use nuTRLan (parallel)
int main( int argn, char **argv ) {
    static const int nrow=1000, lohi=-1, ned=5, maxlan=100, mev=10;
    int lwrk=maxlan*(maxlan+10);
    // local variable declaration
    double eval[mev], evec[mev*nrow], exact[mev];
    double res[lwrk], wrk[lwrk];
    trl_info info;
    int i, check;
    // initialize MPI
    if( MPI_Init(&argn,&argv) != MPI_SUCCESS ) {
	printf( "Failed to initialize MPI.\r\n" );
    }
    // initialize info -- tell TRLAN to compute NEV smallest eigenvalues
    // of the diag_op
    trl_init_info( &info, nrow, maxlan, lohi, ned, 1.4901e-8, 7, 5000, -1 );
    trl_set_iguess( &info, 0, 1, 0, NULL );
    // the Lanczos recurrence is set to start with [1,1,...,1]^T
    memset( eval, 0, mev*sizeof(double) );
    for( i=0; i<nrow; i++ ) {
	evec[i] = 1.0;
    }
    //info.verbose =  8;
    trl_set_debug( &info, 0, "LOG" );
    // call TRLAN to compute the eigenvalues
    trlan(diag_op, &info, nrow, mev, eval, evec, nrow, lwrk, res );
    trl_print_info( &info, 3*nrow );
    //trl_terse_info( &info, 0 );
    for( i=0; i<mev; i++ ) {
	exact[i] = (i+1);
    }
    if( info.nec > 0 ) {
	i = info.nec;
    } else {
	i = mev - 1;
    }
    trl_check_ritz(diag_op, &info, nrow, i, evec, nrow, eval, &check, res,
		   exact, wrk, lwrk);
    if( info.verbose > 1 ) {
	trl_rayleigh_quotients( diag_op, &info, i, evec, nrow, res, NULL );
	trl_check_ritz(diag_op, &info, nrow, i, evec, nrow, res, &check,
		       &res[i], exact, wrk, lwrk);
    }
    MPI_Finalize();
}
