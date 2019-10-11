
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
//#include "f2c.h"
#include "mpi.h"

#include "trlan.h"
#include "trl_map.h"

////
// a simple matrix-vector multiplications routine
// defines a diagonal matrix with values (1, 4, 9, 16, 25, ....)
////
#ifdef TRL_FOTRAN_COMPATIBLE
void diag_op(int *pnrow, int *pncol, double *xin, int *pldx,
	     double *yout, int *pldy ) {
    // local variables
    int i, j, ioff, joff, doff, nrow, ncol, ldx, ldy;

    nrow = *pnrow;
    ncol = *pncol;
    ldx  = *pldx;
    ldy  = *pldy;

    MPI_Comm_rank(MPI_COMM_WORLD, &i);
    doff = nrow*i;
    for( j=0; j<ncol; j++ ) {
	ioff = j*ldx;
	joff = j*ldy;
	for( i=0; i<nrow; i++ ) {
	    //yout[joff+i] = (doff+i+1)*(doff+i+1)*xin[ioff+i];
	    yout[joff+i] = (doff+i+1)*xin[ioff+i];
	}
    }
}
#else
void diag_op(const int nrow, const int ncol, const double *xin, const int ldx,
	     double *yout, const int ldy, void *mvparam ) {
    // local variables
    int i, j, ioff, joff, doff;

    MPI_Comm_rank(MPI_COMM_WORLD, &i);
    doff = nrow*i;
    for( j=0; j<ncol; j++ ) {
	ioff = j*ldx;
	joff = j*ldy;
	for( i=0; i<nrow; i++ ) {
	    //yout[joff+i] = (doff+i+1)*(doff+i+1)*xin[ioff+i];
	    yout[joff+i] = (doff+i+1)*xin[ioff+i];
	}
    }
}
#endif
//// a really simple example of how to use TRLAN
int main( int argn, char **argv ) {
    // This set of parameters tell TRLAN to compute 5 (NED) smallest
    //(LOHI=-1) of a 100 x 100 (NROW) matrix by using a Krylov subspace
    // basis size of 30 (MAXLAN).
    // MEV defines the size of the arrays used to store computed
    // solutions, eval, evec.
    static const int nrow=200, mev=30;
    // local variable declaration
    int lohi, ned, maxlan, restart, lwrk;
    double eval[mev], evec[mev*nrow], exact[mev];
    double *res, *wrk;
    trl_info info;
    int i, j, k, n, fp, check, nbas, ii, nmis;
    char name2[133], name[150], file[150];
    int tmp1, tmp2, nlen;
    FILE *ferr;
    // initialize info -- tell TRLAN to compute NEV smallest eigenvalues
    // of the diag_op
    if( MPI_Init(&argn,&argv) != MPI_SUCCESS ) {
	printf( "Failed to initialize MPI.\r\n" );
    }
    ferr=fopen( "error.txt","w" );
    fclose(ferr);
    for( lohi=-1; lohi<=-1; lohi++ ) {
	for( ned=10; ned<=10; ned+=10 ) {
	    for( maxlan=200; maxlan<=200; maxlan+=100 ) {
		//nbas = max(1, maxlan-mev+1);
		//ii = nbas * ((nrow+3)/4)*4;;
		//nmis = maxlan*(maxlan+10);
		//lwrk=ii;
		//lwrk=0;
		lwrk=maxlan*(maxlan+10);
		if( lwrk > 0 ) {
		    res = (double*)malloc(lwrk*sizeof(double));
		    wrk = (double*)malloc(lwrk*sizeof(double));
		}
		for( restart=1; restart<=1; restart++ ) {
		    trl_init_info( &info, nrow, maxlan, lohi, ned, 1.4901/100000000, restart, 2000000, -1 );
		    if( info.my_pe <= 0 ) {
			ferr=fopen( "error.txt","a" );
			fprintf( ferr," \n ~~ iteration (ned=%d, maxlan=%d, restart=%d lohi=%d) ~~\n",ned,maxlan,restart,lohi );
			printf( " \r\n ~~ iteration (ned=%d, maxlan=%d, restart=%d lohi=%d) ~~\r\n",ned,maxlan,restart,lohi );
			fclose(ferr);
		    }
		    trl_set_iguess( &info, 0, 1, 0, NULL );
		    //trl_set_iguess( &info, 0, 2, 1, "CHECK" );
		    //trl_set_checkpoint( &info, 1, "CHECK" );
		    // the Lanczos recurrence is set to start with [1,1,...,1]^T
		    memset( eval, 0, mev*sizeof(double) );
		    //evec[0]=1.0;
		    //for( i=1; i<nrow; i++ ) {
		    //   evec[i] = 0.0;
		    //}
		    for( i=0; i<nrow; i++ ) {
			evec[i] = 1.0;
		    }
		    if( sprintf( file,"LOG_Nrow%dNed%dLan%dRes%dLoHi%d",nrow,ned,maxlan,restart,lohi ) < 0 ) {
			printf( "error writing file name\r\n" );
		    } else {
			trl_set_debug( &info, 0, file );
			// call TRLAN to compute the eigenvalues
			trlan(diag_op, NULL, &info, nrow, mev, eval, evec, nrow, lwrk, res );
			trl_print_info(&info, 3*nrow );
			if( info.nec > 0 ) {
			    i = info.nec;
			} else {
			    i = mev;
			}
			if( lohi == -1 ) {
			    for( j=0; j<i; j++ ) {
				exact[j] = (j+1)*(j+1);
			    }
			} else if( lohi == 1 ) {
			    for( j=0; j<i; j++ ) {
				exact[j] = (info.ntot-i+j+1)*(info.ntot-i+j+1);
			    }
			}
			lwrk = min( lwrk,i+nrow );
			lwrk = 4*i;
			if( lohi != 0 ) {
			    trl_check_ritz( diag_op, &info, nrow, i, evec, nrow, eval, &check, res, exact, wrk, lwrk );
			} else {
			    trl_check_ritz( diag_op, &info, nrow, i, evec, nrow, eval, &check, res, NULL, wrk, lwrk );
			}

			if( info.my_pe <= 0 ) {
			    ferr = fopen( "error.txt","a" );
			    if( info.stat != 0 ) {
				fprintf( ferr," ** TRLAN FAILED WITH %d **\n **",info.stat );
			    }
			    if( check < 0 ) {
				fprintf( ferr," -------- ERROR (%d converged) --------\n ",info.nec );
				printf( " -------- ERROR (%d converged) --------\r\n ",info.nec );
			    } else if( info.nec < info.ned ) {
				fprintf( ferr," ~~~~~~~~~ only %d converged ~~~~~~~~\n",info.nec );
				printf( " ~~~~~~~~~ only %d converged ~~~~~~~~\r\n",info.nec );
			    } else { 
				printf( " ********* SUCCESS (%d converged) **********\r\n",info.nec );
				fprintf( ferr," ********* SUCCESS (%d converged) **********\n",info.nec );
			    }
			    fclose(ferr);
			}
			if( info.nec == 0 ) info.nec = min(info.ned, mev-1);
		    }
		}
		if( lwrk > 0 ) {
		    free(res);
		    free(wrk);
		}
	    }
	}
    }
    MPI_Finalize();
    return 0;
}

