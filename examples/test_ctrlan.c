
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <math.h>

#include "trl_map.h"
#include "trlan.h"
////
// a simple matrix-vector multiplications routine
// defines a diagonal matrix with values (1, 4, 9, 16, 25, ....)
//
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
	    //yout[j*nrow+i] = (i+1)*xin[j*nrow+i];
	    yout[j*nrow+i] = (i+1)*(i+1)*xin[j*nrow+i];
	    //yout[j*nrow+i] = pow( (double)(i+1), 3.0 )*xin[j*nrow+i];
	}
    }
}
#else
void diag_op(const int nrow, const int ncol, const double *xin, const int ldx,
	     double *yout, const int ldy, void *mvparam) {
    // 
    // ..
    // .. local variables ..
    int i, j;
    //
    // ..
    // .. executable statements ..
    for( j=0; j<ncol; j++ ) {
	for( i=0; i<nrow; i++ ) {
	    //d = (double)(i+1);
	    //yout[j*nrow+i] = d*d*d*xin[j*nrow+i];
	    //yout[j*nrow+i] = (i+1)*xin[j*nrow+i];
	    yout[j*nrow+i] = (i+1)*(i+1)*xin[j*nrow+i];
	    //yout[j*nrow+i] = pow( (double)(i+1), 3.0 )*xin[j*nrow+i];
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
    static const int mev=100, nrow=10000; //1897; 
    // local variable declaration
    int lohi, ned, maxlan, restart, lwrk, lwrk2;
    double eval[mev], evec[mev*nrow], exact[mev];
    double *res, *wrk;
    trl_info info;
    int i, j, k, check, nbas, ii, nmis, ifact, rfact;
    char name2[133], name[150], file[150];
    int tmp1, tmp2, nlen;
    clock_t t1, t2;
    FILE *ferr, *fp;
    // initialize info -- tell TRLAN to compute NEV smallest eigenvalues
    // of the diag_op

    ferr=fopen( "error.txt","w" );
    fclose(ferr);
    for( lohi=-1; lohi<=-1; lohi++ ) {
	for( ned=100; ned<=100; ned+=10 ) {
	    for( maxlan=200; maxlan<=1000; maxlan+=100 ) {
		//nbas = max(1, maxlan-mev+1);
		//ii = nbas * ((nrow+3)/4)*4;;
		//nmis = maxlan*(maxlan+10);
		//lwrk=ii;
		//lwrk=0;
		//lwrk=maxlan*(maxlan+10);
		lwrk=1000*1100 + (nrow+4)*(maxlan-mev+1);
		if( lwrk > 0 ) {
		    res = (double*)malloc(lwrk*sizeof(double));
		    wrk = (double*)malloc(lwrk*sizeof(double));
		}
		for( restart=1; restart<=7; restart+=1 ) {
		    for( rfact=4; rfact<=4; rfact+=1 ) {
			trl_init_info( &info, nrow, maxlan, lohi, ned,
				       1.4901e-8, restart, 500000, -1 );
			trl_set_iguess( &info, 0, 1, 0, NULL );

			//
			// Manually setting some info parameters.....
			info.rfact = ((double)rfact)/10.0;
			//info.rfact = 1.8;
			//trl_set_restart( &info, 1.8 );
			info.mgap = ((double)ifact)/10.0;
			info.ref = 100.0;

			printf( "info.rfact=%e info.mgap=%e info.rfact=%e\n",
				info.rfact,info.mgap,info.rfact );
			ferr=fopen( "error.txt","a" );
			fprintf( ferr, " \n ~~ iteration (ned=%d, maxlan=%d, "
				 "restart=%d lohi=%d) ~~\n", ned, maxlan,
				 restart,lohi );
			printf( " \r\n ~~ iteration (ned=%d, maxlan=%d, "
				"restart=%d lohi=%d) ~~\r\n", ned, maxlan,
				restart,lohi );
			fclose(ferr);

			// start with [1,1,...,1]^T
			memset( eval, 0, mev*sizeof(double) );
			//for( i=0; i<nrow; i++ ) {
			//   evec[i] = 0.0;
			//}
			//evec[0] = 1.0;
			for( i=0; i<nrow; i++ ) {
			    evec[i] = 1.0;
			}
			if( sprintf( file,
				     "LOG_Nrow%dNed%dLan%dRes%dLoHi%dRfact%d",
				     nrow,ned,maxlan,restart,lohi,rfact ) < 0 )  {
			    printf( "error writing file name\r\n" );
			} else {
			    trl_set_debug( &info, 0, file );
			    // call TRLAN to compute the eigenvalues
			    t1 = clock();
			    trlan(diag_op, &info, nrow, mev, eval, evec, nrow,
				  lwrk, res );
			    t2 = clock();
			    printf( "TRLan: %d secs\n",
				    (int)((t2-t1)/CLOCKS_PER_SEC) );
			    trl_print_info(&info, 3*nrow );
			    i = info.nec;
			    if( i > mev ) {
				i = mev;
			    }
			    if( lohi == -1 ) {
				for( j=0; j<i; j++ ) {
				    exact[j] = (j+1)*(j+1);
				}
			    } else if( lohi == 1 ) {
				for( j=0; j<i; j++ ) {
				    exact[j] = (nrow-i+j+1)*(nrow-i+j+1);
				}
			    }
			    lwrk2 = min( lwrk,i+nrow );
			    lwrk2 = 4*i;
			    if( lohi != 0 ) {
				trl_check_ritz( diag_op, &info, nrow, i, evec,
						nrow, eval, &check, res, exact,
						wrk, lwrk2 );
			    } else {
				trl_check_ritz( diag_op, &info, nrow, i, evec,
						nrow, eval, &check, res, NULL,
						wrk, lwrk2 );
			    }

			    ferr = fopen( "error.txt","a" );
			    if( info.stat != 0 ) {
				fprintf( ferr,
					 " ** TRLAN FAILED WITH %d **\n **",
					 info.stat );
			    }
			    if( check < 0 ) {
				fprintf( ferr,
					 " -------- ERROR (%d converged) --------\n ",
					 info.nec );
				printf( " -------- ERROR (%d converged) --------\r\n ",
					info.nec );
			    } else if( info.nec < info.ned ) {
				fprintf( ferr,
					 " ~~~~~~~~~ only %d converged ~~~~~~~~\n",
					 info.nec );
				printf( " ~~~~~~~~~ only %d converged ~~~~~~~~\r\n",
					info.nec );
			    } else { 
				printf( " ********* SUCCESS (%d converged) **********\r\n",
					info.nec );
				fprintf( ferr,
					 " ********* SUCCESS (%d converged) **********\n",
					 info.nec );
			    }
			    fclose(ferr);
			    if( info.nec == 0 ) info.nec = min(info.ned, mev-1);
			}
		    }
		}
		if( lwrk > 0 ) {
		    free(res);
		    free(wrk);
		}
	    }
	}
    }
}
