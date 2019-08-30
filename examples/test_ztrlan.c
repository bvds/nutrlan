//
// When CLAPACk is complied with -DNO_BLAS_WRAP, then "blaswrap.h" is not required here.
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>

#include "trl_map.h"
#include "ztrlan.h"


#define max_rank 500
#define max_nrow 10000

int rank_k;
double rank_alpha[max_rank];
trl_dcomplex rank_v[max_rank*max_nrow];
trl_dcomplex rank_t[max_rank];

////
//
void rankk_init( int nrow, int k ) {
    static int c__1 = 1;
    int i, j;
    trl_dcomplex d__1;
    //
    rank_k = k;
    for( i=0; i<k; i++ ) {
	rank_alpha[i] = (i+1)*(i+1);
	for( j=0; j<nrow; j++ ) {
	    rank_v[i*nrow+j].r = drand48();
	    rank_v[i*nrow+j].i = drand48();
	}

	for( j=0; j<i; j++ ) {
	    trl_zdotc( &d__1, nrow, &(rank_v[j*nrow]), c__1, &(rank_v[i*nrow]), c__1 );
	    d__1.r = -d__1.r;
	    d__1.i = -d__1.i;
	    trl_zaxpy( nrow, d__1, &(rank_v[j*nrow]), c__1, &(rank_v[i*nrow]), c__1 );
	}
	for( j=0; j<i; j++ ) {
	    trl_zdotc( &d__1, nrow, &(rank_v[j*nrow]), c__1, &(rank_v[i*nrow]), c__1 );
	    d__1.r = -d__1.r;
	    d__1.i = -d__1.i;
	    trl_zaxpy( nrow, d__1, &(rank_v[j*nrow]), c__1, &(rank_v[i*nrow]), c__1 );
	}
	trl_zdotc( &d__1, nrow, &(rank_v[i*nrow]), c__1, &(rank_v[i*nrow]), c__1 );
	d__1.r = 1.0/sqrt(d__1.r);
	d__1.i = 0.0;
	trl_zscal( nrow, d__1, &(rank_v[i*nrow]), c__1 );
	//for( j=0; j<nrow; j++ ) {
	//   printf( "V(%d,%d)=%e+%e*i;\r\n",j+1,i+1,rank_v[i*nrow+j].r,rank_v[i*nrow+j].i );
	//}
    }
}
#ifdef TRL_FOTRAN_COMPATIBLE
// a simple matrix-vector multiplications routine
// defines a diagonal matrix with values (1, 4, 9, 16, 25, ....)
//
void zdiag_op(int *pnrow, int *pncol, trl_dcomplex *xin, int *pldx, trl_dcomplex *yout, int *pldy) {
    //
    // ..
    // .. local variables ..
    int i, j, d, ioff, joff, nrow, ncol, ldx, ldy;
    //
    // ..
    // .. executable statements ..
    nrow = *pnrow;
    ncol = *pncol;
    ldx  = *pldx;
    ldy  = *pldy;
    for( j=0; j<ncol; j++ ) {
	for( i=0; i<nrow; i++ ) {
	    d = (i+1);
	    yout[j*ldy+i].r = d*xin[j*ldx+i].r;
	    yout[j*ldy+i].i = d*xin[j*ldx+i].i;
	}
    }
}
//
void rankk_op(int *pnrow, int *pncol, trl_dcomplex *xin, int *pldx, trl_dcomplex *yout, int *pldy ) {
    int i, j, k, nrow, ncol, ldx, ldy;
    nrow = *pnrow;
    ncol = *pncol;
    ldx  = *pldx;
    ldy  = *pldy;

    for( k=0; k<ncol; k++ ) {
	for( i=0; i<rank_k; i++ ) {
	    rank_t[i].r = 0.0;
	    rank_t[i].i = 0.0;
	    // V^H * xin
	    for( j=0; j<nrow; j++ ) {
		rank_t[i].r += (rank_v[i*nrow+j].r * xin[k*ldx+j].r + rank_v[i*nrow+j].i * xin[k*ldx+j].i );
		rank_t[i].i += (rank_v[i*nrow+j].r * xin[k*ldx+j].i - rank_v[i*nrow+j].i * xin[k*ldx+j].r );
	    }
	    // diag(alpha)^2 * (V^H*xin)
	    rank_t[i].r *= rank_alpha[i];
	    rank_t[i].i *= rank_alpha[i];
	}
	for( i=0; i<nrow; i++ ) {
	    yout[k*ldy+i].r = 0.0;
	    yout[k*ldy+i].i = 0.0;
	}
	for( i=0; i<rank_k; i++ ) {
	    // V * (diag(alpha)^2 * V^H *xin)
	    for( j=0; j<nrow; j++ ) {
		yout[k*ldy+j].r += ( rank_v[i*nrow+j].r * rank_t[i].r - rank_v[i*nrow+j].i * rank_t[i].i );
		yout[k*ldy+j].i += ( rank_v[i*nrow+j].r * rank_t[i].i + rank_v[i*nrow+j].i * rank_t[i].r );
	    }
	}
    }
}
#else
// a simple matrix-vector multiplications routine
// defines a diagonal matrix with values (1, 4, 9, 16, 25, ....)
//
void zdiag_op(const int nrow, const int ncol,
	      const trl_dcomplex *xin, const int ldx,
	      trl_dcomplex *yout, const int ldy, void *mvparam) {
    //
    // ..
    // .. local variables ..
    int i, j;
    double d;
    //
    // ..
    // .. executable statements ..
    for( j=0; j<ncol; j++ ) {
	for( i=0; i<nrow; i++ ) {
	    d = (i+1.0)*(i+1.0);
	    yout[j*ldy+i].r = d*xin[j*ldx+i].r;
	    yout[j*ldy+i].i = d*xin[j*ldx+i].i;
	}
    }
}
//
void rankk_op(const int nrow, const int ncol,
	      const trl_dcomplex *xin, const int ldx,
	      trl_dcomplex *yout, const int ldy, void *mvparam) {
    int i, j, k;

    for( k=0; k<ncol; k++ ) {
	for( i=0; i<rank_k; i++ ) {
	    rank_t[i].r = 0.0;
	    rank_t[i].i = 0.0;
	    // V^H * xin
	    for( j=0; j<nrow; j++ ) {
		rank_t[i].r += (rank_v[i*nrow+j].r * xin[k*ldx+j].r + rank_v[i*nrow+j].i * xin[k*ldx+j].i );
		rank_t[i].i += (rank_v[i*nrow+j].r * xin[k*ldx+j].i - rank_v[i*nrow+j].i * xin[k*ldx+j].r );
	    }
	    // diag(alpha)^2 * (V^H*xin)
	    rank_t[i].r *= rank_alpha[i];
	    rank_t[i].i *= rank_alpha[i];
	}
	for( i=0; i<nrow; i++ ) {
	    yout[k*ldy+i].r = 0.0;
	    yout[k*ldy+i].i = 0.0;
	}
	for( i=0; i<rank_k; i++ ) {
	    // V * (diag(alpha)^2 * V^H *xin)
	    for( j=0; j<nrow; j++ ) {
		yout[k*ldy+j].r += ( rank_v[i*nrow+j].r * rank_t[i].r - rank_v[i*nrow+j].i * rank_t[i].i );
		yout[k*ldy+j].i += ( rank_v[i*nrow+j].r * rank_t[i].i + rank_v[i*nrow+j].i * rank_t[i].r );
	    }
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
    static const int nrow=1000, mev=100, c__1=1;
    // local variable declaration
    int lohi, ned, maxlan, restart, lwrk, ldwrk;
    double eval[mev], exact[mev];
    double *res, *dwrk, d__1;
    trl_dcomplex evec[mev*nrow];
    trl_dcomplex *wrk, z__1;
    trl_info info;
    int i, j, k, check, nbas, ii, nmis;
    char name2[133], name[150], file[150];
    int tmp1, tmp2, nlen;
    FILE *ferr;
    // initialize info -- tell TRLAN to compute NEV smallest eigenvalues
    // of the diag_op

    k=10;
    //rankk_init_( nrow, k );
    ferr=fopen( "error.txt","w" );
    fclose(ferr);
    for( lohi=-1; lohi<=-1; lohi++ ) {
	for( ned=100; ned<=100; ned+=10 ) {
	    for( maxlan=200; maxlan<=200; maxlan+=100 ) {
		//nbas = max(1, maxlan-mev+1);
		//ii = nbas * ((nrow+3)/4)*4;;
		//nmis = maxlan*(maxlan+10);
		//lwrk=ii;
		//lwrk=0;
		lwrk  = 10*nrow;
		wrk   = (trl_dcomplex*)malloc( 2*lwrk*sizeof(trl_dcomplex));
		ldwrk = 10*maxlan;
		dwrk  = (double*)malloc( ldwrk*sizeof(double) );
		for( restart=7; restart<=7; restart++ ) {
		    ferr=fopen( "error.txt","a" );
		    fprintf( ferr," \n ~~ iteration (ned=%d, maxlan=%d, restart=%d lohi=%d) ~~\n",ned,maxlan,restart,lohi );
		    printf( " \r\n ~~ iteration (ned=%d, maxlan=%d, restart=%d lohi=%d) ~~\r\n",ned,maxlan,restart,lohi );
		    fclose(ferr);

		    trl_init_info( &info, nrow, maxlan, lohi, ned, 1.4901/100000000, restart, 500000, -1 );
		    trl_set_iguess( &info, 0, 1, 0, NULL );
		    //trl_set_iguess( &info, 0, 2, 1, "CHECK" );
		    //trl_set_checkpoint( &info, 1, "CHECK" );
		    // the Lanczos recurrence is set to start with [1,1,...,1]^T
		    //info.rfact = 1.5;
		    trl_set_restart( &info, 1.5 );
		    memset( eval, 0, mev*sizeof(double) );
		    for( i=0; i<nrow; i++ ) {
			evec[i].r = 1.0;
			evec[i].i = 0.0;
			// to test pertubation set evec to be the first eigenvector
			//evec[i].r = 0.0;
		    }
		    evec[0].r = 1.0;
		    if( sprintf( file,"LOG_Nrow%dNed%dLan%dRes%dLoHi%d",nrow,ned,maxlan,restart,lohi ) < 0 ) {
			printf( "error writing file name\r\n" );
		    } else {
			trl_set_debug( &info, 0, file );
			// call TRLAN to compute the eigenvalues
			ztrlan(zdiag_op, &info, nrow, mev, eval, evec, nrow, wrk, lwrk, dwrk, ldwrk );
			//ztrlan(rankk_op_, &info, nrow, mev, eval, evec, nrow, wrk, lwrk, dwrk, ldwrk );
		    }
		    trl_print_info(&info, 3*nrow );
		    i = info.nec;
		    res = (double*)malloc(i*sizeof(double));
		    for( j=0; j<i; j++ ) {
			res[j]=wrk[j].r;
		    }

		    // exact eigenvalues for rank-k matrix
		    if( lohi == -1 ) {
			ii = 0;
			for( j=0; j<i; j++ ) {
			    if( eval[j] < 0.5 ) {
				exact[j] = 0.0;
			    } else {
				ii++;
				exact[j] = ii*ii;
			    }
			}
		    } else if( lohi == 1 ) {
			for( j=0; j<i; j++ ) {
			    if( j < i-k ) {
				exact[j] = 0.0;
			    } else {
				exact[j] = (j-i+k+1)*(j-i+k+1);
			    }
			}
		    }

		    // exact eigenvalues for diagonal matrix
		    if( lohi == -1 ) {
			for( j=0; j<i; j++ ) {
			    exact[j] = (j+1)*(j+1);
			}
		    } else if( lohi == 1 ) {
			for( j=0; j<i; j++ ) {
			    exact[j] = (nrow-i+j+1)*(nrow-i+j+1);
			}
		    }

		    if( lohi != 0 ) {
			ztrl_check_ritz( zdiag_op, &info, nrow, i, evec, nrow, eval, &check, res, exact, wrk, lwrk );
			//ztrl_check_ritz( rankk_op, &info, nrow, i, evec, nrow, eval, &check, res, exact, wrk, lwrk );
		    } else {
			ztrl_check_ritz( zdiag_op, &info, nrow, i, evec, nrow, eval, &check, res, NULL, wrk, lwrk );
			//ztrl_check_ritz( rankk_op, &info, nrow, i, evec, nrow, eval, &check, res, NULL, wrk, lwrk );
		    }
		    //
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

		    /*
		    // checking eigenvectors for rank-k matrix.
		    // notice, the eigenvectors can be different by a constand complex with a unit-magnitude.
		    ferr = fopen( "error.txt","a" );
		    if( lohi == -1 ) {
			ii = -1;
			for( i=0; i<info.nec; i++ ) {
			    if( eval[i] > 0.5 ) {
				ii ++;
				d__1 = evec[i*nrow].r * evec[i*nrow].r + evec[i*nrow].i * evec[i*nrow].i;
				z__1.r = ( rank_v[ii*nrow].r * evec[i*nrow].r + rank_v[ii*nrow].i * evec[i*nrow].i ) / d__1;
				z__1.i = ( rank_v[ii*nrow].i * evec[i*nrow].r - rank_v[ii*nrow].r * evec[i*nrow].i ) / d__1;
				fprintf( ferr,"%d, %e: cons scal: %e + i * %e (%e) ",i+1, eval[i],z__1.r, z__1.i, sqrt(z__1.r*z__1.r + z__1.i*z__1.i) );
				zscal_( &nrow, &z__1, &(evec[i*nrow]), &c__1 );
				z__1.r = -1.0;
				z__1.i = 0.0;
				zaxpy_( &nrow, &z__1, &(rank_v[ii*nrow]), &c__1, &(evec[i*nrow]), &c__1 );
				//zdotc_( &z__1, &nrow, &(evec[i*nrow]), &c__1, &(evec[i*nrow]), &c__1 );
				fprintf( ferr,"diff evec: %16.15e\r\n",sqrt(z__1.r) );
			    } else {
				fprintf( ferr,"%d, %e: skipping check.\r\n",i+1,eval[i] );
			    }
			}
		    }
		    fclose(ferr);
		    */
		    free(res);
		}
		free(wrk);
		free(dwrk);
	    }
	}
    }
}
