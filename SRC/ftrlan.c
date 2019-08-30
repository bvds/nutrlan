/*
 * The subroutines included in this file provide the fortran
 * interface to the subroutines in trlan.h
 */
#include <stdio.h>
#include "ftrlan.h"

void trlan_(trl_matvec op, trl_info * info, int nrow, int mev, double *eval,
	    double *evec, int lde, int lwrk, double *wrk ) 
{
  trlan( op, info, nrow, mev, eval, evec, lde, lwrk, wrk );
}
/*
 */
void
ztrlan_(ztrl_matvec op, trl_info * info, int nrow, int mev, double *eval,
        trl_dcomplex * evec, int lde, trl_dcomplex * misc, int nmis,
        double *dwrk, int ldwrk)
{
  ztrlan( op, info, nrow, mev, eval, evec, lde, misc, nmis, dwrk, ldwrk );
}
/*
 */
void trl_init_info_(trl_info * info, int nrow, int mxlan, int lohi,
                    int ned, double tol, int restart, int maxmv,
                    int mpicom)
{
  trl_init_info( info, nrow, mxlan, lohi, ned, tol, restart, maxmv, mpicom );
}
/*
 */
void trl_set_debug_(trl_info * info, int msglvl, char *filename)
{
  trl_set_debug( info, msglvl, filename);
}
/*
 */
void trl_set_iguess_(trl_info * info, int nec, int iguess, int nopts,
                     char *cpf )
{
  trl_set_iguess( info, nec, iguess, nopts, cpf );
}
/*
 */
void trl_set_checkpoint_(trl_info * info, int cpflag, char *file)
{
  trl_set_checkpoint( info, cpflag, file );
}
/*
 */
void trl_set_restart_(trl_info * info, double rfact)
{
  trl_set_restart( info, rfact);
}
/*
 */
void trl_print_info_(trl_info * info, int mvflop)
{
  trl_print_info( info, mvflop);
}
/*
 */
void trl_terse_info_(trl_info * info, FILE * iou)
{
  trl_terse_info( info, iou);
}
/*
 */
void    
trl_check_ritz_(trl_matvec op, trl_info * info, int nrow, int ncol,
                double *rvec, int ldrvec, double *alpha, int *check,
                double *beta, double *eval, double *wrk, int lwrk)
{
  trl_check_ritz( op, info, nrow, ncol, rvec, ldrvec, alpha, check, 
                   beta, eval, wrk, lwrk);
}
/*
 */
void
trl_rayleigh_quotients_(trl_matvec op, trl_info * info, int ncol, double *evec,
                        int lde, double *eres, double *base)
{
  trl_rayleigh_quotients( op, info, ncol, evec, lde, eres, base);
}
/*
 */
void ztrl_check_ritz_(ztrl_matvec op, trl_info * info, int nrow, int ncol,
		      trl_dcomplex * rvec, int ldrvec, double *alpha,
		      int *check, double *beta, double *eval, 
		      trl_dcomplex * wrk, int lwrk)
{
  ztrl_check_ritz( op, info, nrow, ncol, rvec, ldrvec, alpha, check, 
                   beta, eval, wrk, lwrk );
}

