/*
 * The subroutines included in this file provide the fortran
 * interface to the subroutines in trlan.h
 */
#ifndef __FTRLAN_H
#define __FTRLAN_H
#ifdef TRL_FORTRAN_COMPATIBLE
#include <stdio.h>
#include "trl_map.h"
#include "trlan.h"
#include "ztrlan.h"

void trlan_(trl_matvec op,
            trl_info * info, int nrow, int mev, double *eval,
            double *evec, int lde, int lwrk, double *wrk );

void ztrlan_(ztrl_matvec op,
	     trl_info * info, int nrow, int mev, double *eval,
	     trl_dcomplex * evec, int lde, trl_dcomplex * misc, int nmis,
	     double *dwrk, int ldwrk);

void trl_init_info_(trl_info * info, int nrow, int mxlan, int lohi,
                    int ned, double tol, int restart, int maxmv,
                    void *mpicomp);

void trl_set_debug_(trl_info * info, int msglvl, char *filename);

void trl_set_iguess_(trl_info * info, int nec, int iguess, int nopts,
                     char *cpf );

void trl_set_checkpoint_(trl_info * info, int cpflag, char *file);

void trl_set_restart_(trl_info * info, double rfact);

void trl_print_info_(trl_info * info, int mvflop);

void trl_terse_info_(trl_info * info, FILE * iou);

void trl_check_ritz_(trl_matvec op,
		     trl_info * info, int nrow, int ncol, double *rvec,
		     int ldrvec, double *alpha, int *check, double *beta,
		     double *eval, double *wrk, int lwrk);

void trl_rayleigh_quotients_(trl_matvec op, trl_info * info, int ncol,
			     double *evec, int lde, double *eres, double *base);

void ztrl_check_ritz_(ztrl_matvec op, trl_info * info, int nrow, int ncol,
		      trl_dcomplex * rvec, int ldrvec, double *alpha,
		      int *check, double *beta, double *eval,
		      trl_dcomplex * wrk, int lwrk);
#else
#error TRLan fortran interface not built.  Rebuild TRLan with TRL_FOTRAN_COMPATIBLE defined
#endif
#endif
