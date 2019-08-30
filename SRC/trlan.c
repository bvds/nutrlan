/*
// TRLan routine
// Lawrence Berkeley National Lab.
//
*/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdarg.h>

#include "trl_map.h"
#include "trlan_i.h"
#include "trlcore_i.h"
#include "trlaux_i.h"
#include "trl_comm_i.h"

void trlan(trl_matvec op, trl_info * info, int nrow, int mev, double *eval,
	   double *evec, int lde, int lwrk, double *wrk )
{
/*
// Purpose: Top (user) level routines
// ========
// A thick-restart Lanczos routine for computing eigenvalues and
// eigenvectors of a real symmetric operator/matrix (A).
// -- only accept one input vector, the input vector is expected
//    to be stored in the (nec+1)st column of EVEC.
// -- it extends the Lanczos basis one vector at a time.
// -- orthogonality among the Lanczos vectors are maintained using
//    full re-orthogonalization when necessary.
//
// Requirements:
// =============
// 1) User supplies OP with the specified interface.
// 2) If (info->nec>0), evec(1:nrow, 1:info->nec) must contain valid
//    eigenvectors and eval(1:nec) must be the corresponding eigenvalues.
//    These eigenpairs are assumed to have zero residual norms.
// 3) lde >= nrow.
// 4) The arrays evec and eval must be large enough to store the
//    solutions, i.e., mev >= info->ned and mev >= info->nec.
// 5) The array wrk may be of arbitrary size.  Internally, the workspace
//    size is
//        nrow*max(0,info->ned-size(evec,2))+maxlan*(maxlan+10)
//
//    If wrk is smaller than this, trlan routine will allocate additional
//    workspace to accommodate.
//
// Arguments:
// ==========
// op      (input) function pointer
//         On entry, points to a function that comptues op(X) == A*X,
//         when given a set of vectors X.
//         The operator that defines the eigenvalue problem is expected to
//         have the following interface
//          void op(nrow, ncol, xin, ldx, yout, ldy)
//            nrow  (input) integer
//                   On entry, specifies the number of rows in xin and yout.
//            ncol  (input) integer
//                   On entry, specifies, the number of columns in Xin and
//                   Yout.
//            xin   (input) double precision vector of length (ldx*ncol)
//                   On entry, contatins the input vector to be multiplied.
//                   It consists of ncol column vectors with each column
//                   stored in consecutive order.
//            ldx   (input) integer
//                   On entry, specifies the leading dimension of the array
//                   Xin, i.e., the i-th column vector starts with element
//                   (i-1)*ldx+1 and ends with element (i-1)*ldx+nrow in Xin.
//            yout  (output) double precision vector of length (ldy*ncol)
//                   On exit, contains the result array, i.e., it stores the
//                   result of matrix-vector multiplications.
//            ldy   (input) integer
//                   On entry, specifies the leading dimension of the array
//                   yout, i.e., the i-th column vector starts with element
//                   (i-1)*ldy+1 in Yout and ends with element (i-1)*ldy+nrow.
//
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the information
//          about the eigenvalue problem and the progress of TRLAN
//
// nrow    (input) integer
//          On entry, specifies the number of rows that is on this processor.
//
// mev     (input) integer
//          On entry, specifies the number of Ritz pairs, that can be stored in
//          eval and evec.
//
// eval    (output) double precision vector of length (mev)
//          On exist, stores the eigenvalues.
//
// evec    (output) double precision vector of lenvth (nrow*mev)
//          On exit, stores the eigenvectors.
//
// lde     (input) integer
//          On entry, specifies the leading dimension of the array evec, i.e.,
//
// lwrk    (optional) integer
//          On entry, specifies, the size of WRK.  When both WRK and LWRK are
//          present, then LWRK should correctly indicate the size of WRK. If WRK
//          is present by not LWRK, the size of WRK is assumed to be MEV which is
//          only enough to store the residual norms on exit.  If WRK is not
//          present, LWRK is not used even if it is present.
//          (lde >= nrow).
//
// wrk     (optional) workspace
//          If it is provided and there is enough space, the residual norm of
//          the converged eigenpairs will be stored at wrk(1:info->nec) on exit.
//
*/
/*
// ..
// .. local scalars ..
*/
    clock_t clk1;
    int ii, nbas, nmis, ibas, imis, ldb, lwrk0;
    double *base, *misc;
/*
// ..
// .. executables statements ..
*/
    imis = -1;			/* if this routine allocated misc, imis will be 0 */
    ibas = -1;			/* if this routine allocated base, ibas will be 0 */
    clk1 = clock();
    info->clk_tot = clk1;
    if (info->ned > mev) {
	printf
	    ("info->ned (%d) is larger than mev (%d) reducing info->ned to %d\n",
	     info->ned, mev, mev);
	info->ned = mev;
    }
    /* there is nothing to do if there is no more eigenvalue to compute */
    if (info->ned <= info->nec || info->ned <= 0)
	goto end;
    /* determine the workspace size and leading dimensions
    if (nopts > 0) {
	va_list argptr;
	va_start(argptr, nopts);
	lwrk = va_arg(argptr, int);
	wrk = va_arg(argptr, double *);
	va_end(argptr);
    } else {
	lwrk = 0;
	wrk = NULL;
    }
    */
    lwrk0 = lwrk;
    info->stat = 0;
    ldb = ((nrow + 3) / 4) * 4;
    if ((ldb % 4096) == 0)
	ldb = ldb + 4;
    trl_clear_counter(info, nrow, mev, lde);
    if (info->stat != 0)
	goto end;
    /*
       // Internally, the workspace is broken into two parts
       // one to store (maxlan+1) Lanczos vectors, and the other to
       // store all others (size maxlan*(maxlan+ned+14))
       // The next If-block decides how the two arrays are mapped.
     */
    nbas = max(1, info->maxlan - mev + 1);
    ii = nbas * ldb;
    nmis = info->maxlan * (info->maxlan + 10);
    if (lwrk0 >= min(ii, nmis)) {
	/* use wrk either as base or misc or both depending its size    */
	if (lwrk0 >= ii + nmis) {
	    /* WRK is large enough for both arrays */
	    base = wrk;
	    misc = &wrk[ii];
	    nmis = lwrk0 - ii;
	    //printf( "\n\n ******** Large enough workspace ********* \n\n" );
	} else if (lwrk0 >= max(ii, nmis)) {
	    /* associate the larger one of base and misc to WRK */
	    if (ii >= nmis) {
		base = wrk;
		misc = (double *) malloc(nmis * sizeof(double));
		if (misc == NULL)
		    info->stat = -4;
		imis = 0;
	    } else {
		misc = wrk;
		nmis = lwrk0;
		base = (double *) malloc(ii * sizeof(double));
		if (base == NULL)
		    info->stat = -5;
		ibas = 0;
	    }
	} else if (ii <= nmis) {
	    /* base is smaller, associate base with WRK */
	    base = wrk;
	    misc = (double *) malloc(nmis * sizeof(double));
	    if (misc == NULL)
		info->stat = -4;
	    imis = 0;
	} else {
	    /* misc is smaller, associate misc with WRK */
	    misc = wrk;
	    nmis = lwrk0;
	    base = (double *) malloc(ii * sizeof(double));
	    if (base == NULL)
		info->stat = -5;
	    ibas = 0;
	}
    } else {
	/* have to allocate both base and misc */
	base = (double *) malloc(ii * sizeof(double));
	if (base == NULL)
	    info->stat = -5;
	ibas = 0;
	misc = (double *) malloc(nmis * sizeof(double));
	if (misc == NULL)
	    info->stat = -4;
	imis = 0;
    }
    memset(base, 0, ii * sizeof(double));
    memset(misc, 0, nmis * sizeof(double));
    /* make sure every process is successful so far */
    ii = trl_sync_flag(info->mpicom, info->stat);
    info->stat = ii;
    if (ii != 0)
	goto end;
    /* open log and checkpoint files */
    trl_open_logfile(info);
    /* trl_open_cptfile(info); */
    if (info->verbose > 0) {
	trl_time_stamp(info->log_fp);
	trl_print_setup(info, nbas * ldb, nmis, lwrk0);
    }

    /* call trlanczos to do the real work  */
    //printf( "calling trlanczso (%d)\n",info->cpflag );
    trlanczos(op, info, nrow, mev, eval, evec, lde, base, ldb, nbas,
	      misc, nmis);
#ifdef DEBUG
    printf( "DEBUG ** out of trlanczos (locked=%d) **\n",info->locked );
#endif

    /* close log and checkpoint files */
    trl_close_logfile(info);
    if (lwrk0 >= mev) {
	memcpy(wrk, misc, mev * sizeof(double));
    } else {
	memcpy(wrk, misc, lwrk0 * sizeof(double));
    }

    /* DONE, reclaim the space allocated */
  end:
    if (imis == 0)
	free(misc);
    if (ibas == 0)
	free(base);
    clk1 = clock();
    if (clk1 < info->clk_tot) {
        info->tick_t +=
            (info->clk_max -
             info->clk_tot) / (double) (info->clk_rate);
        info->tick_t +=
            (info->clk_max + clk1) / (double) (info->clk_rate);
        info->clk_tot = 0;
    } else if (info->clk_tot < 0 && clk1 >= 0) {
        info->tick_t -= info->clk_tot / (double) (info->clk_rate);
        info->tick_t += clk1 / (double) (info->clk_rate);
        info->clk_tot = 0;
    } else {
        info->tick_t  += (clk1 - info->clk_tot) / (double) (info->clk_rate);
        info->clk_tot  = 0;
    }
/*
    if (clk1 >= info->clk_tot) {
	//printf( "clk_tot = %d-%d = ",ii,info->clk_tot );
	info->clk_tot = ii - info->clk_tot;
	//printf( "%d\n",info->clk_tot );
    } else if (clk1 < 0) {
	//assuming ( info->clk_tot > 0 ), otherwise ?? wrap-around twice ??
	info->tick_t +=
	    (info->clk_max - info->clk_tot) / (double) (info->clk_rate);
	info->clk_tot = info->clk_max + ii;
    } else {
	//assuming ( info->clk_tot < 0 ),
	info->tick_t +=
	    (info->clk_max + info->clk_tot) / (double) (info->clk_rate);
	info->clk_tot = ii;
    }
*/
    return;
/*
// .. end of trlan_ ..
*/
}

//
void trl_set_restart(trl_info * info, double rfact)
{
//
// Purpose
// =======
// Set the (minimum) basis size for the restart 7 and 8 schemes, i.e., the (minimum) basis
// size if rfact * (number of Ritz vectors kept)
//
// Arguments:
// ==========
// info    (input/output) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
// rfact   (input) double precision
//          On entry, specify the (minimum) basis size.
//
    info->rfact = rfact;
}

//
void trl_set_debug(trl_info * info, int msglvl, char *filename)
{
/*
// Purpose:
// ========
// Set information related to debugging, the initialization routine trl_init_info sets the
// parameters so that no debug information is printed.
//
// Arguments:
// ==========
// info    (input/output) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
// msglvl  (input) integer
//          On entry, specifies the new message/verbose level:
//           msglvl <  0         : nothing is printed
//           msglvl = 1, .., 10  : the higher the level, the more debug information is
//                                 printed.
//
// file    (input) character string
//          On entry, specifies the leading part of the debug file name. Each processor will
//          generate its own debug file with the file name formed by appending the processor
//          number to the string FILE. The actual name is composed by the routine
//          TRL_PE_FILENAME in trlaux
//
// ..
// .. executable statements ..
*/
    info->verbose = msglvl;
    if (filename != NULL) {
	strcpy(info->log_file, filename);
	if (msglvl >= 0 && info->my_pe == 0) {
	    printf("TRLan will write diagnostic messages to files with "
		   "prefix %s.\n", info->log_file);
	}
    }
/*
// .. end of trl_set_debug_ ..
*/
}

void trl_set_checkpoint(trl_info * info, int cpflag, char *file)
{
/*
// Purpose:
// ========
// Set up the information related to check-pointing
//
// Arguments:
// ==========
// info    (input/output) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
// cpflag  (input) integer
//          On entry, spcifies roughly how many titmes checkpoints are written.
//
// file    (input) character string
//          On entry, specifies the leading part of the checkpoint file name. Each processor
//          will generate its own debug file with the file name formed by appending the
//          processor number to the string FILE. The actual name is composed by the routine
//          TRL_PE_FILENAME in trlaux
//
// ..
// .. executable statements ..
*/
    info->cpflag = cpflag;
    if (file != NULL) {
	strcpy(info->cpfile, file);
	if (info->verbose >= 0 && info->my_pe == 0) {
	    printf("TRLan will write checkpoint to files with prefix %s.\n",
		   info->cpfile);
	}
    }
/*
// .. end of trl_set_checkpoint_ ..
*/
}

void trl_set_iguess(trl_info * info, int nec, int iguess, int nopts,
		     char *cpf )
{
/*
// Purpose:
// ========
// Set up parameters related to initial guesses of the Lanczos iterations, i.e., the number of
// eigenvector already converged (initially assumed to be zero) and whether the user has
// provided initial guess vector (initially assumed no).  It can also tell TRLan to read check
// point file that is different from the default name.
//
// Arguments:
// ==========
// info    (input/output) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
// nec     (input) integer
//          On entry, specifies the number of eigenvalues, that have been converged.
//
// iguess  (input) integer
//          On entry, specifies one of the following options:
//           iguess <= 0, user did not provide initial guess, use [1,1,..,1].
//           iguess =  1, user has supplied initial guess, will only use the first one.
//           iguess >  1, restart with previous check-point file.
//
// nopts   (optional)
//          If provided, it specifies the name of checkpoints file.
//
// ..
// .. executable statements ..
*/
    /*
    char cpf[STRING_LEN];
    //printf( "nec=%d iguess=%d nopts=%d\n",nec,iguess,nopts );
    va_list argptr;
    va_start(argptr, nopts);
    if (nopts > 0) {
	strcpy(info->oldcpf, va_arg(argptr, char *));
    } else {
	strcpy(info->oldcpf, "");
    }
    va_end(argptr);
    */
    /* assign nec and iguess flags to info */
    info->nec = nec;
    info->guess = iguess;
    if (strlen(info->oldcpf) > 0 && info->guess > 1) {
	/* check to make sure the files exist */
	trl_pe_filename(STRING_LEN, cpf, info->oldcpf, info->my_pe,
			 info->npes);
	if ((info->cpt_fp = fopen(cpf, "r")) != NULL) {
	    if (fclose(info->cpt_fp) != 0) {
		info->stat = -9;
	    }
	} else {
	    info->stat = -8;
	}
	info->stat = trl_sync_flag(info->mpicom, info->stat);
    } else {
	info->stat = 0;
    }
/*
// .. end of trl_set_iguess_ ..
*/
}

void trl_print_info(trl_info * info, int mvflop)
{
/*
// Purpose:
// ========
// Provides an uniform way of printing information stored in TRL_INFO_T.  It needs to be
// called by all PEs. In parallel environment, when writing to standard outputd device, only
// PE0 will actually write its local summary information. Note that this function must be
// called on all PEs at the same time ***
//
// Arguments:
// ==========
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
// mvflop  (input) integer
//          On entry, specifies the number of floating-operations per MATVEC per PE.  This
//          information has to be supplied by user, otherwise related entries are left blank
//          in the print out.
//
// ..
// .. local arrays ..
*/
    double tmp1[12], tmp2[12];
/*
// ..
// .. local variables ..
*/
    int i;
    double t_tot, t_op, t_orth, t_res, t_in, t_out, rinv, r_tot, r_op,
	r_orth, r_res, r_in, r_out;
/*
// ..
// .. executable statements ..
*/
    if (info->clk_rate > 0) {
	rinv = 1.0 / (double) (info->clk_rate);
    } else {
	/* get clock rate */
	rinv = 1.0 / CLOCKS_PER_SEC;
    }
    t_op   = info->tick_o + info->clk_op   * rinv;
    t_tot  = info->tick_t + info->clk_tot  * rinv;
    t_res  = info->tick_r + info->clk_res  * rinv;
    t_orth = info->tick_h + info->clk_orth * rinv;
    t_in   = info->clk_in * rinv;
    t_out  = info->clk_out * rinv;
    //printf( "tick_t=%e clk_tot=%d\n",info->tick_t,info->clk_tot );
    if (t_op != 0 && mvflop != 0) {
	if (mvflop > 0) {
	    r_op = ((double)mvflop) * ((double)info->matvec);
	} else {
	    r_op = 0;
	}
    } else {
	r_op = 0;
    }
    if (t_orth != 0) {
	r_orth = info->rflp_h + info->flop_h;
    } else {
	r_orth = 0;
    }
    if (t_res != 0) {
	r_res = info->rflp_r + info->flop_r;
    } else {
	r_res = 0;
    }
    if (r_op > 0) {
	r_tot =
	    ((double) mvflop) * ((double) info->matvec) + info->rflp +
	    ((double) info->flop);
    } else {
	r_tot = 0;
    }
    if (info->clk_in > 0) {
	r_in = 8.0 * info->wrds_in;
    } else {
	r_in = 0;
    }
    if (info->clk_out > 0) {
	r_out = 8.0 * info->wrds_out;
    } else {
	r_out = 0;
    }
    tmp2[0] = t_tot;
    tmp2[1] = t_op;
    tmp2[2] = t_orth;
    tmp2[3] = t_res;
    tmp2[4] = t_in;
    tmp2[5] = t_out;
    tmp2[6] = r_tot;
    tmp2[7] = r_op;
    tmp2[8] = r_orth;
    tmp2[9] = r_res;
    tmp2[10] = r_in;
    tmp2[11] = r_out;
    //printf( "print info\n" );
    //printf( "calling g_sum\n" );
    trl_g_sum(info->mpicom, 12, tmp2, tmp1);
    if (info->log_fp == NULL) {
	trl_reopen_logfile(info);
    }
    trl_time_stamp(info->log_fp);
    //printf( "printing\n" );
    if (info->npes > 1) {
	fprintf(info->log_fp,
		"TRLAN execution summary (exit status = %d) on PE %d\n",
		info->stat, info->my_pe);
    } else {
	fprintf(info->log_fp,
		"TRLAN execution summary (exit status =%d)\n", info->stat);
    }
    if (info->lohi > 0) {
	fprintf(info->log_fp,
		"Number of LARGEST eigenpairs      %10d (computed) %11d (wanted)\n",
		info->nec, info->ned);
    } else if (info->lohi < 0) {
	fprintf(info->log_fp,
		"Number of SMALLEST eigenpairs    %10d (computed) %11d (wanted)\n",
		info->nec, info->ned);
    } else {
	fprintf(info->log_fp,
		"Number of EXTREME eigenpairs     %10d (computed) %11d (wanted)\n",
		info->nec, info->ned);
    }
    fprintf(info->log_fp,
	    "Times the operator is applied:   %10d (MAX: %16d )\n",
	    info->matvec, info->maxmv);
    fprintf(info->log_fp,
	    "Problem size:                    %10d (PE: %4d) %11d (Global)\n",
	    info->nloc, info->my_pe, info->ntot);
    fprintf(info->log_fp,
	    "Convergence tolerance:           %10.3e (rel) %16.3e (abs)\n",
	    info->tol, info->tol * info->anrm);
    fprintf(info->log_fp, "Maximum basis size:              %10d\n",
	    info->maxlan);
    fprintf(info->log_fp, "Restarting scheme:               %10d\n",
	    info->restart);
    fprintf(info->log_fp, "Number of re-orthogonalizations: %10d\n",
	    info->north);
    fprintf(info->log_fp, "Number of (re)start loops:       %10d\n",
	    info->nloop);
    if (info->nrand > 0) {
	fprintf(info->log_fp, "Number of random vectors used:   %10d\n",
		info->nrand);
    }
    if (info->npes > 1) {
	fprintf(info->log_fp, "Number of MPI processes:         %10d\n",
		info->npes);
    }
    fprintf(info->log_fp, "Number of eigenpairs locked:     %10d\n",
	    info->locked);
    if (t_op > 0) {
	fprintf(info->log_fp,
		"OP(MATVEC):            %12.4e sec, %12.4e FLOP/S (%11.4e FLOP)\n",
		t_op, r_op / t_op, r_op);
    } else {
	fprintf(info->log_fp, "time in OP:            %12.4e sec\n", t_op);
    }
    if (t_orth > 0) {
	fprintf(info->log_fp,
		"Re-Orthogonalization:: %12.4e sec, %12.4e FLOP/S (%11.4e FLOP)\n",
		t_orth, r_orth / t_orth, r_orth);
    } else {
	fprintf(info->log_fp, "time in orth:          %12.4e sec\n",
		t_orth);
    }
    if (t_res > 0) {
	fprintf(info->log_fp,
		"Restarting::           %12.4e sec, %12.4e FLOP/S (%11.4e FLOP)\n",
		t_res, r_res / t_res, r_res);
    } else {
	fprintf(info->log_fp, "time in restarting:    %12.4e sec\n",
		t_res);
    }
    if (t_tot > 0) {
	fprintf(info->log_fp,
		"TRLAN on this PE:      %12.4e sec, %12.4e FLOP/S (%11.4e FLOP)\n",
		t_tot, r_tot / t_tot, r_tot);
    } else {
	fprintf(info->log_fp, "total time in TRLAN:   %12.4e sec\n",
		t_tot);
    }
    /*
       //if( info->verbose > 0 && info->log_fp != info->log_fp) {
       //  fprintf( info->log_fp, "Debug infomation written to files %s ####\n",info->log_file );
       //}
     */
    if (info->guess > 1 && info->wrds_in > 0) {
	if (strlen(info->oldcpf) <= 0) {
	    fprintf(info->log_fp,
		    "TRLAN restarted with checkpoint files %s ####\n",
		    info->oldcpf);
	} else {
	    fprintf(info->log_fp,
		    "TRLAN restarted with checkpoint files %s ####\n",
		    info->cpfile);
	}
	fprintf(info->log_fp,
		"Bytes read   %12.5e, Time(sec): %12.5e, Rate(B/s): %12.5e\n",
		r_in, t_in, r_in / t_in);
    }
    if (info->clk_out > 0 && info->wrds_out > 0) {
	fprintf(info->log_fp, "Checkpoint files are %s ####\n",
		info->cpfile);
	fprintf(info->log_fp,
		"Bytes read   %12.5e, Time(sec): %12.5e, Rate(B/s): %12.5e\n",
		r_out, t_out, r_out / t_out);
    }
    if (info->npes > 1) {
	/* write global performance information */
	rinv = 1.0 / info->npes;
	for (i = 0; i < 12; i++) {
	    tmp1[i] = tmp1[i] * rinv;
	}
	for (i = 0; i < 6; i++) {
	    if (tmp1[i] > 0) {
		tmp1[i + 6] = tmp1[i + 6] / tmp1[i];
	    } else {
		tmp1[i + 6] = 0.0;
	    }
	}
	if (tmp1[4] == tmp1[5] && tmp1[4] == 0) {
	    fprintf(info->log_fp,
		    " -- Global summary -- \n" );
            fprintf(info->log_fp,
                    "                       Overall,\t\t  MATVEC,\t  Re-orth,\t  Restart,\n");
	    fprintf(info->log_fp,
		    "Time(ave)             %11.4e,\t %11.4e,\t %11.4e,\t %11.4e\n",
		    tmp1[0], tmp1[1], tmp1[2], tmp1[3]);
	    fprintf(info->log_fp,
		    "Rate(tot)             %11.4e,\t %11.4e,\t %11.4e,\t %11.4e\n",
		    tmp1[6], tmp1[7], tmp1[8], tmp1[9]);
	} else {
	    fprintf(info->log_fp,
		    " -- Global summary -- \n" );
            fprintf(info->log_fp,
                    "                       Overall,\t\t  MATVEC,\t  Re-orth,\t  Restart,\t  Read,\t  Write\n");
	    fprintf(info->log_fp,
		    "Time(ave)             %11.4e,\t %11.4e,\t %11.4e,\t %11.4e,\t %11.4e,\t %11.4e\n",
		    tmp1[0], tmp1[1], tmp1[2], tmp1[3], tmp1[4], tmp1[5]);
	    fprintf(info->log_fp,
		    "Rate(tot)             %11.4e,\t %11.4e,\t %11.4e,\t %11.4e,\t %11.4e,\t %11.4e\n",
		    tmp1[6], tmp1[7], tmp1[8], tmp1[9], tmp1[10],
		    tmp1[11]);
	}
    }
    trl_close_logfile(info);
    return;
/*
// .. end of trl_print_info ..
*/
}

void trl_terse_info(trl_info * info, FILE * iou)
{
/*
// Purpose:
// ========
// It is a more compact version of trl_print_info, i.e., this is a local routine, indivadual
// PE can call it without regard of whether other PEs do the same and the output may be
// written to a different I/O unit number than the log_io
//
// Arguments:
// ==========
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
// iou     (input) pointer to a file
//          On entry, points to the file, here the information is outputed.
//
// ..
// .. local scalars ..
*/
    int rate;
    double t_tot, t_op, t_orth, t_res;
/*
// ..
// .. executable statements ..
*/
    if (iou == NULL) {
	if (info->log_fp == NULL) {
	    iou = stdout;
	} else {
	    iou = info->log_fp;
	}
    }
    if (info->clk_rate > 0) {
	t_op = info->tick_o + info->clk_op / (double) (info->clk_rate);
	t_tot = info->tick_t + info->clk_tot / (double) (info->clk_rate);
	t_res = info->tick_r + info->clk_res / (double) (info->clk_rate);
	t_orth = info->tick_h + info->clk_orth / (double) (info->clk_rate);
    } else {
	/* get clock rate */
	rate = CLOCKS_PER_SEC;
	t_op = info->tick_o + info->clk_op / (double) (rate);
	t_tot = info->tick_t + info->clk_tot / (double) (rate);
	t_res = info->tick_r + info->clk_res / (double) (rate);
	t_orth = info->tick_h + info->clk_orth / (double) (rate);
    }
    if (info->lohi > 0) {
	fprintf(iou,
		"MAXLAN:%10d, Restart:%10d,   NED: + %7d,      NEC:%10d\n",
		info->maxlan, info->restart, info->ned, info->nec);
    } else if (info->lohi < 0) {
	fprintf(iou,
		"MAXLAN:%10d, Restart:%10d,   NED: - %7d,      NEC:%10d\n",
		info->maxlan, info->restart, info->ned, info->nec);
    } else {
	fprintf(iou,
		"MAXLAN:%10d, Restart:%10d,   NED: 0 %7d,      NEC:%10d\n",
		info->maxlan, info->restart, info->ned, info->nec);
    }
    fprintf(iou,
	    "MATVEC:%10d,  Reorth:%10d, Nloop:   %7d,  Nlocked:%10d\n",
	    info->matvec, info->north, info->nloop, info->locked);
    if (t_tot > 0.001 && max(t_tot, max(t_op, max(t_res, t_orth))) < 1000) {
	fprintf(iou,
		"Ttotal:%10.6f,    T_op:%10.6f, Torth:%10.6f,   Tstart:%10.6f\n",
		t_tot, t_op, t_orth, t_res);
    } else {
	fprintf(iou,
		"Ttotal:%10.3e,    T_op:%10.3e, Torth:%10.3e,   Tstart:%10.3e\n",
		t_tot, t_op, t_orth, t_res);
    }
/*
// .. end of trl_terse_info_ ..
*/
}

void
trl_check_ritz(trl_matvec op, trl_info * info, int nrow, int ncol, double *rvec,
	       int ldrvec, double *alpha, int *check, double *beta,
	       double *eval, double *wrk, int lwrk)
{
/*
// Purpose:
// ========
// Performs a standard check on the computed Ritz pairs.
//
// Arguments:
// ==========
// op       (input) function pointer
//           On entry, points to the matrix-vector multiplication routine.
//
// info     (input) pointer to the structure trl_info_
//           On entry, points to the data structure to store the information
//           about the eigenvalue problem and the progress of TRLAN
//
// nrow     (input) integer
//           On entry, specifies the problem size.
//
// ncol     (input) integer
//           On entry, specifies the number of Ritz values computed.
//
// rvec     (input) double precision array of dimension (nrow,ncol)
//           On entry, specifies the array storing the Ritz vectors.
//
// alpha    (input) double precision array of dimension (ncol)
//           On entry, contains he Ritz values computed.
//
// beta     (optional) double precision array of dimension (ncol)
//           If provided, contaions the residual norms returned from a Lanczos routine.
//
// eval     (optional) double precision array of dimension (ncol)
//           If provided, contains the actual eigenvalues.
//
// lwrk     (optional) integer
//           If provided, specifies the size of workspace provided.
//
// wrk      (optional) double precision array of size(lwrk)
//           If provided, double precidion workspace.
//
// ..
// .. local parameters ..
*/
    long c__1 = 1;
    int i__1 = 1;
    double d__1;
/*
// ..
// .. local variables ..
// aq -- store the result of matrix-vector multiplication, size nrow
// rq -- store the Rayleigh-Quotient and the residual norms
// gsumwrk -- workspace left over for trl_g_sum to use dimension of the input arrays
*/
    double *aq, *rq, *gsumwrk, *res, *err;
    FILE *fp;
    int i, aqi, rqi, gsumwrki, icheck;
    double gapl, gapr;
/*
// ..
// .. executable statements ..
*/
    if (ncol <= 0)
	return;			/* there is nothing to do */

/* figure out whether it is necessary to allocate new workspace */
    *check = 0;
    aqi = 0;
    rqi = 0;
    gsumwrki = 0;
    if (lwrk > nrow + (4 * ncol)) {
	aq = &wrk[0];
	rq = &wrk[nrow];
	gsumwrk = &wrk[nrow + (3 * ncol)];
    } else if (lwrk >= (nrow + ncol)) {
	aq = &wrk[0];
	gsumwrk = &wrk[nrow];
	rq = (double *) malloc(3 * ncol * sizeof(double));
	if (rq == NULL) {
	    fprintf(stderr,
		    "TRL_CHECK_RITZ: Failed to allocated workspace RQ.\n");
	    exit(0);
	}
	rqi = 1;
    } else if (lwrk >= (4 * ncol)) {
	rq = &wrk[0];
	gsumwrk = &wrk[3 * ncol];
	aq = (double *) malloc(nrow * sizeof(double));
	if (aq == NULL) {
	    fprintf(stderr,
		    "TRL_CHECK_RITZ: Failed to allocated workspace AQ.\n");
	    exit(0);
	}
	aqi = 1;
    } else if (lwrk >= ncol) {
	gsumwrk = wrk;
	aq = (double *) malloc(nrow * sizeof(double));
	if (aq == NULL) {
	    fprintf(stderr,
		    "TRL_CHECK_RITZ: Failed to allocated workspace AQ.\n");
	    exit(0);
	}
	aqi = 1;
	rq = (double *) malloc(3 * ncol * sizeof(double));
	if (rq == NULL) {
	    fprintf(stderr,
		    "TRL_CHECK_RITZ: Failed to allocated workspace RQ.\n");
	    free(aq);
	    exit(0);
	}
	rqi = 1;
    } else {
	/* WRK not provided -- allocate space for AQ and RQ,  */
	/* gsumwrk points to the last third of RQ             */
	aq = (double *) malloc(nrow * sizeof(double));
	if (aq == NULL) {
	    printf("TRL_CHECK_RITZ: Failed to allocated workspace AQ.\n");
	    return;
	}
	aqi = 1;
	rq = (double *) malloc((3 * ncol) * sizeof(double));
	if (rq == NULL) {
	    printf("TRL_CHECK_RITZ: Failed to allocated workspace RQ.\n");
	    free(aq);
	    return;
	}
	rqi = 1;
	gsumwrk = (double *) malloc(ncol * sizeof(double));
	if (gsumwrk == NULL) {
	    printf
		("TRL_CHECK_RITZ: Failed to allocate workspace GSUMWRK.\n");
	    free(aq);
	    free(rq);
	    return;
	}
	gsumwrki = 1;
    }
    memset(aq, 0, nrow * sizeof(double));
    memset(rq, 0, 2 * ncol * sizeof(double));
    memset(gsumwrk, 0, ncol * sizeof(double));
    /* go through each Ritz pair one at a time, compute Rayleigh  */
    /* quotient and the corresponding residual norm               */
    res = &rq[ncol];
    for (i = 0; i < ncol; i++) {
#ifdef TRL_FORTRAN_COMPATIBLE
	op(&nrow, &i__1, &rvec[i * ldrvec], &ldrvec, aq, &nrow);
#else
	op(nrow, i__1, rvec + i*ldrvec, ldrvec, aq, nrow, info->mvparam);
#endif
	/* Rayleigh quotient -- assuming rvec(:,i) has unit norm */
	rq[i] = trl_ddot(nrow, &rvec[i * ldrvec], c__1, aq, c__1);
	trl_g_sum(info->mpicom, 1, &rq[i], gsumwrk);
	d__1 = -rq[i]; /* indent separated =- into = - */
	trl_daxpy(nrow, d__1, &rvec[i * ldrvec], c__1, aq, c__1);
	res[i] = trl_ddot(nrow, aq, c__1, aq, c__1);
    }
    trl_g_sum(info->mpicom, ncol, res, gsumwrk);
    for (i = 0; i < ncol; i++) {
	res[i] = sqrt(res[i]);
    }
    /* compute the error estimate based on computed residual norms */
    /*  and the Ritz values                                        */
    err = &rq[2 * ncol];
    gapl = alpha[ncol - 1] - alpha[0];
    for (i = 0; i < ncol - 1; i++) {
	gapr = alpha[i + 1] - alpha[i];
	gapl = min(gapl, gapr);
	if (res[i] >= gapl) {
	    err[i] = res[i];
	} else {
	    err[i] = res[i] * res[i] / gapl;
	}
	gapl = gapr;
    }
    if (res[ncol - 1] >= gapl) {
	err[ncol - 1] = res[ncol - 1];
    } else {
	err[ncol - 1] = res[ncol - 1] * res[ncol - 1] / gapl;
    }
    /* if writing to stdout, only PE 0 does it */
    fp = info->log_fp;
    if (fp == NULL) {
	trl_reopen_logfile(info);
	fp = info->log_fp;
    }
    if (fp != stdout || info->my_pe <= 0) {
	if (info->stat != 0) {
#ifdef __DEBUG_OUT
	    ferr = fopen("error.txt", "a");
	    fprintf(ferr, " ** exit error (%d) ** \n", info->stat);
	    fclose(ferr);
#endif
	    *check = -4;
	}
	/* print out the information */
	fprintf(fp, "TRL_CHECK_RITZ: \n");
	fprintf(fp,
		"           Ritz value       res norm   res diff  est error  diff w rq  act. error\n");
#ifdef __DEBUG_OUT
	ferr = fopen("error.txt", "a");
#endif
	if (beta != NULL && eval != NULL) {
	    for (i = 0; i < ncol; i++) {
		icheck = 0;
		fprintf(fp,
			"%21.14f    %11.3e%11.3e%11.3e%11.3e %11.3e%11.3e\n",
			alpha[i], res[i], beta[i] - res[i], err[i],
			rq[i] - alpha[i], eval[i] - alpha[i], eval[i]);
		/* check the accuracy of results.. */
		if (fabs(beta[i] - res[i]) > 0.00001) {
#ifdef __DEBUG_OUT
		    fprintf(ferr,
			    " res diff[%d] = |beta-res| = %5.3e - %5.3e = %5.3e > 0.00001\n",
			    i, beta[i], res[i], fabs(beta[i] - res[i]));
#endif
		    *check = *check - 1;
		    icheck++;
		}

		if (fabs(rq[i] - alpha[i]) > nrow * nrow * info->tol) {
#ifdef __DEBUG_OUT
		    fprintf(ferr,
			    " diff rq[%d] = |rq-alpha| = %5.3e - %5.3e = %5.3e > nrow*nor*tau = %5.3e\n",
			    i, rq[i], alpha[i], fabs(rq[i] - alpha[i]),
			    nrow * nrow * info->tol);
#endif
		    *check = *check - 1;
		    icheck++;
		}

		if (fabs(eval[i] - alpha[i]) > 10 * nrow * nrow * info->tol ||
		    fabs(eval[i] - alpha[i]) > 10 * err[i]) {
#ifdef __DEBUG_OUT
		    fprintf(ferr,
			    " act. error[%d] = |exact-alpha| = %5.3e - %5.3e = %5.3e > 10*nrow*nrow*tau =%5.3e or 10*est err = %5.3e\n",
			    i, eval[i], alpha[i], fabs(eval[i] - alpha[i]),
			    10 * nrow * nrow * info->tol, 10 * err[i]);
#endif
		    *check = *check - 1;
		    icheck++;
		}
	    }

	} else if (beta != NULL) {
	    for (i = 0; i < ncol; i++) {
		fprintf(fp, "%21.14f    %11.3e%11.3e%11.3e%11.3e\n",
			alpha[i], res[i], beta[i] - res[i], err[i],
			rq[i] - alpha[i]);
		/* check the accuracy of results.. */
		if (fabs(beta[i] - res[i]) > 0.00001) {
#ifdef __DEBUG_OUT
		    fprintf(ferr,
			    " res diff[%d] = |beta-res| = %5.3e - %5.3e = %5.3e > 0.00001\n",
			    i, beta[i], res[i], fabs(beta[i] - res[i]));
#endif
		    *check = *check - 1;
		    icheck++;
		}

		if (fabs(rq[i] - alpha[i]) > nrow * nrow * info->tol) {
#ifdef __DEBUG_OUT
		    fprintf(ferr,
			    " diff rq[%d] = |rq-alpha| = %5.3e - %5.3e = %5.3e > nrow*nor*tau = %5.3e\n",
			    i, rq[i], alpha[i], fabs(rq[i] - alpha[i]),
			    nrow * nrow * info->tol);
#endif
		    *check = *check - 1;
		    icheck++;
		}
	    }
	} else if (eval != NULL) {
	    for (i = 0; i < ncol; i++) {
		fprintf(fp,
			"%21.14f     %11.3e           %11.3e%11.3e%11.3e%11.3e\n",
			alpha[i], res[i], err[i], rq[i] - alpha[i],
			eval[i] - alpha[i], eval[i]);
	    }
	} else {
	    for (i = 0; i < ncol; i++) {
		fprintf(fp, "%21.14f    %11.5e           %11.3e%11.3e\n",
			alpha[i], res[i], err[i], rq[i] - alpha[i]);
	    }
	}
#ifdef __DEBUG_OUT
	fclose(ferr);
#endif
    }
    if (info->nec < info->ned)
	*check = 1;
    if (rqi > 0) {
	free(rq);
    }
    if (aqi > 0) {
	free(aq);
    }
    if (gsumwrki > 0) {
	free(gsumwrk);
    }
    trl_close_logfile(info);
/*
// .. end of trl_check_ritz_
*/
}

void
trl_rayleigh_quotients(trl_matvec op, trl_info * info, int ncol, double *evec,
		       int lde, double *eres, double *base)
{
/*
// Purpose:
// ========
// Compute Rayleigh quotients, when it is given a set of Ritz vectors and Ritz values,
// normalize the Ritz vectors, and compute their Rayleigh quotients to replace the Ritz values.
//
// Arguments:
// ==========
// op       (input) function pointer
//           On entry, points to the matrix-vector multiplication routine.
//
// info     (input) pointer to the structure trl_info_
//           On entry, points to the data structure to store the current information about
//           the eigenvalue problem and the progress of TRLAN.
//
// evec     (input) double precision array of dimension (nloc,ncol)
//           On entry, stores the portion of eigenvectors on this PE.
//
// base     (workspace)
//           The workspace used to store results of MATVEC
//
// eres     (output) double precision array of dimension (ncol)
//           On exist, store new Ritz values and new residual norms, i.e., if there are NEV
//           Ritz pairs, eres(1:NEV) stores the new Rayleigh quotient and eres(nev+1:nev+nev)
//           stores the new residual norms.
//
// base     (optional) double precision array od dimension (nloc)
//           If provided, double precision workspace.
//
// ..
// .. local parameters ..
*/
    int c__1 = 1;
    int i__1 = 1;
    double d__1;
/*
// ..
// .. local variables ..
*/
    int i, nrow;
    double wrk[4], *avec;
/*
// ..
// .. executable statements ..
*/
    nrow = info->nloc;
    if (ncol <= 0)
	return;
    if (base != NULL) {
	avec = base;
    } else {
	avec = (double *) malloc(nrow * sizeof(double));
    }
    memset(avec, 0, nrow * sizeof(double));
    if (info->verbose >= 0) {
	FILE *fp = info->log_fp;
	if (fp == NULL) {
	    trl_reopen_logfile(info);
	}
	fp = info->log_fp;
	fprintf(fp,
		"TRLAN computing Rayleigh Quotients for %d Ritz pairs\n",
		ncol);
    }
    /* loop through each vector to normalize the vector, compute Rayleigh  */
    /* quotient and compute residual norm of the new Ritz pairs            */
    for (i = 0; i < ncol; i++) {
	wrk[0] =
	    trl_ddot(nrow, evec + i * lde, c__1, evec + i * lde, c__1);
#ifdef TRL_FORTRAN_COMPATIBLE
	op(&nrow, &i__1, evec + i * lde, &lde, avec, &nrow);
#else
	op(nrow, i__1, evec + i * lde, lde, avec, nrow, info->mvparam);
#endif
	wrk[1] = trl_ddot(nrow, evec+ i * lde, c__1, avec, c__1);
	trl_g_sum(info->mpicom, 2, wrk, &wrk[2]);
	info->matvec = info->matvec + 1;
	info->flop = info->flop + 4 * nrow;
	if (wrk[0] > 0.0) {
	    eres[i] = wrk[1] / wrk[0];
	    d__1 = -eres[i];
	    trl_daxpy(nrow, d__1, evec + i * lde, c__1, avec, c__1);
	    wrk[1] = trl_ddot(nrow, avec, c__1, avec, c__1);
	    trl_g_sum(info->mpicom, 1, &wrk[1], &wrk[2]);
	    wrk[0] = 1.0 / sqrt(wrk[0]);
	    eres[ncol + i] = wrk[0] * sqrt(wrk[1]);
	    d__1 = wrk[0];
	    trl_dscal(nrow, d__1, evec + i * lde, c__1);
	    info->flop = info->flop + 6 * nrow + 3;
	} else {
	    eres[i] = -DBL_MAX;
	    eres[ncol + i] = -DBL_MAX;
	}
    }
    if (base == NULL)
	free(avec);
    trl_close_logfile(info);
/*
// .. end of trl_rayleigh_quotients_ ..
//
*/
}

void trl_clear_counter(trl_info * info, int nrow, int mev, int lde)
{
/*
// Purpose:
// ========
// Clears the counters inside info and performs a minimal check on the input parameters.
//
// Arguments:
// ==========
// info    (input/output) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//          On exit, the information is cleared.
//
// nrow    (input) integer
//          On entry, specifies the number of rows that is on this proccesor.
//
// mev     (input) integer
//          On entry, specifies the number of Ritz pairs, that can be stored in
//          eval and evec.
//
// lde     (input) integer
//          On entry, specifies the leading dimension of the array evec, i.e.,
//          (lde >= nrow).
//
// ..
// .. local scalars ..
*/
    int ntmp;
/*
// ..
// .. executable statements ..
*/
    info->stat = 0;
    if (nrow != info->nloc || nrow > info->ntot) {
	printf("TRLAN: info not setup for this problem.\n");
	printf("       Please reset info before calling TRLAN.\n");
	info->stat = -1;
    }
    if (info->nec < 0)
	info->nec = 0;
    if (lde < nrow) {
	printf("TRLAN: leading dimension of EVEC to small.\n");
	info->stat = -2;
    }
    if (info->tol >= 1.0) {
	info->tol = sqrt(DBL_EPSILON);
    } else if (info->tol <= DBL_MIN) {
	info->tol = DBL_EPSILON;
    }
    if (info->ned + info->ned >= info->ntot) {
	printf
	    ("TRLAN: info->ned (%d) is large relative to the matrix dimension (%d)\n",
	     info->ned, info->ntot);
	printf
	    (" **    It is more appropriate to use LAPACK dsyev/ssyev.\n");
	if (info->ned > info->ntot) {
	    info->ned = min(info->ntot - 1, info->maxlan - 3);
	    printf("TRLAN: ** reduced ned to %d **\n", info->ned);
	}
    }
    if (mev < info->ned) {
	printf
	    ("TRLAN: array EVAL and EVEC can not hold wanted number of eigenpairs.\n");
	info->stat = -3;
    }
    if (info->ntot < 10) {
	printf
	    ("TRLAN is not designed to work with such a small matrix(%dx%d).  Please use LAPACK or EISPACK instead.\n",
	     info->ntot, info->ntot);
	info->stat = -4;
    }
    info->nrand = info->stat;
    info->stat = trl_sync_flag(info->mpicom, info->nrand);

    /* decide what is a good maximum basis size to use */
    if (info->maxlan < info->ned + 3) {
	info->maxlan =
	    info->ned + min(info->ned,
			    20) + (int) (log((double) info->ntot));
	info->maxlan = min(info->maxlan, info->ntot);
	printf("TRLAN: ** reset maxlan to %d! **\n", info->maxlan);
    }
    if (info->maxlan < mev) {
	ntmp = min(info->ntot, max(100 + info->ned, 10 * info->ned));
	if (mev < ntmp) {
	    info->maxlan = mev;
	} else {
	    info->maxlan = ntmp;
	}
    }
    if (info->maxlan < 5) {
	printf
	    ("TRLAN must have at least 5 basis vectors, it is currently %d.\n",
	     info->maxlan);
	info->stat = -5;
    }

    /* clear regular counters */
    info->tmv = -1;
    info->klan = min(info->maxlan, info->ntot);
    if (info->restart >= 7) {
	info->klan =
	    min(info->maxlan, max(100, min(info->klan, 2 * (info->ned))));
    }
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
    info->tick_t = 0.0;
    info->clk_op = 0;
    info->tick_o = 0.0;
    info->clk_orth = 0;
    info->tick_h = 0.0;
    info->clk_res = 0;
    info->tick_r = 0.0;
    info->clk_in = 0;
    info->clk_out = 0;
    info->wrds_in = 0;
    info->wrds_out = 0;
    info->avgm = 0.0;
    return;
/*
// .. end of trl_clear_counter ..
*/
}

void trl_print_setup(trl_info * info, int lbas, int lmis, int lwrk)
{
/*
// Purpose:
// ========
// Print the definition of the eigenvalue problme.
//
// Arguments:
// ==========
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
// lbas    (input) integer
//          On entry, specifies the size of workspace required to store Lanczos basis, i.e.,
//          nrow*(maxlan-mev).
//
// lmis    (input) integer
//          On entry, specifies the size of miscellenious workspace required to solve the
//          eigen problem.
//
// lwrk    (input) integer
//          On entry, specifies the size of workspace provided by the user.
//
// ..
// .. executable statements ..
// print the problem parameters
*/
    if (info->lohi > 0) {
	fprintf(info->log_fp,
		"TRLAN is to compute %6d largest eigenpair(s).\n",
		info->ned);
    } else if (info->lohi < 0) {
	fprintf(info->log_fp,
		"TRLAN is to compute %6d smallest eigenpair(s).\n",
		info->ned);
    } else {
	fprintf(info->log_fp,
		"TRLAN is to compute %6d first converged eigenpair(s).\n",
		info->ned);
    }
    fprintf(info->log_fp,
	    "Problem dimension: %9d (PE:%4d) %12d (Global)\n", info->nloc,
	    info->my_pe, info->ntot);
    fprintf(info->log_fp, "Maximum basis size:                   %10d\n",
	    info->maxlan);
    fprintf(info->log_fp, "Dynamic restarting scheme:            %10d\n",
	    info->restart);
    fprintf(info->log_fp, "Maximum applications of the operator: %10d\n",
	    info->maxmv);
    fprintf(info->log_fp, "Relative convergence tolerance: %10e\n",
	    info->tol);
    /* initial guess */
    if (info->guess == 1) {
	fprintf(info->log_fp, "User provided the starting vector.\n");
    } else if (info->guess == 0) {
	fprintf(info->log_fp, "TRLAN uses [1,1,...] as starting vctor.\n");
    } else if (info->guess < 0) {
	fprintf(info->log_fp,
		"TRLAN generates a random starting vector.\n");
    } else if (info->oldcpf == 0 || strlen(info->oldcpf) == 0) {
	fprintf(info->log_fp,
		"Restarting with existing checkpoint files %s ####\n",
		info->oldcpf);
    } else {
	fprintf(info->log_fp,
		"Restarting with existing checkpoint files %s ####\n",
		info->cpfile);
    }
    if (info->cpflag > 0) {
	fprintf(info->log_fp,
		"TLRAN will write about %d sets of checkpointing files %s ####.\n",
		info->cpflag, info->cpfile);
    }
    /* print the workspace size parameters */
    fprintf(info->log_fp, "(required) array BASE size is %d\n", lbas);
    fprintf(info->log_fp, "(required) array MISC size is %d\n", lmis);
    if (lwrk > 0) {
	fprintf(info->log_fp,
		"Caller has supplied a work array with %d elements.\n",
		lwrk);
    } else {
	fprintf(info->log_fp, "Caller did not supply work array.\n");
    }
}

void
trl_ritz_projection(trl_matvec op, trl_info * info, int mev, double *evec,
		    int lde, double *eres, double *wrk, int lwrk, double *base)
{
/*
// Purpose
// =======
// A separate Rayleigh-Ritz projection routine
// Given a set of approximately orthonormal vectors (V), this routine
// performs the following operations
//  (1) V'*V ==> G
//  (2) R'*R :=  G
//  (3) V'*A*V => H1, inv(R')*H1*inv(R) => H
//  (4) Y*D*Y' := H
//  (5) V*inv(R)*Y => V, diag(D) => lambda,
//      r(i) = ||A*V(:,i)-lambda(i)*V(:,i)||
//
// Arguments:
// ==========
// op       (input) function pointer
//           On entry, points to the matrix-vector multiplication routine.
//
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the information
//          about the eigenvalue problem and the progress of TRLAN
//
// evec    (output) double precision vector of lenvth (nrow*mev)
//          On exit, stores the eigenvectors.
//
// lde     (input) integer
//          On entry, specifies the leading dimension of the array evec, i.e.,
//
// eres    (output) double precision vector of length (2*nev)
//          the array to store new Ritz values and residual norms
//
// base    (workspace) double precision vector (nrow)
//          Workspace to store the result of matrix-vector operation.
//
// wrk     (workspace) double precision vector of length (lwrk)
//          Workspace to store projection matrix, etc.
//
// local variables
*/
    extern int dsyev_();
    extern void dpotrf_();
    extern void dtrtrs_();

    char trans = 'T', notrans = 'N', upper = 'U', job = 'V';
    double one = 1.0, zero = 0.0;
    int i__1 = 1;
    int i, j, ierr, nev, nsqr, nrow, iuau, irvv, lwrk2;
    double d__1;
    double *rvv, *uau, *wrk2, *avec;
/*
// ..
// .. executable statements ..
*/
    nrow = info->nloc;
    if (info->nec > 0) {
	nev = info->nec + 1;
    } else {
	nev = min(info->ned, mev - 1);
	if (info->lohi != 0)
	    nev++;
    }
    nsqr = nev * nev;
    if (lwrk < 0) {
	lwrk = 0;
    }
    if (base != NULL) {
	avec = base;
    } else if (mev > nev) {
	avec = evec + (mev - 1) * lde;
    } else {
	avec = (double *) malloc(sizeof(double) * nrow);
    }
    if (info->verbose >= 0) {
	if (info->log_fp == NULL) {
	    trl_reopen_logfile(info);
	}
	fprintf(info->log_fp,
		"TRLAN performing a Rayleigh-Ritz project for %d vectors.",
		nev);
    }
    /* memory allocation -- need 3*nev*nev elements, will allocate them     */
    /* in two consecutive blocks, uau(nev*nev), rvv(2*nev*nev)              */
    /* in actual use, rvv is further split in two until the last operation  */
    iuau = nsqr;
    irvv = nsqr + nsqr;
    if (lwrk >= iuau + irvv) {
	uau = wrk;
	rvv = &wrk[nsqr];
	wrk2 = &wrk[nsqr + nsqr];
	lwrk2 = lwrk - nsqr - nsqr;
    } else if (lwrk >= irvv) {
	rvv = wrk;
	wrk2 = &wrk[nsqr];
	lwrk2 = lwrk - nsqr;
	uau = (double *) malloc(nsqr * sizeof(double));
	if (uau == NULL) {
	    info->stat = -231;
	    goto end;
	}
    } else if (lwrk >= iuau) {
	uau = wrk;
	rvv = (double *) malloc((nsqr + nsqr) * sizeof(double));
	if (rvv == NULL) {
	    info->stat = -232;
	    goto end;
	}
	wrk2 = &rvv[nsqr];
	lwrk2 = nsqr;
    } else {
	uau = (double *) malloc(nsqr * sizeof(double));
	if (uau == NULL) {
	    info->stat = -231;
	    goto end;
	}
	rvv = (double *) malloc((nsqr + nsqr) * sizeof(double));
	if (rvv == NULL) {
	    info->stat = -232;
	    goto end;
	}
	wrk2 = &rvv[nsqr];
	lwrk2 = nsqr;
    }
    /* step (1) : V'*V ==> G */

    trl_dgemm(&trans, &notrans, nev, nev, nrow, one, evec, lde, evec, lde,
	      zero, rvv, nev);
    trl_g_sum(info->mpicom, nsqr, rvv, wrk2);

    /* step (2) : Choleskey factorization of G */
    dpotrf_(&upper, &nev, rvv, &nev, &ierr);
    if (ierr != 0) {
	info->stat = -234;
	goto end;
    }
    /* step (3) : compute H_1 = V'*A*V                              */
    /* use the first nrow elements of avec to store the results of  */
    /* matrix-vector multiplication                                 */
    memset(wrk2, 0, lwrk2 * sizeof(double));
    for (i = 1; i <= nev; i++) {
#ifdef TRL_FORTRAN_COMPATIBLE
	op(&nrow, &i__1, evec + (i-1)*lde, &lde, avec, &nrow);
#else
	op(nrow, i__1, evec + (i-1)*lde, lde, avec, nrow, info->mvparam);
#endif
	trl_dgemv(&trans, nrow, i, one, evec, lde, avec, i__1, zero,
		  &wrk2[(i - 1) * nev], i__1);
    }
    trl_g_sum(info->mpicom, nsqr, wrk2, uau);
    for (i = 1; i < nev; i++) {
	for (j = 0; j < i; j++) {
	    wrk2[i + j * nev] = wrk2[(i - 1) * nev + j];
	}
    }
    /* compute solution of R^T H_2 = H_1 */
    dtrtrs_(&upper, &trans, &notrans, &nev, &nev, rvv, &nev, wrk2, &nev,
	    &ierr);
    if (ierr != 0) {
	info->stat = -235;
	goto end;
    }
    /* compute solution of R^T H = H_2^T */
    for (i = 1; i < nev; i++) {
	for (j = 0; j < nev; j++) {
	    uau[i + j * nev] = wrk2[(i - 1) * nev + j];
	}
    }
    dtrtrs_(&upper, &trans, &notrans, &nev, &nev, rvv, &nev, uau, &nev,
	    &ierr);
    if (ierr != 0) {
	info->stat = -236;
	goto end;
    }
    /* solve the small symmetric eigenvalue problem */
    dsyev_(&job, &upper, &nev, uau, &nev, eres, wrk2, &nsqr, &ierr);
    if (ierr != 0) {
	info->stat = -237;
	goto end;
    }
    /* solve R Y = Y to prepare for multiplying with V */
    dtrtrs_(&upper, &notrans, &notrans, &nev, &nev, rvv, &nev, uau, &nev,
	    &ierr);
    if (ierr != 0) {
	info->stat = -238;
	goto end;
    }
    /* call trl_ritz_vector to do the final multiplication */
    if (lwrk >= 3 * nsqr) {
	wrk2 = &wrk[nsqr];
    } else if (lwrk >= nsqr + nsqr) {
	wrk2 = wrk;
    } else {
	wrk2 = rvv;
    }
    i = lwrk2;
    trl_ritz_vectors(nrow, 0, nev, uau, nev, evec, lde, nev, avec, nrow,
		     0, wrk2, i);
    /* compute the residual norms */
    for (i = 0; i < nev; i++) {
#ifdef TRL_FORTRAN_COMPATIBLE
	op(&nrow, &i__1, evec + i * lde, &lde, avec, &nrow);
#else
	op(nrow, i__1, evec + i * lde, lde, avec, nrow, info->mvparam);
#endif
	d__1 = eres[i];
	trl_daxpy(nrow, d__1, evec + i * lde, i__1, avec, i__1);
	eres[nev + i] = trl_ddot(nrow, avec, i__1, avec, i__1);
    }
    trl_g_sum(info->mpicom, nev, &eres[nev], avec);
    for (i = nev; i < nev + nev; i++) {
	if (eres[i] > 0.0) {
	    eres[i] = sqrt(eres[i]);
	} else {
	    eres[i] = -DBL_MAX;
	}
    }
    if (info->lohi < 0) {
	for (i = nev - 1; i < nev + nev - 2; i++) {
	    eres[i] = eres[i + 1];
	}
    } else if (info->lohi > 0) {
	for (i = 0; i < nev - 1; i++) {
	    eres[i] = eres[i + 1];
	    memcpy(evec + i * lde, evec + (i + 1) * lde, nrow);
	}
	for (i = nev - 1; i < nev + nev - 2; i++) {
	    eres[i] = eres[i + 2];
	}
    }
  end:
    if (lwrk < iuau) {
	free(uau);
	free(rvv);
    } else if (lwrk < irvv) {
	free(rvv);
    } else if (lwrk < iuau + irvv) {
	free(uau);
    }
}

