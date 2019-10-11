/*
// ZTRLan routine 
// Lawrence Berkeley National Lab.
//
//
// This file contains most auxillary routines which are not extensively
// used in during the normal operations of the TRLan, e.g., printing,
// error checking, etc.
*/
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>

#include "trl_map.h"
#include "trlan.h"
#include "trlan_i.h"
#include "trlaux_i.h"
#include "trl_comm_i.h"


int close_file(FILE * fp, int err1, int err2)
{
/*
// Purpose
// =======
// Closes the file handler.
//
// Arguments
// =========
// fp      (input) pointer to file
//          On entry, points to the file to be closed.
//
// err1    (input) integer
//          On entry, specifies the return value when closing file failed.
//
// err2    (input) integer
//          On entry, specifies the return value when closing file succeeded.
//
// ..
// .. executable statements ..
*/
    if (fclose(fp) != 0) {
	return err2;
    }
    return err1;
/*
// .. end of close_file ..
*/
}

void trl_open_logfile(trl_info * info)
{
/*
// Purpose:
// ========
// Opens log file.
//
// Arguments:
// ==========
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the information
//          about the eigenvalue problem and the progress of TRLAN.
//
// ..
// .. local arrays ..
*/
    char filename[STRING_LEN];
/*
// ..
// .. executable statements ..
*/
    if (info->log_file != 0 && strlen(info->log_file) > 0) {
	trl_pe_filename(STRING_LEN, filename, info->log_file, info->my_pe,
			 info->npes);
	info->log_fp = fopen(filename, "w");
    } else {
	info->log_fp = stdout;
    }
/*
// .. end of trl_open_logfilie_ ..
*/
}

void trl_reopen_logfile(trl_info * info)
{
/*
// Purpose:
// ========
// Opens log file.
//
// Arguments:
// ==========
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the information
//          about the eigenvalue problem and the progress of TRLAN.
//
// ..
// .. local arrays ..
*/
    char filename[STRING_LEN];
/*
// ..
// .. executable statements ..
*/
    if (info->log_file != 0 && strlen(info->log_file) > 0) {
	trl_pe_filename(STRING_LEN, filename, info->log_file, info->my_pe,
			 info->npes);
	info->log_fp = fopen(filename, "a");
    } else {
	info->log_fp = stdout;
    }
/*
// .. end of trl_open_logfilie_ ..
*/
}

void trl_close_logfile(trl_info * info)
{
/*
// Purpose:
// ========
// Closes log file.
//
// Arguments:
// ==========
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the information
//          about the eigenvalue problem and the progress of TRLAN.
//
// ..
// .. executable statements ..
*/
    if (info->log_fp != NULL && info->log_fp != stdout) {
	fclose(info->log_fp);
    }
    info->log_fp = NULL;
/*
// .. end of trl_close_logfile ..
*/
}

void trl_open_cptfile(trl_info * info)
{
/*
// Purpose:
// ========
// Opens check points file.
//
// Arguments:
// ==========
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the information
//          about the eigenvalue problem and the progress of TRLAN.
//
// ..
// .. local arrays ..
*/
    char filename[STRING_LEN];
/*
// ..
// .. executable statements ..
*/
    if (info->cpfile != 0 && strlen(info->cpfile) > 0) {
	trl_pe_filename(STRING_LEN, filename, info->cpfile, info->my_pe,
			 info->npes);
	info->cpt_fp = fopen(filename, "w");
    } else {
	info->cpt_fp = stdout;
    }
/*
// .. end of trl_open_cptfile ..
*/
}

void trl_close_cptfile(trl_info * info)
{
/*
// Purpose:
// ========
// Closes check point file.
//
// Arguments:
// ==========
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the information
//          about the eigenvalue problem and the progress of TRLAN.
//
// ..
// .. executable statements ..
*/
    if (info->cpt_fp != stdout) {
	fclose(info->cpt_fp);
    }
    info->cpt_fp = NULL;
/*
// .. end of trl_close_cptfile ..
*/
}

void trl_print_int(trl_info * info, char *title, int size_array,
		    int *array, int inc)
{
/*
// Purpose:
// ========
// Print an integer array for debugging.
//
// Arguments:
// ==========
// info        (input) pointer to the structure trl_info_
//              On entry, points to the data structure to store the information
//              about the eigenvalue problem and the progress of TRLAN
//
// title       (input) character string
//              On entry, specifies the title of the information to be printed.
//
// size_array  (input) integer
//              On entry specifies, the number of integers to be printed.
//
// array       (input) integer array of length ((size_array-1)*inc+1)
//              On entry, contains the integer to be printed.
//
// inc         (input) integer
//              On entry, specifies how the index to array should be incremented.
//
// ..
// .. local scalars ..
*/
    int i;
/*
// ..
// .. executable statements ..
*/
    fprintf(info->log_fp, "PE %d : %s", info->my_pe, title);
    if (size_array > 2) {
	fprintf(info->log_fp, "\n");
    }
    for (i = 0; i < size_array; i += inc) {
	fprintf(info->log_fp, "%10d", array[i]);
	if ((i % 8) == 7)
	    fprintf(info->log_fp, "\n");
    }
    if (((size_array - 1) % 8) != 7)
	fprintf(info->log_fp, "\n");
/*
// .. end of trl_print_int ..
//
*/
}

void trl_print_real(trl_info * info, char *title, int size_array,
		     double *array, int inc)
{
/*
// Purpose
// =======
// Print a double precision array for debugging.
//
// Arguments:
// ==========
// info        (input) pointer to the structure trl_info_
//              On entry, points to the data structure to store the information
//              about the eigenvalue problem and the progress of TRLAN
//
// title       (input) character string
//              On entry, specifies the title of the information to be printed.
//
// size_array  (input) integer
//              On entry specifies, the number of doubles to be printed.
//
// array       (input) double array of length ((size_array-1)*inc+1)
//              On entry, contains the doubles to be printed.
//
// inc         (input) integer
//              On entry, specifies how the index to array should be incremented.
//
// ..
// .. local scalars ..
*/
    int i;
/*
// ..
// .. executable statements ..
*/
    fprintf(info->log_fp, "PE %d : %s", info->my_pe, title);
    if (size_array > 1) {
	fprintf(info->log_fp, "\n");
    }
    for (i = 0; i < size_array; i += inc) {
	fprintf(info->log_fp, " %10.7e", array[i]);
	if ((i % 8) == 7)
	    fprintf(info->log_fp, "\n");
    }
    if (((size_array - 1) % 8) != 7)
	fprintf(info->log_fp, "\n");
/*
// .. end of trl_print_real ..
*/
}

void trl_print_progress(trl_info * info)
{
/*
// Purpose
// =======
// Print the current progress of eigenvalue solution
//
// Arguments:
// ==========
// info        (input) pointer to the structure trl_info_
//              On entry, points to the data structure to store the information
//              about the eigenvalue problem and the progress of TRLAN
//
// ..
// .. executable statements ..
*/
    fprintf(info->log_fp, "MATVEC: %10d,    Nloop: %10d,     Nec: %10d\n",
	    info->matvec, info->nloop, info->nec);
    fprintf(info->log_fp, "Reorth: %10d,    Nrand: %10d,    Ierr: %10d\n",
	    info->north, info->nrand, info->stat);
    fprintf(info->log_fp,
	    "Target: %10.3e,   ResNrm: %10.3e,    CFact: %10.3e\n",
	    info->trgt, info->tres, info->crat);
/*
// .. end of trl_print_progress ..
*/
}

void trl_check_orth(trl_info * info, int nrow, double *v1, int ld1,
		     int j1, double *v2, int ld2, int j2, double *wrk,
		     int lwrk)
{
/*
// Purpose:
// ========
// Check orthogonality of the basis.
//
// Arguments:
// ==========
// info     (input) pointer to the structure trl_info_
//           On entry, points to the data structure to store the information 
//           about the eigenvalue problem and the progress of TRLAN.
//
// nrow     (input) integer
//           On entry, specifies the problem size, i.e., the number of rows in 
//           v1 and v2.
//
// v1       (input) double precision array of diimension (ld1,j1)
//           On entry, contains the first part of the basis.
//
// ld1      (input) integer
//           On entry, specifies the leading diimension of v1.
//
// j1       (input) integer
//           On entry, specifies the last column index of v1, containing the 
//           basis.
//
// v2       (input) double precision array of dimension (ld2,j2)
//           On entry, contains the second part of the basis.
//
// ld2      (input) integer
//           On entry, specifies the leading dimension of v2.
//
// j2       (input) integer
//           On entry, specifies the last column index of v2, containing the 
//           basis.
//
// wrk      (workspace) double precision array of length (lwrk)
//
// lwrk     (input) integer
//           On entry, specifies the size of workspace.
//
// ..
// .. local parameters ..
*/
    double one = 1.0, zero = 0.0;
    long c__1 = 1;
/*
// ..
// .. local variables
*/
    int i, j, k, jnd;
    double nrmfro, nrminf;
/*
// ..
// .. executable statements ..
*/
    jnd = j1 + j2;
    nrmfro = zero;
    nrminf = zero;
    if (jnd <= 0)
	return;
    if (lwrk < (jnd + jnd)) {
	fprintf(info->log_fp, "TRL_CHECK_ORTH: workspace too small.\n");
	return;
    }
    fprintf(info->log_fp,
	    "TRL_CHECK_ORTH: check orthogonality of %d basis vectors.\n",
	    jnd);
    /*
       // check orthognality of the basis vectors
     */
    for (i = 0; i < j1; i++) {
	trl_g_dot_(info->mpicomp, nrow, v1, ld1, i + 1, v2, ld2, 0,
		   v1 + i * ld1, wrk);
	wrk[i] = wrk[i] - one;
	if (info->verbose > 7) {
	    fprintf(info->log_fp, "Orthogonality level of v(%d) ..\n",
		    i + 1);
	    for (j = 0; j <= i; j++) {
		fprintf(info->log_fp, " %10.3e", wrk[j]);
		if ((j % 8) == 7)
		    fprintf(info->log_fp, "\n");
	    }
	    if ((i % 8) != 7)
		fprintf(info->log_fp, "\n");
	}
	nrmfro =
	    nrmfro + 2 * trl_ddot(i, wrk, c__1, wrk,
				  c__1) + wrk[i] * wrk[i];
	if (i == 0) {
	    wrk[i + 1] = fabs(wrk[i]);
	} else {
	    wrk[i + 1] = max(wrk[i], wrk[i - 1]);
	}
	nrminf = max(nrminf, wrk[i + 1]);
    }
    for (i = 0; i < j2; i++) {
	j = j1 + i;
	trl_g_dot_(info->mpicomp, nrow, v1, ld1, j1, v2, ld2, i + 1,
		   v2 + i * ld2, wrk);
	wrk[j] = wrk[j] - one;
	if (info->verbose > 7) {
	    fprintf(info->log_fp, "Orthogonality level of v(%d) ..\n",
		    j + 1);
	    for (k = 0; k <= j; k++) {
		fprintf(info->log_fp, " %10.3e", wrk[k]);
		if ((k % 8) == 7)
		    fprintf(info->log_fp, "\n");
	    }
	    if ((j % 8) != 7)
		fprintf(info->log_fp, "\n");
	}
	nrmfro =
	    nrmfro + 2 * trl_ddot(j, wrk, c__1, wrk,
				  c__1) + wrk[j] * wrk[j];
	nrminf = max(nrminf, fabs(wrk[j]));
    }
    fprintf(info->log_fp,
	    "Frobenius norm of orthogonality level %10i %4i  %14.5e\n",
	    info->matvec, jnd, sqrt(nrmfro));
    fprintf(info->log_fp,
	    "Maximum abs. value of orthogonality level is  %14.5e\n",
	    nrminf);
/*
// .. end of trl_check_orth ..
*/
}

void
trl_check_recurrence(trl_matvec op,
		      trl_info * info, int nrow, int ncol, double *v1,
		      int ld1, int m1, double *v2, int ld2, int m2,
		      int kept, double *alpha, double *beta, double *wrk,
		      int lwrk)
{
/*
// Purpose
// =======
// Check Lanczos recurrence relation for debug purpose.
//
// Arguments:
// ==========
// op       (input) function pointer
//           On entry, points to the matrix-vector multiplication routine.
//
// info     (input) pointer to the structure trl_info_
//           On entry, points to the data structure to store the information about 
//           the eigenvalue problem and the progress of TRLAN.
//
// nrow     (input) integer
//           On entry, specifies the problem size, i.e., the number of ros in v1 
//           and v2.
//
// ncol     (input) integer
//           On entry, specifies the maximum number of eigenvalues that can be 
//           stored.
//
// v1       (input) double precision array of dimension (ld1,m1)
//           On entry, contains the first part of basis.
//
// ld1      (input) integer
//           On entry, specifies the leading dimension of v1.
//
// m1       (input) integer
//           On entry, specifies the last column index of v1 that contains a 
//           base vector.
//
// v2       (input) double precision array of dimension (ld2,m2)
//           On entry, contains the second part of basis.
//
// ld2      (input) integer
//           On entry, specifies the leading dimension of v2.
//
// m2       (input) integer
//           On entry, specifies the last column index of v2 that contains a 
//           base vector.
//
// kept     (input) integer
//           On entry, specifies the number of basis kept at the last restart.
//
// alpha    (input) integer
//           On entry, contains the alpha values computed so far.
//
// beta     (input) integer
//           On entry, contains the beta values computed so far.
//
// wrk      (workspace) double precision vector of length (lwrk)
//
// lwrk     (input) integer
//           On entry, specifies the size of workspace.
//
// ..
// .. local parameters ..
*/
    long c__1 = 1;
    int i__1 = 1;
    double zero = 0.0, one = 1.0;
/*
// ..
// .. local variables ..
*/
    int i, ii, j, j1, j2, jnd, mv1;
    char title[STRING_LEN];
    double d__1;
    double *aq, *qkp1, *cs, *alf, *bet;
/*
// ..
// .. executable statements ..
*/
    mv1 = m1;
    if (m2 > 0) {
	j2 = m2 - 1;
	j1 = m1;
    } else {
	j2 = 0;
	j1 = m1 - 1;
    }
    jnd = j1 + j2;
    if (lwrk < jnd * 4 + max(jnd * 4, nrow)) {
	fprintf(info->log_fp,
		"TRL_CHECK_RECURRENCE: not enough workspace.\n");
	return;
    }
    if (lwrk >= jnd * 4 + nrow) {
	aq = wrk + (lwrk - nrow);
    } else {
	aq = (double *) malloc(nrow * sizeof(double));
	if (aq == NULL) {
	    fprintf(info->log_fp,
		    "TRL_CHECK_RECURRENCE: failed to allcoate workspace.\n");
	    return;
	}
    }
    memset(wrk, 0, 4 * jnd * sizeof(double));
    cs = &wrk[jnd];
    alf = &wrk[2 * jnd];
    bet = &wrk[3 * jnd];
    /*
       // first type of relation
       // A q_i = Alpha_i q_i + Beta_i q_{k+1}
     */
    if (kept < ncol) {
	qkp1 = v1 + kept * ld1;
    } else {
	qkp1 = v2 + (kept - j1) * ld2;
    }
    for (i = 0; i < min(j1, kept); i++) {
#ifdef TRL_FORTRAN_COMPATIBLE
	op(&nrow, &i__1, v1 + i * ld1, &ld1, aq, &nrow);
#else
	op(nrow, i__1, v1 + i * ld1, ld1, aq, nrow, info->mvparam);
#endif
	for (ii = 0; ii < nrow; ii++) {
	    alf[i] += aq[ii] * v1[i * nrow + ii];
	    aq[ii] -= alpha[i] * v1[i * nrow + ii];
	    bet[i] += aq[ii] * aq[ii];
	    cs[i] += aq[ii] * qkp1[ii];
	    aq[ii] -= beta[i] * qkp1[ii];
	    wrk[i] += aq[ii] * aq[ii];
	}
    }
    for (i = 0; i < (kept - j1); i++) {
	j = i + j1;
#ifdef TRL_FORTRAN_COMPATIBLE
	op(&nrow, &i__1, v2 + i * ld2, &ld2, aq, &nrow);
#else
	op(nrow, i__1, v2 + i * ld2, ld2, aq, nrow, info->mvparam);
#endif
	for (ii = 0; ii < nrow; ii++) {
	    alf[j] += aq[ii] * v2[i * nrow + ii];
	    aq[ii] -= alpha[j] * v2[i * nrow + ii];
	    bet[j] += aq[ii] * aq[ii];
	    cs[j] += aq[ii] * qkp1[ii];
	    aq[ii] -= beta[j] * qkp1[ii];
	    wrk[j] += aq[ii] * aq[ii];
	}
    }
    /*
       // the (k+1)st base vector need to orthogonalize against all previous
       // vectors
     */
    if (jnd > kept) {
#ifdef TRL_FORTRAN_COMPATIBLE
	op(&nrow, &i__1, qkp1, &nrow, aq, &nrow);
#else
	op(nrow, i__1, qkp1, nrow, aq, nrow, info->mvparam);
#endif
	alf[kept] = trl_ddot(nrow, aq, c__1, qkp1, c__1);
	d__1 = -alpha[kept];
	trl_daxpy(nrow, d__1, qkp1, c__1, aq, c__1);
	for (i = 0; i < min(j1, kept); i++) {
	    d__1 = -beta[i];
	    trl_daxpy(nrow, d__1, v1 + i * ld1, c__1, aq, c__1);
	}
	for (i = 0; i < kept - j1; i++) {
	    j = j1 + i;
	    d__1 = -beta[j];
	    trl_daxpy(nrow, d__1, v2+ i * ld2, c__1, aq, c__1);
	}
	bet[kept] = trl_ddot(nrow, aq, c__1, aq, c__1);
	if (kept + 2 <= j1) {
	    cs[kept] =
		trl_ddot(nrow, aq, c__1, v1 + (kept + 1) * ld1, c__1);
	    d__1 = -beta[kept];
	    trl_daxpy(nrow, d__1, v1 + (kept + 1) * ld1, c__1, aq, c__1);
	} else {
	    cs[kept] =
		trl_ddot(nrow, aq, c__1, v2 + (kept + 1 - j1) * ld2,
			 c__1);
	    d__1 = -beta[kept];
	    trl_daxpy(nrow, d__1, v2 + (kept + 1 - j1) * ld2, c__1, aq,
		      c__1);
	}
	wrk[kept] = trl_ddot(nrow, aq, c__1, aq, c__1);
    }
    /*
       // the third kind of relation -- normal three term recurrence
       // depending the fact that if the lower-bound of loop is less than
       // upper bound, the look should not be executed
     */
    for (i = kept + 1; i < j1; i++) {
#ifdef TRL_FORTRAN_COMPATIBLE
	op(&nrow, &i__1, v1 + ld1 * i, &ld1, aq, &nrow);
#else
	op(nrow, i__1, v1 + ld1 * i, ld1, aq, nrow, info->mvparam);
#endif
	if (i < (mv1 - 1)) {
	    for (ii = 0; ii < nrow; ii++) {
		alf[i] += aq[ii] * v1[i * ld1 + ii];
		aq[ii] -=
		    (alpha[i] * v1[i * ld1 + ii] +
		     beta[i - 1] * v1[(i - 1) * nrow + ii]);
		bet[i] += aq[ii] * aq[ii];
		cs[i] += aq[ii] * v1[(i + 1) * ld1 + ii];
		aq[ii] -= beta[i] * v1[(i + 1) * ld1 + ii];
		wrk[i] += aq[ii] * aq[ii];
	    }
	} else {
	    for (ii = 0; ii < nrow; ii++) {
		alf[i] += aq[ii] * v1[i * ld1 + ii];
		aq[ii] -=
		    (alpha[i] * v1[i * ld1 + ii] +
		     beta[i - 1] * v1[(i - 1) * ld1 + ii]);
		bet[i] += aq[ii] * aq[ii];
		cs[i] += aq[ii] * v2[ii];
		aq[ii] -= beta[i] * v2[ii];
		wrk[i] += aq[ii] * aq[ii];
	    }
	}
    }
    for (i = max(0, kept - j1 + 1); i < j2; i++) {
	j = i + j1;
#ifdef TRL_FORTRAN_COMPATIBLE
	op(&nrow, &i__1, v2 + i * ld2, &ld2, aq, &nrow);
#else
	op(nrow, i__1, v2 + i * ld2, ld2, aq, nrow, info->mvparam);
#endif
	if (i > 0) {
	    for (ii = 0; ii < nrow; ii++) {
		alf[j] += aq[ii] * v2[i * ld2 + ii];
		aq[ii] -=
		    (beta[j - 1] * v2[(i - 1) * ld2 + ii] +
		     alpha[j] * v2[i * ld2 + ii]);
		bet[j] += aq[ii] * aq[ii];
		cs[j] += aq[ii] * v2[(i + 1) * ld2 + ii];
		aq[ii] -= beta[j] * v2[(i + 1) * ld2 + ii];
		wrk[j] += aq[ii] * aq[ii];
	    }
	} else {
	    for (ii = 0; ii < nrow; ii++) {
		alf[j] += aq[ii] * v2[ii];
		aq[ii] -=
		    (beta[j - 1] * v1[(j1 - 1) * ld1 + ii] +
		     alpha[j] * v2[ii]);
		bet[j] += aq[ii] * aq[ii];
		cs[j] += aq[ii] * v2[ld2 + ii];
		aq[ii] -= beta[j] * v2[ld2 + ii];
		wrk[j] += aq[ii] * aq[ii];
	    }
	}
    }

    trl_g_sum(info->mpicomp, jnd * 4, wrk, &wrk[jnd * 4]);
    aq[0] = zero;
    for (ii = 0; ii < jnd; ii++) {
	aq[0] += wrk[ii];
    }
    aq[0] = sqrt(aq[0]);
    for (ii = 0; ii < jnd; ii++) {
	wrk[ii] = sqrt(wrk[ii]);
    }
    for (ii = 0; ii < jnd; ii++) {
	if (bet[ii] > zero) {
	    if (beta[ii] < zero) {
		bet[ii] = -sqrt(bet[ii]);
	    } else {
		bet[ii] = sqrt(bet[ii]);
	    }
	    cs[ii] = cs[ii] / bet[ii];
	} else {
	    bet[ii] = zero;
	}
    }
    strcpy(title, "Alpha computed by TRLAN ..");
    trl_print_real(info, title, jnd, alpha, 1);
    strcpy(title, "Alpha computed explicitly in TRL_CHECK_RECURRENCE ..");
    trl_print_real(info, title, jnd, alf, 1);
    strcpy(title, "Differences in alpha ..");
    d__1 = -one;
    trl_daxpy(jnd, d__1, alpha, c__1, alf, c__1);
    trl_print_real(info, title, jnd, alf, 1);
    strcpy(title, "Beta computed by TRLAN ..");
    trl_print_real(info, title, jnd, beta, 1);
    strcpy(title, "Beta computed explicitly in TRL_CHECK_RECURRENCE ..");
    trl_print_real(info, title, jnd, bet, 1);
    strcpy(title, "Differences in beta ..");
    d__1 = -one;
    trl_daxpy(jnd, d__1, beta, c__1, bet, c__1);
    trl_print_real(info, title, jnd, bet, 1);
    strcpy(title, "Error in Lanczos recurrence (overall) =");
    trl_print_real(info, title, 1, aq, 1);
    if (info->verbose > 7) {
	strcpy(title,
	       "|| A q_i - alpha_i q_i - beta_{i-1} q_{i-1} - beta_i q_{i+1} ||..");
	trl_print_real(info, title, jnd, wrk, 1);
	strcpy(title,
	       "(A q_i - alpha_i q_i - beta_{i-1} q_{i-1})*q_{i+1}/beta_i ..");
	trl_print_real(info, title, jnd, cs, 1);
	strcpy(title, "Sine of the angles ..");
	for (ii = 0; ii < jnd; ii++) {
	    cs[ii] = cs[ii] * cs[ii];
	    if (cs[ii] < one) {
		cs[ii] = sqrt(one - cs[ii]);
	    } else {
		cs[ii] = -one;
	    }
	}
	trl_print_real(info, title, jnd, cs, 1);
    }
    if (lwrk < jnd * 4 + nrow)
	free(aq);
/*
// .. end of trl_check_recurrence ..
*/
}

int trl_write_checkpoint(char *filename, int nrow, double *alpha,
			  double *beta, double *evec, int lde, int me,
			  double *base, int ldb, int nb)
{
/*
// Purpose
// =======
// Write a check-point file.
//
// Arguments
// =========
// filename   (input) character string
//             On entry, specifies the name of checkpoint file.
//
// nrow       (input) integer
//             On entry, specifies the problem size.
//
// alpha      (input) double precision array of length (me+nb-1)
//             On entry, contains the alpha values computed so far.
//
// beta       (input) double precision array of length (me+ne-1)
//             On entry, contains the beta values computed so far.
//
// evec       (input) double precision array of dimensioni (lde,me)
//             On entry, contains the first part of basis vectors.
//
// lde        (input) integer
//             On entry, specifies the leading dimension of evec.
//
// me         (input) integer
//             On entry, specifies the last column index of evec, that contains a base vector.
//
// base       (input) double precision array of dimension (ldb,nb)
//             On entry, contains the second part of basis.
//
// ldb        (input) integer
//             On entry, specifies the leading dimension of base.
//
// nb         (input) integer
//             On entry, specifies the last column index of base, that contains a base vector.
//
// ..
// .. local variables ..
*/
    int jnd, i, j;
    FILE *io_fp;
/*
// ..
// .. executable statements ..
*/
    jnd = me + nb - 1;
    io_fp = fopen(filename, "w");
    if (io_fp == NULL) {
	printf("TRL_WRITE_CHECKPOINT: failed to open file: %s.\n",
	       filename);
	return -221;
    }
    if (fwrite(&nrow, sizeof(nrow), 1, io_fp) < 1) {
	return close_file(io_fp, -223, -222);
    }
    if (fwrite(&jnd, sizeof(jnd), 1, io_fp) < 1) {
	return close_file(io_fp, -223, -222);
    }

    for (i = 0; i < jnd; i++) {
	if (fwrite(&alpha[i], sizeof(double), 1, io_fp) < 1) {
	    return close_file(io_fp, -223, -222);
	}
    }
    for (i = 0; i < jnd; i++) {
	if (fwrite(&beta[i], sizeof(double), 1, io_fp) < 1) {
	    return close_file(io_fp, -223, -222);
	}
    }
    for (i = 0; i < me; i++) {
	for (j = 0; j < nrow; j++) {
	    if (fwrite
		(&evec[i * lde + j], sizeof(double), 1,
		 io_fp) < 1) {
		return close_file(io_fp, -223, -222);
	    }
	}
    }
    for (i = 0; i < nb; i++) {
	for (j = 0; j < nrow; j++) {
	    if (fwrite
		(&base[i * ldb + j], sizeof(double), 1,
		 io_fp) < 1) {
		return close_file(io_fp, -223, -222);
	    }
	}
    }
    return close_file(io_fp, 0, -223);
/*
// .. end of trl_write_checkpoint ..
*/
}

int trl_read_checkpoint(char *filename, int nrow, double *evec, int lde,
			 int mev, int *j1, double *base, int ldb, int nbas,
			 int *j2, int nalpha, double *alpha, int nbeta,
			 double *beta)
{
/*
// Purpose
// =======
// Read check-point file
//
// Arguments:
// ==========
// filename   (input) character string
//             On entry, specifies the name of checkpoint file.
//
// nrow       (input) integer
//             On entry, specifies the problem size.
//
// evec       (output) double precision array of dimensioni (lde,j1)
//             On exit, contains the first part of basis vectors stored in checkpoint.
//
// lde         (input) integer
//             On entry, specifies the leading dimension of evec.
//
// mev        (input) integer
//             On entry, specifies the number of eigenvalues converged.
//
// j1         (input) integer
//             On entry, specifies the last column index of evec, that contains a base vector.
//
// base       (output) double precision array of dimension (ldb,nbas)
//             On exit, contains the second part of basis stored in the checkpoint.
//
// ldb        (input) integer
//             On entry, specifies the leading dimension of base.
//
// nbas       (input) integer
//             On entry, specifies the number of columns in base.
//
// j2         (input) integer
//             On entry, specifies the last column index of base, that contains a base vector.
//
// nalpha     (input) integer
//             On entry, specifies the size of alpha
//
// alpha      (output) double precision array of length (nalpha)
//             On exit, contains the alpha values stored in checkpoint.
//
// nbeta      (input) integer
//             On entry, specifies the size of beta.
//
// beta       (output) double precision array of length (nbeta)
//             On exit, contains the beta values stored in checkpoint.
//
// ..
// .. local variables ..
*/
    int i, j;
    FILE *io_fp;
/*
// ..
// .. executable statements ..
*/
    if (lde < nrow || ldb < nrow) {
	printf("TRL_READ_CHECKPOINT: leading dimensions too small.\n");
	return -211;
    }
    /* open file */
    io_fp = fopen(filename, "r");
    if (io_fp == NULL) {
	printf
	    ("TRL_READ_CHECKPOINT: failed to open check-point file %s.\n",
	     filename);
	return -212;
    }
    /* read size information */
    if (fread(j1, sizeof(*j1), 1, io_fp) <= 0) {
	return close_file(io_fp, -215, -216);
    }
    if (fread(j2, sizeof(*j2), 1, io_fp) <= 0) {
	return close_file(io_fp, -215, -216);
    }
    if (*j1 != nrow) {
	printf("TRL_READ_CHECKPOINT: Nrow mismatch.\n");
	return -213;
    }
    if (*j2 > min(nalpha, min(nbeta, mev + nbas - 1))) {
	printf("TRL_READ_CHECKPOINT: MAXLAN too small.");
	return -214;
    }
    /* can continue read all data */
    for (i = 0; i < *j2; i++) {
	if (fread(&alpha[i], sizeof(double), 1, io_fp) <= 0) {
	    return close_file(io_fp, -215, -216);
	}
    }
    for (i = 0; i < *j2; i++) {
	if (fread(&beta[i], sizeof(double), 1, io_fp) <= 0) {
	    return close_file(io_fp, -215, -216);
	}
    }
    *j1 = min(mev, *j2);
    *j2 = *j2 - *j1;
    if (*j1 < mev) {
	for (i = 0; i <= *j1; i++) {
	    for (j = 0; j < nrow; j++) {
		if (fread
		    (&evec[i * lde + j], sizeof(double), 1,
		     io_fp) <= 0) {
		    return close_file(io_fp, -215, -216);
		}
	    }
	}
    } else {
	for (i = 0; i < *j1; i++) {
	    for (j = 0; j < nrow; j++) {
		if (fread
		    (&evec[i * nrow + j], sizeof(double), 1,
		     io_fp) <= 0) {
		    return close_file(io_fp, -215, -216);
		}
	    }
	}
	for (i = 0; i < *j2; i++) {
	    for (j = 0; j < nrow; j++) {
		if (fread
		    (&base[i * ldb + j], sizeof(double), 1, io_fp) <= 0) {
		    return close_file(io_fp, -215, -216);
		}
	    }
	}
    }
    return close_file(io_fp, 0, -215);
/*
// .. end of trl_read_checkpoint ..
*/
}

int indchar(char *a, char b)
{
/*
// Purpose:
// ========
// Used in trl_pe_filename to look for the first occurence of a character b in
// a string a.
//
// Arguments:
// ==========
// a       (input) character string
//          On entry, contains the string to search for a character b.
//
// b       (input) chararcter
//          On entry, specifies the character to look for.
//
// ..
// .. local scalars ..
*/
    char *t = strchr(a, b);
/*
// ..
// .. executable statements ..
*/
    if (t != NULL) {
	return (1 + (t - a));
    } else {
	return strlen(a) + 1;
    }
}

void trl_pe_filename(int nlen, char *filename, char *base, int my_rank, int npe)
{
/*
// Purpose
// =======
// Generates file name from a base name and the PE number.
//
// Arguments:
// ==========
// nlen         (input) integer
//               On entry, specifies the size of filiename.
//
// filename     (output) character string of length <= nlen.
//               On exit, specifies the file name.
//
// base         (input) character string
//               On entry, specifies the leading part of the file name.
//
// my_rank      (input) integer
//               On entry, specifies the PE number.
//
// npe          (input) integer
//               On entry, specifies the number of processors.
//
// ..
// .. local variable ..
*/
    int lead, ndig, len, off;
    char *format;
/*
// ..
// .. executable statements ..
*/
    ndig = 1;
    lead = npe;
    while (lead > 9) {
	lead /= 10;
	++ndig;
    }
    len = indchar(base, ' ') - 1;
    if (nlen < len + ndig + 1) {
	fprintf(stderr,
		"error: not enough space for filename (%d+%d chars).\n",
		len, ndig);
	exit(0);
    }
    memset(filename, 0, nlen * sizeof(char));
    strncpy(filename, base, len);
    off = 1 + ndig % 10;
    off = 5 + 2 * off;
    format = (char *) malloc(off * sizeof(char));
    sprintf(format, "%%s%%0%d.%dd", ndig, ndig);
    sprintf(filename, format, filename, my_rank);
    free(format);
    return;
/*
// .. end of trl_pe_filename ..
*/
}

void trl_time_stamp(FILE * fp)
{
/*
// Purpose
// =======
// Print the current date and time.
//
// Arguments:
// ==========
// fp         (input) pointer to a file
//             On entry, points to the file to output.
//
// ..
// .. local variables ..
*/
    time_t clck;
/*
// ..
// .. executable statements ..
*/
    clck = time(NULL);
    fprintf(fp, "                                                  %s",
	    asctime(localtime(&clck)));
/*
// .. end of trl_time_stamp ..
*/
}
