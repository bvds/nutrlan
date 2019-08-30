/**
   @file Internal header file for the symetric version of the nuTRLan code.
 */
#ifndef	 __TRLANI_H
#define __TRLANI_H
#include <trlan.h> /* trlan_info */
#ifdef BLASWRAP
/* When CLAPACk is complied with -DNO_BLAS_WRAP, then "blaswrap.h" is not required here. */
#include "blaswrap.h"
#endif
/* The maximum number of string allowed for a information titile. */
#define STRING_LEN 132


void trl_clear_counter(trl_info * info, int nrow, int mev, int lde);
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
////
*/
void trl_print_setup(trl_info * info, int lbas, int lmis, int lwrk);
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
//          On entry, specifies the size of miscellenious workspace requred to solve the
//          eigen problem.
//
// lwrk    (input) integer
//          On entry, specifies the size of workspace provided by the user.
//
*/

#endif
