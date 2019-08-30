/*
// ZTRLan routine 
// Lawrence Berkeley National Lab.
//
*/
#include <math.h>

#include "trl_map.h"
#include "dsort2_i.h"

void dsort2(int N, double *ARRAY1, double *ARRAY2)
{
/*
// Purpose
// =======
// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing 
// order of ARRAY1.
//
// Arguments
// =========
// N       (input) INTEGER
//          On entry, specifies the size of the arrays.
//
// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ..
// .. local scalars ..
*/
    int IGAP, I, J;
    double TEMP;
/*
// ..
// .. executable statements ..
*/
    IGAP = N / 2;
    while (IGAP > 0) {
	for (I = IGAP; I < N; I++) {
	    J = I - IGAP;
	    while (J >= 0) {
		if (ARRAY1[J] > ARRAY1[J + IGAP]) {
		    TEMP = ARRAY1[J];
		    ARRAY1[J] = ARRAY1[J + IGAP];
		    ARRAY1[J + IGAP] = TEMP;
		    TEMP = ARRAY2[J];
		    ARRAY2[J] = ARRAY2[J + IGAP];
		    ARRAY2[J + IGAP] = TEMP;
		    J = J - IGAP;
		} else {
		    break;
		}
	    }
	}
	IGAP = IGAP / 2;
    }
/*
// .. end of dsort2_c_
*/
}

void dsort2a(int N, double *ARRAY1, double *ARRAY2)
{
/*
// Purpose
// =======
// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing 
// order of abs(ARRAY1).
//
// Arguments
// =========
// N       (input) INTEGER
//          On entry, specifies the size of the arrays.
//
// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ..
// .. local scalars ..
*/
    int IGAP, I, J;
    double TEMP;
/*
// ..
// .. executable statements ..
*/
    IGAP = N / 2;
    while (IGAP > 0) {
	for (I = IGAP; I < N; I++) {
	    J = I - IGAP;
	    while (J >= 0) {
		if (fabs(ARRAY1[J]) > fabs(ARRAY1[J + IGAP])) {
		    TEMP = ARRAY1[J];
		    ARRAY1[J] = ARRAY1[J + IGAP];
		    ARRAY1[J + IGAP] = TEMP;
		    TEMP = ARRAY2[J];
		    ARRAY2[J] = ARRAY2[J + IGAP];
		    ARRAY2[J + IGAP] = TEMP;
		    J = J - IGAP;
		} else {
		    break;
		}
	    }
	}
	IGAP = IGAP / 2;
    }
/*
// .. end of dsort2ac_
*/
}

void dsort2s(int N, double s, double *ARRAY1, double *ARRAY2)
{
/*
// Purpose
// =======
// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing 
// order of abs(ARRAY1-s).
//
// Arguments
// =========
// N       (input) INTEGER
//          On entry, specifies the size of the arrays.
//
// s       (input) DOUBLE PRECISION
//          On entry, specifies the reference value s.
//          
// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ..
// .. local scalars ..
*/
    int IGAP, I, J;
    double TEMP;
/*
// ..
// .. executable statements ..
*/
    IGAP = N / 2;
    while (IGAP > 0) {
	for (I = IGAP; I < N; I++) {
	    J = I - IGAP;
	    while (J >= 0) {
		if (fabs(ARRAY1[J] - s) > fabs(ARRAY1[J + IGAP] - s)) {
		    TEMP = ARRAY1[J];
		    ARRAY1[J] = ARRAY1[J + IGAP];
		    ARRAY1[J + IGAP] = TEMP;
		    TEMP = ARRAY2[J];
		    ARRAY2[J] = ARRAY2[J + IGAP];
		    ARRAY2[J + IGAP] = TEMP;
		    J = J - IGAP;
		} else {
		    break;
		}
	    }
	}
	IGAP = IGAP / 2;
    }
/*
// .. end of dsort2ac_
*/
}

void dsort2su_(int N, double s, double *ARRAY1, double *ARRAY2)
{
/*
// Purpose
// =======
// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing order 
// of ARRAY1-s if ARRAY1-s is non-negative. Negative ARRAY1-s are ordered after
// those with non-negative ARRAY1-s and in the increasing order of ARRAY1.
// 
//
// Arguments
// =========
// N       (input) INTEGER
//          On entry, specifies the size of the arrays.
//
// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ..
// .. local scalars ..
*/
    int IGAP, I, J;
    double TEMP, v1, v2, d1, d2, maxd;
/*
// ..
// .. executable statements ..
*/
    IGAP = N / 2;
    maxd = fabs(ARRAY1[0]);
    for (I = 1; I < N; I++) {
	if (maxd < fabs(ARRAY1[I])) {
	    maxd = fabs(ARRAY1[I]);
	}
    }
    while (IGAP > 0) {
	for (I = IGAP; I < N; I++) {
	    J = I - IGAP;
	    while (J >= 0) {
		v1 = fabs(ARRAY1[J]);
		v2 = fabs(ARRAY1[J + IGAP]);
		d1 = v1 - s;
		d2 = v2 - s;
		if (d1 < 0.0) {
		    d1 = maxd + v1;
		}
		if (d2 < 0.0) {
		    d2 = maxd + v2;
		}
		if (d1 > d2) {
		    TEMP = ARRAY1[J];
		    ARRAY1[J] = ARRAY1[J + IGAP];
		    ARRAY1[J + IGAP] = TEMP;
		    TEMP = ARRAY2[J];
		    ARRAY2[J] = ARRAY2[J + IGAP];
		    ARRAY2[J + IGAP] = TEMP;
		    J = J - IGAP;
		} else {
		    break;
		}
	    }
	}
	IGAP = IGAP / 2;
    }
/*
// .. end of dsort2ac_
*/
}

void dsort2sd(int N, double s, double *ARRAY1, double *ARRAY2)
{
/*
// Purpose
// =======
// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing order
// of abs(ARRAY1-s) if ARRAY1-s is non-positive. Positive ARRAY1-s are ordered 
// after those with non-positive ARRAY1-s and in the increasing order of -ARRAY1.
//
// Arguments
// =========
// N       (input) INTEGER
//          On entry, specifies the size of the arrays.
//
// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ..
// .. local scalars ..
*/
    int IGAP, I, J;
    double TEMP, v1, v2, d1, d2, maxd;
/*
// ..
// .. executable statements ..
*/
    IGAP = N / 2;
    maxd = fabs(ARRAY1[0]);
    for (I = 1; I < N; I++) {
	if (maxd < fabs(ARRAY1[I])) {
	    maxd = fabs(ARRAY1[I]);
	}
    }
    maxd = maxd + 1.0;
    while (IGAP > 0) {
	for (I = IGAP; I < N; I++) {
	    J = I - IGAP;
	    while (J >= 0) {
		v1 = fabs(ARRAY1[J]);
		v2 = fabs(ARRAY1[J + IGAP]);
		d1 = s - v1;
		d2 = s - v2;
		if (d1 < 0.0) {
		    d1 = s + (maxd - v1);
		}
		if (d2 < 0.0) {
		    d2 = s + (maxd - v2);
		}
		if (d1 > d2) {
		    TEMP = ARRAY1[J];
		    ARRAY1[J] = ARRAY1[J + IGAP];
		    ARRAY1[J + IGAP] = TEMP;
		    TEMP = ARRAY2[J];
		    ARRAY2[J] = ARRAY2[J + IGAP];
		    ARRAY2[J + IGAP] = TEMP;
		    J = J - IGAP;
		} else {
		    break;
		}
	    }
	}
	IGAP = IGAP / 2;
    }
/*
// .. end of dsort2ac_
*/
}
