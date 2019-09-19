/*
// ZTRLan routine 
// Lawrence Berkeley National Lab.
//
*/

#include <math.h>
#include <stdio.h>

#include "trl_map.h"
#include "trlan.h"
#include "trlan_i.h"
#include "dstqrb_i.h"

static inline double d_sign(doublereal * a, doublereal * b)
{
    double x;
    x = (*a >= 0 ? *a : -*a);
    return (*b >= 0 ? x : -x);
}

int dstqrb_(integer * n, doublereal * d__, doublereal * e,
	    doublereal * z__, doublereal * work, integer * info)
{
/*

    Purpose
    =======
    DSTQRB computes all eigenvalues and the last component of the eigenvectors of a
    symmetric tridiagonal matrix using the implicit QL or QR method.
    This is mainly a modification of the CLAPACK subroutine dsteqr.c

    Arguments
    =========

    N       (input) INTEGER
            The order of the matrix.  N >= 0.

    D       (input/output) DOUBLE PRECISION array, dimension (N)
            On entry, the diagonal elements of the tridiagonal matrix.
            On exit, if INFO = 0, the eigenvalues in ascending order.

    E       (input/output) DOUBLE PRECISION array, dimension (N-1)
            On entry, the (n-1) subdiagonal elements of the tridiagonal
            matrix.
            On exit, E has been destroyed.

    Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
            On entry, if  COMPZ = 'V', then Z contains the orthogonal
            matrix used in the reduction to tridiagonal form.
            On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
            orthonormal eigenvectors of the original symmetric matrix,
            and if COMPZ = 'I', Z contains the orthonormal eigenvectors
            of the symmetric tridiagonal matrix.
            If COMPZ = 'N', then Z is not referenced.

    WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
            If COMPZ = 'N', then WORK is not referenced.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
            > 0:  the algorithm has failed to find all the eigenvalues in
                  a total of 30*N iterations; if INFO = i, then i
                  elements of E have not converged to zero; on exit, D
                  and E contain the elements of a symmetric tridiagonal
                  matrix which is orthogonally similar to the original
                  matrix.

    =====================================================================
*/
    /* Table of constant values */
    doublereal c_b10 = 1.;
    integer c__0 = 0;
    integer c__1 = 1;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;
    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);
    extern /* Subroutine */ int dlae2_(doublereal *, doublereal *, doublereal
				       *, doublereal *, doublereal *);
    doublereal b, c__, f, g;
    integer i__, j, k, l, m;
    doublereal p, r__, s;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dlasr_(char *, char *, char *, integer *,
				       integer *, doublereal *,
				       doublereal *, doublereal *,
				       integer *);
    doublereal anorm;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *,
				       doublereal *, integer *);
    integer l1;
    extern /* Subroutine */ int dlaev2_(doublereal *, doublereal *,
					doublereal *, doublereal *,
					doublereal *, doublereal *,
					doublereal *);
    integer lendm1, lendp1;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    integer ii;
    extern doublereal dlamch_(char *);
    integer mm, iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *,
					doublereal *, doublereal *,
					integer *, integer *, doublereal *,
					integer *, integer *),
	dlaset_(char *, integer *, integer *, doublereal *, doublereal *,
		doublereal *, integer *);
    doublereal safmin;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *,
					doublereal *, doublereal *,
					doublereal *);
    doublereal safmax;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern doublereal dlanst_(char *, integer *, doublereal *,
			      doublereal *);
    extern /* Subroutine */ int dlasrt_(char *, integer *, doublereal *,
					integer *);
    /* Local variables */
    integer lend, jtot;
    integer lendsv;
    doublereal ssfmin;
    integer nmaxit;
    doublereal ssfmax;
    integer lm1, mm1, nm1;
    doublereal rt1, rt2, eps;
    integer lsv;
    doublereal tst, eps2;

    --d__;
    --e;
    --z__;
    /* z_dim1 = *ldz;             */
    /* z_offset = 1 + z_dim1 * 1; */
    /* z__ -= z_offset;           */
    --work;

    /* Function Body */
    *info = 0;
    /* Taken out for TRLan
       if (lsame_(compz, "N")) {
       icompz = 0;
       } else if (lsame_(compz, "V")) {
       icompz = 1;
       } else if (lsame_(compz, "I")) {
       icompz = 2;
       } else {
       icompz = -1;
       }
       if (icompz < 0) {
       *info = -1;
       } else if (*n < 0) {
       *info = -2;
       } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
       *info = -6;
       }
       if (*info != 0) {
       i__1 = -(*info);
       xerbla_("DSTEQR", &i__1);
       return 0;
       }
    */
    /*  icompz = 2; */

    /*	Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (*n == 1) {
	z__[1] = 1;
	return 0;
    }

    /*	Determine the unit roundoff and over/underflow thresholds. */

    eps = dlamch_("E");
    /*	Computing 2nd power */
    d__1 = eps;
    eps2 = d__1 * d__1;
    safmin = dlamch_("S");
    safmax = 1. / safmin;
    ssfmax = sqrt(safmax) / 3.;
    ssfmin = sqrt(safmin) / eps2;

    /*	Compute the eigenvalues and eigenvectors of the tridiagonal
	matrix. */
    /* Taken out for TRLan
       if (icompz == 2) {
       dlaset_("Full", n, n, &c_b9, &c_b10, &z__[z_offset], ldz);
       }
    */
    for (j = 1; j < *n; j++) {
	z__[j] = 0.0;
    }
    z__[*n] = 1.0;
    nmaxit = *n * 30;
    jtot = 0;

    /*	Determine where the matrix splits and choose QL or QR iteration
	for each block, according to whether top or bottom diagonal
	element is smaller. */

    l1 = 1;
    nm1 = *n - 1;

 L10:
    if (l1 > *n) {
	goto L160;
    }
    if (l1 > 1) {
	e[l1 - 1] = 0.;
    }
    if (l1 <= nm1) {
	i__1 = nm1;
	for (m = l1; m <= i__1; ++m) {
	    tst = (d__1 = e[m], fabs(d__1));
	    if (tst == 0.) {
		goto L30;
	    }
	    if (tst <=
		sqrt((d__1 = d__[m], fabs(d__1))) * sqrt((d__2 =
							  d__[m + 1],
							  fabs(d__2))) *
		eps) {
		e[m] = 0.;
		goto L30;
	    }
	    /* L20: */
	}
    }
    m = *n;

 L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

    /*	Scale submatrix in rows and columns L to LEND */

    i__1 = lend - l + 1;
    anorm = dlanst_("I", &i__1, &d__[l], &e[l]);
    iscale = 0;
    if (anorm == 0.) {
	goto L10;
    }
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l],
		n, info);
	i__1 = lend - l;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n,
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l],
		n, info);
	i__1 = lend - l;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n,
		info);
    }

    /*	Choose between QL and QR iteration */

    if ((d__1 = d__[lend], fabs(d__1)) < (d__2 = d__[l], fabs(d__2))) {
	lend = lsv;
	l = lendsv;
    }

    if (lend > l) {

	/*	QL Iteration

		Look for small subdiagonal element. */

    L40:

	if (l != lend) {
	    lendm1 = lend - 1;
	    i__1 = lendm1;
	    for (m = l; m <= i__1; ++m) {
		/*		Computing 2nd power */
		d__2 = (d__1 = e[m], fabs(d__1));
		tst = d__2 * d__2;
		if (tst <=
		    eps2 * (d__1 = d__[m], fabs(d__1)) * (d__2 =
							  d__[m + 1],
							  fabs(d__2)) +
		    safmin) {
		    goto L60;
		}
		/* L50: */
	    }
	}

	m = lend;

    L60:
	if (m < lend) {
	    e[m] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L80;
	}

	/*	If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
		to compute its eigensystem. */

	if (m == l + 1) {
	    dlaev2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
	    work[l] = c__;
	    work[*n - 1 + l] = s;
	    /* Taken out for TRLan
	       dlasr_("R", "V", "B", n, &c__2, &work[l], &work[*n - 1 + l], &
	       z___ref(1, l), ldz);
	    */
	    tst = z__[l + 1];
	    z__[l + 1] = c__ * tst - s * z__[l];
	    z__[l] = s * tst + c__ * z__[l];
	    d__[l] = rt1;
	    d__[l + 1] = rt2;
	    e[l] = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L40;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

	/*	Form shift. */

	g = (d__[l + 1] - p) / (e[l] * 2.);
	r__ = dlapy2_(&g, &c_b10);
	g = d__[m] - p + e[l] / (g + d_sign(&r__, &g));

	s = 1.;
	c__ = 1.;
	p = 0.;

	/*	Inner loop */

	mm1 = m - 1;
	i__1 = l;
	for (i__ = mm1; i__ >= i__1; --i__) {
	    f = s * e[i__];
	    b = c__ * e[i__];
	    dlartg_(&g, &f, &c__, &s, &r__);
	    if (i__ != m - 1) {
		e[i__ + 1] = r__;
	    }
	    g = d__[i__ + 1] - p;
	    r__ = (d__[i__] - g) * s + c__ * 2. * b;
	    p = s * r__;
	    d__[i__ + 1] = g + p;
	    g = c__ * r__ - b;

	    /*		If eigenvectors are desired, then save rotations. */

	    work[i__] = c__;
	    work[*n - 1 + i__] = -s;

	    /* L70: */
	}

	/*	If eigenvectors are desired, then apply saved rotations. */

	mm = m - l + 1;
	/* Taken out for TRLan
	   dlasr_("R", "V", "B", n, &mm, &work[l], &work[*n - 1 + l], &
	   z___ref(1, l), ldz);
	*/
	dlasr_("R", "V", "B", &c__1, &mm, &work[l], &work[*n - 1 + l],
	       &z__[l], &c__1);

	d__[l] -= p;
	e[l] = g;
	goto L40;

	/* 	Eigenvalue found. */

    L80:
	d__[l] = p;

	++l;
	if (l <= lend) {
	    goto L40;
	}
	goto L140;

    } else {

	/*	QR Iteration

		Look for small superdiagonal element. */

    L90:
	if (l != lend) {
	    lendp1 = lend + 1;
	    i__1 = lendp1;
	    for (m = l; m >= i__1; --m) {
		/*			Computing 2nd power */
		d__2 = (d__1 = e[m - 1], fabs(d__1));
		tst = d__2 * d__2;
		if (tst <=
		    eps2 * (d__1 = d__[m], fabs(d__1)) * (d__2 =
							  d__[m - 1],
							  fabs(d__2)) +
		    safmin) {
		    goto L110;
		}
		/* L100: */
	    }
	}

	m = lend;

    L110:
	if (m > lend) {
	    e[m - 1] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L130;
	}

	/*	If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
		to compute its eigensystem. */

	if (m == l - 1) {
	    dlaev2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s);
	    /* Taken out for TRLan
	       work[m] = c__;
	       work[*n - 1 + m] = s;
	       dlasr_("R", "V", "F", n, &c__2, &work[m], &work[*n - 1 + m], &
	       z___ref(1, l - 1), ldz);
	    */
	    tst = z__[l];
	    z__[l] = c__ * tst - s * z__[l - 1];
	    z__[l - 1] = s * tst + c__ * z__[l - 1];

	    d__[l - 1] = rt1;
	    d__[l] = rt2;
	    e[l - 1] = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L90;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

	/*	Form shift. */

	g = (d__[l - 1] - p) / (e[l - 1] * 2.);
	r__ = dlapy2_(&g, &c_b10);
	g = d__[m] - p + e[l - 1] / (g + d_sign(&r__, &g));

	s = 1.;
	c__ = 1.;
	p = 0.;

	/*	Inner loop */

	lm1 = l - 1;
	i__1 = lm1;
	for (i__ = m; i__ <= i__1; ++i__) {
	    f = s * e[i__];
	    b = c__ * e[i__];
	    dlartg_(&g, &f, &c__, &s, &r__);
	    if (i__ != m) {
		e[i__ - 1] = r__;
	    }
	    g = d__[i__] - p;
	    r__ = (d__[i__ + 1] - g) * s + c__ * 2. * b;
	    p = s * r__;
	    d__[i__] = g + p;
	    g = c__ * r__ - b;

	    /*		If eigenvectors are desired, then save rotations. */

	    work[i__] = c__;
	    work[*n - 1 + i__] = s;

	    /* L120: */
	}

	/*	If eigenvectors are desired, then apply saved rotations. */

	mm = l - m + 1;
	/*
	  dlasr_("R", "V", "F", n, &mm, &work[m], &work[*n - 1 + m], &
	  z___ref(1, m), ldz);
	*/
	dlasr_("R", "V", "F", &c__1, &mm, &work[m], &work[*n - 1 + m],
	       &z__[m], &c__1);

	d__[l] -= p;
	e[lm1] = g;
	goto L90;

	/*        Eigenvalue found. */

    L130:
	d__[l] = p;

	--l;
	if (l >= lend) {
	    goto L90;
	}
	goto L140;

    }

    /*     Undo scaling if necessary */

 L140:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1,
		&d__[lsv], n, info);
	i__1 = lendsv - lsv;
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv],
		n, info);
    } else if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1,
		&d__[lsv], n, info);
	i__1 = lendsv - lsv;
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv],
		n, info);
    }

    /*     Check for no convergence to an eigenvalue after a total
	   of N*MAXIT iterations. */

    if (jtot < nmaxit) {
	goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (e[i__] != 0.) {
	    ++(*info);
	}
	/* L150: */
    }
    goto L190;

    /*     Order eigenvalues and eigenvectors. */

 L160:

    /*        Use Selection Sort to minimize swaps of eigenvectors */

    i__1 = *n;
    for (ii = 2; ii <= i__1; ++ii) {
	i__ = ii - 1;
	k = i__;
	p = d__[i__];
	i__2 = *n;
	for (j = ii; j <= i__2; ++j) {
	    if (d__[j] < p) {
		k = j;
		p = d__[j];
	    }
	    /* L170: */
	}
	if (k != i__) {
	    d__[k] = d__[i__];
	    d__[i__] = p;
	    /* Taken out for TRLan
	       dswap_(n, &z___ref(1, i__), &c__1, &z___ref(1, k), &c__1);
	    */
	    p = z__[k];
	    z__[k] = z__[i__];
	    z__[i__] = p;
	}
	/* L180: */
    }

 L190:
    return 0;
}				/* dstqrb_ */
