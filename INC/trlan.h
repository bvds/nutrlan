#ifndef __TRLAN_H
#define __TRLAN_H
/**
   @file The public header file for C version of the nuTRLan code.
*/

#include <time.h>	/* clock_t */

#ifdef TRL_FORTRAN_COMPATIBLE
/**
   Prototype matrix-vector multiplication function.  Useful if you have a
   matrix-vector multiplication function defined in Fortran.

   @arg pnrow Pointer to integer value nrow, the number of rows locally on
   this processor.

   @arg pncol Pointer to integer value ncol, the number of columns in
   vectors x and y.

   @arg x Pointer to the elements of input vectors x.  It is assumed to
   have [*pldx * *pncol] elements, with the ith column starting at element
   *pldx * i.

   @arg pldx Pointer to integer value ldx, the leading dimension of x.
   Note that *pldx must not be less than *pnrow.

   @arg y Pointer to the elements of output vectors y.  It is assumed to
   have [*pldy * *pncol] elements, with the ith column starting at element
   *pldy * i.

   @arg pldy Pointer to integer value ldy, the leadying dimension of y.
   Note that *pldy must not be less than *pnrow.
*/
typedef void (*trl_matvec) (int *pnrow, int *pncol, double *x, int *pldx,
			    double *y, int *pldy);
#else
/**
   Prototype matrix-vector multiplication function.  Defined to be easier
   to use with C/C++ code.

   @arg nrow The number of rows locally on this processor.

   @arg ncol The number of columns in vectors x and y.

   @arg x Pointer to the elements of input vectors x.  It is assumed to
   have [ldx * ncol] elements, with the ith column starting at element
   ldx * i.

   @arg ldx The leading dimension of x, the ith column of x starts at
   element ldx * i.
   Note that ldx must not be less than nrow.

   @arg y Pointer to the elements of output vectors y.  It is assumed to
   have [ldy * ncol] elements, with the ith column starting at element
   ldy * i.

   @arg ldy The leadying dimension of y.  The ith column of y starts at
   position ldy * i.
   Note that ldy must not be less than nrow.

   @arg mvparam The extra parameter to be passed to the matrix-vector
   multiplication through trl_info.  This parameter mvparam is used
   exclusively in the matrix-vector multiplication function and nowhere
   else in TRLan.
*/
typedef void (*trl_matvec) (const int nrow, const int ncol,
			    const double *x, const int ldx,
			    double *y, const int ldy, void *mvparam);
#endif

/**
   The data structure to store the current information about
   the eigenvalue problem and the progress of TRLAN.
*/
typedef struct strct_trl_info {
    int stat;			/* status  (error code) of TRLAN */

    /* specification of the eigenvalue problem */
    /** Which end of spectrum to compute.
	- lohi < 0 --> the smallest eigenvalues
	- lohi = 0 --> whichever converge first
	- lohi > 0 --> the largest eigenvalues
    */
    int lohi;
    /** Number of eigenpairs wanted.  */
    int ned;
    /** Number of eigenpairs converged.  if the user has nec correct
    eigenvectors on input, then they are expected to be stored at the
    beginning of the eigenvector array.  */
    int nec;
    /* Convergence tolerance.  An eigenpair is declared converged if its
       residual norm is less than tol*||OP||.  */
    double tol;

    /* specification of resource allowed to use by TRLAN */
    /** The MPI communicator.  */
    int mpicom;
    /** The maximum basis size to be used.  */
    int maxlan;
    /** The actual basis size currently. This value may be smaller than
    maxlan.  It is set during restarting.  */
    int klan;
    /** The maximum number of MATVEC allowed. Note one MATVEC == one
	application of the operator on one vector.  */
    int maxmv;

    /** The restarting scheme to use.  */
    int restart;
    /** The number of eigenvalues locked.  */
    int locked;
    /** Option for handling initial guesses:
       - <= 0, user did not provide initial guess, use    
           [1,1,..,1].
       -  = 1, user has supplied initial guess, will only 
           use the first one.
       -  > 1, restart with previous check-point file.      */
    int guess;

    /* some information about the progress and resouce comsumption    */
    /** The number of MATVEC used by TRLAN.                    */
    int matvec;
    /** The number of restart of the Lanczos iterations.       */
    int nloop;
    /** The number of full orthogonalization invoked.          */
    int north;
    /** The number of times a random element is introduced. Random elements
	are introduced when an invariant subspace is found, but the number
	of converged eigen-pairs is less than desired.  */
    int nrand;
    /** Floating-point operations count (EXCLUDING MATVEC). */
    int flop;
    /** FLOPS used in re-orthogonalization.                */
    int flop_h;
    /** FLOPS used in restarting.                          */
    int flop_r;
    double rflp;
    double rflp_h;
    double rflp_r;

    /* variables to store timing results */
    /** system clock rate (SYSTEM_CLOCK)                */
    clock_t clk_rate;
    /** Maximum counter value                           */
    clock_t clk_max;
    /** Total time spent in TRLAN (in clock ticks)      */
    clock_t clk_tot;
    /** Time in applying the operator (MATVEC)          */
    clock_t clk_op;
    /** Time in re-orthogonalization                    */
    clock_t clk_orth;
    /** Time in restarting the Lanczos iterations       */
    clock_t clk_res;
    /** The sum of clk_tot and tick_t is the actual time */
    double tick_t;
    double tick_o;
    double tick_h;
    double tick_r;
    /** Time spent in reading input data file           */
    int clk_in;
    /** Number of real(8) words read                    */
    int wrds_in;
    /** Time spent in writing output data file          */
    int clk_out;
    /** Number of real(8) words written to file         */
    int wrds_out;

    /** Norm of the operator used.  This is an estimate based on the
	largest absolute value of a Rayleigh quotient.  */
    double anrm;

    /** The PE number of current processor (start with 0).  */
    int my_pe;
    /** number of PEs in the group                      */
    int npes;
    /** Local problem size                              */
    int nloc;
    /** Global problem size                             */
    int ntot;

    /** How much inforation to output during the execution of TRLAN.  By
       default, it only print information related to fatal errors.  If
       verbose > 0, more diagnostic messages are printed.  */
    int verbose;
    /** Variables needed to measure convergence factor (crat).  The
       convergence rate of the restarted Lanczos algorithm is measured by
       the reduction of residual norm per MATVEC.  The residual norm of the
       target is used.  */
    double crat;
    /** The Ritz value that might convege next.  */
    double trgt;
    /** The residual norm of the target.  */
    double tres;
    /** MATVEC used when target and tres were recorded  */
    int tmv;
    double avgm;

    /* Stores some convergence history for restart scheme 9           */
    double cfac;
    double ptres;
    double predicted_crate;
    double mgap;
    double mgamma, gamma0;
    double old_target;
    int target_id;
    int old_locked;
    int k1, k2, k;
    double rfact;

    /* Store "shift" */
    double ref;

    /* Fortran I/O unit number to be used for          */
    /* debug output. Used if verbose > 0.              */
    int log_io;
    FILE *log_fp;
    /** Base of the file names used by TRLAN to store debug info if verbose
     > 0, the filenames are computed by appending 'PE#' to this base.  */
    char log_file[128];

    /** check-pointing parameters.
       when cpflag is greater than 0, the basis vectors will be written
       out roughly 'cpflag' times.  For simplicitly, each PE writes its
       own portion of the basis vectors to a file with cpfile followed by
       the processor number.  The file is written as unformatted fortran
       files with the following content:
       @pre
       nrow, kb(basis size)
       alpha(1:kb)
       beta(1:kb)
       1st basis vector, 2nd basis vector, ..., last basis vector
       the residual vector
       @endpre
     */
    int cpflag, cpio;
    FILE *cpt_fp;
    char cpfile[128], oldcpf[128];
    /**
       This parameter is only used by the matrix-vector multiplication
       function when the conditional macro TRL_FORTRAN_COMPATIBLE is not
       defined.  It is used nowhere else in TRLan.
     */
    void* mvparam;
} trl_info;

/**
   The main user function for computing eigen-pairs.

 A thick-restart Lanczos routine for computing eigenvalues and
 eigenvectors of a real symmetric operator/matrix (A).
 -- It only accept one input vector, the input vector is expected
    to be stored in the (nec+1)st column of EVEC.
 -- It extends the Lanczos basis one vector at a time.
 -- Orthogonality among the Lanczos vectors are maintained using
    full re-orthogonalization when necessary.
 -- User must initialize the trl_info object passed to this function with
 valid information.

 Requirements:
 1) User supplies OP with the specified interface.
 2) If (info%nec>0), evec(1:nrow, 1:info%nec) must contain valid
    eigenvectors and eval(1:nec) must be the corresponding eigenvalues.
    These eigenpairs are assumed to have zero residual norms inside
    TRLAN.
 3) lde >= nrow.
 4) The arrays evec and eval must be large enough to store the
    solutions, i.e., mev >= info%ned and mev >= info%nec.
 5) The array wrk may be of arbitrary size.  Internally, the workspace
    size is
        nrow*max(0,info%ned-size(evec,2))+maxlan*(maxlan+10)

    If wrk is smaller than this, trlan routine will allocate additional
    workspace to accommodate.

@arg op    (input) function pointer to the matrix-vector multiplication routine.
         It points to a function that comptues op(X) == A*X,
         when given a set of vectors X.

@arg info    (input) [pointer to the structure trl_info]
          It points to the data structure to store the information
          about the eigenvalue problem and the progress of TRLAN

@arg nrow    (input) [integer]
          It specifies the number of rows that is on this processor.

@arg mev     (input) [integer]
          It specifies the number of Ritz pairs, that can be stored in
          eval and evec.

@arg eval    (output) [double precision vector of length (mev)]
          On exist, stores the eigenvalues.

@arg evec    (output) [double precision vector of length (nrow*mev)]
          On exit, stores the eigenvectors.

@arg lde     (input) [integer]
          On entry, specifies the leading dimension of the array evec.  The
          ith vector starts at element lde*i.

@arg lwrk (input) [integer] It specifies the size of WRK.  It should
          correctly indicate the size of wrk.

@arg wrk (workspace) If is enough space, the residual norm of the converged
          eigenpairs will be stored at wrk(1:info%nec) on exit.

*/
void trlan(trl_matvec op,
	   trl_info * info, int nrow, int mev, double *eval,
	   double *evec, int lde, int lwrk, double *wrk );
/**
 Check the validity of the computed Ritz pairs.

@arg op       (input) [function pointer]
           It points to the matrix-vector multiplication routine.

@arg info     (input) [pointer to the structure trl_info]
           It points to the data structure to store the information
           about the eigenvalue problem and the progress of TRLAN.

@arg nrow     (input) [integer]
           It specifies the local problem size, the number of rows
           in this processor.

@arg ncol     (input) [integer]
           It specifies the number of Ritz values computed.

@arg rvec     (input) [double precision array of dimension (ldrvec*ncol)]
           It specifies the array storing the Ritz vectors.

@arg alpha    (input) [double precision array of dimension (ncol)]
           On entry, contains the Ritz values computed.

@arg beta     (input) [double precision array of dimension (ncol)]
           It contaions the residual norms returned from a Lanczos routine.

@arg eval     (input) [double precision array of dimension (ncol)]
           It contains the actual eigenvalues to compute the error in
           Ritz values.

@arg lwrk     (input) [integer]
           It specifies the size of workspace provided.

@arg wrk      (workspace) [double precision array of size(lwrk)]
*/
void trl_check_ritz(trl_matvec op,
		    trl_info * info, int nrow, int ncol, double *rvec,
		    int ldrvec, double *alpha, int *check, double *beta,
		    double *eval, double *wrk, int lwrk);

/**
 Initializes a TRL_INFO variable.  This function must be called before calling
 any other user level routine in TRLAN package.

@arg info    (input) [pointer to the structure trl_info]
          It points to the data structure to store the information
          about the eigenvalue problem and the progress of TRLAN.

@arg nrow    (input) [integer]
          It specifies the local dimension of the problem.

@arg mxlan   (input) [integer]
          It specifies the maximum number of basis vectors to be used.

@arg lohi    (input)  [integer]
          It specifies, which end of the spectrum to compute:
          - lohi < 0, then lower end, the smallest eigenvalues
          - lohi > 0, then high end, the largest eigenvalues
          - lohi = 0, then either lower and or high end, whoever converges first
          are computed.

@arg ned      (input) [integer]
           It specifies the number of wanted eigenvalues and
           eigenvectors.

@arg tol      (input) double precision
           If provided, specifies the tolerance on residual norm. By default,
           tol is set to be sqrt(epsilon).

@arg trestart (input) [integer]
           If it is an integer between 1 and 8, it specifies the
           thick-restarting scheme, otherwise the default restarting scheme
           is 1.

@arg mxmv     (input) [integer]
           If input is a positive number, it specifies the maximum number
           of matrix-vector multiplication allowed.  By default, mxmv is
           set to be info%ntot*info%ned.

@arg mpicom   (input) [integer]
           It specifites the MPI communicator.  If MPI is used, this has to
           be a valid MPI communicator.  In the sequential case, this
           parameter is never really used (even though passed around).
*/
void trl_init_info(trl_info *info, int nrow, int mxlan, int lohi,
		   int ned, double tol, int restart, int maxmv,
		   int mpicom);

/**
 Set the (minimum) basis size for the restart schemes 7 and 8, i.e., the
 (minimum) basis size if rfact * (number of Ritz vectors kept)

@arg info    (input/output) [pointer to the structure trl_info]
          It points to the data structure to store the current
          information about the eigenvalue problem and the progress of
          TRLAN.

@arg rfact   (input) double precision
          On entry, specify the (minimum) basis size.
*/
void trl_set_restart(trl_info * info, double rfact);

/**
 Set information related to debugging.  The initialization routine
 trl_init_info sets the parameters so that no debug information is printed.

@arg info    (input/output) [pointer to the structure trl_info]
          It points to the data structure to store the current information
          about the eigenvalue problem and the progress of TRLAN.

@arg msglvl  (input) [integer]
          It specifies the new message/verbose level:
          - msglvl <  0         : nothing is printed
          - msglvl = 1, .., 10  : the higher the level, the more debug
                                 information is printed.

@arg file    (input) [character string]
          It specifies the leading part of the debug file name.  Each
          processor will generate its own debug file with the file name
          formed by appending the processor number to the string FILE. The
          actual name is composed by the routine TRL_PE_FILENAME in trlaux.

*/
void trl_set_debug(trl_info * info, int msglvl, char *filename);

/**
 Set up the information related to check-pointing

@arg info    (input/output) [pointer to the structure trl_info]
          It points to the data structure to store the current information
          about the eigenvalue problem and the progress of TRLAN.

@arg cpflag  (input) [integer]
          It spcifies roughly how many titmes checkpoints are written.

@arg file    (input) [character string]
          It specifies the leading part of the checkpoint file name. Each
          processor will generate its own debug file with the file name
          formed by appending the processor number to the string FILE. The
          actual name is composed by the routine TRL_PE_FILENAME in trlaux.
*/
void trl_set_checkpoint(trl_info * info, int cpflag, char *file);

/**
   Set up parameters related to initial guesses of the Lanczos iterations.
   It specifies the number of eigenvector already converged (initially
   assumed to be zero) and whether the user has provided initial guess
   vector (initially assumed no).  It can also tell TRLan to read a check
   point file that is different from the default name.

@arg info    (input/output) [pointer to the structure trl_info]
          It points to the data structure to store the current information
          about the eigenvalue problem and the progress of TRLAN.

@arg nec     (input) [integer]
          It specifies the number of eigenvalues, that have been converged.

@arg iguess  (input) [integer]
          It specifies one of the following options:
          - iguess <= 0, user did not provide initial guess, use [1,1,..,1].
          - iguess =  1, user has supplied initial guess, will only use the
            first one.
          - iguess >  1, restart with previous check-point file.

@arg ncps    (input) [integer]
          It specifies the number of time a checkpoint file will be
          written.

@arg cpf  (input) [string]
          It is the prefix of names of the checkpoint files.  Each MPI
          process will generate its own checkpoint file with its processor
          id appended to this prefix.
*/
void trl_set_iguess(trl_info * info, int nec, int iguess, int ncps,
		    char *cpf );

/**
   Provides an uniform way of printing information stored in TRL_INFO_T.
   It needs to be called by all PEs. In parallel environment, when writing
   to standard outputd device, only PE0 will actually write its local
   summary information.  NOTE that this function must be called on all PEs
   at the same time.

@arg info    (input) [pointer to the structure trl_info]
          It points to the data structure to store the current information
          about the eigenvalue problem and the progress of TRLAN.

@arg mvflop  (input) [integer]
          It specifies the number of floating-operations per MATVEC per PE.
          This information has to be supplied by user, otherwise related
          entries are left blank in the print out.
*/
void trl_print_info(trl_info * info, int mvflop);

/**
 It is a more compact version of trl_print_info.  This is a local routine.
 Indivadual PE can call it without regard of whether other PEs do the same
 and the output may be written to a different I/O unit number than the
 log_io.

@arg info    (input) [pointer to the structure trl_info]
          It points to the data structure to store the current information
          about the eigenvalue problem and the progress of TRLAN.

@arg iou     (input) Pointer to a file to write information from this PE.
*/
void trl_terse_info(trl_info * info, FILE * iou);

/**
 A separate Rayleigh-Ritz projection routine.
 Given a set of approximately orthonormal vectors (V), this routine
 performs the following operations
  (1) V'*V ==> G
  (2) R'*R :=  G
  (3) V'*A*V => H1, inv(R')*H1*inv(R) => H
  (4) Y*D*Y' := H
  (5) V*inv(R)*Y => V, diag(D) => lambda,
      r(i) = ||A*V(:,i)-lambda(i)*V(:,i)||

@arg op       (input) [function pointer]
          It points to the matrix-vector multiplication routine.

@arg info    (input) [pointer to the structure trl_info]
          It points to the data structure to store the information
          about the eigenvalue problem and the progress of TRLAN.

@arg mev     (input) [integer]  The number of eigen-pairs to check.

@arg evec    (output) [double precision vector of lenvth (lde*mev)]
          On exit, stores the eigenvectors.

@arg lde     (input) integer
          On entry, specifies the leading dimension of the array evec.  The
          ith vector in evec is assumed to start at element lde*i.

@arg eres    (output) [double precision vector of length (2*nev)]
          The array to store new Ritz values and residual norms.

@arg wrk     (worksapce) [double precision vector of length (lwrk)]
          Workspace to store projection matrix, etc.

@arg lwrk    (input) [integer] The size of array wrk.

@arg base    (workspace) [double precision vector]
          Workspace to store the result of matrix-vector operation.
*/
void
trl_ritz_projection(trl_matvec op, trl_info *info, int mev,
		    double *evec, int lde, double *eres, double *wrk, int lwrk,
		    double *base);

/**
   Compute Rayleigh quotients.  When it is given a set of Ritz vectors and
 Ritz values, normalize the Ritz vectors, and compute their Rayleigh
 quotients to replace the existing Ritz values.

@arg op       (input) [function pointer]
           It is the matrix-vector multiplication routine.

@arg info     (input) [pointer to the structure trl_info]
           It points to the data structure to store the current information
           about the eigenvalue problem and the progress of TRLAN.

@arg evec     (input) [double precision array of dimension (lde*ncol)]
           It stores the portion of eigenvectors on this PE.

@arg lde      (input) [integer]
           The leading dimension of array evec.

@arg base     (workspace)
           The workspace used to store results of MATVEC

@arg eres     (output) [double precision array of dimension (ncol)]
           On exist, store new Ritz values and new residual norms, i.e., if
           there are NEV Ritz pairs, eres(1:NEV) stores the new Rayleigh
           quotient and eres(nev+1:nev+nev) stores the new residual norms.

@arg base     (wordspace) [double precision array of dimension (info->nloc)]
           If provided, double precision workspace.
*/
void
trl_rayleigh_quotients(trl_matvec op, trl_info * info, int ncol, double *evec,
		       int lde, double *eres, double *base);

#endif
