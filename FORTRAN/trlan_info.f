      Module trlan_info
      Type TRLAN_INFO_T
      Integer :: stat    ! status  (error code) of TRLAN

      ! specification of eigenpairs wanted
      Integer :: lohi    ! which end of spectrum to compute
                         ! lohi < 0 --> the smallest eigenvalues
                         ! lohi = 0 --> whichever converge first
                         ! lohi > 0 --> the largest eigenvalues
      Integer :: ned     ! number of eigenpairs wanted
      Integer :: nec     ! number of eigenpairs converged
                         ! if the user has nec correct eigenvectors, then
                         ! they are expected to be stored at the beginning
                         ! of the eigenvector array
      Real(8) :: tol     ! an eigenpair is declared converged if its
                         ! residual norm is less than tol*||OP||

      ! specification of resource allowed to use by TRLAN
      Integer :: mpicom  ! the MPI communicator
      Integer :: maxlan  ! maximum basis size to be used
      Integer :: klan    ! the actual basis size, this value may be smaller
                         ! than maxlan.  It is set when restarting.
      Integer :: maxmv   ! maximum number of MATVEC allowed
                         ! one MATVEC == one application of the operator on
                         ! one vector
 
      Integer :: restart ! index of restarting schemes
      Integer :: locked  ! number of eigenvalue locked
      Integer :: guess   ! initial guess
      ! <= 0, user did not provide initial guess, use [1,1,..,1]
      ! =  1, user has supplied initial guess, will only use the first one
      ! >  1, restart with previous check-point file
 
      ! some information about the progress and resouce comsumption
      Integer :: matvec  ! number of MATVEC used by TRLAN
      Integer :: nloop   ! number of restart of the Lanczos iterations
      Integer :: north   ! number of full orthogonalization invoked
      Integer :: nrand   ! number of times random elements are introduced.
                         ! Random elements are introduced when an invariant
                         ! subspace is found but not all wanted eigenvalues
                         ! are computed.
      Integer :: flop    ! Floating-point operations count (EXCLUDING MATVEC)
      Integer :: flop_h  ! FLOPS used in re-orthogonalization
      Integer :: flop_r  ! FLOPS used in restarting
      Real(8) :: rflp
      Real(8) :: rflp_h
      Real(8) :: rflp_r
 
      ! variables to store timing results
      Integer :: clk_rate! system clock rate (SYSTEM_CLOCK)
      Integer :: clk_max ! maximum counter value
      Integer :: clk_tot ! total time spent in TRLAN (in clock ticks)
      Integer :: clk_op  ! time in applying the operator (MATVEC)
      Integer :: clk_orth! time in re-orthogonalization
      Integer :: clk_res ! time in restarting the Lanczos iterations
      Real(8) :: tick_t  ! the sum of clk_tot and tick_t is the actual time
      Real(8) :: tick_o
      Real(8) :: tick_h
      Real(8) :: tick_r
      Integer :: clk_in  ! time spent in reading input data file
      Integer :: wrds_in ! number of real(8) words read
      Integer :: clk_out ! time spent in writing output data file
      Integer :: wrds_out! number of real(8) words written to file
 
      Real(8) :: anrm    ! norm of the operator used
 
      Integer :: my_pe   ! the PE number of current processor (start with 0)
      Integer :: npes    ! number of PEs in the group
      Integer :: nloc    ! local problem size
      Integer :: ntot    ! global problem size
 
      ! how much inforation to output during the execution of TRLAN
      Integer :: verbose ! default only print information related to
                         ! fatal errors
                         ! if verbose > 0, more diagnostic messages
                         ! are printed to the following files
      ! variables needed to measure convergence factor (crat)
      ! convergence rate of the restarted Lanczos algorithm is measure by
      ! the reduction of residual norm per MATVEC.  The residual norm of
      ! the target is used.
      Real(8) :: crat
      Real(8) :: trgt    ! the Ritz value that might convege next
      Real(8) :: tres    ! residual norm of the target
      Integer :: tmv     ! MATVEC used when target and tres were recorded
      Real(8) :: avgm

      ! Stores some convergence history for restart scheme 9
      Real(8) cfac
      Real(8) ptres
      Real(8) predicted_crate
      Real(8) mgap
      Real(8) mgamma, gamma0
      Real(8) old_target
      Integer target_id
      Integer old_locked
      Integer k1, k2, k
      Real(8) rfact
 
      ! Store "shift" */
      Real(8) ref;
 
      Integer :: log_io  ! Fortran I/O unit number to be used for
                         ! debug output. Used if verbose > 0.
      Integer :: lop_fp
      Character(128) :: log_file
      ! base of the file names used by TRLAN to store debug information if
      ! verbose > 0.  The filenames are computed by appending 'PE#' to
      ! this basis.
 
      ! check-pointing parameters
      ! when cpflag is greater than 0, the basis vectors will be written
      ! out roughly 'cpflag' times.  For simplicitly, each PE writes its
      ! own portion of the basis vectors to a file with cpfile followed by
      ! the processor number.  The file is written as unformatted fortran
      ! files with the following content:
      ! nrow, kb(basis size)
      ! alpha(1:kb)
      ! beta(1:kb)
      ! 1st basis vector, 2nd basis vector, ..., last basis vector
      ! the residual vector
      Integer :: cpflag, cpio
      Integer :: cpg_fp
      Character(128) :: cpfile, oldcpf
 
 
      Integer :: dummy   ! a dummy variable to fill space
      End Type TRLAN_INFO_T
      End Module trlan_info
