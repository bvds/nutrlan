
      Program simple
      Use trlan_info
      Implicit None
      include 'mpif.h'
      
      Integer, Parameter :: nrow=1000, lohi=-1, mev=100
      Integer ned, maxlan, restart, lwrk, ldwrk, ierr
      Real(8) :: eval(mev), exact(mev)
      Real(8), DIMENSION(:), ALLOCATABLE :: dwrk, res
      Complex*16 :: evec(nrow, mev)
      Complex*16, DIMENSION(:), ALLOCATABLE :: wrk
      Type(trlan_info_t) :: info
      Integer :: i, ii, check
      External diag_op

      Call Mpi_init(ierr)      
      DO ned = 10, 10, 2
      
       DO maxlan = 200, 200, 10
        lwrk = maxlan*(maxlan+10)
        ldwrk = 10*maxlan
        ALLOCATE( wrk(lwrk) )
        ALLOCATE( dwrk(ldwrk) )
        
        DO restart = 7, 7
         WRITE(6,*) 'ned=',ned,'maxlan=',maxlan
         Call trl_init_info_(info,%VAL(nrow),%VAL(maxlan),
     +       %VAL(lohi),%VAL(ned),%VAL(1.4901D-8),
     +       %VAL(restart),%VAL(2000000))
         WRITE(6,*) 'calling iguess'
         Call trl_set_iguess_(info,%VAL(0),%VAL(1),%VAL(0),
     +                        %VAL(0))
         WRITE(6,*) ' calling set debug'
         Call trl_set_debug_( info, %VAL(0), 'LOG' )
         eval = 0.0D0
         evec(1:nrow,1) = 1.0D0

         Call trl_set_restart_( info,%VAL(1.5),%VAL(0.0) )
         write(6,*) 'calling trlan',nrow,mev,nrow,lwrk,ldwrk,
     +               info%rfact
         Call ztrlan_(diag_op, info, %VAL(nrow), %VAL(mev),
     +      eval, evec, %VAL(nrow),wrk,%VAL(lwrk),dwrk, 
     +      %VAL(ldwrk) )
         WRITE(6,*) '   info->stat',info%stat,eval(1)
         Call trl_print_info_(info, %VAL(3*nrow))
         Do i = 1, mev
          exact(i) = i*i
         End Do
         If (info%nec.Gt.0) Then
          i = info%nec
         Else
         End If
         DO II=1,nrow
          WRITE(10,*) evec(II,1:info%nec)
         ENDDO
         ALLOCATE( res(i) );
         Call ztrl_check_ritz_(diag_op,info,%VAL(nrow),
     +       %VAL(i),evec(:,1:i),%val(nrow),eval(1:i), 
     +       check, res, exact,%val(lwrk),wrk )
         If (info%nec .Eq. 0) info%nec = Min(info%ned, mev-1)
         DEALLOCATE( res )
        ENDDO
        DEALLOCATE( wrk,dwrk )
       ENDDO
      ENDDO
 100  Format('E(', I2, ') = ', 1PG25.17, 1PE16.4)
      Call mpi_finalize(ierr)
      End Program simple

      Subroutine diag_op(nrow, ncol, xin, ldx, yout, ldy)
      Implicit None
      include 'mpif.h'
      
      Integer, Intent(in) :: nrow, ncol, ldx, ldy
      Complex*16, Dimension(ldx*ncol), Intent(in)
     +                                     :: xin
      Complex*16, Dimension(ldy*ncol), Intent(out) 
     +                                     :: yout
      Integer :: i, j, ioff, joff, doff, ierr
      Real(8) :: re, im
      
      Call MPI_Comm_rank(MPI_COMM_WORLD, i, ierr)
      doff = nrow*i;
      Do j = 1, ncol
        ioff = (j-1)*ldx
        joff = (j-1)*ldy
        Do i = 1, nrow
           re = (i+doff)*real(xin(ioff+i))
           im = (i+doff)*aimag(xin(ioff+i))
           yout(joff+i) = dcmplx(re,im)
        End Do
      End Do
      End Subroutine diag_op
