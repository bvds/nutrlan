      Program simple
      Use trlan_info
      Implicit None
      Integer, Parameter :: nrow=1000, lohi=-1, mev=100
      Integer ned, maxlan, restart, lwrk, ldwrk, ierr
      Real(8) :: eval(mev), exact(mev)
      Real(8), DIMENSION(:), ALLOCATABLE :: dwrk, res
      Complex*16 :: evec(nrow, mev)
      Complex*16, DIMENSION(:), ALLOCATABLE :: wrk
      Type(trlan_info_t) :: info
      Integer :: i, ii, check
      External diag_op
      
      DO ned = 100, 100, 2
      
       DO maxlan = 200, 200, 10
        lwrk = maxlan*(maxlan+10)
        ldwrk = 10*maxlan
        ALLOCATE( wrk(lwrk) )
        ALLOCATE( dwrk(ldwrk) )
        
        DO restart = 7, 7
         WRITE(6,*) 'ned=',ned,'maxlan=',maxlan
         Call trl_init_info(info,%VAL(nrow),%VAL(maxlan),
     +       %VAL(lohi),%VAL(ned),%VAL(1.4901D-8),
     +       %VAL(restart),%VAL(2000000))
         Call trl_set_iguess(info,%VAL(0),%VAL(1),%VAL(0),
     +                        %VAL(0))
         Call trl_set_debug( info, %VAL(0), 'LOG' )
         eval(1:nrow) = 0.0D0
         evec(1:nrow,1) = 1.0D0
         Call ztrlan(diag_op, info, %VAL(nrow), %VAL(mev),
     +      eval, evec, %VAL(nrow), wrk, %VAL(lwrk), dwrk, 
     +      %VAL(ldwrk) )
         WRITE(6,*) '   info->stat',info%stat
         Call trl_print_info(info, %VAL(3*nrow))
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
         Call ztrl_check_ritz(diag_op,info,%VAL(nrow),
     +       %VAL(i),evec(:,1:i),%val(nrow),eval(1:i), 
     +       check, res, exact,%val(lwrk),wrk )
         If (info%nec .Eq. 0) info%nec = Min(info%ned, mev-1)
         DEALLOCATE( res )
        ENDDO
        DEALLOCATE( wrk,dwrk )
       ENDDO
      ENDDO
 100  Format('E(', I2, ') = ', 1PG25.17, 1PE16.4)
      End Program simple

      Subroutine diag_op(nrow, ncol, xin, ldx, yout, ldy)
      Implicit None
      Integer, Intent(in) :: nrow, ncol, ldx, ldy
      Complex*16, Dimension(ldx*ncol), Intent(in)
     +                                     :: xin
      Complex*16, Dimension(ldy*ncol), Intent(out) 
     +                                     :: yout
      Integer :: i, j, ioff, joff, ierr
      Real(8) :: re, im
      
      Do j = 1, ncol
        ioff = (j-1)*ldx
        joff = (j-1)*ldy
        Do i = 1, nrow
           re = (i*i)*real(xin(ioff+i))
           im = (i*i)*aimag(xin(ioff+i))
           yout(joff+i) = dcmplx(re,im)
        End Do
      End Do
      End Subroutine diag_op
