      Program simple
      Use trlan_info
      Implicit None
      Integer, Parameter :: nrow=1000, lohi=-1, mev=100
      Integer ned, maxlan, restart, lwrk
      Real(8) :: eval(mev), evec(nrow, mev), exact(mev)
      Real(8), DIMENSION(:), ALLOCATABLE :: res, wrk
      Type(trlan_info_t) :: info
      Integer :: i, ii, check
      External diag_op, user_ortho
      DO ned = 10, 10, 2
      DO maxlan = 200, 200, 10
      lwrk = maxlan*(maxlan+10)
      ALLOCATE( res(lwrk) )
      ALLOCATE( wrk(lwrk) )
      DO restart = 1, 1
      WRITE(6,*) 'ned=',ned,'maxlan=',maxlan
      Call trl_init_info(info,%VAL(nrow),%VAL(maxlan),%VAL(lohi),
     +       %VAL(ned),%VAL(1.4901D-8),%VAL(restart),%VAL(2000000))
      Call trl_set_iguess(info,%VAL(0),%VAL(1),%VAL(0),%VAL(0))
      eval(1:nrow) = 0.0D0
      evec(1:nrow,1) = 1.0D0
      Call trlan(diag_op, user_ortho, info, %VAL(nrow), %VAL(mev), eval, evec, 
     +           %VAL(nrow), %VAL(lwrk), res )
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
      Call trl_check_ritz(diag_op, info, %VAL(nrow),
     +       %VAL(i),evec(:,1:i),%val(nrow),eval(1:i), 
     +       check, res, exact, %val(lwrk), wrk )
      If (info%nec .Eq. 0) info%nec = Min(info%ned, mev-1)
      ENDDO
      DEALLOCATE( res,wrk )
      ENDDO
      ENDDO
 100  Format('E(', I2, ') = ', 1PG25.17, 1PE16.4)
      End Program simple

      Subroutine diag_op(nrow, ncol, xin, ldx, yout, ldy)
      Implicit None
      Integer, Intent(in) :: nrow, ncol, ldx, ldy
      Real(8), Dimension(ldx*ncol), Intent(in) :: xin
      Real(8), Dimension(ldy*ncol), Intent(out) :: yout
      Integer :: i, j, ioff, joff
      Do j = 1, ncol
        ioff = (j-1)*ldx
        joff = (j-1)*ldy
        Do i = 1, nrow
           yout(joff+i) = i*xin(ioff+i)
        End Do
      End Do
      End Subroutine diag_op

      Subroutine user_ortho(nrow, x, y)
      IMPLICIT NONE
      INTEGER, intent(in)  :: nrow
      Real(8), intent(input) :: x
      Real(8), intent(output) :: y
      INTEGER i
      Do i = 1, nrow
         y(i) = x(i)
      End Do
      End Subroutine user_ortrho
