
      module trlan_interface
      integer key
      integer trl_mrn, trl_ilocal, trl_ifso, trl_imthd
      double precision trl_Eref
      complex*16,allocatable,dimension(:,:) :: trl_wrk
      complex*16,allocatable,dimension(:,:) :: trl_ugh
      
      contains
     
      subroutine trlan_allocate( np,mr_n,mg_nx,ilocal,imthd,if_so )
      integer mr_n, ilocal, imthd, if_so, lwrk
      
      trl_mrn    = mr_n
      trl_mgnx   = mg_nx
      trl_ilocal = ilocal
      trl_ifso   = if_so
      trl_imthd  = imthd

      write(6,*) 'info',if_so,'mrn,mgnx',mr_n,mg_nx
      if( if_so.eq.1 ) then
         allocate( trl_wrk(2*mr_n,np) )
         allocate( trl_ugh(2*mg_nx,np) )
      else
         allocate( trl_wrk(mr_n+1,np) )
         allocate( trl_ugh(mg_nx,np) )
      endif
      
      end subroutine

      subroutine trlan_setkey( key2 )
      integer key2
      key = key2
      end subroutine

      subroutine trlan_deallocate()
      deallocate( trl_wrk,trl_ugh )
      end subroutine

      end module trlan_interface

