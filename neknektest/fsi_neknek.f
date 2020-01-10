!----------------------------------------------------------------------
!
!     NekNek routines for FSI
!     Author: Prabal Negi
!     Comments: Mostly slightly modified routines from multimesh.f      
! 
!---------------------------------------------------------------------- 
!====================================================================== 
      subroutine fsi_neknek_velex

      implicit none            

      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      include 'CTIMER'
      include 'STRUCT'

      integer lt,lxyz
      parameter (lt=lx1*ly1*lz1*lelt,lxyz=lx1*ly1*lz1)

      real pm1,wk1,wk2
      common /scrcg/ pm1(lt),wk1(lxyz),wk2(lxyz)

      real fieldout(nmaxl_nn,nfldmax_nn)
      real field(lx1*ly1*lz1*lelt)
      integer nv,nt,i,j,k,n,ie,ix,iy,iz,idx,ifld

      real etime,tsync

      if (nio.eq.0) write(6,98) 
     $   ' Multidomain data exchange ... ', nfld_neknek
 98   format(12x,a,i3)

      etime0 = dnekclock_sync()
      call neknekgsync()
      etime1 = dnekclock()

      call mappr(pm1,pr,wk1,wk2)  ! Map pressure to pm1 
      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

c     Interpolate using findpts_eval
      if (fsi_iffluid) then
        call field_eval(fieldout(1,1),1,vx)
        call field_eval(fieldout(1,2),1,vy)
        if (ldim.eq.3) call field_eval(fieldout(1,ldim),1,vz)
      elseif (fsi_ifstruct) then
!        call field_eval(fieldout(1,1),1,velx)
!        call field_eval(fieldout(1,2),1,vely)
!        if (ldim.eq.3) call field_eval(fieldout(1,ldim),1,velz)

        call field_eval(fieldout(1,1),1,vx)
        call field_eval(fieldout(1,2),1,vy)
        if (ldim.eq.3) call field_eval(fieldout(1,ldim),1,vz)
           
      endif 

      call field_eval(fieldout(1,ldim+1),1,pm1)
      if (nfld_neknek.gt.ldim+1) then 
        do i=ldim+2,nfld_neknek
          call field_eval(fieldout(1,i),1,t(1,1,1,1,i-ldim-1))
        enddo
      endif
         
c     Now we can transfer this information to valint array from which
c     the information will go to the boundary points
      do i=1,npoints_nn
        idx = iList(1,i)
        do ifld=1,nfld_neknek
          valint(idx,1,1,1,ifld)=fieldout(i,ifld)
        enddo
      enddo

      call nekgsync()
      etime = dnekclock() - etime1
      tsync = etime1 - etime0

      if (nio.eq.0) write(6,99) istep,
     $              '  done :: Multidomain data exchange',
     $              etime, etime+tsync
 99   format(i11,a,1p2e13.4)

      return
      end subroutine fsi_neknek_velex
c--------------------------------------------------------------------------


