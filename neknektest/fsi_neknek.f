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
     $              '  done :: Multidomain data exchange (Vel)',
     $              etime, etime+tsync
 99   format(i11,a,1p2e13.4)

      return
      end subroutine fsi_neknek_velex
c--------------------------------------------------------------------------

      subroutine fsi_neknek_stressex

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
        call field_eval(fieldout(1,1),1,struct_ssx)
        call field_eval(fieldout(1,2),1,struct_ssy)
        if (ldim.eq.3) call field_eval(fieldout(1,ldim),1,struct_ssz)
      elseif (fsi_ifstruct) then
        call field_eval(fieldout(1,1),1,struct_ssx)
        call field_eval(fieldout(1,2),1,struct_ssy)
        if (ldim.eq.3) call field_eval(fieldout(1,ldim),1,struct_ssz)
      endif 
         
c     Now we can transfer this information to valint array from which
c     the information will go to the boundary points
      do i=1,npoints_nn
        idx = iList(1,i)
        struct_ssx(idx,1,1,1)=fieldout(i,1)
        struct_ssy(idx,1,1,1)=fieldout(i,2)
        if (if3d) struct_ssz(idx,1,1,1)=fieldout(i,3)
      enddo

      call nekgsync()
      etime = dnekclock() - etime1
      tsync = etime1 - etime0

      if (nio.eq.0) write(6,99) istep,
     $              '  done :: Multidomain data exchange (Stress)',
     $              etime, etime+tsync
 99   format(i11,a,1p2e13.4)

      return
      end subroutine fsi_neknek_stressex
c--------------------------------------------------------------------------

      subroutine fluid_forces(fx,fy,fz,ux,uy,uz,pr2,scale)

c     Compute fluid forces

      implicit none

      INCLUDE 'SIZE'  
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'   ! ifaxis
      INCLUDE 'TSTEP'   ! ifield
      INCLUDE 'SOLN'    ! vdiff
      INCLUDE 'NEKNEK'
!      INCLUDE 'OBJDATA' !
      INCLUDE 'TOPOL'

!      INCLUDE 'TOTAL' 

      real flow_rate,base_flow,domain_length,xsec,scale_vf
      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)

      real x0(3),w1(0:maxobj)
      logical ifdout,iftout

      real sij,pm1,xm0,ym0,zm0
      common /scrns/         sij (lx1,ly1,lz1,3*ldim-3,lelv)
      common /scrcg/         pm1 (lx1,ly1,lz1,lelv)
      common /scrsf/         xm0(lx1,ly1,lz1,lelt)
     $,                      ym0(lx1,ly1,lz1,lelt)
     $,                      zm0(lx1,ly1,lz1,lelt)

      integer lr
      parameter (lr=lx1*ly1*lz1)

      real ur,us,ut,vr,vs,vt,wr,ws,wt
      common /scruz/         ur(lr),us(lr),ut(lr)
     $                     , vr(lr),vs(lr),vt(lr)
     $                     , wr(lr),ws(lr),wt(lr)

!     Input      
      real ux(lx1,ly1,lz1,lelt),uy(lx1,ly1,lz1,lelt)
      real uz(lx1,ly1,lz1,lelt),pr2(lx2,ly2,lz2,lelt)

!     Output      
      real fx(lx1,ly1,lz1,lelt),fy(lx1,ly1,lz1,lelt)
      real fz(lx1,ly1,lz1,lelt)

      real dpdx_mean,dpdy_mean,dpdz_mean
      real scale

      integer e,f,nfaces,i
      integer n,nij
      integer pf,js1,jf1,jskip1,js2,jf2,jskip2
      integer j1,j2

      real s11,s21,s31,s12,s22,s32,s13,s23,s33
      real v,a,n1,n2,n3
      real dg(3,2)

      n = lx1*ly1*lz1*nelv

      call mappr(pm1,pr2,xm0,ym0) ! map pressure onto Mesh 1

!     Add mean_pressure_gradient.X to p:

      if (param(55).ne.0) then
         dpdx_mean = -scale_vf(1)
         dpdy_mean = -scale_vf(2)
         dpdz_mean = -scale_vf(3)
      endif

      call add2s2(pm1,xm1,dpdx_mean,n)  ! Doesn't work if object is cut by 
      call add2s2(pm1,ym1,dpdy_mean,n)  ! periodicboundary.  In this case,
      call add2s2(pm1,zm1,dpdz_mean,n)  ! set ._mean=0 and compensate in

!    Compute sij

      nij = 3
      if (if3d.or.ifaxis) nij=6
      call comp_sij(sij,nij,ux,uy,uz,ur,us,ut,vr,vs,vt,wr,ws,wt)

!     Fill up viscous array w/ default

      if (istep.lt.1) call cfill(vdiff,param(2),n)

      call opzero(fx,fy,fz)

      nfaces = 2*ldim
      ifield = 1
      do e  = 1,nelt
      do f  = 1,nfaces
        if (intflag(f,e).eq.1) then
          call dsset(lx1,ly1,lz1)    ! set up counters
          pf     = eface1(f)         ! convert from preproc. notation
          js1    = skpdat(1,pf)
          jf1    = skpdat(2,pf)
          jskip1 = skpdat(3,pf)
          js2    = skpdat(4,pf)
          jf2    = skpdat(5,pf)
          jskip2 = skpdat(6,pf)

          if (if3d.or.ifaxis) then
            i = 0
            a = 0
            do j2=js2,jf2,jskip2
            do j1=js1,jf1,jskip1
              i = i+1
              n1 = unx(i,1,f,e)*area(i,1,f,e)
              n2 = uny(i,1,f,e)*area(i,1,f,e)
              n3 = unz(i,1,f,e)*area(i,1,f,e)
              a  = a +          area(i,1,f,e)

              v  = vdiff(j1,j2,1,e,ifield)

              s11 = sij(j1,j2,1,1,e)
              s21 = sij(j1,j2,1,4,e)
              s31 = sij(j1,j2,1,6,e)

              s12 = sij(j1,j2,1,4,e)
              s22 = sij(j1,j2,1,2,e)
              s32 = sij(j1,j2,1,5,e)

              s13 = sij(j1,j2,1,6,e)
              s23 = sij(j1,j2,1,5,e)
              s33 = sij(j1,j2,1,3,e)

              dg(1,1) = pm1(j1,j2,1,e)*n1     ! pressure drag
              dg(2,1) = pm1(j1,j2,1,e)*n2
              dg(3,1) = pm1(j1,j2,1,e)*n3

              dg(1,2) = -v*(s11*n1 + s12*n2 + s13*n3) ! viscous drag
              dg(2,2) = -v*(s21*n1 + s22*n2 + s23*n3)
              dg(3,2) = -v*(s31*n1 + s32*n2 + s33*n3)

              call cmult(dg,scale,6)

!             f = pressure force + viscous force
              fx(j1,j2,1,e) = dg(1,1) + dg(1,2)
              fy(j1,j2,1,e) = dg(2,1) + dg(2,2)
              fz(j1,j2,1,e) = dg(3,1) + dg(3,2)

            enddo       ! j2
            enddo       ! j1
          else ! 2D
            i = 0
            a = 0
            do j2=js2,jf2,jskip2
            do j1=js1,jf1,jskip1
              i = i+1
              n1 = unx(i,1,f,e)*area(i,1,f,e)
              n2 = uny(i,1,f,e)*area(i,1,f,e)
              a  = a +          area(i,1,f,e)
              v  = vdiff(j1,j2,1,e,ifield)

              s11 = sij(j1,j2,1,1,e)
              s12 = sij(j1,j2,1,3,e)
              s21 = sij(j1,j2,1,3,e)
              s22 = sij(j1,j2,1,2,e)

              dg(1,1) = pm1(j1,j2,1,e)*n1     ! pressure drag
              dg(2,1) = pm1(j1,j2,1,e)*n2

              dg(1,2) = -v*(s11*n1 + s12*n2) ! viscous drag
              dg(2,2) = -v*(s21*n1 + s22*n2)

              call cmult(dg,scale,6)

!             f = pressure force + viscous force
              fx(j1,j2,1,e) = dg(1,1) + dg(1,2)
              fy(j1,j2,1,e) = dg(2,1) + dg(2,2)
            enddo       ! j2
            enddo       ! j1
          endif         ! if3d
        endif           ! if intflag(f,e).eq.1
      enddo             ! f=1,nfaces
      enddo             ! e=1,nelv

      return
      end subroutine fluid_forces
c-----------------------------------------------------------------------

