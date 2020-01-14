!====================================================================== 
!
!     NekNek routines for FSI
!     Author: Prabal Negi
!     Comments:  
! 
!---------------------------------------------------------------------- 
!====================================================================== 

      subroutine fsi_coupling()

      implicit none
      include 'SIZE'
      include 'NEKNEK'        ! igeom
      include 'GEOM'          ! ifgeom
      include 'INPUT'         ! ifrich
      include 'CTIMER'
      include 'STRUCT'
      include 'FSI'           
      include 'MVGEOM'        ! wx
      include 'TSTEP'
      include 'PARALLEL'      ! iglobalcomm
      include 'SOLN'

      include 'FSI_DEBUG'

      integer i,itr
      logical ifconverged
      real scale

      real            rv1 (lx1,ly1,lz1,lelv)
     $ ,              rv2 (lx1,ly1,lz1,lelv)
     $ ,              rv3 (lx1,ly1,lz1,lelv)

      real            prt (lx2,ly2,lz2,lelv)

      integer maxloops
      real enorm
      real beta
      real alph
      integer nt

      if (istep.eq.0) then
        call opzero(ext_vx,ext_vy,ext_vz)
        call opzero(res_vx(1,1,1,1,1),res_vy(1,1,1,1,1),
     $              res_vz(1,1,1,1,1))
        call opzero(res_vx(1,1,1,1,2),res_vy(1,1,1,1,2),
     $              res_vz(1,1,1,1,2)) 
      endif

!     Do nothing if there's no FSI      
      if ((.not.fsi_ifstruct).and.(.not.fsi_iffluid)) return 


!     Just initialization
      if (istep.eq.0) then

!       I'll control the neknek processes myself 
        ifneknekc = .false.

        nt = nx1*ny1*nz1*nelv
        do i=1,nt
          irmask(i,1,1,1)=imask(i,1,1,1)+0.
        enddo              

        if (fsi_iffluid) then            
!          call opzero(wx,wy,wz)
          call chk_neknek
             
        else
          call opzero(velx,vely,velz)
          call opzero(accx,accy,accz)

!         Correction for interpolation            
          call chk_neknek
          continue
        endif  

        return
      endif        

      if (istep.eq.1) then
        maxloops = 10 
      else        
        maxloops = 10
      endif        

      scale = 1.0
      ifconverged = .false.
      itr = 0

      call neknekgsync()      ! sync fluid and structural sessions

      if (fsi_iffluid) then

!       Calculate fluid stresses                
        call fluid_forces(struct_ssx,struct_ssy,struct_ssz,
     $                    vx,vy,vz,pr,scale)

!       debugging        
!        call outpost(struct_ssx,struct_ssy,struct_ssz,pr,t,'ssf')

!       Send fluid stresses to structure
        call fsi_neknek_stressex(struct_ssx,struct_ssy,struct_ssz)

!        ifconverged = .true.  ! for now no iterations

!       Have to figure out how to do this        
!!       Get broadcasted convergence
!        call setnekcomm(iglobalcomm)
!        call lbcast(ifconverged)
!        call setnekcomm(intracomm)

        call opzero(sc_vx,sc_vy,sc_vz)
        call rzero(sc_pr,nx2*ny2*nz2*nelv)

!       call check_fsi_convergence()         
        do while (.not.ifconverged.and.(itr.lt.maxloops))

          itr = itr + 1

          if (nio.eq.0) write(6,*) 'Iteration cycle:', itr

!         Get new interface velocity          
          call fsi_neknek_velex(vx,vy,vz)
!         Apply neknek mask            
          call opsub2(valint(1,1,1,1,1),valint(1,1,1,1,2),
     $                valint(1,1,1,1,3),vx,vy,vz)
          call opcolv(valint(1,1,1,1,1),valint(1,1,1,1,2),
     $                valint(1,1,1,1,3),irmask)

!          call outpost(valint(1,1,1,1,1),valint(1,1,1,1,2),
!     $                 valint(1,1,1,1,3),pr,t,'val')

!         Stokes correction step          
          call stokes_solve()

!         Calculate fluid stresses                
          call fluid_forces(stokes_ssx,stokes_ssy,stokes_ssz,
     $                    sc_vx,sc_vy,sc_vz,sc_pr,scale)

          call opcopy(rv1,rv2,rv3,struct_ssx,struct_ssy,struct_ssz)
          call opadd2(rv1,rv2,rv3,stokes_ssx,stokes_ssy,stokes_ssz)

!         Send (total) fluid stresses to structure
          call fsi_neknek_stressex(rv1,rv2,rv3)

!!         debugging 
!          call outpost(sc_vx,sc_vy,sc_vz,sc_pr,t,'stk')
!          call outpost(vx,vy,vz,pr,t,'std')
!          call opcopy(rv1,rv2,rv3,vx,vy,vz)
!          call add2(prt,sc_pr,nx2*ny2*nz2*nelv)
!          call opadd2(rv1,rv2,rv3,sc_vx,sc_vy,sc_vz) 
!          call outpost(rv1,rv2,rv3,pr,t,'ful')

        enddo

!       Update velocities        
        call opadd2(vx,vy,vz,sc_vx,sc_vy,sc_vz)
        call add2(pr,sc_pr,nx2*ny2*nz2*nelv)

      else

        do while (.not.ifconverged.and.(itr.le.maxloops))

          itr = itr + 1

          if (nio.eq.0) write(6,*) 'Iteration cycle:', itr

!         Send extrapolated velocities to fluid
          if (itr.gt.1) then
!           This contains extrapolated velocities
            call fsi_neknek_velex(rv1,rv2,rv3)
          endif            

!         Need previous fields for update            
          if (itr.eq.1) then
            call opcopy(vxlag,vylag,vzlag,vx,vy,vz)
            call opcopy(velxlag,velylag,velzlag,velx,vely,velz)
            call opcopy(axlag,aylag,azlag,accx,accy,accz)

!           This should be extrapolated velocity. But just keeping this
!           for now since we are doing zeroth order extrapolation 
            call opcopy(rv1,rv2,rv3,velx,vely,velz)
          else
            call opcopy(vx,vy,vz,vxlag,vylag,vzlag)
            call opcopy(velx,vely,velz,velxlag,velylag,velzlag)
            call opcopy(accx,accy,accz,axlag,aylag,azlag)
          endif            

!         Get fluid stresses 
          call fsi_neknek_stressex(rv1,rv2,rv3) ! just dummy variables

!         Solve structural equation
          call plan_s
          
!         Check Get extrapolated velocities
          call extrap_int_vel(enorm,rv1,rv2,rv3,itr)

!          call opcopy(velx,vely,velz,rv1,rv2,rv3)
          call outpost(rv1,rv2,rv3,pr,t,'cor')

!         if they match          
!           ifconverged=.true.
!           Broadcast flag to fluid

!         Figure out how to bcast from structural part            
!!         For now doing just one iteration            
!          ifconverged = .true.
!          call setnekcomm(iglobalcomm)
!          call lbcast(ifconverged)
!          call setnekcomm(intracomm)

        enddo 

 
      endif

      
      if (nio.eq.0) write(6,*) 'Outside Iteration loop'

!     Send extrapolated velocities to fluid
      if (fsi_ifstruct) then 
        call fsi_neknek_velex(velx,vely,velz)
      else        
        call fsi_neknek_velex(vx,vy,vz)
        call opcopy(ext_vx,ext_vy,ext_vz,valint(1,1,1,1,1),
     $            valint(1,1,1,1,2),valint(1,1,1,1,3))
!       Apply neknek mask            
        call opcolv(ext_vx,ext_vy,ext_vz,irmask)
       
      endif        


!      call neknekgsync()      ! sync fluid and structural sessions

      return            
      end subroutine fsi_coupling
!---------------------------------------------------------------------- 

      subroutine fsi_neknek_velex(ux,uy,uz)

      implicit none            

      include 'SIZE'
!      include 'TOTAL'
      include 'TSTEP'
      include 'SOLN'
      include 'NEKNEK'
      include 'CTIMER'
      include 'STRUCT'

      include 'GEOM'    ! debugging

      integer lt,lxyz
      parameter (lt=lx1*ly1*lz1*lelt,lxyz=lx1*ly1*lz1)

      real ux(lt),uy(lt),uz(lt)

      real pm1,wk1,wk2
      common /scrcg/ pm1(lt),wk1(lxyz),wk2(lxyz)

      real fieldout(nmaxl_nn,nfldmax_nn)
      real field(lx1*ly1*lz1*lelt)
      integer nv,nt,i,j,k,n,ie,ix,iy,iz,idx,ifld

      real etime,tsync

!      if (nio.eq.0) write(6,98) 
!     $   ' Multidomain data exchange ... ', nfld_neknek
! 98   format(12x,a,i3)

      etime0 = dnekclock_sync()
      call neknekgsync()
      etime1 = dnekclock()

!      call mappr(pm1,pr,wk1,wk2)  ! Map pressure to pm1 
      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

c     Interpolate using findpts_eval
      call field_eval(fieldout(1,1),1,ux)
      call field_eval(fieldout(1,2),1,uy)
      if (ldim.eq.3) call field_eval(fieldout(1,ldim),1,uz)

!      call field_eval(fieldout(1,ldim+1),1,pm1)
!      if (nfld_neknek.gt.ldim+1) then 
!        do i=ldim+2,nfld_neknek
!          call field_eval(fieldout(1,i),1,t(1,1,1,1,i-ldim-1))
!        enddo
!      endif
         
c     Now we can transfer this information to valint array from which
c     the information will go to the boundary points
      do i=1,npoints_nn
        idx = iList(1,i)
        do ifld=1,nfld_neknek
          valint(idx,1,1,1,ifld)=fieldout(i,ifld)
        enddo

!       debugging        
!        write(6,*) xm1(idx,1,1,1),ym1(idx,1,1,1),
!     $             fieldout(i,1),fieldout(i,2)

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

      subroutine fsi_neknek_stressex(sx,sy,sz)

      implicit none            

      include 'SIZE'
!      include 'TOTAL'
      include 'TSTEP'
      include 'NEKNEK'
      include 'CTIMER'
      include 'STRUCT'

      integer lt,lxyz
      parameter (lt=lx1*ly1*lz1*lelt,lxyz=lx1*ly1*lz1)

      real sx(lt),sy(lt),sz(lt)

      real pm1,wk1,wk2
      common /scrcg/ pm1(lt),wk1(lxyz),wk2(lxyz)

      real fieldout(nmaxl_nn,nfldmax_nn)
      real field(lx1*ly1*lz1*lelt)
      integer nv,nt,i,j,k,n,ie,ix,iy,iz,idx,ifld

      real etime,tsync

!      if (nio.eq.0) write(6,98) 
!     $   ' Multidomain data exchange ... ', nfld_neknek
! 98   format(12x,a,i3)

      etime0 = dnekclock_sync()
      call neknekgsync()
      etime1 = dnekclock()

      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

c     Interpolate using findpts_eval
      call field_eval(fieldout(1,1),1,sx)
      call field_eval(fieldout(1,2),1,sy)
      if (ldim.eq.3) call field_eval(fieldout(1,ldim),1,sz)

      if (fsi_ifstruct) then
        do i=1,npoints_nn
          idx = iList(1,i)
          struct_ssx(idx,1,1,1)=fieldout(i,1)
          struct_ssy(idx,1,1,1)=fieldout(i,2)
          if (ldim.eq.3) struct_ssz(idx,1,1,1)=fieldout(i,3)
        enddo
      endif        

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

      real tmp_mult(lx1,ly1,lz1,lelt)     ! since vertices/edges are on multiple faces
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

      call opzero(fx,fy,fz)
      call rzero(tmp_mult,n)

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
      call rzero(tmp_mult(1,1,1,e),nx1*ny1*nz1)
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
              n1 = unx(i,1,f,e) !*area(i,1,f,e)
              n2 = uny(i,1,f,e) !*area(i,1,f,e)
              n3 = unz(i,1,f,e) !*area(i,1,f,e)
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
              fx(j1,j2,1,e) = fx(j1,j2,1,e) + dg(1,1) + dg(1,2)
              fy(j1,j2,1,e) = fy(j1,j2,1,e) + dg(2,1) + dg(2,2)
              fz(j1,j2,1,e) = fz(j1,j2,1,e) + dg(3,1) + dg(3,2)

              tmp_mult(j1,j2,1,e) = tmp_mult(j1,j2,1,e) + 1.

            enddo       ! j2
            enddo       ! j1
          else ! 2D
            i = 0
            a = 0
            do j2=js2,jf2,jskip2
            do j1=js1,jf1,jskip1
              i = i+1
              n1 = unx(i,1,f,e) !*area(i,1,f,e)
              n2 = uny(i,1,f,e) !*area(i,1,f,e)
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
              fx(j1,j2,1,e) = fx(j1,j2,1,e) + dg(1,1) + dg(1,2)
              fy(j1,j2,1,e) = fy(j1,j2,1,e) + dg(2,1) + dg(2,2)

              tmp_mult(j1,j2,1,e) = tmp_mult(j1,j2,1,e) + 1.
            enddo       ! j2
            enddo       ! j1
          endif         ! if3d

!         This is if one element has multiple faces on the interface          
          do j2=js2,jf2,jskip2
          do j1=js1,jf1,jskip1
            fx(j1,j2,1,e) = fx(j1,j2,1,e)/tmp_mult(j1,j2,1,e) 
            fy(j1,j2,1,e) = fy(j1,j2,1,e)/tmp_mult(j1,j2,1,e)
            if (if3d.or.ifaxis) then
              fx(j1,j2,1,e) = fx(j1,j2,1,e)/tmp_mult(j1,j2,1,e)
            endif
          enddo   ! j1
          enddo   ! j2

        endif           ! if intflag(f,e).eq.1
      enddo             ! f=1,nfaces
      enddo             ! e=1,nelv

!      do e  = 1,nelt
!      do f  = 1,nfaces
!        if (intflag(f,e).eq.1) then
!          call dsset(lx1,ly1,lz1)    ! set up counters
!          pf     = eface1(f)         ! convert from preproc. notation
!          js1    = skpdat(1,pf)
!          jf1    = skpdat(2,pf)
!          jskip1 = skpdat(3,pf)
!          js2    = skpdat(4,pf)
!          jf2    = skpdat(5,pf)
!          jskip2 = skpdat(6,pf)
!
!          do j2=js2,jf2,jskip2
!          do j1=js1,jf1,jskip1
!            fx(j1,j2,1,e) = fx(j1,j2,1,e)/tmp_mult(j1,j2,1,e) 
!            fy(j1,j2,1,e) = fy(j1,j2,1,e)/tmp_mult(j1,j2,1,e)
!            if (if3d.or.ifaxis) then
!              fx(j1,j2,1,e) = fx(j1,j2,1,e)/tmp_mult(j1,j2,1,e)
!            endif
!          enddo   ! j1
!          enddo   ! j2
!        endif     ! intflag.eq.1          
!      enddo       ! f=1,nfaces
!      enddo       ! e=1,nelt

      return
      end subroutine fluid_forces
c-----------------------------------------------------------------------

      subroutine extrap_int_vel(enorm,rv1,rv2,rv3,itr)

      implicit none

      include 'SIZE'
      include 'STRUCT'
      include 'FSI'
      include 'NEKNEK'
      include 'TOPOL'
      include 'GEOM'

      integer itr

!     Contains previous velocity prediction      
      real rv1(lx1,ly1,lz1,lelv)
      real rv2(lx1,ly1,lz1,lelv)
      real rv3(lx1,ly1,lz1,lelv)

      integer i,e,f,nfaces
      integer pf,js1,js2,jf1,jf2,jskip1,jskip2
      integer j1,j2
      real enorm,a,ta
      real e2,ta2
      real glsum
      real resid,resid_diff,resid_abs
      real rs_num,rs_den


      nfaces=2*ndim

      ta = 0.
      enorm = 0.
      resid_diff = 0.
      resid_abs  = 0.
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

          i = 0
          do j2=js2,jf2,jskip2
          do j1=js1,jf1,jskip1
            i = i +1
            a =  area(i,1,f,e)
            ta = ta + a
            
!           Save new residual x            
            resid = velx(j1,j2,1,e)-rv1(j1,j2,1,e)
            res_vx(j1,j2,1,e,2) = res_vx(j1,j2,1,e,1) ! copy old residual vector
            res_vx(j1,j2,1,e,1) = resid
            enorm = enorm + a*(resid**2)

!           Numerator for \omega_i+1        
            resid      = res_vx(j1,j2,1,e,1)-res_vx(j1,j2,1,e,2)
            resid_diff = resid_diff + res_vx(j1,j2,1,e,2)*resid*a
!           Denominator for \omega_i+1
            resid_abs  = resid_abs + a*resid**2

!           Save new residual y           
            resid = vely(j1,j2,1,e)-rv2(j1,j2,1,e)
            res_vy(j1,j2,1,e,2) = res_vy(j1,j2,1,e,1) ! copy old residual vector
            res_vy(j1,j2,1,e,1) = resid
            enorm = enorm + a*(resid**2)

!           Numerator for \omega_i+1      
            resid      = res_vy(j1,j2,1,e,1)-res_vy(j1,j2,1,e,2)
            resid_diff = resid_diff + res_vy(j1,j2,1,e,2)*resid*a 
!           Denominator for \omega_i+1
            resid_abs  = resid_abs + a*resid**2            

            if (ldim.eq.3) then
!             Save new residual z                      
              resid = velz(j1,j2,1,e)-rv3(j1,j2,1,e)
              res_vz(j1,j2,1,e,2) = res_vz(j1,j2,1,e,1) ! copy old residual vector
              res_vz(j1,j2,1,e,1) = resid
              enorm = enorm + a*(resid**2)

!             Numerator for \omega_i+1 
              resid      = res_vy(j1,j2,1,e,1)-res_vy(j1,j2,1,e,2)
              resid_diff = resid_diff + res_vy(j1,j2,1,e,2)*resid*a 
!             Denominator for \omega_i+1
              resid_abs  = resid_abs + a*resid**2
            endif          
          enddo   ! j1
          enddo   ! j2
        endif     ! intflag.eq.1          
      enddo       ! f=1,nfaces
      enddo       ! e=1,nelt

      e2  = glsum(enorm,1)
      ta2 = glsum(ta,1)
     
      enorm = sqrt(e2/ta2)

      rs_num = glsum(resid_diff,1)
      rs_den = glsum(resid_abs,1)

      if (itr.lt.2) then
        fsi_omega = 0.001
      else
        fsi_omega = -fsi_omega*rs_num/rs_den
      endif            
      
      write(6,*) 'Omega:', fsi_omega
      write(6,*) 'Enorm:', enorm

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

          i = 0
          do j2=js2,jf2,jskip2
          do j1=js1,jf1,jskip1
            i = i +1
            a =  area(i,1,f,e)
            ta = ta + a
            
!           New predicted velocity
            rv1(j1,j2,1,e) = rv1(j1,j2,1,e) + 
     $                              fsi_omega*res_vx(j1,j2,1,e,1)

            rv2(j1,j2,1,e) = rv2(j1,j2,1,e) + 
     $                              fsi_omega*res_vy(j1,j2,1,e,1)

            if (ldim.eq.3) then
              rv3(j1,j2,1,e) = rv3(j1,j2,1,e) + 
     $                              fsi_omega*res_vz(j1,j2,1,e,1)
            endif 
          enddo   ! j1
          enddo   ! j2
        endif     ! intflag.eq.1          
      enddo       ! f=1,nfaces
      enddo       ! e=1,nelt


      return
      end subroutine extrap_int_vel
!----------------------------------------------------------------------

      subroutine chk_neknek

      implicit none
      
      include 'SIZE'
      include 'NEKNEK'
      include 'WZ'

      integer i,j,k
      integer ix,iy,iz
      real xmin,ymin,zmin,d
      real xdist(1000),ydist(1000),zdist(1000)

      do i=1,npoints_nn
        j=(i-1)*ldim+1
!        write(6,*) rst(j),rst(j+1),rcode(i),elid(i)

        xmin=10.
        ymin=10.
        if (ldim.eq.3) zmin=10.
        ix=1
        iy=1
        iz=1
        do k=1,lx1
          d = abs(rst(j)-zgm1(k,1))
          if (d.lt.xmin) then
            xmin = d
            ix   = k
          endif

          d = abs(rst(j+1)-zgm1(k,2))
          if (d.lt.ymin) then
            ymin = d
            iy   = k
          endif

          if (ldim.eq.3) then
            d = abs(rst(j+2)-zgm1(k,3))
            if (d.lt.zmin) then
              ymin = d
              iz   = k
            endif
          endif              
        enddo
        rst(j)=zgm1(ix,1)
        rst(j+1)=zgm1(iy,2)
        if (ldim.eq.3) rst(j+2)=zgm1(iz,3)

        xdist(i)=xmin
        ydist(i)=ymin
        if (ldim.eq.3) zdist(i)=zmin

!       Reset rcode            
        rcode(i)=0
      enddo        

!      write(6,*) 'Corrected'
!      do i=1,npoints_nn
!        j=(i-1)*ldim+1
!        write(6,*) rst(j),rst(j+1),rcode(i),elid(i)
!      enddo        

      return
      end subroutine chk_neknek
!---------------------------------------------------------------------- 
