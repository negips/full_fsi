!======================================================================
!     Routines for coupled fluid structure interaction problems
!     with linear structural equations      
!     Author: Prabal S. Negi
!
!======================================================================       

      subroutine fsi_advance()

      implicit none
      include 'SIZE'
      include 'NEKNEK'        ! igeom
      include 'GEOM'          ! ifgeom
      include 'INPUT'         ! ifrich
      include 'CTIMER'
      include 'STRUCT'

      integer ntot
      logical ifconverged

      ifconverged = .false.

      if (fsi_iffluid) then

!       call check_fsi_convergence()         
        do while (.not.ifconverged)            
!         Stokes correction step            
!          call stokes_solve()

!         call check_fsi_convergence(ifconverged)

        enddo
      else            
!       First Structural solve
        do igeom=1,ngeom

          if (ifgeom) then
             if (.not.ifrich) call gengeom (igeom)
             call geneig  (igeom)
          endif
        
!          call struct(igeom) 
        enddo
      endif


      return            
      end subroutine fsi_advance
!---------------------------------------------------------------------- 
      subroutine plan_s 

      implicit none

      include 'SIZE'
      include 'STRUCT'
      include 'INPUT'
      include 'NEKNEK'        ! igeom
      include 'GEOM'          ! ifgeom
      include 'TSTEP'         ! ifield,DT

      include 'SOLN'          ! just testing

!      real resv1,resv2,resv3,resx1,resx2,resx3,h2

      real            resv1 (lx1,ly1,lz1,lelv)
     $ ,              resv2 (lx1,ly1,lz1,lelv)
     $ ,              resv3 (lx1,ly1,lz1,lelv)
     $ ,              resx1 (lx1,ly1,lz1,lelv)
     $ ,              resx2 (lx1,ly1,lz1,lelv)
     $ ,              resx3 (lx1,ly1,lz1,lelv)
     $ ,              h2    (lx1,ly1,lz1,lelv)

      real            dv1   (lx1,ly1,lz1,lelv)
     $ ,              dv2   (lx1,ly1,lz1,lelv)
     $ ,              dv3   (lx1,ly1,lz1,lelv)


!      real h1
      real            h1    (lx1,ly1,lz1,lelv)

      real const


!     For now.
      ifgeom = .false.
      ifield = 1

!     First Structural solve
      call opzero(resv1,resv2,resv3)
      call opzero(resx1,resx2,resx3)

      do igeom=1,ngeom

        if (igeom.eq.1) then
          call struct_makef
        else

!!         testing              
!          call opcopy(ts1,ts2,ts3,bfx,bfy,bfz)
!          call opcopy(ts4,ts5,ts6,vx,vy,vz)
         
!          call outpost(ts1,ts2,ts3,pr,t,'dbg') 
!          call outpost(ts4,ts5,ts6,pr,t,'dbg')

!         initial residual for correction 
          call struct_cresvif(resv1,resv2,resv3)

!         Solve for increment          
          call solve_elasticity(resv1,resv2,resv3)

!         debugging            
          call opcopy(ts1,ts2,ts3,resv1,resv2,resv3)

          call update_fields(resv1,resv2,resv3) 

!!         New displacements          
!!          call opcopy(vx,vy,vz,resv1,resv2,resv3)
!          call opadd2(vx,vy,vz,resv1,resv2,resv3)
!
!!         Un+1 - Un
!          call opcopy(resv1,resv2,resv3,vx,vy,vz)
!          call sub2(resv1,resv2,resv3,vxlag,vylag,vzlag)
!
!!         New velocity            
!!         Vn+1 = 2/DT*(Un+1 - Un) - Vn            
!          const=2./DT
!          call opcmult(resv1,resv2,resv3,const)
!          call opcopy(resx1,resx2,resx3,velx,vely,velz)     ! save  
!          call opsub2(velx,vely,velz,resv1,resv2,resv3)
!          call opcmult(velx,vely,velz,-1.)
!          
!!         New Acceleration            
!!         An+1 = -[An - (4/DT^2)*(Un+1 - Un) + 4/DT*Vn]
!          const=2./DT
!          call opcmult(resv1,resv2,resv3,const) ! already has 2/DT factor
!          call opsub2(accx,accy,accz,resv1,resv2,resv3)
!
!          const=4./DT
!          call opadd2cm(accx,accy,accz,resx1,resx2,resx3,const)
!          call opcmult(accx,accy,accz,-1.)

        endif
      enddo

      ifield = 2


      return
      end subroutine plan_s
!----------------------------------------------------------------------
      subroutine struct_makef

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'CTIMER'

      etime1 = dnekclock()

      call struct_makeuf                  ! body forcing in the structure
                                          ! Contains mass matrix

      if (iftran) then
        call make_accf                    ! terms due to acceleration
        call make_velf                    ! terms due to velocity
        call make_dxf                     ! terms due to displacement
      endif        

      tmakf=tmakf+(dnekclock()-etime1)

      return
      end subroutine struct_makef            
!----------------------------------------------------------------------

      subroutine struct_makeuf

      implicit none

      include 'SIZE'
      include 'INPUT'   ! iftran
      include 'SOLN'    ! bfx,bfy,bfz
      include 'STRUCT'  ! struct_bfdx, ...
      include 'MASS'    ! BM1
      include 'TSTEP'   ! time

c
      time = time-dt
      call struct_nekuf(bfx,bfy,bfz)
      call opcolv(bfx,bfy,bfz,bm1)

      time = time+dt
C
      return
      end subroutine struct_makeuf
!----------------------------------------------------------------------       

      subroutine struct_nekuf (f1,f2,f3)

      implicit none            

      include 'SIZE'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'CTIMER'

      real f1 (lx1,ly1,lz1,lelv)
      real f2 (lx1,ly1,lz1,lelv)
      real f3 (lx1,ly1,lz1,lelv)

      integer i,j,k,iel,ielg

! #ifdef TIMER
!       etime1=dnekclock_sync()
! #endif

      call oprzero (f1,f2,f3)
      do 100 iel=1,nelv
         ielg = lglel(iel)
         do 100 k=1,lz1
         do 100 j=1,ly1
         do 100 i=1,lx1
            if (optlevel.le.2) call nekasgn (i,j,k,iel)
            call struct_userf   (i,j,k,ielg)     ! define this
            f1(i,j,k,iel) = ffx
            f2(i,j,k,iel) = ffy
            f3(i,j,k,iel) = ffz
 100  continue

! #ifdef TIMER
!       tusfq=tusfq+(dnekclock()-etime1)
! #endif

      return
      end subroutine struct_nekuf 
!---------------------------------------------------------------------- 
      subroutine struct_userf(ix,iy,iz,ieg)

!     User defined body forcing for structural equations
      implicit none

      include 'SIZE'
      include 'NEKUSE'

      integer ix,iy,iz,ieg

      ffx = 0.
      ffy = 0.
      ffz = 0.


      return
      end subroutine struct_userf
!----------------------------------------------------------------------
      subroutine make_accf

!     Add acceleration terms to rhs            

      implicit none            

      include 'SIZE'
      include 'SOLN'    ! vtrans
      include 'MASS'    ! bm1
      include 'STRUCT'

      real ta1,ta2,ta3,tb1,tb2,tb3,h2
      common /scrns/ ta1(lx1,ly1,lz1,lelv)
     $ ,             ta2(lx1,ly1,lz1,lelv)
     $ ,             ta3(lx1,ly1,lz1,lelv)
     $ ,             tb1(lx1,ly1,lz1,lelv)
     $ ,             tb2(lx1,ly1,lz1,lelv)
     $ ,             tb3(lx1,ly1,lz1,lelv)
     $ ,             h2 (lx1,ly1,lz1,lelv)

      integer nt
      real const
      integer ifld
      
      ifld = 1

      nt = lx1*ly1*lz1*nelv
      const = 1.
      call cmult2(h2,vtrans(1,1,1,1,ifld),const,nt)
      call col2(h2,bm1,nt)

      call opcolv3(ta1,ta2,ta3,accx,accy,accz,h2)

      call opadd2 (bfx,bfy,bfz,ta1,ta2,ta3)

      return
      end subroutine make_accf
!----------------------------------------------------------------------
      subroutine make_velf

!     Add vel terms to rhs            

      implicit none            

      include 'SIZE'
      include 'SOLN'    ! vtrans,bfx,...
      include 'MASS'    ! bm1
      include 'STRUCT'
      include 'TSTEP'   ! DT

      real ta1,ta2,ta3,tb1,tb2,tb3,h2
      common /scrns/ ta1(lx1,ly1,lz1,lelv)
     $ ,             ta2(lx1,ly1,lz1,lelv)
     $ ,             ta3(lx1,ly1,lz1,lelv)
     $ ,             tb1(lx1,ly1,lz1,lelv)
     $ ,             tb2(lx1,ly1,lz1,lelv)
     $ ,             tb3(lx1,ly1,lz1,lelv)
     $ ,             h2 (lx1,ly1,lz1,lelv)

      integer nt
      real const
      real damp_const
      integer ifld
      
      ifld = 1

      nt = lx1*ly1*lz1*nelv
      const = 4./DT
      call cmult2(h2,vtrans(1,1,1,1,ifld),const,nt)
      call col2(h2,bm1,nt)

      call opcolv3(ta1,ta2,ta3,velx,vely,velz,h2)

      call opadd2 (bfx,bfy,bfz,ta1,ta2,ta3)

!     Damping terms here. For now set to zero
      damp_const = 0.
      const = 1.*damp_const
      call copy(h2,bm1,nt)
      call cmult(h2,const,nt)

      call opcolv3(ta1,ta2,ta3,velx,vely,velz,h2)

      call opadd2 (bfx,bfy,bfz,ta1,ta2,ta3)

      return
      end subroutine make_velf
!----------------------------------------------------------------------
      subroutine make_dxf

!     Add displacement terms to rhs            

      implicit none            

      include 'SIZE'
      include 'SOLN'    ! vtrans
      include 'MASS'    ! bm1
      include 'STRUCT'
      include 'TSTEP'   ! DT

      real ta1,ta2,ta3,tb1,tb2,tb3,h2
      common /scrns/ ta1(lx1,ly1,lz1,lelv)
     $ ,             ta2(lx1,ly1,lz1,lelv)
     $ ,             ta3(lx1,ly1,lz1,lelv)
     $ ,             tb1(lx1,ly1,lz1,lelv)
     $ ,             tb2(lx1,ly1,lz1,lelv)
     $ ,             tb3(lx1,ly1,lz1,lelv)
     $ ,             h2 (lx1,ly1,lz1,lelv)

      integer nt
      real const
      real damp_const
      integer ifld
      
      ifld = 1

      nt = lx1*ly1*lz1*nelv
      const = 4./(DT**2)
      call cmult2(h2,vtrans(1,1,1,1,ifld),const,nt)
      call col2(h2,bm1,nt)

      call opcolv3(ta1,ta2,ta3,vx,vy,vz,h2)

      call opadd2 (bfx,bfy,bfz,ta1,ta2,ta3)

!     Damping terms here. For now set to zero
      damp_const = 0.
      const = 2./DT*damp_const
      call copy(h2,bm1,nt)
      call cmult(h2,const,nt)

      call opcolv3(ta1,ta2,ta3,vx,vy,vz,h2)

      call opadd2 (bfx,bfy,bfz,ta1,ta2,ta3)

      return
      end subroutine make_dxf
!----------------------------------------------------------------------

      subroutine struct_bcneutr

      implicit none

      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'

      real trx,try,trz,stc
      common /scrsf/ trx(lx1,ly1,lz1)
     $             , try(lx1,ly1,lz1)
     $             , trz(lx1,ly1,lz1)
      common /ctmp0/ stc(lx1,ly1,lz1)
      real sigst(lx1,ly1)

      logical ifalgn,ifnorx,ifnory,ifnorz
      
      character cb*3
      common  /nekcb/ cb

      integer iel,ifc,ifld,nface,nxy1,nxyz1
      real bc1,bc2,bc3,bc4

      integer j

      ifld  = 1
      nface = 2*ldim
      nxy1  = lx1*ly1
      nxyz1 = lx1*ly1*lz1
C
      do 104 iel=1,nelv
      do 104 ifc=1,nface
c
         cb  = cbc (ifc,iel,ifld)
         bc1 = bc(1,ifc,iel,ifld)
         bc2 = bc(2,ifc,iel,ifld)
         bc3 = bc(3,ifc,iel,ifld)
         bc4 = bc(4,ifc,iel,ifld)

         call rzero3 (trx,try,trz,nxyz1)

C
C        Prescribed tractions and shear tractions
C
         if (cb.eq.'S  ' .or. cb.eq.'SL ' .or.
     $       cb.eq.'SH ' .or. cb.eq.'SHL' ) then
             call trcon (trx,try,trz,bc1,bc2,bc3,iel,ifc)
             if (ifqinp(ifc,iel)) call globrot (trx,try,trz,iel,ifc)
             goto 120
         endif

         if (cb.eq.'s  ' .or. cb.eq.'sl ' .or.
     $       cb.eq.'sh ' .or. cb.eq.'shl' ) then
             call struct_faceiv (cb,trx,try,trz,iel,ifc,
     $                           lx1,ly1,lz1)
             call faccvs (trx,try,trz,area(1,1,ifc,iel),ifc)
             if (ifqinp(ifc,iel)) call globrot (trx,try,trz,iel,ifc)
             goto 120
         endif

!!        Surface-tension
!         if (cb.eq.'MS ' .or. cb.eq.'MSI' .or.
!     $       cb.eq.'MM ' .or. cb.eq.'MM ' .or.
!     $       cb.eq.'MS ' .or. cb.eq.'MSI') then
!             if (cb.eq.'MS '.or.cb.eq.'MM ') then
!                bcn = -bc1
!                call trcon   (trx,try,trz,bcn,bc2,bc3,iel,ifc)
!                call globrot (trx,try,trz,iel,ifc)
!             endif
!c            if (cb.eq.'ms '.or.cb.eq.'mm ') then
!             if (cb.eq.'ms '.or.cb.eq.'msi') then
!                call struct_faceiv  (cb,trx,try,trz,iel,ifc,lx1,ly1,lz1)
!                call faccvs  (trx,try,trz,area(1,1,ifc,iel),ifc)
!                call globrot (trx,try,trz,iel,ifc)
!             endif
!             if (cb(1:1).eq.'M') then
!                call cfill  (sigst,bc4,nxy1)
!             else
!                call faceis (cb,stc,iel,ifc,lx1,ly1,lz1)
!                call facexs (sigst,stc,ifc,0)
!             endif
!             if (ifaxis) then
!                call trstax (trx,try,sigst,iel,ifc)
!             elseif (ldim.eq.2) then
!                call trst2d (trx,try,sigst,iel,ifc)
!             else
!                call trst3d (trx,try,trz,sigst,iel,ifc)
!             endif
!         endif
C

  120    call add2 (bfx(1,1,1,iel),trx,nxyz1)
         call add2 (bfy(1,1,1,iel),try,nxyz1)
         if (ldim.eq.3) call add2 (bfz(1,1,1,iel),trz,nxyz1)

  104 continue

      return
      end subroutine struct_bcneutr
c-----------------------------------------------------------------------

      subroutine struct_faceiv (cb,v1,v2,v3,iel,iface,nx,ny,nz)

c     Assign fortran function boundary conditions to 
c     face IFACE of element IEL for vector (V1,V2,V3).

      implicit none

      INCLUDE 'SIZE'
      INCLUDE 'NEKUSE'
      INCLUDE 'PARALLEL'

      integer nx,ny,nz,iel,ieg,iface

      real v1,v2,v3
      dimension v1(nx,ny,nz),v2(nx,ny,nz),v3(nx,ny,nz)
      character cb*3
c
      character*1 cb1(3)
c
      common  /nekcb/ cb3
      character*3 cb3

      integer ix,iy,iz
      integer kx1,ky1,kz1,kx2,ky2,kz2

      cb3 = cb

      call chcopy(cb1,cb,3)

      ieg = lglel(iel)
      call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)

      if (cb.eq.'v  ' .or. cb.eq.'ws ' .or. cb.eq.'mv '.or. 
     $    cb.eq.'mvn') then
c
         do 105 iz=kz1,kz2
         do 105 iy=ky1,ky2
         do 105 ix=kx1,kx2
            if (optlevel.le.2) call nekasgn (ix,iy,iz,iel)
            call struct_userbc  (ix,iy,iz,iface,ieg)
            v1(ix,iy,iz) = ux
            v2(ix,iy,iz) = uy
            v3(ix,iy,iz) = uz
  105    continue
         return

      elseif (cb.eq.'s  ' .or. cb.eq.'sh ') then
         do 106 iz=kz1,kz2
         do 106 iy=ky1,ky2
         do 106 ix=kx1,kx2
            if (optlevel.le.2) call nekasgn (ix,iy,iz,iel)
            call struct_userbc  (ix,iy,iz,iface,ieg)
            v1(ix,iy,iz) = trx
            v2(ix,iy,iz) = try
            v3(ix,iy,iz) = trz
  106    continue
         return

      elseif (cb.eq.'sl ' .or. cb.eq.'shl') then

         do 107 iz=kz1,kz2
         do 107 iy=ky1,ky2
         do 107 ix=kx1,kx2
            if (optlevel.le.2) call nekasgn (ix,iy,iz,iel)
            call struct_userbc  (ix,iy,iz,iface,ieg)
            v1(ix,iy,iz) = trn
            v2(ix,iy,iz) = tr1
            v3(ix,iy,iz) = tr2
  107    continue
C
      elseif (cb.eq.'ms ') then
c
         do 108 iz=kz1,kz2
         do 108 iy=ky1,ky2
         do 108 ix=kx1,kx2
            if (optlevel.le.2) call nekasgn (ix,iy,iz,iel)
            call struct_userbc(ix,iy,iz,iface,ieg)
            v1(ix,iy,iz) = -pa
            v2(ix,iy,iz) = tr1
            v3(ix,iy,iz) = tr2
  108    continue

      endif

      return
      end subroutine struct_faceiv
c-----------------------------------------------------------------------
      subroutine struct_userbc (ix,iy,iz,iside,ieg)

      implicit none
       
      include 'SIZE'
      include 'NEKUSE'          ! UX, UY, UZ, TEMP, X, Y, PA

      integer ix,iy,iz,iside,ieg
      real amp

      ux = 0.0
      uy = 0.
      uz = 0.

!     traction forces on the structure
      amp = 1.0e-2
      trx = amp*exp(-((y-0.8)/0.025)**2 - ((x-3.0)/0.1)**2) !*cos(x)
      try = 0.
      trz = 0.
      
      return
      end
c -----------------------------------------------------------------------
      subroutine struct_cresvif (resv1,resv2,resv3)
C
C     Compute startresidual/right-hand-side in the velocity solver
C

      implicit none

      include 'SIZE'
      include 'STRUCT'
      include 'NEKNEK'        ! igeom
      include 'INPUT'         ! iftran
      include 'SOLN'
      include 'MASS'          ! BM1
      include 'TSTEP'         ! ifield

!      include 'TOTAL'
      real           resv1 (lx1,ly1,lz1,lelv)
      real           resv2 (lx1,ly1,lz1,lelv)
      real           resv3 (lx1,ly1,lz1,lelv)

      logical ifmsk,ifdss


!     Save fields to lag arrays
      if (igeom.eq.2) then

!       This is actually displacement            
        call lagvel

        call struct_bcdirvc(vx,vy,vz,v1mask,v2mask,v3mask)
        call struct_bcneutr       ! add traction to rhs 

        call opcopy(resv1,resv2,resv3,vx,vy,vz)

        ifmsk = .false.
        ifdss = .false.

        call struct_Ax(resv1,resv2,resv3,ifdss,ifmsk)

!       solve for increment 
        call opsub2(resv1,resv2,resv3,bfx,bfy,bfz)

        call opcmult(resv1,resv2,resv3,-1.)

!       debugging 
!        call opcopy(ts4,ts5,ts6,resv1,resv2,resv3)

      endif

      return
      end subroutine struct_cresvif
c-----------------------------------------------------------------------
      subroutine struct_bcdirvc(v1,v2,v3,mask1,mask2,mask3)

C     Apply Dirichlet boundary conditions to surface of vector (V1,V2,V3).
C     Use IFIELD as a guide to which boundary conditions are to be applied.

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'TOPOL'
      include 'CTIMER'
      common /scruz/ tmp1(lx1,ly1,lz1,lelv)
     $             , tmp2(lx1,ly1,lz1,lelv)
     $             , tmp3(lx1,ly1,lz1,lelv)
      common /scrmg/ tmq1(lx1,ly1,lz1,lelv)
     $             , tmq2(lx1,ly1,lz1,lelv)
     $             , tmq3(lx1,ly1,lz1,lelv)
c
      real v1(lx1,ly1,lz1,lelv),v2(lx1,ly1,lz1,lelv)
     $    ,v3(lx1,ly1,lz1,lelv)
      real mask1(lx1,ly1,lz1,lelv),mask2(lx1,ly1,lz1,lelv)
     $    ,mask3(lx1,ly1,lz1,lelv)
c
      common  /nekcb/ cb
      character cb*3
      character*1 cb1(3)
      equivalence (cb1,cb)
c
      logical ifonbc

c
      ifonbc = .false.
c
      if (icalld.eq.0) then
         tusbc=0.0
         nusbc=0
         icalld=icalld+1
      endif
      nusbc=nusbc+1
      etime1=dnekclock()

      nfaces=2*ldim
      nxyz  =lx1*ly1*lz1
      nel   =nelfld(ifield)
      ntot  =nxyz*nel
c
      call rzero(tmp1,ntot)
      call rzero(tmp2,ntot)
      if (if3d) call rzero(tmp3,ntot)
c
c     velocity boundary conditions
c
c     write(6,*) 'bcdirv: ifield',ifield
      do 2100 isweep=1,2
         do 2000 ie=1,nel
         do 2000 iface=1,nfaces
            cb  = cbc(iface,ie,ifield)
            bc1 = bc(1,iface,ie,ifield)
            bc2 = bc(2,iface,ie,ifield)
            bc3 = bc(3,iface,ie,ifield)

            if (cb.eq.'V  ' .or. cb.eq.'VL '  .or.
     $          cb.eq.'WS ' .or. cb.eq.'WSL') then
               call facev (tmp1,ie,iface,bc1,lx1,ly1,lz1)
               call facev (tmp2,ie,iface,bc2,lx1,ly1,lz1)
               if (if3d) call facev (tmp3,ie,iface,bc3,lx1,ly1,lz1)
               if ( ifqinp(iface,ie) )
     $         call globrot (tmp1(1,1,1,ie),tmp2(1,1,1,ie),
     $                       tmp3(1,1,1,ie),ie,iface)
            endif

            IF (CB.EQ.'v  ' .OR. CB.EQ.'vl ' .OR. 
     $          CB.EQ.'ws ' .OR. CB.EQ.'wsl' .OR.
     $          CB.EQ.'mv ' .OR. CB.EQ.'mvn' .OR.
     $          cb1(1).eq.'d'.or.cb1(2).eq.'d'.or.cb1(3).eq.'d') then

                call struct_faceiv (cb,tmp1(1,1,1,ie),tmp2(1,1,1,ie),
     $                       tmp3(1,1,1,ie),ie,iface,lx1,ly1,lz1)

                if ( ifqinp(iface,ie) )
     $          call globrot (tmp1(1,1,1,ie),tmp2(1,1,1,ie),
     $                        tmp3(1,1,1,ie),ie,iface)
            endif

 2000    continue
         do 2010 ie=1,nel
         do 2010 iface=1,nfaces
            if (cbc(iface,ie,ifield).eq.'W  ') then
               call facev (tmp1,ie,iface,0.0,lx1,ly1,lz1)
               call facev (tmp2,ie,iface,0.0,lx1,ly1,lz1)
               if (if3d) call facev (tmp3,ie,iface,0.0,lx1,ly1,lz1)
            endif
 2010    continue
C
C        Take care of Neumann-Dirichlet shared edges...
C
         if (isweep.eq.1) then
            call opdsop(tmp1,tmp2,tmp3,'MXA')
         else
            call opdsop(tmp1,tmp2,tmp3,'MNA')
         endif
 2100 continue
C
C     Copy temporary array to velocity arrays.
C
      if (.not.ifstrs ) then
         call col2(v1,mask1,ntot)
         call col2(v2,mask2,ntot)
         if (if3d) call col2(v3,mask3,ntot)
      else
         call rmask (v1,v2,v3,nelv)
      endif

      call add2(v1,tmp1,ntot)
      call add2(v2,tmp2,ntot)
      if (if3d) call add2(v3,tmp3,ntot)

!      if (ifneknekc) call fix_surface_flux

      tusbc=tusbc+(dnekclock()-etime1)

      return
      end subroutine 
c-----------------------------------------------------------------------
      subroutine struct_elast3d_e(w1,w2,w3,u1,u2,u3,e,g,lambda)

!     Taken from NekExamples            
!     Apply elasticity operator to u:   w = Eu

      implicit none

      include 'SIZE'
      include 'GEOM'    ! jacmi,rxm1, etc.
      include 'INPUT'   ! if3d
      include 'MASS'    ! bm1
      include 'SOLN'    ! vtrans
      include 'TSTEP'   ! dt
      include 'WZ'      ! w3m1

      real w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      integer e,i

!      common /elastr/ lambda,g,gamma_n,beta_n,rho,damp
      real lambda,g
      real a,b,ba

      integer lxyz
      parameter (lxyz=lx1*ly1*lz1)

      real ur,us,ut,vr,vs,vt,wr,ws,wt
      common /ctmp1/ ur(lxyz),us(lxyz),ut(lxyz)
     $             , vr(lxyz),vs(lxyz),vt(lxyz)
     $             , wr(lxyz),ws(lxyz),wt(lxyz)

      real u11,u21,u31,u12,u22,u32,u13,u23,u33
      real u1d,u2d,u3d,div

      real w

      call gradl_rst(ur,us,ut,u1,nx1,if3d) ! Grad on GLL
      call gradl_rst(vr,vs,vt,u2,nx1,if3d)
      call gradl_rst(wr,ws,wt,u3,nx1,if3d)

      a  = g
      b  = lambda 
      ba = b/a

      do i=1,lxyz

c        uij := jac*( du_i / dx_j )

         u11=ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e)+ut(i)*txm1(i,1,1,e)
         u21=vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e)+vt(i)*txm1(i,1,1,e)
         u31=wr(i)*rxm1(i,1,1,e)+ws(i)*sxm1(i,1,1,e)+wt(i)*txm1(i,1,1,e)
         u12=ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e)+ut(i)*tym1(i,1,1,e)
         u22=vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e)+vt(i)*tym1(i,1,1,e)
         u32=wr(i)*rym1(i,1,1,e)+ws(i)*sym1(i,1,1,e)+wt(i)*tym1(i,1,1,e)
         u13=ur(i)*rzm1(i,1,1,e)+us(i)*szm1(i,1,1,e)+ut(i)*tzm1(i,1,1,e)
         u23=vr(i)*rzm1(i,1,1,e)+vs(i)*szm1(i,1,1,e)+vt(i)*tzm1(i,1,1,e)
         u33=wr(i)*rzm1(i,1,1,e)+ws(i)*szm1(i,1,1,e)+wt(i)*tzm1(i,1,1,e)

         div = u11 + u22 + u33

         u1d = 2*u11 + ba*div
         u2d = 2*u22 + ba*div
         u3d = 2*u33 + ba*div

         w   = a*w3m1(i,1,1)*jacmi(i,e)  ! note, ry has jac in it.

         u12=u12+u21
         u13=u13+u31
         u23=u23+u32

         ur(i)=(u1d*rxm1(i,1,1,e)+u12*rym1(i,1,1,e)+u13*rzm1(i,1,1,e))*w
         us(i)=(u1d*sxm1(i,1,1,e)+u12*sym1(i,1,1,e)+u13*szm1(i,1,1,e))*w
         ut(i)=(u1d*txm1(i,1,1,e)+u12*tym1(i,1,1,e)+u13*tzm1(i,1,1,e))*w
         vr(i)=(u12*rxm1(i,1,1,e)+u2d*rym1(i,1,1,e)+u23*rzm1(i,1,1,e))*w
         vs(i)=(u12*sxm1(i,1,1,e)+u2d*sym1(i,1,1,e)+u23*szm1(i,1,1,e))*w
         vt(i)=(u12*txm1(i,1,1,e)+u2d*tym1(i,1,1,e)+u23*tzm1(i,1,1,e))*w
         wr(i)=(u13*rxm1(i,1,1,e)+u23*rym1(i,1,1,e)+u3d*rzm1(i,1,1,e))*w
         ws(i)=(u13*sxm1(i,1,1,e)+u23*sym1(i,1,1,e)+u3d*szm1(i,1,1,e))*w
         wt(i)=(u13*txm1(i,1,1,e)+u23*tym1(i,1,1,e)+u3d*tzm1(i,1,1,e))*w

      enddo

2     format(2i5,4e15.6,a7)

      call gradl_rst_t(w1,ur,us,ut,nx1,if3d)
      call gradl_rst_t(w2,vr,vs,vt,nx1,if3d)
      call gradl_rst_t(w3,wr,ws,wt,nx1,if3d)

      return
      end subroutine struct_elast3d_e
c-----------------------------------------------------------------------
      subroutine struct_elast2d_e(w1,w2,w3,u1,u2,u3,e,lambda,g)

!     Taken from NekExamples
!     Apply elasticity operator to u:   w = Eu

      implicit none

      include 'SIZE'
      include 'GEOM'    ! jacmi,rxm1, etc.
      include 'INPUT'   ! if3d
      include 'MASS'    ! bm1
      include 'SOLN'    ! vtrans
      include 'TSTEP'   ! dt
      include 'WZ'      ! w3m1


      real w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      integer e,i

!      common /elastr/ lambda,g,gamma_n,beta_n,rho,damp
      real lambda,g
      real a,b,ba

      integer lxyz
      parameter (lxyz=lx1*ly1*lz1)

      real ur,us,ut,vr,vs,vt,wr,ws,wt
      common /ctmp1/ ur(lxyz),us(lxyz),ut(lxyz)
     $             , vr(lxyz),vs(lxyz),vt(lxyz)
     $             , wr(lxyz),ws(lxyz),wt(lxyz)

      real u11,u21,u12,u22
      real u1d,u2d,div

      real w


      call gradl_rst(ur,us,ut,u1,nx1,if3d) ! Grad on GLL
      call gradl_rst(vr,vs,vt,u2,nx1,if3d)

      a =          g
      b =      lambda 
      ba = b/a

      do i=1,lxyz

c        uij := jac*( du_i / dx_j )

         u11=ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e)
         u21=vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e)
         u12=ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e)
         u22=vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e)

         div = u11 + u22 

         u1d = 2*u11 + ba*div
         u2d = 2*u22 + ba*div

         w   = a*w3m1(i,1,1)*jacmi(i,e)  ! note, ry has jac in it.

         u12=u12+u21

         ur(i)=(u1d*rxm1(i,1,1,e)+u12*rym1(i,1,1,e))*w
         us(i)=(u1d*sxm1(i,1,1,e)+u12*sym1(i,1,1,e))*w
         vr(i)=(u12*rxm1(i,1,1,e)+u2d*rym1(i,1,1,e))*w
         vs(i)=(u12*sxm1(i,1,1,e)+u2d*sym1(i,1,1,e))*w

      enddo

      call gradl_rst_t(w1,ur,us,ut,nx1,if3d)
      call gradl_rst_t(w2,vr,vs,vt,nx1,if3d)

      return
      end subroutine struct_elast2d_e
C---------------------------------------------------------------------------
      subroutine struct_elast(w1,w2,w3,u1,u2,u3,lambda,g,ifmsk,ifdss)

      implicit none      

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'    ! v1mask,...

      real w1(1),w2(1),w3(1),u1(1),u2(1),u3(1),lambda,g
      logical ifmsk,ifdss

      integer k,e,nxyz
      real vsum

      vsum = 0

      nxyz = nx1*ny1*nz1
      k = 1
      do e=1,nelv

         if (if3d) then
            call struct_elast3d_e(w1(k),w2(k),w3(k),u1(k),u2(k),u3(k),e,
     $                            lambda,g)
         else
            call struct_elast2d_e(w1(k),w2(k),w3(k),u1(k),u2(k),u3(k),e,
     $                            lambda,g)
         end if

         if (ifmsk) then
            call col2 (w1(k),v1mask(1,1,1,e),nxyz)
            call col2 (w2(k),v2mask(1,1,1,e),nxyz)
            if (if3d) call col2 (w3(k),v3mask(1,1,1,e),nxyz)
         endif

         k = k+nxyz

      enddo

      if (ifdss) call opdssum(w1,w2,w3)

      return
      end subroutine struct_elast
c-----------------------------------------------------------------------

      subroutine solve_elasticity(rv1,rv2,rv3)

!     Solve the structural system using GMRES

      implicit none

      include 'SIZE'
      include 'INPUT'   ! iftran
      include 'SOLN'
      include 'MASS'    ! bm1,v1mask,...
      include 'TSTEP'
      include 'STRUCT'

      integer i,j,k,ipass,npass,miter
      real betav,betax,beta,resid0,resid01

      integer nt

      logical ifdss     ! if dssum in elast
      logical ifmsk     ! Not sure about this yet
      integer nk

      real rv1(lx1*ly1*lz1*lelv)
      real rv2(lx1*ly1*lz1*lelv)
      real rv3(lx1*ly1*lz1*lelv)

      real w1(lx1*ly1*lz1*lelv)
      real w2(lx1*ly1*lz1*lelv)
      real w3(lx1*ly1*lz1*lelv)

      real w4(lx1*ly1*lz1*lelv)
      real w5(lx1*ly1*lz1*lelv)
      real w6(lx1*ly1*lz1*lelv)

      real resv(lx1*ly1*lz1*lelv,3)       ! initial residuals       

      real solv(lx1*ly1*lz1*lelv,3)       ! solution       

      real Axv(lx1*ly1*lz1*lelv,3)       ! intermediate Ax (vel)       
     
      real h2(lx1*ly1*lz1*lelv)
      real bmm(lx1*ly1*lz1*lelv)

      real op_glsc2_wt

      real young,g,nu,lambda
      character*32 outfmt

      real rad,cs,sn          ! Given's rotations
      real r1(2),r2(2)
      real rhs(struct_nkryl+1)
      real lsq_resid(struct_nkryl)      ! residual after each iteration
      real rel_res

      real tol
      integer ic,ncycl    ! no of gmres cycles
      integer ikryl       ! current krylov space size

      logical ifconv
      real red_soln(struct_nkryl)
      real red_soln2(struct_nkryl+1)
     
      real s
      real const


!     Taking relative tolerance      
      tol = 1.0e-6
      npass = 1

      nk = struct_nkryl
      call rzero(struct_hessen,nk*(nk+1))
      call rzero(struct_R,nk*(nk+1))
      call rzero(struct_GivensC,nk+1)
      call rzero(struct_GivensS,nk+1)
      call rzero(rhs,nk+1)
      call rzero(lsq_resid,nk+1)

      miter  = nk 
      ncycl  = 100

!     Continuous edges      
      call opdssum(rv1,rv2,rv3)
      call opcol2(rv1,rv2,rv3,v1mask,v2mask,v3mask)

      call opcopy(resv(1,1),resv(1,2),resv(1,3),rv1,rv2,rv3)


!     Use this for inner products
      nt=nx1*ny1*nz1*nelv 
      call copy(bmm,bm1,nt)
      call col2(bmm,vmult,nt)

      ifconv = .false.        ! if converged

      ic = 0
      do while ((ic.lt.ncycl).and.(.not.ifconv))

        ic=ic+1
!       Norm      
        betax  = op_glsc2_wt(rv1,rv2,rv3,rv1,rv2,rv3,bmm)

        beta   = sqrt(betax)

        if (ic.eq.1) resid0 = beta
        resid01 = beta

        rhs(1) = beta

!       normalize
        call opcmult(rv1,rv2,rv3,1./beta)

!       save v part of first vector.
        call opcopy(struct_krylv(1,1,1),struct_krylv(1,2,1),
     $                struct_krylv(1,3,1),rv1,rv2,rv3)

        ikryl = 1
        do i=1,miter

          ifmsk = .true.
          ifdss = .true.

!         y = Ax
          call struct_Ax(rv1,rv2,rv3,ifdss,ifmsk)

!         Grahm-Schmidt 
          do ipass=1,npass
            call opcopy(w1,w2,w3,rv1,rv2,rv3)
            do j=1,i
            
              beta  = op_glsc2_wt(struct_krylv(1,1,j),
     $                 struct_krylv(1,2,j),struct_krylv(1,3,j),
     $                  w1,w2,w3,bmm)
                 
              struct_hessen(j,i) = struct_hessen(j,i)+beta
              struct_R(j,i) = struct_R(j,i)+beta

              call opadd2cm(rv1,rv2,rv3,struct_krylv(1,1,j),
     $                   struct_krylv(1,2,j),struct_krylv(1,3,j),-beta)


            enddo  ! j
          enddo  ! ipass

!         Residual after ortho-normalization 
          betax = op_glsc2_wt(rv1,rv2,rv3,rv1,rv2,rv3,bmm)
          
          beta  = sqrt(betax)
          struct_hessen(i+1,i) = beta
          struct_R(i+1,i)      = beta

!         normalize v
          call opcmult(rv1,rv2,rv3,1./beta)

!         save v part of the vector        
          call opcopy(struct_krylv(1,1,i+1),struct_krylv(1,2,i+1),
     $                struct_krylv(1,3,i+1),rv1,rv2,rv3)

!         Increase krylov size            
          ikryl = ikryl + 1

!         QR decomposition of Hessenberg matrix using Given's rotations
!         Do previous rotations        
          do j=1,i-1
            r1(1) = struct_R(j,i)
            r1(2) = struct_R(j+1,i)

!           Apply rotation  
            r2(1) = struct_GivensC(j)*r1(1) + struct_GivensS(j)*r1(2)
            r2(2) = -struct_GivensS(j)*r1(1) + struct_GivensC(j)*r1(2)

            struct_R(j,i)=r2(1)
            struct_R(j+1,i)=r2(2)
          enddo         

!         Rotation for the ith column        
          r1(1) = struct_R(i,i)
          r1(2) = struct_R(i+1,i)
          rad   = sqrt(r1(1)**2 + r1(2)**2)
          cs    = r1(1)/rad 
          sn    = r1(2)/rad

          struct_GivensC(i)=cs 
          struct_GivensS(i)=sn

!         Apply rotations 
          r2(1) =  struct_GivensC(i)*r1(1) + struct_GivensS(i)*r1(2)
          r2(2) = -struct_GivensS(i)*r1(1) + struct_GivensC(i)*r1(2)

          struct_R(i,i)=r2(1)
          struct_R(i+1,i)=r2(2)

!         Rotate rhs
          r1(1) =  rhs(i)
          r1(2) =  rhs(i+1)
          r2(1) =  struct_GivensC(i)*r1(1) + struct_GivensS(i)*r1(2)
          r2(2) = -struct_GivensS(i)*r1(1) + struct_GivensC(i)*r1(2)
          rhs(i)   = r2(1)
          rhs(i+1) = r2(2)

          lsq_resid(i)=rhs(i+1)

          rel_res = abs(rhs(i+1)/resid0)

          if (rel_res.lt.tol) then
            ifconv = .true.
            exit ! move out of loop 
          endif              

        enddo       ! i=1,miter 

!!       All below is debugging
!!--------------------------------------------------       
!!       for writing the hessenberg matrix      
!        call blank(outfmt,32)
!        write(outfmt,'(A7,I2,A13)') '(A3,2x,',miter,'(E18.8E2,2x))'
!
!!        write(6,*) outfmt
!
!        if (nid.eq.0) then
!          do i=1,miter+1
!            write(6,outfmt), 'Hes',(struct_hessen(i,j),j=1,miter)
!          enddo  
!        endif
!!       Write out Upper triangular R
!        if (nid.eq.0) then
!          do i=1,miter+1
!            write(6,outfmt), 'QRR',(struct_R(i,j),j=1,miter)
!          enddo  
!        endif
!
        call blank(outfmt,32)
        write(outfmt,'(A7,I2,A13)') '(A3,2x,',miter+1,'(E18.8E2,2x))'
!!     
!        write(6,outfmt), 'rhs',(rhs(j),j=1,miter+1)
!        write(6,outfmt), 'lsq',(lsq_resid(j),j=1,miter)

        write(6,'(A18,2x,2(I3,2x),3(E13.5E2,2x))')
     $      'ic,ikryl,rr,r01,r0',ic,ikryl,rel_res,resid01,resid0

!        call check_ortho(ikryl,bmm)

!--------------------------------------------------
        

!       Generate solution
!       Back substitution
        do j=ikryl-1,1,-1           ! row
          s = 0.
          do k=ikryl-1,j+1,-1       ! column
            s = s + struct_R(j,k)*red_soln(k)
          enddo
          red_soln(j)=(rhs(j)-s)/struct_R(j,j)
        enddo          

!       debugging           
!        write(6,outfmt), 'sol',(red_soln(j),j=1,ikryl-1)

!       Recombine Krylov space for solution
!       soln = Vn*yn
!        call opzero(solv(1,1),solv(1,2),solv(1,3))
        do i=1,ikryl-1
          call opadd2cm(solv(1,1),solv(1,2),solv(1,3),
     $                  struct_krylv(1,1,i),
     $                  struct_krylv(1,2,i),
     $                  struct_krylv(1,3,i),red_soln(i))
        enddo


        if (.not.ifconv) then
!         restart loop
!              _        
!         yn+1 = H*yn
!          call mxm(struct_hessen,ikryl,red_soln,
!     $                    ikryl-1,red_soln2,1)
!         Somehow this is bugging out              
          do i=1,ikryl
            red_soln2(i)=0.
            do j=1,ikryl-1
              red_soln2(i)=red_soln2(i)+struct_hessen(i,j)*red_soln(j)
            enddo
          enddo

!         debugging           
!          write(6,outfmt), 'sl2',(red_soln2(j),j=1,ikryl)

          call opzero(Axv(1,1),Axv(1,2),Axv(1,3))
!         Recombine Krylov vectors        
!         Ax = Vn+1*yn+1
          do i=1,ikryl
            call opadd2cm(Axv(1,1),Axv(1,2),Axv(1,3),
     $                    struct_krylv(1,1,i),
     $                    struct_krylv(1,2,i),
     $                    struct_krylv(1,3,i),red_soln2(i))
          enddo

!         Remove part solution
!         b = (b - Ax)
          call opsub2(resv(1,1),resv(1,2),resv(1,3),
     $                Axv(1,1),Axv(1,2),Axv(1,3))          

          call opcopy(rv1,rv2,rv3,resv(1,1),resv(1,2),resv(1,3))

          call rzero(struct_hessen,nk*(nk+1))
          call rzero(struct_R,nk*(nk+1))

          call rzero(rhs,nk+1)
          call rzero(lsq_resid,nk+1)


        else
          ifconv = .true.
          exit
        endif

      enddo       ! ic=1,ncycl 

      call opcopy(rv1,rv2,rv3,solv(1,1),solv(1,2),solv(1,3))


      return
      end subroutine solve_elasticity            

c------------------------------------------------------------------------

      subroutine struct_Ax(rv1,rv2,rv3,ifdss,ifmsk)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'SOLN'          ! v1mask,...
      include 'TSTEP'         ! dt

      real rv1(lx1*ly1*lz1*lelv)
      real rv2(lx1*ly1*lz1*lelv)
      real rv3(lx1*ly1*lz1*lelv)

      real w1(lx1*ly1*lz1*lelv)
      real w2(lx1*ly1*lz1*lelv)
      real w3(lx1*ly1*lz1*lelv)

      real w4(lx1*ly1*lz1*lelv)
      real w5(lx1*ly1*lz1*lelv)
      real w6(lx1*ly1*lz1*lelv)

      real h2(lx1*ly1*lz1*lelv)

      logical ifdss,ifmsk
      logical ifdss1,ifmsk1

      real young,g,nu,lambda
      real damp_const
      real const

      integer nt,ifld

!     move somewhere else      
      young=100.
      nu=0.3

      lambda = young*nu/( (1+nu)*(1-2*nu) ) ! for 3D, and 2D plane strain

      g      = .5*young/(1+nu)
!----------------------------------------     

!     Apply elasticity operator. Assuming transient simulation 
!     y = Ax

      ifld = 1

      nt = lx1*ly1*lz1*nelv
      const = 4./(DT**2)
      call cmult2(h2,vtrans(1,1,1,1,ifld),const,nt)
      call col2(h2,bm1,nt)

      call opcolv3(w1,w2,w3,rv1,rv2,rv3,h2)

!!     damping term
!      damp_const = 0.
!      const = damp_const*2./DT
!      call col2(h2,bm1,nt)
!      call cmult(h2,const,nt)
!
!      call opcolv3(w4,w5,w6,rv1,rv2,rv3,h2)
!
!      call add2(w1,w2,w3,w4,w5,w6)

      ifdss1=.false.
      ifmsk1=.false.

!     Elasticity operator           
      call struct_elast(w4,w5,w6,rv1,rv2,rv3,lambda,g,ifdss1,ifmsk1)

!     w1+w4,... 
      call opadd2(w1,w2,w3,w4,w5,w6)

      call opcopy(rv1,rv2,rv3,w1,w2,w3)

!!     Make continuous
      if (ifdss) call opdssum(rv1,rv2,rv3)
      if (ifmsk) call opcol2(rv1,rv2,rv3,v1mask,v2mask,v3mask)

      return
      end subroutine struct_Ax

!---------------------------------------------------------------------- 


      subroutine struct_seth2(h2)

      implicit none

      include 'SIZE'
      include 'TSTEP'   ! DT
      include 'MASS'    ! bm1
      include 'SOLN'    ! vtrans

      real h2(lx1*ly1*lz1*lelv)
      integer nt
      real const

      nt = lx1*ly1*lz1*nelv
      const = 1./(DT**2)
      call cmult2(h2,vtrans(1,1,1,1,1),const,nt)
      call col2(h2,bm1,nt)

      return
      end subroutine struct_seth2
!---------------------------------------------------------------------- 

      subroutine struct_lagdisp

!     Keep old velocity field(s)

      implicit none            

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      integer nt,ilag

      nt = lx1*ly1*lz1*nelv

      do 100 ilag=3-1,2,-1
         call copy (vxlag (1,1,1,1,ilag),vxlag (1,1,1,1,ilag-1),nt)
         call copy (vylag (1,1,1,1,ilag),vylag (1,1,1,1,ilag-1),nt)
         if (ldim.eq.3)
     $   call copy (vzlag (1,1,1,1,ilag),vzlag (1,1,1,1,ilag-1),nt)
 100  continue
c
      call opcopy (vxlag,vylag,vzlag,vx,vy,vz)
c
      return
      end

!---------------------------------------------------------------------- 

      subroutine update_fields(delx,dely,delz)

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'STRUCT'
      include 'TSTEP'

      real resv1,resv2,resv3,resx1,resx2,resx3,h2
      common /scrns/  resv1 (lx1,ly1,lz1,lelv)
     $ ,              resv2 (lx1,ly1,lz1,lelv)
     $ ,              resv3 (lx1,ly1,lz1,lelv)
     $ ,              resx1 (lx1,ly1,lz1,lelv)
     $ ,              resx2 (lx1,ly1,lz1,lelv)
     $ ,              resx3 (lx1,ly1,lz1,lelv)
     $ ,              h2    (lx1,ly1,lz1,lelv)

      real delx(1),dely(2),delz(3)

      real const


!     New displacements          
      call opadd2(vx,vy,vz,delx,dely,delz)

!     Un+1 - Un
      call opcopy(resv1,resv2,resv3,vx,vy,vz)
      call sub2(resv1,resv2,resv3,vxlag,vylag,vzlag)

!     New Acceleration            
!     An+1 = -[An - (4/DT^2)*(Un+1 - Un) + 4/DT*Vn]
      const = 4./DT
      call opadd2cm(accx,accy,accz,velx,vely,velz,const)

      const = -4./(DT**2)
      call opadd2cm(accx,accy,accz,resv1,resv2,resv3,const)
      call opcmult(accx,accy,accz,-1.)


!     New velocity            
!     Vn+1 = -[Vn - 2/DT*(Un+1 - Un)]
      const=-2./DT
      call opadd2cm(velx,vely,velz,resv1,resv2,resv3,const)
      call opcmult(velx,vely,velz,-1.)


      return
      end subroutine update_fields
!---------------------------------------------------------------------- 




