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

      real resv1,resv2,resv3,dv1,dv2,dv3,h1,h2

      common /scrns/  resv1 (lx1,ly1,lz1,lelv)
     $ ,              resv2 (lx1,ly1,lz1,lelv)
     $ ,              resv3 (lx1,ly1,lz1,lelv)
     $ ,              dv1   (lx1,ly1,lz1,lelv)
     $ ,              dv2   (lx1,ly1,lz1,lelv)
     $ ,              dv3   (lx1,ly1,lz1,lelv)

      real            resx4 (lx1,ly1,lz1,lelv)
     $ ,              resx5 (lx1,ly1,lz1,lelv)
     $ ,              resx6 (lx1,ly1,lz1,lelv)

      common /scrvh/  h1    (lx1,ly1,lz1,lelv)
     $ ,              h2    (lx1,ly1,lz1,lelv)


!     For now.
      ifgeom = .false.


!     First Structural solve
      do igeom=1,ngeom

        if (igeom.eq.1) then
          call struct_makef
        else

!          call struct_sethlm
          call struct_cresvif(resv1,resv2,resv3,resx1,resx2,resx3)

          call solve_elasticity                          

!          call struct(igeom)
        endif
      enddo



      return
      end subroutine plan_s
!----------------------------------------------------------------------
      subroutine struct_sethlm

      implicit none
c     Set the variable property arrays H1 and H2
c     in the Helmholtz equation.
c     INTLOC =      integration type

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      include 'STRUCT'

      real dtbd
      integer nel,ntot1

      nel   = nelv ! nelfld(ifield)
      ntot1 = lx1*ly1*lz1*nel

      if (iftran) then
         dtbd = bd(1)/dt

!        Equation 1         
         call rzero  (e1hv1,ntot1)
         call cmult2 (e1hv2,vtrans(1,1,1,1,ifield),dtbd,ntot1)

         call copy  (e1hx1,vdiff(1,1,1,1,ifield),ntot1)
         call rzero (e1hx2,ntot1)

!        Equation 2
!        Adding density here because I add density into the body 
!        forcing for this equation. Which probably should just be
!        removed.
         call cmult2 (e2hv1,vtrans(1,1,1,1,ifield),-1.,ntot1)
         call rzero  (e2hv2,ntot1) 

         call rzero  (e1hx1,ntot1)
         call cmult2 (e1hx2,vtrans(1,1,1,1,ifield),dtbd,ntot1)

      else

!        Equation 1
         call rzero (e1hv1,ntot1)
         call rzero (e1hv2,ntot1)

         call copy  (e1hx1,vdiff (1,1,1,1,ifield),ntot1)
         call rzero (e1hx2,ntot1)

!        Equation 2         
         call rzero (e2hx1,ntot1)
         call rzero (e2hx2,ntot1)

         call rzero (e2hv1,ntot1)
         call rzero (e2hv2,ntot1)

      endif


      return 
      end subroutine struct_sethlm
!----------------------------------------------------------------------       
      subroutine struct_makef

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'CTIMER'

      etime1 = dnekclock()

      call struct_makeuf                  ! body forcing in the structure
                                          ! Contains mass matrix

      if (iftran) call struct_makeabf     ! extrapolate forcing terms
      if (iftran) call struct_makebdf     ! Backward difference for time

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

      integer seqn     ! forcing equation  .
                       ! seqn=1        ==> V - div(sigma) = f1
                       !                   .
                       ! seqn=2        ==> X - V          = f2

c
      time = time-dt
      seqn = 1
      call struct_nekuf(bfx,bfy,bfz,seqn)
      call opcolv(bfx,bfy,bfz,bm1)

!     Displacement time derivative equation      
      if (iftran) then
        seqn = 2
        call struct_nekuf(struct_bfdx,struct_bfdy,struct_bfdz,seqn)
        call opcolv(struct_bfdx,struct_bfdy,struct_bfdz,bm1)
      endif        

      time = time+dt
C
      return
      end subroutine struct_makeuf
!----------------------------------------------------------------------       

      subroutine struct_nekuf (f1,f2,f3,seqn)

      implicit none            

      include 'SIZE'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'CTIMER'

      real f1 (lx1,ly1,lz1,lelv)
      real f2 (lx1,ly1,lz1,lelv)
      real f3 (lx1,ly1,lz1,lelv)
      integer seqn

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
            call struct_userf   (i,j,k,ielg,seqn)     ! define this
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
      subroutine struct_userf(ix,iy,iz,ieg,seqn)

!     User defined body forcing for structural equations
!     seqn == 1 ==> Main eqn 1
!     seqn == 2 ==> Position derivative eqn            
      implicit none

      include 'SIZE'
      include 'NEKUSE'

      integer ix,iy,iz,ieg,seqn      


      return
      end subroutine struct_userf
!----------------------------------------------------------------------

      subroutine struct_makeabf

C     Sum up contributions to kth order extrapolation scheme.

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      include 'STRUCT'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1,ab0,ab1,ab2

      ntot1 = lx1*ly1*lz1*nelv

      ab0 = ab(1)
      ab1 = ab(2)
      ab2 = ab(3)
      call add3s2 (ta1,abx1,abx2,ab1,ab2,ntot1)
      call add3s2 (ta2,aby1,aby2,ab1,ab2,ntot1)
      call copy   (abx2,abx1,ntot1)
      call copy   (aby2,aby1,ntot1)
      call copy   (abx1,bfx,ntot1)
      call copy   (aby1,bfy,ntot1)
      call add2s1 (bfx,ta1,ab0,ntot1)
      call add2s1 (bfy,ta2,ab0,ntot1)
      call col2   (bfx,vtrans,ntot1)          ! multiply by density
      call col2   (bfy,vtrans,ntot1)
      if (ldim.eq.3) then
         call add3s2 (ta3,struct_abz1,struct_abz2,ab1,ab2,ntot1)
         call copy   (struct_abz2,struct_abz1,ntot1)
         call copy   (struct_abz1,bfz,ntot1)
         call add2s1 (bfz,ta3,ab0,ntot1)
         call col2   (bfz,vtrans,ntot1)
      endif

!     also for second equation
      call add3s2 (ta1,struct_abx1,struct_abx2,ab1,ab2,ntot1)
      call add3s2 (ta2,struct_aby1,struct_aby2,ab1,ab2,ntot1)
      call copy   (struct_abx2,struct_abx1,ntot1)
      call copy   (struct_aby2,struct_aby1,ntot1)
      call copy   (struct_abx1,struct_bfdx,ntot1)
      call copy   (struct_aby1,struct_bfdy,ntot1)
      call add2s1 (struct_bfdx,ta1,ab0,ntot1)
      call add2s1 (struct_bfdy,ta2,ab0,ntot1)
      call col2   (struct_bfdx,vtrans,ntot1)   ! multiply by density
      call col2   (struct_bfdy,vtrans,ntot1)

      if (ldim.eq.3) then
         call add3s2 (ta3,struct_abz1,struct_abz2,ab1,ab2,ntot1)
         call copy   (struct_abz2,struct_abz1,ntot1)
         call copy   (struct_abz1,struct_bfdz,ntot1)
         call add2s1 (struct_bfdz,ta3,ab0,ntot1)
         call col2   (struct_bfdz,vtrans,ntot1)
      endif

      return
      end subroutine struct_makeabf
!----------------------------------------------------------------------
      subroutine struct_makebdf

c     Add contributions to F from lagged BD terms.

      implicit none            

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'STRUCT'

      real ta1,ta2,ta3,tb1,tb2,tb3,h2
      common /scrns/ ta1(lx1,ly1,lz1,lelv)
     $ ,             ta2(lx1,ly1,lz1,lelv)
     $ ,             ta3(lx1,ly1,lz1,lelv)
     $ ,             tb1(lx1,ly1,lz1,lelv)
     $ ,             tb2(lx1,ly1,lz1,lelv)
     $ ,             tb3(lx1,ly1,lz1,lelv)
     $ ,             h2 (lx1,ly1,lz1,lelv)

      integer ilag,ntot1
      real const

      ntot1 = lx1*ly1*lz1*nelv
      const = 1./dt

!     Equation 1
!     Need to think about this H2 
      call cmult2(h2,vtrans(1,1,1,1,ifield),const,ntot1)

      call opcolv3c (tb1,tb2,tb3,vx,vy,vz,bm1,bd(2))

      do 101 ilag=2,nbd
         if (ifgeom) then
            call opcolv3c(ta1,ta2,ta3,vxlag (1,1,1,1,ilag-1),
     $                                vylag (1,1,1,1,ilag-1),
     $                                vzlag (1,1,1,1,ilag-1),
     $                                bm1lag(1,1,1,1,ilag-1),bd(ilag+1))
         else
            call opcolv3c(ta1,ta2,ta3,vxlag (1,1,1,1,ilag-1),
     $                                vylag (1,1,1,1,ilag-1),
     $                                vzlag (1,1,1,1,ilag-1),
     $                                bm1                   ,bd(ilag+1))
         endif
         call opadd2  (tb1,tb2,tb3,ta1,ta2,ta3)
 101  continue
      call opadd2col (bfx,bfy,bfz,tb1,tb2,tb3,h2)


!     Equation 2
      call cmult2(h2,vtrans(1,1,1,1,ifield),const,ntot1)
      call opcolv3c (tb1,tb2,tb3,dispx,dispy,dispz,bm1,bd(2))

      do 102 ilag=2,nbd
         if (ifgeom) then
            call opcolv3c(ta1,ta2,ta3,dxlag (1,1,1,1,ilag-1),
     $                                dylag (1,1,1,1,ilag-1),
     $                                dzlag (1,1,1,1,ilag-1),
     $                                bm1lag(1,1,1,1,ilag-1),bd(ilag+1))
         else
            call opcolv3c(ta1,ta2,ta3,dxlag (1,1,1,1,ilag-1),
     $                                dylag (1,1,1,1,ilag-1),
     $                                dzlag (1,1,1,1,ilag-1),
     $                                bm1                   ,bd(ilag+1))
         endif
         call opadd2  (tb1,tb2,tb3,ta1,ta2,ta3)
 102  continue
      call opadd2col(struct_bfdx,struct_bfdy,struct_bfdz,
     $               tb1,tb2,tb3,h2)
     

      return
      end subroutine struct_makebdf            
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

      integer seqn

      seqn = 1

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
     $                           lx1,ly1,lz1,seqn)
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
c
  104 continue
c
      return
      end subroutine struct_bcneutr
c-----------------------------------------------------------------------

      subroutine struct_faceiv (cb,v1,v2,v3,iel,iface,nx,ny,nz,seqn)

c     Assign fortran function boundary conditions to 
c     face IFACE of element IEL for vector (V1,V2,V3).

      implicit none

      INCLUDE 'SIZE'
      INCLUDE 'NEKUSE'
      INCLUDE 'PARALLEL'

      integer nx,ny,nz,iel,ieg,iface,seqn

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
            call userbc  (ix,iy,iz,iface,ieg,seqn)
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
            call struct_userbc  (ix,iy,iz,iface,ieg,seqn)
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
            call struct_userbc  (ix,iy,iz,iface,ieg,seqn)
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
            call struct_userbc(ix,iy,iz,iface,ieg,seqn)
            v1(ix,iy,iz) = -pa
            v2(ix,iy,iz) = tr1
            v3(ix,iy,iz) = tr2
  108    continue

      endif

      return
      end subroutine struct_faceiv
c-----------------------------------------------------------------------
      subroutine struct_userbc (ix,iy,iz,iside,ieg,seqn)

      implicit none
       
      include 'SIZE'
      include 'NEKUSE'          ! UX, UY, UZ, TEMP, X, Y, PA

      integer ix,iy,iz,iside,ieg
      integer seqn
      real amp

      ux = 0.0
      uy = 0.
      uz = 0.

!     traction forces on the structure
      amp = 1.0
      trx = amp !*cos(x)
      try = 0.
      trz = 0.
      
      return
      end
c -----------------------------------------------------------------------
      subroutine struct_cresvif (resv1,resv2,resv3,
     $                           resx1,resx2,resx3)
C
C     Compute startresidual/right-hand-side in the velocity solver
C

      implicit none

      include 'SIZE'
      include 'STRUCT'
      include 'NEKNEK'        ! igeom
      include 'INPUT'         ! iftran
      include 'SOLN'

!      include 'TOTAL'
      real           resv1 (lx1,ly1,lz1,lelv)
      real           resv2 (lx1,ly1,lz1,lelv)
      real           resv3 (lx1,ly1,lz1,lelv)
      real           resx1 (lx1,ly1,lz1,lelv)
      real           resx2 (lx1,ly1,lz1,lelv)
      real           resx3 (lx1,ly1,lz1,lelv)

      real           h1    (lx1,ly1,lz1,1)
      real           h2    (lx1,ly1,lz1,1)

      real w1,w2,w3
      common /scruz/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)

      integer seqn

!      common /cgeom/ igeom


!     Save fields to lag arrays
      if (igeom.eq.2) call lagvel
      seqn = 1
      call struct_bcdirvc(vx,vy,vz,v1mask,v2mask,v3mask,seqn)
      call struct_bcneutr       ! add traction to rhs 

      call opcopy(resv1,resv2,resv3,bfx,bfy,bfz)

!     displacements      
      if (iftran.and.igeom.eq.2) call lagdx 

!     prabal. Need to create new masks for dx... 
      seqn = 2
      call struct_bcdirvc(dispx,dispy,dispz,v1mask,v2mask,v3mask,seqn)

      call opcopy(resx1,resx2,resx3,
     $            struct_bfdx,struct_bfdy,struct_bfdz)


      RETURN
      END
c-----------------------------------------------------------------------
      subroutine struct_bcdirvc(v1,v2,v3,mask1,mask2,mask3,seqn)

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

      integer seqn
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
     $                       tmp3(1,1,1,ie),ie,iface,lx1,ly1,lz1,seqn)

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
            call elast3d_e(w1(k),w2(k),w3(k),u1(k),u2(k),u3(k),e)
         else
            call elast2d_e(w1(k),w2(k),w3(k),u1(k),u2(k),u3(k),e)
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

      subroutine lagdx

!     Keep old velocity field(s) 
      
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      include 'STRUCT'

      integer ntot1,ilag

      ntot1 = lx1*ly1*lz1*nelv

!      do 100 ilag=nbdinp-1,2,-1
      do 100 ilag=3-1,2,-1
         call copy (dxlag(1,1,1,1,ilag),dxlag(1,1,1,1,ilag-1),ntot1)
         call copy (dylag(1,1,1,1,ilag),dylag(1,1,1,1,ilag-1),ntot1)
         if (ldim.eq.3)
     $   call copy (dzlag(1,1,1,1,ilag),dzlag(1,1,1,1,ilag-1),ntot1)
 100  continue

      call opcopy (vxlag,vylag,vzlag,dispx,dispy,dispz)

      return
      end subroutine lagdx
!---------------------------------------------------------------------- 
      subroutine solve_elasticity

!     Solve the structural system using GMRES

      implicit none

      include 'SIZE'
      include 'INPUT'   ! iftran
      include 'SOLN'
      include 'MASS'    ! bm1
      include 'TSTEP'

      integer i,j,miter
      real betav,betax,beta,resid0

      integer nt

      logical ifdss     ! if dssum in elast
      logical ifmask    ! Not sure about this yet

      miter = 2

!     Continuous edges      
      call opdssum(resv1,resv2,resv3)
      call opcolv(resv1,resv2,resv3,vmult)

      betav  = glsc3_wt(resv1,resv2,resv3,resv1,resv2,resv3,bm1)
      beta   = sqrt(betav)
      resid0 = beta
      if (iftran) then
!       Continuous edges      
        call opdssum(resv4,resv5,resv6)
        call opcolv(resv4,resv5,resv6,vmult)
           
        betax  = glsc3_wt(resv4,resv5,resv6,resv4,resv5,resv6,bm1)
        beta   = sqrt(betav + betax)
        resid0 = beta
!       normalize x        
        call opcmult(resv4,resv5,resv6,1./beta)
!       save x part of first vector        
        call opcopy(struct_krylx(1,1,1),struct_krylx(1,2,1),
     $              struct_krylx(1,3,1),resv4,resv5,resv6)
      endif

!     save v      
      call opcmult(resv1,resv2,resv3,1./beta)

!     save v part of first vector. This is the x part if its not a
!     transient simulation 
      call opcopy(struct_krylv(1,1,1),struct_krylv(1,2,1),
     $              struct_krylv(1,3,1),resv1,resv2,resv3)

      ifmsk = .false. ! ?
      ifdss = .true.

      nt = lx1*ly1*lz1*nelv
      call copy(htmp,bm1,nt)
      call cmult(htmp,bd(1)/DT,nt)

      do i=1,miter

!        Apply elasticity operator      
         if (iftran) then

!          bd(1)M/DT*v               
           call opcopy(w4,w5,w6,resv1,resv2,resv3)
           call opcolv(w4,w5,w6,htmp)
           call opcolv(w4,w5,w6,resv1,resv2,resv3)

!          Elasticity operator           
           call elast(w1,w2,w3,resv4,resv5,resv6,lambda,g,ifmsk,ifdss)
!          w1+w4,... 
           call opadd2(w1,w2,w3,w4,w5,w6)

!          2nd equation
           call opcolv3c(w4,w5,w6,resv1,resv2,resv3,bm1,-1.)
           call opadd2col(w4,w5,w6,resv4,resv5,resv6,htmp)


         else  
           call elast(w1,w2,w3,resv1,resv2,resv3,lambda,g,ifmsk,ifdss)
         endif


      enddo         



      return
      end subroutine solve_elasticity            

c------------------------------------------------------------------------




