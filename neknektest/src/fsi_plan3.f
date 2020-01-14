!======================================================================
!     Routines for fluid structure interaction
!     Author: Prabal S. Negi
!     Description: Solver routines for the Stokes' solve.
!
!======================================================================       

c-----------------------------------------------------------------------
      subroutine stokes_solve()

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      include 'FSI'
      include 'FSI_DEBUG'
C

      real resv1,resv2,resv3,dv1,dv2,dv3

      common /scrns/  resv1 (lx1,ly1,lz1,lelv)
     $ ,              resv2 (lx1,ly1,lz1,lelv)
     $ ,              resv3 (lx1,ly1,lz1,lelv)
     $ ,              dv1   (lx1,ly1,lz1,lelv)
     $ ,              dv2   (lx1,ly1,lz1,lelv)
     $ ,              dv3   (lx1,ly1,lz1,lelv)

      real h1,h2
      common /scrvh/  h1    (lx1,ly1,lz1,lelv)
     $ ,              h2    (lx1,ly1,lz1,lelv)

      logical fsi_hlm

      integer icalld
      save icalld
      data icalld /0/

      real dttol
      parameter (dttol=1.0E-12)  ! Tolerance to check if DT has changed

      integer intype
      real mstep
      real tolv_old,tolp_old
      real tolv_new,tolp_new
      real tolv_max,tolp_max


      tolv_old = param(22)
      tolp_old = param(21)


!      if (nio.eq.0) write(6,11) 'Solving Stokes'' step', fsi_timea,
!     $       fsi_timei
!   11 format(A21,1x,2(E12.5E2,1x))

      call opzero(sc_vx,sc_vy,sc_vz)
      call rzero(sc_pr,nx2*ny2*nz2*nelt)
      call opzero(resv1,resv2,resv3)

      intype = -1
      call sethlm  (h1,h2,intype)
      call fsi_cresvif (resv1,resv2,resv3,h1,h2)


      mstep = 0 !abs(param(94))
!      if (param(94).ne.0. .and. istep.ge.mstep) then
!        Projections. No projection of AMP right now. 
!        call ophinvpr(dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxh)
!      else
!       Just solve Helmholz equation right now
      fsi_hlm=.true. 
      call ophinv  (dv1,dv2,dv3,resv1,resv2,resv3,
     $                  h1,h2,tolhv,nmxh)

      call opadd2  (sc_vx,sc_vy,sc_vz,dv1,dv2,dv3)

      call fsi_incomprn(sc_vx,sc_vy,sc_vz,sc_pr)

         
      return
      end subroutine stokes_solve
C---------------------------------------------------------------------

      subroutine fsi_cresvif (resv1,resv2,resv3,h1,h2)

      implicit none            
      
      include 'SIZE'
      include 'NEKNEK'
      include 'FSI'
      include 'SOLN'    ! v1mask,...

      include 'FSI_DEBUG'

      real           resv1 (lx1,ly1,lz1,lelv)
      real           resv2 (lx1,ly1,lz1,lelv)
      real           resv3 (lx1,ly1,lz1,lelv)
      real           h1    (lx1,ly1,lz1,lelv)
      real           h2    (lx1,ly1,lz1,lelv)

      real w1,w2,w3      
      common /scruz/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)


!     Set Bdry conditions 
      call fsi_bcdirvc (sc_vx,sc_vy,sc_vz,
     $       v1mask,v2mask,v3mask)

!     call bcneutr

!     There's no RHS forcing
      call ophx (w1,w2,w3,sc_vx,sc_vy,sc_vz,h1,h2)

 
      call opzero  (resv1,resv2,resv3) 
      call opsub2  (resv1,resv2,resv3,w1,w2,w3)


      return
      end subroutine fsi_cresvif

!----------------------------------------------------------------------
      subroutine fsi_incomprn (ux,uy,uz,up)
c
c     Project U onto the closest incompressible field
c
c     Input:  U     := (ux,uy,uz)
c
c     Output: updated values of U, iproj, proj; and
c             up    := pressure currection req'd to impose div U = 0
c
c
c     Dependencies: ifield ==> which "density" (vtrans) is used.
c
c     Notes  1.  up is _not_ scaled by bd(1)/dt.  This should be done
c                external to incompr().
c
c            2.  up accounts _only_ for the perturbation pressure,
c                not the current pressure derived from extrapolation.
c

      implicit none

      include 'SIZE'
!      include 'TOTAL'
      include 'TSTEP'         ! DT
      include 'SOLN'
      include 'CTIMER'
      include 'FSI_DEBUG'

      real ux(1),uy(1),uz(1),up(1)

      real w1,w2,w3,dv1,dv2,dv3,dp
      common /scrns/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             dp    (lx2,ly2,lz2,lelv)

      real h1,h2      
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)

      real h2inv
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

!      parameter(nset = 1 + lbelv/lelv)
!      common /orthov/ pset(lx2*ly2*lz2*lelv*mxprev,nset)
!      common /orthbi/ nprv(2)

      logical ifprjp
      integer ntot1,ntot2,nprev,intype
      real bdti,scaledt,scaledi,dtb

!     If Project out previous pressure solutions?
      if (istep.ge.10) then
        ifprjp=.true.
      else
        nprev=0
        ifprjp=.false.
      endif

      npres  = npres+1
      etime1 = dnekclock()

      ntot1  = nx1*ny1*nz1*nelv
      ntot2  = nx2*ny2*nz2*nelv
      intype = 1

      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call invers2 (h2inv,h2,ntot1)

      call opdiv   (dp,ux,uy,uz)

      bdti = -bd(1)/dt
      call cmult   (dp,bdti,ntot2)

!      call add2col2(dp,bm2,usrdiv,ntot2) ! User-defined divergence.

      call ortho   (dp)

      if (ifprjp)  call fsi_setrhsp(dp,h1,h2,h2inv)
                   scaledt = dt/bd(1)
                   scaledi = 1./scaledt
                   call cmult(dp,scaledt,ntot2)        ! scale for tol

                   call esolver  (dp,h1,h2,h2inv,intype)
                   call cmult(dp,scaledi,ntot2)
      if (ifprjp)  call fsi_gensolnp(dp,h1,h2,h2inv)

      call add2(up,dp,ntot2)

      call opgradt  (w1 ,w2 ,w3 ,dp)
      call opbinv   (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
      dtb  = dt/bd(1)
      call opadd2cm (ux ,uy ,uz ,dv1,dv2,dv3, dtb )

      tpres=tpres+(dnekclock()-etime1)

      return
      end subroutine fsi_incomprn
c-----------------------------------------------------------------------

      subroutine fsi_bcdirvc(v1,v2,v3,mask1,mask2,mask3)
C
C     Apply Dirichlet boundary conditions to surface of vector (V1,V2,V3).
C     Use IFIELD as a guide to which boundary conditions are to be applied.

      implicit none

      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
      INCLUDE 'TOPOL'
      INCLUDE 'CTIMER'

      real tmp1,tmp2,tmp3
      common /scruz/ tmp1(lx1,ly1,lz1,lelv)
     $             , tmp2(lx1,ly1,lz1,lelv)
     $             , tmp3(lx1,ly1,lz1,lelv)

      real tmq1,tmq2,tmq3
      common /scrmg/ tmq1(lx1,ly1,lz1,lelv)
     $             , tmq2(lx1,ly1,lz1,lelv)
     $             , tmq3(lx1,ly1,lz1,lelv)

      real v1(nx1,ny1,nz1,lelv),v2(nx1,ny1,nz1,lelv)
     $    ,v3(nx1,ny1,nz1,lelv)
      real mask1(nx1,ny1,nz1,lelv),mask2(nx1,ny1,nz1,lelv)
     $    ,mask3(nx1,ny1,nz1,lelv)

      common  /nekcb/ cb
      character cb*3
      character*1 cb1(3)
      equivalence (cb1,cb)

      logical ifonbc

      real val                ! value at the faces
      integer nfaces,nxyz,nel,ntot,isweep,ie,iface

      real bc1,bc2,bc3

      etime1=dnekclock()

      nfaces=2*ndim
      nxyz  =nx1*ny1*nz1
      nel   =nelfld(ifield)
      ntot  =nxyz*nel

      call rzero(tmp1,ntot)
      call rzero(tmp2,ntot)
      if (if3d) call rzero(tmp3,ntot)

!     Velocity boundary conditions

      do 2100 isweep=1,2
         do 2000 ie=1,nel
         do 2000 iface=1,nfaces
            cb  = cbc(iface,ie,ifield)
            bc1 = bc(1,iface,ie,ifield)
            bc2 = bc(2,iface,ie,ifield)
            bc3 = bc(3,iface,ie,ifield)

            if (cb.eq.'ON ' .or. cb.eq.'on ') then   ! 5/21/01 pff
                ifonbc =.true.
                call fsi_faceiv ('v  ',tmp1(1,1,1,ie),tmp2(1,1,1,ie),
     $                       tmp3(1,1,1,ie),ie,iface,nx1,ny1,nz1)
            else 
                call fsi_faceiv (cb,tmp1(1,1,1,ie),tmp2(1,1,1,ie),
     $                       tmp3(1,1,1,ie),ie,iface,nx1,ny1,nz1)

                if ( ifqinp(iface,ie) )
     $          call globrot (tmp1(1,1,1,ie),tmp2(1,1,1,ie),
     $                        tmp3(1,1,1,ie),ie,iface)
            endif

 2000    continue

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
      if ( .not.ifstrs ) then
         call col2(v1,mask1,ntot)
         call col2(v2,mask2,ntot)
         if (if3d) call col2(v3,mask3,ntot)
         if (ifonbc) then
            call antimsk1(tmp1,mask1,ntot)
            call antimsk1(tmp2,mask2,ntot)
            if (if3d) call antimsk1(tmp3,mask3,ntot)
         endif
      else
         if (ifmodel) then
             call copy (tmq1,tmp1,ntot)
             call copy (tmq2,tmp2,ntot)
             if (ndim.eq.3) call copy (tmq3,tmp3,ntot)
             call amask (tmp1,tmp2,tmp3,tmq1,tmq2,tmq3,nelv)
         endif
         call rmask (v1,v2,v3,nelv)
      endif

      call add2(v1,tmp1,ntot)
      call add2(v2,tmp2,ntot)
      if (if3d) call add2(v3,tmp3,ntot)

      tusbc=tusbc+(dnekclock()-etime1)

      return
      end
c-----------------------------------------------------------------------

      subroutine fsi_faceiv (cb,v1,v2,v3,iel,iface,nx,ny,nz)
C
C     Assign fortran function boundary conditions to 
C     face IFACE of element IEL for vector (V1,V2,V3).
C
      INCLUDE 'SIZE'
      INCLUDE 'NEKUSE'
      INCLUDE 'PARALLEL'
C
      dimension v1(nx,ny,nz),v2(nx,ny,nz),v3(nx,ny,nz)
      character cb*3
c
      character*1 cb1(3)
c
      common  /nekcb/ cb3
      character*3 cb3
      cb3 = cb
c
      call chcopy(cb1,cb,3)
c
      ieg = lglel(iel)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
C
      if (cb.eq.'mv '.or.cb.eq.'mvn'.or.cb.eq.'v  ') then
c
         do 100 iz=kz1,kz2
         do 100 iy=ky1,ky2
         do 100 ix=kx1,kx2
            call nekasgn (ix,iy,iz,iel)
            call fsi_userbc  (ix,iy,iz,iface,ieg)
            v1(ix,iy,iz) = ux
            v2(ix,iy,iz) = uy
            v3(ix,iy,iz) = uz
  100    continue
         return

!      else 
!         do iz=kz1,kz2
!         do iy=ky1,ky2
!         do ix=kx1,kx2
!            v1(ix,iy,iz) = 0.
!            v2(ix,iy,iz) = 0.
!            v3(ix,iy,iz) = 0.
!         enddo
!         enddo
!         enddo
         return
      endif

C
      return
      end subroutine fsi_faceiv
c-----------------------------------------------------------------------

      subroutine fsi_userbc (ix,iy,iz,iside,ieg)

      implicit none

      include 'SIZE'
      include 'NEKUSE'
      include 'PARALLEL'
      include 'FSI'
      include 'TSTEP'
      include 'INPUT'
      include 'STRUCT'
      include 'SOLN'
      include 'NEKNEK'

      integer ix,iy,iz,ieg,iside,iel
      real dx,dy

      iel=gllel(ieg)

!     Correction of velocity is the boundary condition
!     Assuming valint has been filled up
      if (imask(ix,iy,iz,iel).eq.1) then 
        ux = valint(ix,iy,iz,iel,1)
        uy = valint(ix,iy,iz,iel,2)
        if (if3d) uz = valint(ix,iy,iz,iel,3)
      else
        ux = 0.
        uy = 0.
        uz = 0.
      endif        

!     apply boundary conditions obtained from structural equations      

      return
      end

!---------------------------------------------------------------------- 

      subroutine fsi_hmhzsf (name,u1,u2,u3,r1,r2,r3,h1,h2,
     $                   rmask1,rmask2,rmask3,rmult,
     $                   tol,maxit,matmod)
C-----------------------------------------------------------------------
C
C     Compute solution to coupled Helmholtz equations 
C     (stress formulation)
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
!      include 'SOLN'   ! For outpost diagnostic call
      include 'TSTEP'
      include 'ORTHOSTRS'

      DIMENSION U1(LX1,LY1,LZ1,LELV)
     $        , U2(LX1,LY1,LZ1,LELV)
     $        , U3(LX1,LY1,LZ1,LELV)
     $        , R1(LX1,LY1,LZ1,LELV)
     $        , R2(LX1,LY1,LZ1,LELV)
     $        , R3(LX1,LY1,LZ1,LELV)
     $        , H1(LX1,LY1,LZ1,LELV)
     $        , H2(LX1,LY1,LZ1,LELV)
     $        , RMASK1(LX1,LY1,LZ1,LELV)
     $        , RMASK2(LX1,LY1,LZ1,LELV)
     $        , RMASK3(LX1,LY1,LZ1,LELV)
     $        , RMULT (LX1,LY1,LZ1,LELV)
      CHARACTER NAME*4

      common /cpfjunk/ y(lx1*ly1*lz1*lelt,3)
      common /cpfjun2/ v(lx1*ly1*lz1*lelt,3)

      nel = nelfld(ifield)
      vol = volfld(ifield)
      n   = nx1*ny1*nz1*nel

      call rmask   (r1,r2,r3,nel)
      call opdssum (r1,r2,r3)
      call rzero3  (u1,u2,u3,n)

c     call set_up_h1_crs_strs(h1,h2,ifield,matmod)
      
      if (imesh.eq.1) then
         call chktcgs (r1,r2,r3,rmask1,rmask2,rmask3,rmult,binvm1
     $                ,vol,tol,nel)

!         if (matmod.lt.0) then
!          napprox(1) = 0
!          iproj      = 0 !param(94)
!          if (iproj.gt.0.and.istep.gt.iproj) napprox(1)=param(93)
!          napprox(1)=min(napprox(1),istep/3)
!          call strs_project_a(r1,r2,r3,h1,h2,rmult,ifield,ierr,matmod)

c         call opcopy(y(1,1),y(1,2),y(1,3),x(1),x(1+n),x(1+2*n))

!         endif

         call cggosf  (u1,u2,u3,r1,r2,r3,h1,h2,rmult,binvm1
     $                ,vol,tol,maxit,matmod)


!         if (matmod.lt.0)
!     $    call strs_project_b(u1,u2,u3,h1,h2,rmult,ifield,ierr)

      else
         call chktcgs (r1,r2,r3,rmask1,rmask2,rmask3,rmult,bintm1
     $                ,vol,tol,nel)
         call cggosf  (u1,u2,u3,r1,r2,r3,h1,h2,rmult,bintm1
     $                ,vol,tol,maxit,matmod)
      endif

      return
      end
!---------------------------------------------------------------------- 
      subroutine fsi_setrhsp(p,h1,h2,h2inv)
C
C     Project soln onto best fit in the "E" norm.
C
      implicit none
      
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
      include 'FSI'

      real p    (lx2,ly2,lz2,lelv)
      real h1   (lx1,ly1,lz1,lelv)
      real h2   (lx1,ly1,lz1,lelv)
      real h2inv(lx1,ly1,lz1,lelv)
!      real pset (lx2,ly2,lz2,lelv,sc_mproj)

      integer ltot2
      parameter (ltot2=lx2*ly2*lz2*lelv)

      integer i,intetype,ntot2,ierr

      real vlsc2                ! function

!      real pbar,pnew,alpha,work
!      common /fsi_orthox/ pbar(ltot2),pnew(ltot2)
!      common /fsi_orthos/ alpha(sc_mproj),work(sc_mproj)

      if (Nprev.eq.0) return

      ntot2  = nx2*ny2*nz2*nelv

      ierr = 0
      call fsi_updrhse(p,h1,h2,h2inv,ierr) ! update rhs's if E-matrix has changed
      if (ierr.ne.0) then
        if (nio.eq.0) then
          write(6,*) 'Error in fsi_updrhse. Setting Nprev=0' 
        endif  
        Nprev=0           ! Doesn't happen w/ new formulation
      endif  

      do i=1,Nprev  ! Perform Gram-Schmidt for previous soln's.
         alpha(i) = vlsc2(p,sc_prproj(1,1,1,1,i),ntot2)
      enddo
      call gop(alpha,work,'+  ',Nprev)

      call rzero(pbar,ntot2)
      do i=1,Nprev
         call add2s2(pbar,sc_prproj(1,1,1,1,i),alpha(i),ntot2)
      enddo
C
      intetype = 1
      call cdabdtp(pnew,pbar,h1,h2,h2inv,intetype)
      call sub2   (p,pnew,ntot2)

      return
      end subroutine fsi_setrhsp
c-----------------------------------------------------------------------
      subroutine fsi_gensolnp(p,h1,h2,h2inv)
C
C     Reconstruct the solution to the original problem by adding back
C     the previous solutions
C
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FSI'

      real p    (lx2,ly2,lz2,lelv)
      real h1   (lx1,ly1,lz1,lelv)
      real h2   (lx1,ly1,lz1,lelv)
      real h2inv(lx1,ly1,lz1,lelv)
!      real pset (lx2*ly2*lz2*lelv,sc_mproj)

      integer ltot2
      parameter (ltot2=lx2*ly2*lz2*lelv)

!      real pbar,pnew,alpha,work
!      common /fsi_orthox/ pbar(ltot2),pnew(ltot2)
!      common /fsi_orthos/ alpha(sc_mproj),work(sc_mproj)

      integer ierr,ntot2,mprv

      integer isave,isave_freq
      save isave
      data isave /0/

      ierr = 0

      mprv=sc_mproj

      isave_freq=1
      isave=mod(isave,isave_freq)

      ntot2=nx2*ny2*nz2*nelv

      if (Nprev.lt.mprv) then
        if (isave.eq.0) then
!          if (nid.eq.0) write(6,*) 'Saving new AMP projection'
          Nprev = Nprev + 1
          call copy  (sc_prproj(1,1,1,1,Nprev),p,ntot2)          ! Save current solution
          call econjp(sc_prproj,Nprev,h1,h2,h2inv,ierr)          ! Orthonormalize set
        endif  
        call add2  (p,pbar,ntot2)                               ! Reconstruct solution.

        if (ierr.eq.1) then
          if (nid.eq.0) write(6,*) 'Error in AMP projection', ierr
          isave = 0
          Nprev = 1
          call copy  (sc_prproj(1,1,1,1,Nprev),p,ntot2)        ! Save current solution
          call econjp(sc_prproj,Nprev,h1,h2,h2inv,ierr)        ! and orthonormalize.
        endif
      else                                                      ! (uses pnew).
!        if (nid.eq.0) write(6,*) 'Resetting AMP projections'
        Nprev = 1
        call add2  (p,pbar,ntot2)                               ! Reconstruct solution.
        call copy  (sc_prproj(1,1,1,1,Nprev),p,ntot2)          ! Save current solution
        call econjp(sc_prproj,Nprev,h1,h2,h2inv,ierr)          ! and orthonormalize.
      endif
      isave = isave + 1

      return
      end subroutine fsi_gensolnp
c-----------------------------------------------------------------------
      subroutine fsi_updrhse(p,h1,h2,h2inv,ierr)
C
C     Update rhs's if E-matrix has changed
C
C
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'FSI'
C
      parameter (ltot2=lx2*ly2*lz2*lelv)
!      common /fsi_orthox/ pbar(ltot2),pnew(ltot2)
!      common /fsi_orthos/ alpha(mxprev), work(mxprev), alphan, dtlast
!      common /fsi_orthoi/ nprev,mprev
      common /orthol/ ifnewe
!      real alpha,work
      logical ifnewe

      real             p    (lx2,ly2,lz2,lelv)
      real             h1   (lx1,ly1,lz1,lelv)
      real             h2   (lx1,ly1,lz1,lelv)
      real             h2inv(lx1,ly1,lz1,lelv)

      integer icalld
      save    icalld
      data    icalld/0/

      ntot2=nx2*ny2*nz2*nelv
c
c     first, we have to decide if the e matrix has changed.
c
      if (icalld.eq.0) then
         icalld=1
         dtlast=dt
      endif
c
      ifnewe=.false.
      if (ifmvbd) then
         ifnewe=.true.
         call invers2(bm2inv,bm2,ntot2)
      elseif (dtlast.ne.dt) then
         ifnewe=.true.
         dtlast=dt
      endif
      if (ifnewe.and.nio.eq.0) write(6,*) 'reorthogo AMP:',nprev
     
C     Next, we reconstruct a new rhs set.
     
      if (ifnewe) then

         nprevt = nprev
         do 100 iprev=1,nprevt
c           orthogonalize this rhs w.r.t. previous rhs's
            call fsi_econj (iprev,h1,h2,h2inv,ierr)
            if (ierr.eq.1) then
               nprev = 0
               return
            endif
  100    continue
c
      endif
c
      return
      end subroutine fsi_updrhse
!---------------------------------------------------------------------- 

      subroutine fsi_econj(kprev,h1,h2,h2inv,ierr)
C
C     Orthogonalize the rhs wrt previous rhs's for which we already
C     know the soln.
C
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
      include 'FSI'
C
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)
C
      INTEGER LTOT2
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
!      COMMON /FSI_ORTHOX/ Pbar(ltot2),Pnew(ltot2),Pbrr(ltot2)
!      COMMON /FSI_ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
!      COMMON /FSI_ORTHOI/ Nprev,Mprev
!      REAL ALPHA,WORK
      real ALPHAd
      real alpham
      real GLSC2,VLSC2

      integer Kprev,Kprev1,ipass,npass,ntot2
      integer I,ierr
      integer intetype
C
C
      ierr  = 0
      ntot2 = nx2*ny2*nz2*nelv
      intetype=1
c
c     gram schmidt, w re-orthogonalization
c
      npass=1
      if (abs(param(105)).eq.2) npass=2
      do ipass=1,npass
c
         call cdabdtp(pbrr,sc_prproj(1,1,1,1,kprev),h1,h2,h2inv,
     $                  intetype)
c
c        compute part of the norm
         alphad = glsc2(sc_prproj(1,1,1,1,kprev),pbrr,ntot2)
c
c        gram-schmidt
         kprev1=kprev-1
         do 10 i=1,kprev1
            alpha(i) = vlsc2(pbrr,sc_prproj(1,1,1,1,i),ntot2)
   10    continue
         if (kprev1.gt.0) call gop(alpha,work,'+  ',kprev1)
c
         do 20 i=1,kprev1
            alpham = -alpha(i)
            call add2s2(sc_prproj(1,1,1,1,kprev),sc_prproj(1,1,1,1,i),
     $             alpham,ntot2)
            alphad = alphad - alpha(i)**2
   20    continue
      enddo
c
c    .normalize new element in p~
c
      if (alphad.le.0.0) then
         write(6,*) 'ERROR:  alphad .le. 0 in ECONJ',alphad,Kprev
         ierr = 1
         return
      endif
      alphad = 1.0/sqrt(alphad)
      alphan = alphad
      call cmult(sc_prproj (1,1,1,1,kprev),alphan,ntot2)
c
      return
      end subroutine fsi_econj
!----------------------------------------------------------------------

