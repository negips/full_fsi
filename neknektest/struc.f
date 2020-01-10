C-----------------------------------------------------------------------
c
c     user subroutines required by nek5000
c
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'NEKUSE'          ! UDIFF, UTRANS

      UDIFF =0.
      UTRANS=0.

      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'NEKUSE'          ! FF[XYZ]

      FFX = 0.0
      FFY = 0.0
      FFZ = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'NEKUSE'          ! QVOL

      QVOL   = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      implicit none
      include 'SIZE'
      include 'NEKUSE'          ! UX, UY, UZ, TEMP, X, Y, PA
      include 'NEKNEK'
      include 'PARALLEL'
      integer ix,iy,iz,iside,ieg,iel

      ux = 1.0
      uy = 0.
      uz = 0. 
 
      iel = gllel(ieg)
      if (imask(ix,iy,iz,iel).eq.1) then
        ux = valint(ix,iy,iz,iel,1)
        uy = valint(ix,iy,iz,iel,2)
        uz = valint(ix,iy,iz,iel,3)
        if (nfld_neknek.gt.3) temp = valint(ix,iy,iz,iel,ldim+2)
      end if
      
      return
      end
c -----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)

      implicit none

      include 'SIZE'
      include 'NEKUSE'          ! UX, UY, UZ, TEMP, X, Y, Z
      include 'STRUCT'
      integer ix,iy,iz,ieg

      real amp, ran
      

      amp = 1.0
      if (fsi_ifstruct) then
        ux  = 0.1 + 0*cos(x)
      else
        ux = 1.0 + amp*cos(x)
      endif  
     
      uy = amp*sin(y)
      
      uz = 0.

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat
      implicit none
      include 'SIZE'
      include 'INPUT'         ! ngeom
      include 'NEKNEK'        ! nfld_neknek
      include 'TSTEP'         ! ninter
      include 'STRUCT'

      ngeom = 2   ! >2 => internal iterations

      ninter = 1  ! order of interface extrapolation

      nfld_neknek = 3   ! field to interpolate
                        ! 3: u,v,pr (in 2D)


      if (uparam(1).eq.1) then
        fsi_ifstruct = .true.
        fsi_iffluid  = .false.
      else
        fsi_ifstruct = .false.
        fsi_iffluid  = .true.
      endif            

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      implicit none
      include 'SIZE'
      include 'SOLN'            ! vx,vy,vz,pr,t
      include 'GEOM'            ! boundaryID
      include 'INPUT'           ! cbc
      
      
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      implicit none
      include 'SIZE'
      include 'INPUT'           ! param, if3d
      include 'MASS'            ! volvm1      
      include 'NEKNEK'

      

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'
      include 'STRUCT'
      include 'MASS'
      include 'NEKNEK'        ! igeom
      include 'GEOM'

      integer ltmp
      parameter (ltmp=100)
      integer*8 newnum(ltmp)
      integer gsh_tmp
      real tst(ltmp)

      integer len2
      integer i,nt

      real scale

      call outpost(vx,vy,vz,pr,t,'  ')

      call fsi_neknek_velex
      igeom = 2
      ifield = 1
      call bcdirvc(vx,vy,vz,v1mask,v2mask,v3mask)

      call outpost(vx,vy,vz,pr,t,'  ')

      nt = nx1*ny1*nz1*nelt
!      do i=1,nt
!        vy(i,1,1,1)=imask(i,1,1,1)+0.
!      enddo

      call opzero(bfx,bfy,bfz)
      scale = 1.
      if (fsi_iffluid) then
        call fluid_forces(struct_ssx,struct_ssy,struct_ssz,
     $                    vx,vy,vz,pr,scale)
      else
        call opzero(struct_ssx,struct_ssy,struct_ssz)
      endif

      call outpost(struct_ssx,struct_ssy,struct_ssz,pr,t,'  ')

      call fsi_neknek_stressex

      call outpost(struct_ssx,struct_ssy,struct_ssz,pr,t,'  ')
     
      call exitt

      if (istep.gt.0) then

        fsi_iftran = .true. 

        call fsi_coupling

        call exitt

      endif        
      
      return
      end
c -----------------------------------------------------------------------


c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
