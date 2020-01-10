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
      
      amp = 0.0

      ran = 3.e4*(ieg+X*sin(Y)+Z*cos(Y))
     $     + 4.7e2*ix*iy*iz - 1.5e3*ix*iy + .5e5*ix
      ran = 6.e3*sin(ran)
      ran = 3.e3*sin(ran)
      ran = cos(ran)

      if (fsi_ifstruct) then
        ux = 0.1 + ran*amp
      else
        ux = 1.0 + ran*amp
      endif  
     
      ran = (2+ran)*1.e4*(ieg+Y*sin(Z)+X*cos(Z))
     $     + 1.5e3*ix*iy*iz - 2.5e3*ix*iy + 8.9e4*ix
      ran = 2.e3*sin(ran)
      ran = 7.e3*sin(ran)
      ran = cos(ran)
      uy = ran*amp
      
      ran = (4+ran)*5.1e4*(ieg+Z*sin(X)+Y*cos(X))
     $     + 4.6e3*ix*iy*iz - 2.9e4*ix*iy + 3.7e3*ix
      ran = 9.e3*sin(ran)
      ran = 4.e3*sin(ran)
      ran = cos(ran)
      uz = ran*amp

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

      call outpost(vx,vy,vz,pr,t,'  ')

      call fsi_neknek_velex
      igeom = 2
      ifield = 1
      call bcdirvc(vx,vy,vz,v1mask,v2mask,v3mask)

      call outpost(vx,vy,vz,pr,t,'  ')

      nt = nx1*ny1*nz1*nelt
      do i=1,nt
        vy(i,1,1,1)=imask(i,1,1,1)+0.
      enddo  

      call outpost(v1mask,vy,v3mask,pr,t,'  ')

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
