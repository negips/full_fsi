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
      include 'FSI'             ! EXT_VX,...
      include 'STRUCT'          ! FSI_IFFLUID

      integer ix,iy,iz,iside,ieg,iel

      ux = 1.0
      uy = 0.
      uz = 0. 
 
      iel = gllel(ieg)
      if ((fsi_iffluid).and.(imask(ix,iy,iz,iel).eq.1)) then
        ux = ext_vx(ix,iy,iz,iel)
        uy = ext_vy(ix,iy,iz,iel)
        uz = ext_vz(ix,iy,iz,iel)
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
      if (fsi_ifstruct) then
        ux  = 0.0 + 0*cos(x)
      else
        ux = 0.0 + amp*cos(x)
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

      if (if3d) then
        nfld_neknek = 4   ! field to interpolate
                          ! 4: u,v,w,pr (in 3D)
      else
        nfld_neknek = 3   ! field to interpolate
                          ! 3: u,v,pr (in 2D)
      endif

      if (uparam(1).eq.1) then
        fsi_ifstruct = .false.
        fsi_iffluid  = .true.
      elseif (uparam(1).eq.2) then
        fsi_ifstruct = .true.
        fsi_iffluid  = .false.
      else
        fsi_ifstruct = .false.
        fsi_iffluid  = .false.
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

      common /c_is1/ glo_num(1*lx1*ly1*lz1*lelv)
      integer*8 glo_num
      real tmp(lx1*ly1*lz1*lelv)
      real tmp2(lx1*ly1*lz1*lelv)
      real tmp3(lx1*ly1*lz1*lelv)


      if (istep.eq.0) then
        nt = nx1*ny1*nz1*nelv
        do i=1,nt
          tmp(i)=glo_num(i)+0.
          tmp2(i)=v1mask(i,1,1,1)+0.
          tmp3(i)=v2mask(i,1,1,1)+0.
        enddo
        
        call outpost(tmp2,tmp3,tmp,pr,t,'glo') 
      endif        

      fsi_iftran = .true. 

      call fsi_coupling

!      call outpost(ts1,ts2,ts3,pr,t,'ts1')
!      call outpost(ts4,ts5,ts6,pr,t,'ts4')

      if (fsi_ifstruct.and.(mod(istep,iostep).eq.0)
     $      .and.(istep.gt.0)) then
        call opcopy(ts1,ts2,ts3,xm1,ym1,zm1)
        call opadd2(xm1,ym1,zm1,vx,vy,vz)

        call outpost(velx,vely,velz,pr,t,'vel')
        call opcopy(xm1,ym1,zm1,ts1,ts2,ts3)
      endif        


      ifto = .true.
      nt = nx1*ny1*nz1*nelv
      call copy(t,cflf,nt)

      
      return
      end
c -----------------------------------------------------------------------


c automatically added by makenek
      subroutine usrdat0() 

      return
      end

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