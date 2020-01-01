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
       integer ix,iy,iz,iside,ieg

       ux = 1.0
       uy = 0.
       uz = 0. 
       
       return
       end
c -----------------------------------------------------------------------
       subroutine useric (ix,iy,iz,ieg)
       implicit none
      include 'SIZE'
      include 'NEKUSE'          ! UX, UY, UZ, TEMP, X, Y, Z
      integer ix,iy,iz,ieg

      real amp, ran
      
      amp = 0.0

      ran = 3.e4*(ieg+X*sin(Y)+Z*cos(Y))
     $     + 4.7e2*ix*iy*iz - 1.5e3*ix*iy + .5e5*ix
      ran = 6.e3*sin(ran)
      ran = 3.e3*sin(ran)
      ran = cos(ran)
      ux = 1. + ran*amp
      
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

!      ngeom = 20   ! >2 => internal iterations
!
!      ninter = 2  ! order of interface extrapolation
!
!      nfld_neknek = 3   ! field to interpolate
                        ! 3: u,v,pr (in 2D)


      fsi_ifstruct = .true.
      fsi_iffluid  = .false. 

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

      integer ltmp
      parameter (ltmp=100)
      integer*8 newnum(ltmp)
      integer gsh_tmp
      real tst(ltmp)

      integer len2

!      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer i

      call opzero(ts1,ts2,ts3)
      call opzero(ts4,ts5,ts6)
!      call oprone(dispx,dispy,dispz)

      if (istep.eq.1) then
        call plan_s

        call outpost(ts1,ts2,ts3,pr,t,'dbg')
        call outpost(ts4,ts5,ts6,pr,t,'dbg')

        do i=1,struct_nkryl
          call outpost(struct_krylv(1,1,i),struct_krylv(1,2,i),
     $            struct_krylv(1,3,i),pr,t,'slv')
          call outpost(struct_krylx(1,1,i),struct_krylx(1,2,i),
     $            struct_krylx(1,3,i),pr,t,'slx')
        enddo

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
