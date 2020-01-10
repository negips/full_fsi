!---------------------------------------------------------------------- 
!     Debugging routines here
!     Not to clutter the main file
!
!---------------------------------------------------------------------- 

      subroutine check_ortho(ikryl,wts)

      implicit none

      include 'SIZE'
      include 'STRUCT'
      include 'MASS'

      integer ikryl
      integer i,j
      real wts(lx1,ly1,lz1,lelv)

      real vtv(struct_nkryl+1,struct_nkryl+1)

      real op_glsc2_wt
      real beta


      character*32 outfmt

      call blank(outfmt,32)
      write(outfmt,'(A1,I2,A13)') '(',ikryl,'(E12.4E2,2x))'

      open(unit=109,file='vtv.out',status='unknown',form='formatted')

      call opcopy(struct_krylx(1,1,1),struct_krylx(1,2,1),
     $            struct_krylx(1,3,1),struct_krylv(1,1,1),
     $            struct_krylv(1,2,1),struct_krylv(1,3,1))

!      call opcolv(struct_krylx(1,1,1),struct_krylx(1,2,1),
!     $            struct_krylx(1,3,1),bm1)

      do i=1,ikryl
        do j=1,ikryl
            
          beta = op_glsc2_wt(struct_krylv(1,1,i),
     $              struct_krylv(1,2,i),struct_krylv(1,3,i),
     $              struct_krylv(1,1,j),struct_krylv(1,2,j),
     $              struct_krylv(1,3,j),wts)
          vtv(i,j) = beta
        enddo
        write(109,outfmt) (vtv(i,j), j=1,ikryl)
      enddo

      close(109) 

      return
      end subroutine check_ortho
!---------------------------------------------------------------------- 


