!----------------------------------------------------------------------
!
!     NekNek routines for FSI
!     Author: Prabal Negi
!     Comments: Mostly slightly modified routines from multimesh.f      
! 
!---------------------------------------------------------------------- 
!====================================================================== 
      subroutine set_intflag 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      character*3 cb
      character*2 cb2
      equivalence (cb2,cb)
      integer j,e,f

c     Set interpolation flag: points with boundary condition = 'int' 
c     get intflag=1. 
c
c     Boundary conditions are changed back to 'v' or 't'.

      nfaces = 2*ldim
      
      nflag=nelt*nfaces
      call izero(intflag,nflag)

      do j=1,nfield
         nel = nelfld(j)
      do e=1,nel
      do f=1,nfaces
         cb=cbc(f,e,j)
         if (cb2.eq.'in') then
            intflag(f,e)=1
            if (j.ge.2) cbc(f,e,j)='t  '
            if (j.eq.1) cbc(f,e,j)='mv '
c            if (cb.eq.'inp') cbc(f,e,j)='on ' ! Pressure
            if (cb.eq.'inp') cbc(f,e,j)='o  ' ! Pressure
         endif
      enddo
      enddo
      enddo

c     zero out valint
      do i=1,nfld_neknek
        call rzero(valint(1,1,1,1,i),lx1*ly1*lz1*nelt)
      enddo

      return
      end
c------------------------------------------------------------------------

