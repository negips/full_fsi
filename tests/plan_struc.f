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
      include 'CTIMER'
      include 'STRUCT'

      integer ntot
      logical ifconverged

      ifconverged = .false.

      if (fsi_fluid) then

!       call check_fsi_convergence()         
        do while (.not.ifconverged)            
!         Do Stokes correction step            
          call stokes_solve()

!         call check_fsi_convergence(ifconverged)

        enddo
      else            
!       First Structural solve
        do igeom=1:ngeom

          if (ifgeom) then
             if (.not.ifrich) call gengeom (igeom)
             call geneig  (igeom)
          endif
        
      call struct(igeom) 


      return            
      end subroutine fsi_advance
!---------------------------------------------------------------------- 
      subroutine struc_solve(igeom)

      implicit none

      include 'SIZE'
      include 'STRUCT'




      return
      end subroutine struc_solve             
