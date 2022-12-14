c-----------------------------------------------------------------------------
!      subroutine userchk
!      include 'SIZE'
!      include 'TOTAL'
!
!      call steady_elast_solve
!
!      call error_check
!
!      ifxyo = .true.
!      call outpost(vx,vy,vz,pr,t,'   ')
!
!      if (nio.eq.0) write(*,*) 'Done elast solve; nelv',nelv
!      call exitt0
!
!      return
!      end
c-----------------------------------------------------------------------
!      subroutine usrdat2
!      include 'SIZE'
!      include 'TOTAL'
!      common /material/ young,nu
!      common /elastr/ lambda,g,gamma_n,beta_n,rho,damp
!      logical ifPlaneStress
!
!      real nu, lambda
!
!
!      ifPlaneStress=.false.
!
!      young=100.
!      nu=0.3
!
!      if (ifPlaneStress) then
!      lambda = young*nu/( (1+nu)*(1-nu) )   ! for plane stress
!      else
!      lambda = young*nu/( (1+nu)*(1-2*nu) ) ! for 3D, and 2D plane strain
!      end if
!
!      g      = .5*young/(1+nu)
!
!      return
!      end
c-----------------------------------------------------------------------
      subroutine steady_elast_solve
      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (l3=lx1*ly1*lz1*lelt*ldim)
      common /scruz/ w1(lt),w2(lt),w3(lt)
      common /scrns/ x(l3),f(l3),r(l3)
      common /scrcg/ w(l3),p(l3),z(l3)
      common /scrsf/ wt(l3)
      real loadx,loady,loadz

      real SX(LX1,LY1,LZ1,LELT), SY(LX1,LY1,LZ1,LELT), 
     &     SZ(LX1,LY1,LZ1,LELT)


      common /elasti/ ielast   ! Eqn type
      common /elastr/ lambda,g,gamma_n,beta_n,rho,damp

      real lambda

      logical ifmsk,ifdss

      ifield = 1
      ielast = 1   ! Eqn type

      call bcdirvc (vx,vy,vz,v1mask,v2mask,v3mask)


      ifmsk = .false.
      ifdss = .false.
      call elast   (w1,w2,w3,vx,vy,vz,lambda,g,ifmsk,ifdss) ! w/ dssum & mask

      n = nx1*ny1*nz1*nelv
      i0 = 0
      i1 = n
      i2 = n*2

      loadx=0.
      loady=0.
      loadz=0.

      do i=1,n

        f (i0+i) =  bm1(i,1,1,1)*loadx-w1(i)
        f (i1+i) =  bm1(i,1,1,1)*loady-w2(i)
        f (i2+i) =  bm1(i,1,1,1)*loadz-w3(i)
        wt(i0+i) =  vmult(i,1,1,1)
        wt(i1+i) =  vmult(i,1,1,1)
        wt(i2+i) =  vmult(i,1,1,1)
      end do

c        Add stuff to the rhs due to the stress bc ('S')

       call GETSTRESS(SX,SY,SZ)
       
       do i=1,n
       f(i0+i)=(f(i0+i)+SX(i,1,1,1))*v1mask(i,1,1,1)
       f(i1+i)=(f(i1+i)+SY(i,1,1,1))*v2mask(i,1,1,1)
       f(i2+i)=(f(i2+i)+SZ(i,1,1,1))*v3mask(i,1,1,1)
       end do

      call opdssum(f(i0+1),f(i1+1),f(i2+1))

      n3    = 3*n
      niter = 3*n

      call cg(x,f,r,w,p,z,wt,n3,niter,nid)

      do i=1,n
         vx(i,1,1,1) = x(i0+i) + vx(i,1,1,1)
         vy(i,1,1,1) = x(i1+i) + vy(i,1,1,1)
         vz(i,1,1,1) = x(i2+i) + vz(i,1,1,1)
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine elast3d_e(w1,w2,w3,u1,u2,u3,e)

c     Apply elasticity operator to u:   w = Eu

      real w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      integer e

      include 'SIZE'
      include 'GEOM'    ! jacmi,rxm1, etc.
      include 'INPUT'   ! if3d
      include 'MASS'    ! bm1
      include 'SOLN'    ! vtrans
      include 'TSTEP'   ! dt
      include 'WZ'      ! w3m1

      common /elasti/ ielast   ! Eqn type
      common /elastr/ lambda,g,gamma_n,beta_n,rho,damp
      real lambda

      parameter (lxyz=lx1*ly1*lz1)
      common /ctmp1/ ur(lxyz),us(lxyz),ut(lxyz)
     $             , vr(lxyz),vs(lxyz),vt(lxyz)
     $             , wr(lxyz),ws(lxyz),wt(lxyz)

      call gradl_rst(ur,us,ut,u1,nx1,if3d) ! Grad on GLL
      call gradl_rst(vr,vs,vt,u2,nx1,if3d)
      call gradl_rst(wr,ws,wt,u3,nx1,if3d)

      a =          g
      b =      lambda 
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

      if (ielast.eq.2) then
         sd = beta_n*dt*dt
         do i=1,lxyz
            w1(i) = sd*w1(i) + bm1(i,1,1,e)*vtrans(i,1,1,e,1)*u1(i)
            w2(i) = sd*w2(i) + bm1(i,1,1,e)*vtrans(i,1,1,e,1)*u2(i)
            w3(i) = sd*w3(i) + bm1(i,1,1,e)*vtrans(i,1,1,e,1)*u3(i)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine elast2d_e(w1,w2,w3,u1,u2,u3,e)

c     Apply elasticity operator to u:   w = Eu

      real w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      integer e

      include 'SIZE'
      include 'GEOM'    ! jacmi,rxm1, etc.
      include 'INPUT'   ! if3d
      include 'MASS'    ! bm1
      include 'SOLN'    ! vtrans
      include 'TSTEP'   ! dt
      include 'WZ'      ! w3m1

      common /elasti/ ielast   ! Eqn type
      common /elastr/ lambda,g,gamma_n,beta_n,rho,damp
      real lambda

      parameter (lxyz=lx1*ly1*lz1)
      common /ctmp1/ ur(lxyz),us(lxyz),ut(lxyz)
     $             , vr(lxyz),vs(lxyz),vt(lxyz)
     $             , wr(lxyz),ws(lxyz),wt(lxyz)

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

      if (ielast.eq.2) then
         sd = beta_n*dt*dt
         do i=1,lxyz
            w1(i) = sd*w1(i) + bm1(i,1,1,e)*vtrans(i,1,1,e,1)*u1(i)
            w2(i) = sd*w2(i) + bm1(i,1,1,e)*vtrans(i,1,1,e,1)*u2(i)
         enddo
      endif

      return
      end
C---------------------------------------------------------------------------
      subroutine elast(w1,w2,w3,u1,u2,u3,lambda,g,ifmsk,ifdss)
      include 'SIZE'
      include 'TOTAL'

      real w1(1),w2(1),w3(1),u1(1),u2(1),u3(1),lambda,g
      logical ifmsk,ifdss

      integer e

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
      end
c-----------------------------------------------------------------------
      subroutine solveM(z,r,n)

      include 'SIZE'
      include 'TOTAL'

      real z(nx1*ny1*nz1*nelt,3),r(nx1*ny1*nz1*nelt,3)

      m = nx1*ny1*nz1*nelt
      do i=1,m
         z(i,1) = binvm1(i,1,1,1)*r(i,1)*v1mask(i,1,1,1)
         z(i,2) = binvm1(i,1,1,1)*r(i,2)*v2mask(i,1,1,1)
         z(i,3) = binvm1(i,1,1,1)*r(i,3)*v3mask(i,1,1,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cg(x,f,r,w,p,z,wt,n,niter,nid)

c     Solve Ax=f where A is SPD and is invoked by the routine ax()
c
c     Output:  x - vector of length n
c
c     Input:   f - vector of length n
c
c     Work arrays:   r,w,p,z  - vectors of length n
c
c     User-provided routine ax(w,z,n) returns  w := Az,  
c
c     User-provided routine solveM(z,r,n) ) returns  z := M^-1 r,  
c
      real x(n),f(n),r(n),w(n),p(n),z(n),wt(n)

      pap = 0.0
c
c     set machine tolerances
c
      one = 1.
      eps = 1.e-20
      if (one+eps .eq. one) eps = 1.e-14
      if (one+eps .eq. one) eps = 1.e-7
      tol=eps

      rtz1=1.0

      call rzero(x,n)
      call copy (r,f,n)
      rnorm = sqrt(glsc3(r,wt,r,n))
      iter = 0
      if (nid.eq.0) write(6,6) iter,rnorm

      miter = niter
      do 1000 iter=1,miter
         call solveM(z,r,n)  !  Invert preconditioner here

         rtz2=rtz1
         rtz1=glsc3(r,wt,z,n)   ! parallel inner product

         beta = rtz1/rtz2
         if (iter.eq.1) beta=0.0

         call add2s1(p,z,beta,n)

         call axcg(w,p,n)

         pap=glsc3(w,wt,p,n)
         alpha=rtz1/pap
         alphm=-alpha
         call add2s2(x,p,alpha,n)
         call add2s2(r,w,alphm,n)

         rtr = glsc3(r,wt,r,n)
         rnorm = sqrt(rtr)
         if (iter.eq.1) rlim  = rnorm*tol
         if (iter.eq.1) rnm0  = rnorm
         if (rnm0.gt.0) rrel  = rnorm/rnm0
         if (nid.eq.0.and.(mod(iter,50).eq.0.or.iter.lt.10)) 
     $      write(6,6) iter,rnorm,rlim,rnm0,rrel
    6    format('cg:',i6,1p4e12.4)

         if (rnorm.le.rlim.and.iter.gt.1) then
         if (nid.eq.0) 
     $      write(6,6) iter,rnorm,rlim,rnm0,rrel
	 goto 1001
	 end if

 1000 continue
 1001 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine axcg(w,z,n)

      include 'SIZE'
      include 'TOTAL'

      real w(nx1*ny1*nz1*nelt,3),z(nx1*ny1*nz1*nelt,3)

      common /elastr/ lambda,g,gamma_n,beta_n,rho,damp
      real lambda,g

      logical ifmsk,ifdss
      ifmsk = .true.
      ifdss = .true.


c     with dssum & mask
      call elast
     $  (w(1,1),w(1,2),w(1,3),z(1,1),z(1,2),z(1,3),lambda,g,ifmsk,ifdss)

      return
      end

c-----------------------------------------------------------------------
      subroutine GETSTRESS(SX,SY,SZ)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      character*3 cb

      real SX(LX1,LY1,LZ1,LELT), SY(LX1,LY1,LZ1,LELT), 
     &     SZ(LX1,LY1,LZ1,LELT)


      n = nx1*ny1*nz1*nelv
           nface = 2*ldim


	 call opzero(SX,SY,SZ)

	 do 1000 iel=1,nelv

             ieg=lglel(iel)

	 do 1000 iface=1,nface
         cb  = cbc (iface,iel,ifield)
    	 
	 if (cb.eq.'S') then 

               IA=0
  	 CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFACE)  

               DO 100 IZ=KZ1,KZ2
               DO 100 IY=KY1,KY2
               DO 100 IX=KX1,KX2
                  IA = IA + 1

            CALL NEKASGN (IX,IY,IZ,IEL)
            CALL USERBC  (IX,IY,IZ,IFACE,IEG)

c      Rotate stress components from local to global system    

            TRX         = TRN*UNX(IA,1,IFACE,IEL) +
     $                    TR1*T1X(IA,1,IFACE,IEL) +
     $                    TR2*T2X(IA,1,IFACE,IEL)
            TRY         = TRN*UNY(IA,1,IFACE,IEL) +
     $                    TR1*T1Y(IA,1,IFACE,IEL) +
     $                    TR2*T2Y(IA,1,IFACE,IEL)
            TRZ         = TRN*UNZ(IA,1,IFACE,IEL) +
     $                    TR1*T1Z(IA,1,IFACE,IEL) +
     $                    TR2*T2Z(IA,1,IFACE,IEL)


        SX(IX,IY,IZ,IEL) = SX(IX,IY,IZ,IEL) + TRX*AREA(IA,1,IFACE,IEL)
        SY(IX,IY,IZ,IEL) = SY(IX,IY,IZ,IEL) + TRY*AREA(IA,1,IFACE,IEL)
        SZ(IX,IY,IZ,IEL) = SZ(IX,IY,IZ,IEL) + TRZ*AREA(IA,1,IFACE,IEL)
100	continue
	end if

1000    continue
	

	return
	end
c------------------------------------------------------------------------
      subroutine error_check
      include 'SIZE'
      include 'TOTAL'
      parameter (lt=lx1*ly1*lz1*lelt)
      common /exacu/ ue(lt),ve(lt)
      common /exacd/ ud(lt),vd(lt)
      common /exacw/ we(lt),wd(lt)

c     Error check

      n = nx1*ny1*nz1*nelv
      
      call exact  (ue,ve,we,xm1,ym1,zm1,n)

      call sub3   (ud,ue,vx,n)
      call sub3   (vd,ve,vy,n)
      call sub3   (wd,we,vz,n)
      call energy_norm(ud,vd,wd,n,enorm_dif)
      call energy_norm(ue,ve,we,n,enorm_abs)

      enorm=enorm_dif/enorm_abs

      v_dif_l2=sqrt((glsc3(ud,bm1,ud,n)+glsc3(vd,bm1,vd,n)+
     & glsc3(wd,bm1,wd,n))/3/volvm1)
      v_abs_l2=sqrt((glsc3(ue,bm1,ue,n)+glsc3(ve,bm1,ve,n)+
     & glsc3(we,bm1,we,n))/3/volvm1)

      v_l2=v_dif_l2/v_abs_l2

        if (nid.eq.0) then
          write(6,11) istep,time,v_l2,enorm,'error'
        end if

11    format(i5,1p3e14.6,a7)

      return
      end
c----------------------------------------------------------------------------
      subroutine energy_norm(u1,u2,u3,n,enorm)
      include 'SIZE'
      include 'TOTAL' 
      
      real sij (lx1,ly1,lz1,3*ldim-3,lelv)
      real fun (lx1,ly1,lz1,lelv)
      real fun1 (lx1,ly1,lz1,lelv),fun2 (lx1,ly1,lz1,lelv)
      parameter (lr=lx1*ly1*lz1)
      common /scruz/         ur(lr),us(lr),ut(lr)
     $                     , vr(lr),vs(lr),vt(lr)
     $                     , wr(lr),ws(lr),wt(lr)
      common /elastr/ lambda,g,gamma_n,beta_n,rho,damp
      real lambda,g
      real u1(lx1,ly1,lz1,lelt), u2(lx1,ly1,lz1,lelt)
      real u3(lx1,ly1,lz1,lelt)
      integer e

      nij = 3
      if (if3d.or.ifaxis) nij=6
      call comp_sij(sij,nij,u1,u2,u3,ur,us,ut,vr,vs,vt,wr,ws,wt)


      do 2010 e=1,nelv
               DO 100 IZ=1,nz1
               DO 100 IY=1,ny1
               DO 100 IX=1,nx1

         IF (IF3D.OR.IFAXIS) THEN
c
         s11 = sij(ix,iy,iz,1,e)
         s21 = sij(ix,iy,iz,4,e)
         s31 = sij(ix,iy,iz,6,e)
c
         s12 = sij(ix,iy,iz,4,e)
         s22 = sij(ix,iy,iz,2,e)
         s32 = sij(ix,iy,iz,5,e)
c
         s13 = sij(ix,iy,iz,6,e)
         s23 = sij(ix,iy,iz,5,e)
         s33 = sij(ix,iy,iz,3,e)

         div2=0.25*(s11+s22+s33)*(s11+s22+s33)

         strain=0.25*(s11*s11+s22*s22+s33*s33+
     &   2*s12*s12+2*s13*s13+2*s23*s23) 

	 else

            s11 = sij(ix,iy,iz,1,e)
            s12 = sij(ix,iy,iz,3,e)
            s21 = sij(ix,iy,iz,3,e)
            s22 = sij(ix,iy,iz,2,e)

            div2=0.25*(s11+s22)*(s11+s22)

            strain=0.25*(s11*s11+s22*s22+2*s12*s12) 

	 end if    

	    fun1(ix,iy,iz,e)=strain
	    fun2(ix,iy,iz,e)=div2		    

100 	    continue
2010	    continue

	    enorm1=sqrt(glsc2(fun1,bm1,n)/volvm1)
	    enorm2=sqrt(glsc2(fun2,bm1,n)/volvm1)	    	    

	    enorm=2*g*enorm1+lambda*enorm2

12     format(1p3e14.6,a7)

      return
      end

c------------------------------------------------------------------------
      subroutine exact (uu,vv,ww,xx,yy,zz,n)
      common /material/ young,nu
      real nu
      real uu(n),vv(n),ww(n),xx(n),yy(n),zz(n)

      p=100.

      r1=0.5
      r2=1.
   
      do i=1,n

      rr=sqrt(xx(i)*xx(i)+yy(i)*yy(i)+zz(i)*zz(i))

      sigma_r=-p/((r2/r1)**3-1)*(r2**3/rr**3-1)
      sigma_theta=p/((r2/r1)**3-1)*(0.5*r2**3/rr**3+1)
      eps_theta=((1.-nu)*sigma_theta-nu*sigma_r)/young

      ur=rr*eps_theta
     
      uu(i)=ur*xx(i)/rr
      vv(i)=ur*yy(i)/rr
      ww(i)=ur*zz(i)/rr

      end do


      return
      end
!---------------------------------------------------------------------- 




