module boundary_bd

use mesh
use gauss
use gauss_bd
use utilities
use boundary_condition
use mesh_curved

implicit none

 double precision, allocatable :: cosng(:,:),sinng(:,:)
 double precision, allocatable :: normc(:,:,:),&
                           nxx(:,:), nyy(:,:), &
		dnxx(:,:), dnxy(:,:), dnyx(:,:), dnyy(:,:)
 double precision omega,ramp_period
 double precision, allocatable :: jac(:,:), &
       amm_bd(:,:,:), amm_bdo(:,:,:), &
       amm_bdl(:,:,:), agx_bd(:,:,:), agy_bd(:,:,:)
  
  
 contains

!
!***********************************************************************
! open boundary condition
!
! the elevation is prescribed at the boundary
!***********************************************************************
!
subroutine openbc1_read

 implicit none

    integer i,j,k

    if (nobc1.gt.0) then

      allocate(cosng(ngf_bd,nobc1),sinng(ngf_bd,nobc1))

      open(1,file='obc1.dat',status='old')
      read(1,*) omega,ramp_period
      do i=1,nobc1
        do j=1,ngf_bd
	 read(1,*) cosng(j,i),sinng(j,i)
	enddo
      enddo
      close(1)

    endif

    return

end subroutine openbc1_read

!
!***********************************************************************
! open boundary condition
!
! the elevation is prescribed at the boundary
!***********************************************************************
!
subroutine openbc1(fbd,time,fact)

 implicit none

      integer i,j,k,face
      double precision fbd(ngf_bd,*),time,phase,cosa,sina,fact,ans
      double precision pi,per,ramp

      if (nobc1.gt.0) then
	 if (ramp_period.gt.1.d-12) then
	  ramp = min (time,ramp_period)
	  ramp = ramp / ramp_period
!          ramp = tanh(time/ramp_period)
	 else
	  ramp=1.d0
	 endif
	 phase = time * omega
	 cosa = cos(phase) * fact * ramp
	 sina = sin(phase) * fact * ramp
         do i=1,nobc1
            face = pobc1(i)
 	    do j=1,ngf_bd
	       fbd(j,face) = cosng(j,i) * cosa + sinng(j,i) * sina
	    enddo
         enddo
      endif

      return
end subroutine openbc1

!
!***********************************************************************
! open boundary condition
! condition for mass conservation
!
!***********************************************************************
!
subroutine massbc(fbdx,fbdy)

 implicit none

      integer i,j,k,face,flag
      double precision fbdx(ngf_bd,*),fbdy(ngf_bd,*)

      do face=1,nf
       flag = fflag(face)
       if (flag.eq.0) then
        do j=1,ngf
         fbdx(j,face) = fbdx(j,face) + fbdy(j,face)
        enddo
       else if (flag.eq.1) then
        do j=1,ngf_bd
         fbdx(j,face) = 0.d0
        enddo
       else if (flag.eq.5) then
        do j=1,ngf_bd
         fbdx(j,face) = fbdx(j,face) + fbdy(j,face)
        enddo
       endif
      enddo

      return
end subroutine massbc


!******************************************************************
! calcul du vecteur normal pour les elements courbes
!
      subroutine cal_normc
      implicit none
      double precision af,bf,x1,dx2,dy2,dx,dy,ax,ay,bx,by,u
      integer cell,i,j,k,l,ns1,tk,tk1,tk2,face

      allocate (normc(ngf_bd,2,nfbd))
      
      do i=1,nfbd
         face=sbd(i)
	 cell=tseg(1,face)
         tk = pbdl_inv(cell)
         ax = vect(1,cell)
         bx = vect(2,cell)
         ay = vect(3,cell)
         by = vect(4,cell)
         af = coeff(1,tk)
         bf = coeff(2,tk)

         do j =1,ngf_bd
           u = tif_bd(j)
           dx2 = - bf * u + 1.d0
           dy2 = ( af + bf ) * u
           dx = ax * dx2 + bx * dy2
           dy = ay * dx2 + by * dy2
           normc(j,1,i) =  dy*2.d0
           normc(j,2,i) = -dx*2.d0
         enddo
      enddo

      return
      end subroutine cal_normc


!**************************************************************
! calcul de la fonction de transformation
!
      subroutine cal_courb_jac
      implicit none

      double precision jac1,jac2,jac3,jac4,delta,x,y,x1,y1,af,bf
      integer cell,i,j
      
      allocate(jac(ng_bd,nfbd))

      do cell=1,ne
        if (bdl(cell)) then
	
	 i=pbdl_inv(cell)
         af = coeff(1,i)
         bf = coeff(2,i)

         do j=1,ng_bd
            x = xi_bd(j)
            y = yi_bd(j)

            jac1 =-bf * ( x + 0.5d0 + 0.5d0*y ) + 1.d0
            jac2 =-bf * ( x + 1.d0 ) * 0.5d0 
            jac3 = (bf+af) * ( x + 0.5d0 + 0.5d0*y )
            jac4 = (bf+af) * ( x + 1.d0 )  * 0.5d0 + 1.d0


            jac(j,i)   = jac1 * jac4 - jac2 * jac3
         enddo
        endif
      enddo

      return
      end subroutine cal_courb_jac



!*****************************************************************
! calcul des matrices pour les elements courbes
!*****************************************************************


      subroutine cal_mat_courb
      implicit none
      integer i,cell,j,k,info
      double precision amm1(nm,nm),amm2(nm,nm), &
                      jac0(ng_bd), &
                      jacxx(ng_bd), jacxy(ng_bd),&
                      jacyx(ng_bd), jacyy(ng_bd),&
      jac1,jac2,jac3,jac4,delta,x,y,x1,y1,af,bf
      external f2,dfx,dfy
      double precision f2,dfx,dfy
      
      
      allocate(amm_bdo(nm,nm,nfbd))
      allocate(amm_bdl(nm,nm,nfbd))
      allocate(amm_bd (nm,nm,nfbd))
      allocate(agx_bd (nm,nm,nfbd))
      allocate(agy_bd (nm,nm,nfbd))
      
      
      do cell=1,ne
        if (bdl(cell)) then
	
	 i=pbdl_inv(cell)
         if (abs(coeff(1,i))+abs(coeff(2,i)).gt.1.d-6) then

         af = coeff(1,i)
         bf = coeff(2,i)
	 
         do j=1,ng_bd
            x = xi_bd(j)
            y = yi_bd(j)

            jac1 =-bf * ( x + 0.5d0 + 0.5d0*y ) + 1.d0
            jac2 =-bf * ( x + 1.d0 ) * 0.5d0 
            jac3 = (bf+af) * ( x + 0.5d0 + 0.5d0*y )
            jac4 = (bf+af) * ( x + 1.d0 )  * 0.5d0 + 1.d0

            jacxx(j) = jac4
            jacxy(j) = -jac2
            jacyx(j) = -jac3
            jacyy(j) = jac1
            jac0(j) = jac(j,i)
         enddo
	 

         call cal_matc(amm1,f2,jac0)
         do j=1,nm
            do k=1,nm
            amm_bdo(k,j,i) = amm1(k,j)
            enddo
         enddo


         call cal_matl(nm,nc,amm1,amm2)


         call dpotrf('U',nm,amm1,nm,info)
      if (info.ne.0) write(*,*) 'pb avec matrice',i
         do j=1,nm
            do k=1,nm
            amm_bd(k,j,i) = amm1(k,j)
            enddo
         enddo
         call dpotrf('U',nm-nc-1,amm2,nm,info)
      if (info.ne.0) write(*,*) 'pb avec matrice',i
         do j=1,nm
            do k=1,nm
            amm_bdl(k,j,i) = amm2(k,j)
            enddo
         enddo

!--------------------- matrice gradient dans direction 1 ---------

         call cal_matc(amm1,dfx,jacxx)
         call cal_matc(amm2,dfy,jacyx)

         do j=1,nm
            do k=1,nm
            agx_bd(k,j,i) = amm1(k,j) + amm2(k,j)
            enddo
         enddo


!--------------------- matrice gradient dans direction 2 ---------

         call cal_matc(amm1,dfx,jacxy)
         call cal_matc(amm2,dfy,jacyy)

         do j=1,nm
            do k=1,nm
            agy_bd(k,j,i) = amm1(k,j) + amm2(k,j)
            enddo
         enddo

        else
         do j=1,nm
          do k=1,nm
             amm_bdo(k,j,i)=ammo(k,j)
             amm_bdl(k,j,i)=amml(k,j)
             amm_bd(k,j,i) = amm(k,j)
             agx_bd(k,j,i) = agx(k,j)
             agy_bd(k,j,i) = agy(k,j)
          enddo
         enddo
        endif
       endif
      enddo
      
      return
      end subroutine cal_mat_courb



!******************************************************************
! integration dans triangle dans element courbe
!******************************************************************

      subroutine cal_matc(mat,ff,jac0)
      implicit none
      double precision mat(nm,nm),jac0(ng_bd)
      double precision xa,xb,ans,x,ff
      integer kc1,kc2,kc3,kc4
      common /i4_cheb/ kc1,kc2,kc3,kc4
      integer i,j,i1,i2
      external xa,xb,ff

! integration dans triangle

      mat=0.d0

      i1=0
        do kc1=0,nc
        do kc2=0,nc-kc1
          i1=i1+1
          i2=0
          do kc3=0,nc
          do kc4=0,nc-kc3
             i2=i2+1
             call inttric(ff,ans,jac0)
             mat(i1,i2)=ans
          enddo
          enddo
        enddo
        enddo

      return
      end subroutine cal_matc


!**************************************************
! integration des polynomes sur le triangle courbe
!**************************************************

      subroutine inttric(ff,ans,jac0)
      implicit none
      double precision ff,ans,x,y,jac0(ng_bd)
      external ff
      integer i,j

      ans=0.d0
      do j=1,ng_bd
         x = xi_bd(j)
         y = yi_bd(j)
         ans = ans + wi_bd(j) * ff(x,y) * jac0(j)
      enddo

      return
      end subroutine inttric



!**************************************************************
! calcul de la normal et tangente dans le triangle courbe
!

      subroutine cal_courb_nx
      implicit none 
      double precision dnx,dny,dd,x,y,af,bf,ax,bx,ay,by,dx1,dy1
      integer cell,i,j
      
      allocate(nxx(ng_bd,nfbd))
      allocate(nyy(ng_bd,nfbd))

      do cell=1,ne
        if (bdl(cell)) then
	
	 i=pbdl_inv(cell)
         ax = vect(1,cell)
         bx = vect(2,cell)
         ay = vect(3,cell)
         by = vect(4,cell)
         af = coeff(1,i)
         bf = coeff(2,i)

         do j=1,ng_bd
            x = xi_bd(j)
            y = yi_bd(j)

            dnx =-bf * ( x + 0.5d0 + 0.5d0*y ) + 1.d0
            dny = (af+bf) * ( x + 0.5d0 + 0.5d0*y ) 
            dx1 = dnx * ax + dny * bx
            dy1 = dnx * ay + dny * by
            dd  = sqrt ( dx1**2 + dy1**2 )
            nxx(j,i) = dx1 / dd
            nyy(j,i) = dy1 / dd
         enddo

        endif
      enddo

      return
      end subroutine cal_courb_nx


!**************************************************************
! calcul des derivees la normal et tangente a la face du triangle courbe
!
      subroutine cal_courb_nx2
      implicit none
      double precision jacxx,jacxy,jacyx,jacyy,jac1,jac2,jac3,jac4,&
      dnx,dny,dd,x,y,af,bf,ax,bx,ay,by,dx1,dy1,ddxx,ddxy,ddyx,ddyy,&
      dxx1,dxy1,dyx1,dyy1,dxx2,dxy2,dyx2,dyy2,detai,&
      dxx3,dxy3,dyx3,dyy3
      integer cell,i,j
      
      allocate(dnxx(ng_bd,nfbd), dnxy(ng_bd,nfbd))
      allocate(dnyy(ng_bd,nfbd), dnyx(ng_bd,nfbd))

      do cell=1,ne
        if (bdl(cell)) then
	
	 i=pbdl_inv(cell)
	 
         ax = vect(1,cell)
         bx = vect(2,cell)
         ay = vect(3,cell)
         by = vect(4,cell)
         af = coeff(1,i)
         bf = coeff(2,i)
         detai= vect(5,cell)

         do j=1,ng_bd
            x = xi_bd(j)
            y = yi_bd(j)

            jac1 =-bf * ( x + 0.5d0 + 0.5d0*y ) + 1.d0
            jac2 =-bf * ( x + 1.d0 ) * 0.5d0 
            jac3 = (bf+af) * ( x + 0.5d0 + 0.5d0*y )
            jac4 = (bf+af) * ( x + 1.d0 )  * 0.5d0 + 1.d0

            jacxx = jac4
            jacxy = -jac2
            jacyx = -jac3
            jacyy = jac1

            dnx =-bf * ( x + 0.5d0 + 0.5d0*y ) + 1.d0
            dny = (af+bf) * ( x + 0.5d0 + 0.5d0*y ) 

            dx1 = dnx * ax + dny * bx
            dy1 = dnx * ay + dny * by
            dd  = 1.d0/sqrt ( dx1**2 + dy1**2 )

! ddxx = d(dnx)/dxsi1
! ddxy = d(dnx)/dxsi2
! ddyx = d(dny)/dxsi1  !! derivees tangente locale dans le repere local
! ddyy = d(dny)/dxsi2

             ddxx = - bf
             ddxy = - 0.5d0 * bf
             ddyx = af + bf 
             ddyy = 0.5d0 * ( af + bf )

! dxx1 = d(dx1)/dxsi1
! dxy1 = d(dx1)/dxsi2
! dyx1 = d(dy1)/dxsi1 !! derivees tangente non normalisee dans le repere local
! dyy1 = d(dy1)/dxsi2

            dxx1 = ddxx * ax + ddyx * bx
            dxy1 = ddxy * ax + ddyy * bx
            dyx1 = ddxx * ay + ddyx * by
            dyy1 = ddxy * ay + ddyy * by

! dxx2 = d(dx1*dd)/dxsi1
! dxy2 = d(dx1*dd)/dxsi2
! dyx2 = d(dy1*dd)/dxsi1 !! derivees tangente normalisee dans le repere local
! dyy2 = d(dy1*dd)/dxsi2

            dxx2 = dxx1 * dd - dx1 * (dx1*dxx1+dy1*dyx1) * dd**3
            dxy2 = dxy1 * dd - dx1 * (dx1*dxy1+dy1*dyy1) * dd**3 
            dyx2 = dyx1 * dd - dy1 * (dx1*dxx1+dy1*dyx1) * dd**3 
            dyy2 = dyy1 * dd - dy1 * (dx1*dxy1+dy1*dyy1) * dd**3 

! dxx3 = d(dx1*dd)/dxsi11
! dxy3 = d(dx1*dd)/dxsi21
! dyx3 = d(dy1*dd)/dxsi11 !! derivees tangente normalisee 
! dyy3 = d(dy1*dd)/dxsi21 !! dans le repere local non courbe

            dxx3 = dxx2 * jacxx + dxy2 * jacyx
            dxy3 = dxx2 * jacxy + dxy2 * jacyy
            dyx3 = dyx2 * jacxx + dyy2 * jacyx
            dyy3 = dyx2 * jacxy + dyy2 * jacyy

! dnxx = d(dx1*dd)/dx
! dnxy = d(dx1*dd)/dy
! dnyx = d(dy1*dd)/dx !! derivees tangente normalisee dans le repere absolu
! dnyy = d(dy1*dd)/dy

            dnxx(j,i) = detai * ( by * dxx3 - ay * dxy3)
            dnxy(j,i) = detai * (-bx * dxx3 + ax * dxy3)
            dnyx(j,i) = detai * ( by * dyx3 - ay * dyy3)
            dnyy(j,i) = detai * (-bx * dyx3 + ax * dyy3)

! vraies derivees :

            dnxx(j,i) = dnxx(j,i) / jac(j,i)
            dnxy(j,i) = dnxy(j,i) / jac(j,i)
            dnyx(j,i) = dnyx(j,i) / jac(j,i)
            dnyy(j,i) = dnyy(j,i) / jac(j,i)

         enddo
	 
	endif

      enddo

      return
      end subroutine cal_courb_nx2


!******************************************************************

      subroutine coriolis_c(fc,f0,beta)
      implicit none
      double precision fc(nm,*),dia1(nm),dia2(nm)
      double precision b(nm),jac0(ng_bd),mff(ng_bd)
      double precision f0,beta,af,bf, &
      agx,agy,bgx,bgy,xof,yof,x,y,x1,y1,deta,ans,x2,y2
      integer cell,i1,i2,i3,j,k,tk
      integer kc1,kc2
      double precision amm1(nm,nm)

      do j=1,nm
         dia1(j) = 1.d0/amm(j,j)
      enddo


      do cell=1,ne

         i1 = in(1,cell)
         i2 = in(2,cell)
         i3 = in(3,cell)
         agx = (xgr(i2)-xgr(i1))*0.5d0
         agy = (ygr(i2)-ygr(i1))*0.5d0
         bgx = (xgr(i3)-xgr(i1))*0.5d0
         bgy = (ygr(i3)-ygr(i1))*0.5d0
         xof = xgr(i1)
         yof = ygr(i1)

         if (bdl(cell)) then

         tk = pbdl_inv(cell)
         af = coeff(1,tk)
         bf = coeff(2,tk)

         do j=1,ng_bd
            x = xi_bd(j)
            y = yi_bd(j)
            x2 = x    - bf  * 0.5d0 * (x**2 + x + x*y + y)
            y2 = y+ (af+bf) * 0.5d0 * (x**2 + x + x*y + y)
            x1 = agx * (x2+1.d0) + bgx * (y2+1.d0) + xof
            y1 = agy * (x2+1.d0) + bgy * (y2+1.d0) + yof
            mff(j) = beta * y1
         enddo

         amm1(:,:) = amm_bd(:,:,tk)
         do k=1,nm
          dia2(k) = 1.d0/amm1(k,k)
         enddo

         call xytospec_c(b,mff,jac(:,tk) )
         call dtrsm4(nm,nm,amm1,b,dia2)

        else

         do j=1,ng
            x = xi(j)
            y = yi(j)
            x1 = agx * (x+1.d0) + bgx * (y+1.d0) + xof
            y1 = agy * (x+1.d0) + bgy * (y+1.d0) + yof
            mff(j) = beta * y1
         enddo

         call xytospec(b,mff)
         call dtrsm4(nm,nm,amm,b,dia1)

      endif

        do j=1,nm
            fc(j,cell) = b(j)
        enddo
        fc(1,cell) = fc(1,cell) + f0


      enddo

      return
      end subroutine coriolis_c


!******************************************************************

      subroutine calh_c(ff,h)
      implicit none
      double precision h(nm,*),b(nm),jac0(ng_bd),mff(ng_bd)
      double precision af,bf,dia1(nm),dia2(nm),&
      agx,agy,bgx,bgy,xof,yof,x,y,x1,y1,deta,ans,x2,y2
      integer cell,i1,i2,i3,j,k,tk
      double precision ff
      external ff
      double precision a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
      integer kc1,kc2
      double precision amm1(nm,nm)

      do j=1,nm
         dia1(j) = 1.d0/amm(j,j)
      enddo


      do cell=1,ne

         i1 = in(1,cell)
         i2 = in(2,cell)
         i3 = in(3,cell)
         agx = (xgr(i2)-xgr(i1))*0.5d0
         agy = (ygr(i2)-ygr(i1))*0.5d0
         bgx = (xgr(i3)-xgr(i1))*0.5d0
         bgy = (ygr(i3)-ygr(i1))*0.5d0
         xof = xgr(i1)
         yof = ygr(i1)

         if (bdl(cell)) then

         tk = pbdl_inv(cell)
         af = coeff(1,tk)
         bf = coeff(2,tk)
         
         do j=1,ng_bd
            x = xi_bd(j)
            y = yi_bd(j)
            x2 = x    - bf  * 0.5d0 * (x**2 + x + x*y + y)
            y2 = y+ (af+bf) * 0.5d0 * (x**2 + x + x*y + y)
            x1 = agx * (x2+1.d0) + bgx * (y2+1.d0) + xof
            y1 = agy * (x2+1.d0) + bgy * (y2+1.d0) + yof
            mff(j) = ff(x1,y1)
         enddo

         amm1(:,:) = amm_bd(:,:,tk)
         do k=1,nm
          dia2(k) = 1.d0/amm1(k,k)
         enddo

         call xytospec_c(b,mff,jac(:,tk) )
         call dtrsm4(nm,nm,amm1,b,dia2)

        else

         do j=1,ng
            x = xi(j)
            y = yi(j)
            x1 = agx * (x+1.d0) + bgx * (y+1.d0) + xof
            y1 = agy * (x+1.d0) + bgy * (y+1.d0) + yof
            mff(j) = ff(x1,y1)
         enddo

         call xytospec(b,mff)
         call dtrsm4(nm,nm,amm,b,dia1)

      endif

         do j=1,nm
            h(j,cell) = b(j)
         enddo

      enddo

      return
      end subroutine calh_c



!******************************************************************

      subroutine calh_c2(filename,h)
      implicit none
      double precision h(nm,*),b(nm),jac0(ng_bd),mff(ng_bd)
      double precision af,bf,dia1(nm),dia2(nm),&
      agx,agy,bgx,bgy,xof,yof,x,y,x1,y1,deta,ans,x2,y2
      integer cell,i1,i2,i3,j,k,tk
      double precision ff
      external ff
      double precision a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
      integer kc1,kc2
      double precision amm1(nm,nm)
      character*10 filename

      write(*,*) 'filename for interpolation: ',filename
      open(1,file=filename,status='old')

      do j=1,nm
         dia1(j) = 1.d0/amm(j,j)
      enddo


      do cell=1,ne

         if (bdl(cell)) then

         tk = pbdl_inv(cell)
         af = coeff(1,tk)
         bf = coeff(2,tk)
         
         do j=1,ng_bd
	    read(1,*) mff(j)
         enddo

         amm1(:,:) = amm_bd(:,:,tk)
         do k=1,nm
          dia2(k) = 1.d0/amm1(k,k)
         enddo

         call xytospec_c(b,mff,jac(:,tk) )
         call dtrsm4(nm,nm,amm1,b,dia2)

        else

         do j=1,ng
	    read(1,*) mff(j)
         enddo

         call xytospec(b,mff)
         call dtrsm4(nm,nm,amm,b,dia1)

      endif

         do j=1,nm
            h(j,cell) = b(j)
         enddo

      enddo
      
      close(1)

      return
      end subroutine calh_c2
      

!**************************************************************
! transformation de l'espace physique a l'espace spectral
! pour des elements courbes
!
      subroutine xytospec_c(b,mff,jac0)
      implicit none
      double precision b(*),jac0(ng_bd)
      double precision mff(ng_bd),mff1(ng_bd),ans
      integer tk,i,j,k,i1,kc1,kc2

      do j=1,ng_bd
         mff1(j) =  wi_bd(j) * mff(j) * jac0(j)
      enddo

      do i1 =1,nm
         ans=0.d0
         do j=1,ng_bd
            ans = ans + mff1(j) * pi2_bd(j,i1)
         enddo
         b(i1)=ans
      enddo

      return
      end subroutine xytospec_c

!**************************************************
! integration des polynomes sur le triangle courbe
!**************************************************

      subroutine inttri2c(mff,ans,i1)
      implicit none
      double precision mff(ng_bd),ans
      integer i1,j

      ans=0.d0
      do j=1,ng_bd
         ans = ans + wi_bd(j) * mff(j) * pi2_bd(j,i1)
      enddo

      return
      end subroutine inttri2c

!
!*******************************************************************
!
      subroutine spectoxy_c(mff,ak)
      use gauss
      implicit none
      double precision mff(ng_bd),ak(nm)
      integer j,k,kc1,kc2,i1

	     mff = 0.d0

	 i1=0
	 do kc1=0,nc
	  do kc2=0,nc-kc1
	   i1=i1+1
	   do j=1,ng_bd
	     mff(j) = mff(j) + ak(i1) * pi2_bd(j,i1)
	   enddo
	  enddo
	 enddo

	return
	end subroutine spectoxy_c

!******************************************************************

end module boundary_bd
