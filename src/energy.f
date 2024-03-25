c******************************************************************
c
      subroutine energy(ne,u0,v0,h0,g0,hb,ke,pe,mass,
     &                  in,xgr,ygr,ammo)
      implicit none
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision xi(ng_max),yi(ng_max),wi(ng_max)
      common /gauss/ wi,xi,yi
      integer ne
      double precision u0(nnmod,*),v0(nnmod,*),h0(nnmod,*),hb(nnmod,*),
     &                 xgr(*),ygr(*),
     &                 ammo(nnmod,*),
     &                 g0,ke,pe,mass,
     &                 mffu(ng_max),mffv(ng_max),mffh(ng_max),
     &                 mff(ng_max),mffhb(ng_max),
     &                 uk(nnmod),vk(nnmod),hk(nnmod),hbk(nnmod),
     &                 pi2(ng_max,nnmod),
     &                 ans,ax,ay,bx,by,deta
      integer in(3,*)
      integer i,j,k,i1,i2,i3,kc1,kc2
      common /gauss2/ pi2
c
	ke=0.d0
	pe=0.d0
	mass=0.d0
c
	do i=1,ne
c
	   i1 = in(1,i)
	   i2 = in(2,i)
	   i3 = in(3,i)
	   ax = (xgr(i2)-xgr(i1))*0.5d0
	   bx = (xgr(i3)-xgr(i1))*0.5d0
	   ay = (ygr(i2)-ygr(i1))*0.5d0
	   by = (ygr(i3)-ygr(i1))*0.5d0
	   deta = ax * by - ay * bx
c
	   do j=1,nm
	      uk(j) =  u0(j,i)
	      vk(j) =  v0(j,i)
	      hk(j) =  h0(j,i)
	      hbk(j)=  h0(j,i)+hb(j,i)
	   enddo
c
c variation dans la masse totale
c
	   do j=1,nm
	      mass = mass + ammo(j,1) * deta * hk(j)
	   enddo
c
	 do j=1,ng
	     mff (j) = 0.d0
	     mffu(j) = 0.d0
	     mffh(j) = 0.d0
	     mffhb(j)= 0.d0
	     mffv(j) = 0.d0
	 enddo
c
	 i1 = 0
	 do kc1=0,nc
	  do kc2=0,nc-kc1
	   i1=i1+1
	   do j=1,ng
	     mffu(j) = mffu(j) + uk(i1) * pi2(j,i1)
	     mffv(j) = mffv(j) + vk(i1) * pi2(j,i1)
	     mffh(j) = mffh(j) + hk(i1) * pi2(j,i1)
	     mffhb(j)= mffhb(j)+ hbk(i1)* pi2(j,i1)
	   enddo
	  enddo
	 enddo
c
c kinetic energy
c
	 do j=1,ng
	     mff(j) = mffhb(j) * ( mffu(j) * mffu(j) + mffv(j) * mffv(j) )
	 enddo
c
	 call inttri2(mff,ans,1)
	 ke = ke + ans * deta * 0.5d0
c
c
c potential energy
c
	 do j=1,ng
	     mff(j) = mffh(j) * mffh(j)
	 enddo
c
	 call inttri2(mff,ans,1)
	 pe = pe + ans * deta * 0.5d0 * g0
c
	enddo
c
	return
	end	
