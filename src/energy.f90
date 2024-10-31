!******************************************************************
!
      subroutine energy(u0,v0,h0,g0,hb,ke,pe,mass)
      use mesh
      use gauss
      
      implicit none
      double precision u0(nm,*),v0(nm,*),h0(nm,*),hb(nm,*), &
                      g0,ke,pe,mass,&
                      mffu(ng),mffv(ng),mffh(ng),&
                      mff(ng),mffhb(ng),&
                      uk(nm),vk(nm),hk(nm),hbk(nm),&
                      ans,ax,ay,bx,by,deta
      integer cell,j,k,i1,i2,i3,kc1,kc2

	ke=0.d0
	pe=0.d0
	mass=0.d0

	do cell=1,ne

	   i1 = in(1,cell)
	   i2 = in(2,cell)
	   i3 = in(3,cell)
	   ax = (xgr(i2)-xgr(i1))*0.5d0
	   bx = (xgr(i3)-xgr(i1))*0.5d0
	   ay = (ygr(i2)-ygr(i1))*0.5d0
	   by = (ygr(i3)-ygr(i1))*0.5d0
	   deta = ax * by - ay * bx

	   do j=1,nm
	      uk(j) =  u0(j,cell)
	      vk(j) =  v0(j,cell)
	      hk(j) =  h0(j,cell)
	      hbk(j)=  h0(j,cell)+hb(j,cell)
	   enddo
!
! variation dans la masse totale
!
	   do j=1,nm
	      mass = mass + ammo(j,1) * deta * hk(j)
	   enddo

	 do j=1,ng
	     mff (j) = 0.d0
	     mffu(j) = 0.d0
	     mffh(j) = 0.d0
	     mffhb(j)= 0.d0
	     mffv(j) = 0.d0
	 enddo

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
!
! kinetic energy
!
	 do j=1,ng
	     mff(j) = mffhb(j) * ( mffu(j) * mffu(j) + mffv(j) * mffv(j) )
	 enddo

	 call inttri2(mff,ans,1)
	 ke = ke + ans * deta * 0.5d0

!
! potential energy
!
	 do j=1,ng
	     mff(j) = mffh(j) * mffh(j)
	 enddo

	 call inttri2(mff,ans,1)
	 pe = pe + ans * deta * 0.5d0 * g0

	enddo

	return
	end	
