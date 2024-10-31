!******************************************************************
!
      subroutine energy(u0,v0,h0,g0,hb,ke,pe,mass)
      use mesh
      use mesh_curved
      use gauss
      use gauss_bd
      use boundary_bd
      
      implicit none
      double precision u0(nm,*),v0(nm,*),h0(nm,*),hb(nm,*), &
                      g0,ke,pe,mass,&
                      mffu(ng_bd),mffv(ng_bd),mffh(ng_bd),&
                      mff(ng_bd),mffhb(ng_bd),&
                      uk(nm),vk(nm),hk(nm),hbk(nm),&
                      ans,ax,ay,bx,by,deta
      integer cell,j,i1,i2,i3,tk

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
	   

	      tk=pbdl_inv(cell)

! *************** boundary element

	   if (bdl(cell)) then

              call spectoxy_c(mffu ,uk )
              call spectoxy_c(mffv ,vk )
              call spectoxy_c(mffh ,hk )
              call spectoxy_c(mffhb,hbk)

!
! variation dans la masse totale
!
         do j=1,nm
            mass = mass + amm_bdo(j,1,tk) * deta * hk(j)
         enddo

! kinetic energy

	 do j=1,ng_bd
	     mff(j) = mffhb(j) * ( mffu(j) * mffu(j) + mffv(j) * mffv(j) )
	 enddo
         mff = mff * jac(:,tk)
	 call inttri2c(mff,ans,1)
	 ke = ke + ans * deta * 0.5d0


! potential energy

	 do j=1,ng_bd
	     mff(j) = mffh(j) * mffh(j)
	 enddo
         mff = mff * jac(:,tk)
	 call inttri2c(mff,ans,1)
	 pe = pe + ans * deta * 0.5d0 * g0



! ******** interior element

	   else


! variation dans la masse totale

	   do j=1,nm
	      mass = mass + ammo(j,1) * deta * hk(j)
	   enddo

              call spectoxy(mffu ,uk )
              call spectoxy(mffv ,vk )
              call spectoxy(mffh ,hk )
              call spectoxy(mffhb,hbk)


! kinetic energy

	 do j=1,ng
	     mff(j) = mffhb(j) * ( mffu(j) * mffu(j) + mffv(j) * mffv(j) )
	 enddo

	 call inttri2(mff,ans,1)
	 ke = ke + ans * deta * 0.5d0


! potential energy

	 do j=1,ng
	     mff(j) = mffh(j) * mffh(j)
	 enddo

	 call inttri2(mff,ans,1)
	 pe = pe + ans * deta * 0.5d0 * g0
	 
	 endif

	enddo

	return
	end	
