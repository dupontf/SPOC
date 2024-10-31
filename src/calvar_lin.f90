!
!******************************************************************
! linear version
! no viscosity, no beta plane
! curved elements
! open boudary
!******************************************************************
!
subroutine cal_var(u0,v0,h0,u1,v1,h1,mu0,mb0,fc,w0,z0,g0,hb,time)
     
    use mesh
    use gauss
    use utilities
    use boundary_condition
    use gauss_bd
    use boundary_bd
    
      implicit none
      double precision u0(nm,ne),v0(nm,ne),h0(nm,ne), &
                       u1(nm,ne),v1(nm,ne),h1(nm,ne),&
                 w0(nm,ne),z0(nm,ne),hb(nm,ne),&
                 dia1(nm),dia2(nm),&
                 dia3(nm),dia4(nm), dia5(nm), &
                 ak(nm),b(nm),bu(nm),&
                 uk(nm),vk(nm),hk(nm),wk(nm),zk(nm),&
                 uhk(nm),vhk(nm),&
                 uux(nm),vuy(nm),&
                 uvx(nm),vvy(nm),&
                 gx(nm),gy(nm),gn(nm),gt(nm),&
                 gx1(nm),gy1(nm),&
		 amm0(nm,nm),amm1(nm,nm),amm2(nm,nm), &
		 amm3(nm,nm), amm4(nm,nm), &
                 ax,by,ay,bx,deta,detai,dis,dnx,dny,mu,mv,&
                 mu0,mb0,fc(nm,ne),g0,time,pi,tday,hh0,&
                 gxx(nm),gxy(nm),gyx(nm),gyy(nm),&
                 gnx(nm),gtx(nm),gny(nm),gty(nm),&
                 gnn(nm),gtn(nm),gnt(nm),gtt(nm),&
                 qxx(nm,ne),qxy(nm,ne),&
                 qyx(nm,ne),qyy(nm,ne),&
                 fxx(nm,ne),fxy(nm,ne),&
                 fyx(nm,ne),fyy(nm,ne),&
                 fbdxx(ngf_bd,nf),fbdxy(ngf_bd,nf),&
                 fbdyx(ngf_bd,nf),fbdyy(ngf_bd,nf),&
                 fbdx (ngf_bd,nf),fbdy (ngf_bd,nf),&
                 fbd(ngf_bd,nf),&
                 adu(nm,ne),adv(nm,ne),&
                 uh(nm,ne),vh(nm,ne),&
                 w1(nm,ne),z1(nm,ne),&
                 mff(ng_bd), mff1(ng_bd), &
		 mff2(ng_bd), mff3(ng_bd), mff4(ng_bd), &
		 jac0(ng_bd), jac1(ng_bd), jac2(ng_bd)
      integer i,j,k,l,info,i1,i2,i3,ind,ns1,tk
      integer cell,face


!----------- linear stuff definition---------------

      hh0=hb(1,1)
      w1 = w0 / hh0 - mb0 * u0
      z1 = z0 / hh0 - mb0 * v0
      
!----------- matrix inversion---------------

      do i=1,nm
         dia1(i) = 1.d0/amm(i,i)
      enddo
      do i=1,nm-nc-1
         dia2(i) = 1.d0/amml(i,i)
      enddo

!------------------------------------------------
!--- u,v
!------------------------------------------------

	fxx = -g0 * h0
	fxy = 0.d0
	fyx = 0.d0
	fyy = -g0 * h0

	call cal_flu(fbd,fxx)
	call openbc1(fbd,time,-g0)
	
	call norm_flu(fbd,fbdxx,1)
	call norm_flu(fbd,fbdyy,2)
	fbdxy=0.d0
	fbdyx=0.d0
		

	do cell=1,ne

         ax    = vect(1,cell)
         bx    = vect(2,cell)
         ay    = vect(3,cell)
         by    = vect(4,cell)
         detai = vect(5,cell)
	 	 
	 gx=0.d0
	 gy=0.d0
!
!------------------------------------------------
! wind forcing, coriolis forcing, advection
!
        uk = w1(:,cell)
        vk = z1(:,cell)
        gxx = fxx(:,cell)
        gxy = fxy(:,cell)
        gyx = fyx(:,cell)
        gyy = fyy(:,cell)


       tk=pbdl_inv(cell)


!
!------------------------------------------------
! si la cellule est au bord, rotation
!------------------------------------------------
!
	 if (fflagt(cell).eq.1) then
	 
	 
	 amm0=amm_bdo(:,:,tk)
         call cal_b(gx,uk,amm0,1.d0)
         call cal_b(gy,vk,amm0,1.d0)
	 
	 amm0=amm_bd(:,:,tk)
	 amm1=amm_bdl(:,:,tk)
	 amm2=agx_bd(:,:,tk)
	 amm3=agy_bd(:,:,tk)

         do k=1,nm-nc-1
          dia4(k) = 1.d0/amm1(k,k)
         enddo
         do k=1,nm
          dia5(k) = 1.d0/amm0(k,k)
         enddo

            jac0 = jac(:,tk)
            jac1 = nxx(:,tk)
            jac2 = nyy(:,tk)

!------------------------------------------------
! calcul de tenseur pour u
!
         call cal_b(gx,gxx,amm2,-by*detai)
         call cal_b(gx,gxx,amm3, ay*detai)

         call cal_b(gx,gxy,amm2, bx*detai)
         call cal_b(gx,gxy,amm3,-ax*detai)

         call tri_fc(cell,fbdxx,gx,detai)
         call tri_fc(cell,fbdxy,gx,detai)

!------------------------------------------------
! calcul de tenseur pour v
!
         call cal_b(gy,gyx,amm2,-by*detai)
         call cal_b(gy,gyx,amm3,ay*detai)

         call cal_b(gy,gyy,amm2,bx*detai)
         call cal_b(gy,gyy,amm3,-ax*detai)

         call tri_fc(cell,fbdyx,gy,detai)
         call tri_fc(cell,fbdyy,gy,detai)


!------------------------------------------------
! reduction
!
         call dtrsm4(nm,nm,amm0,gx,dia5)
         call dtrsm4(nm,nm,amm0,gy,dia5)
	 
	 call rot_vec_c(gt,gn,gx,gy,jac0,jac1,jac2)
	 
         call dtrsm4(nm,nm     ,amm0,gt,dia5)
         call reduc_gn(gn,bu)
         call dtrsm4(nm,nm-nc-1,amm1,bu,dia4)
         call reduc_inv_gn(gn,bu)

	 call rot_vec_c(gx,gy,gt,gn,jac0,jac1,jac2)

         call dtrsm4(nm,nm,amm0,gx,dia5)
         call dtrsm4(nm,nm,amm0,gy,dia5)


!******* open boudanry condition *******

	 else if (fflagt(cell).eq.5) then

	 amm0=amm_bdo(:,:,tk)
         call cal_b(gx,uk,amm0,1.d0)
         call cal_b(gy,vk,amm0,1.d0)
	 
	 amm0=amm_bd(:,:,tk)
	 amm2=agx_bd(:,:,tk)
	 amm3=agy_bd(:,:,tk)

         do k=1,nm-nc-1
          dia4(k) = 1.d0/amm1(k,k)
         enddo
         do k=1,nm
          dia5(k) = 1.d0/amm0(k,k)
         enddo

            jac0 = jac(:,tk)
            jac1 = nxx(:,tk)
            jac2 = nyy(:,tk)

!------------------------------------------------
! calcul de tenseur pour u
!
         call cal_b(gx,gxx,amm2,-by*detai)
         call cal_b(gx,gxx,amm3, ay*detai)

         call cal_b(gx,gxy,amm2, bx*detai)
         call cal_b(gx,gxy,amm3,-ax*detai)

         call tri_fc(cell,fbdxx,gx,detai)
         call tri_fc(cell,fbdxy,gx,detai)

!------------------------------------------------
! calcul de tenseur pour v
!
         call cal_b(gy,gyx,amm2,-by*detai)
         call cal_b(gy,gyx,amm3,ay*detai)

         call cal_b(gy,gyy,amm2,bx*detai)
         call cal_b(gy,gyy,amm3,-ax*detai)

         call tri_fc(cell,fbdyx,gy,detai)
         call tri_fc(cell,fbdyy,gy,detai)

!--- solution

         call dtrsm4(nm,nm,amm0,gx,dia5)
         call dtrsm4(nm,nm,amm0,gy,dia5)



! ******* interior element *******
	   else

        call cal_b(gx,uk,ammo,1.d0)
        call cal_b(gy,vk,ammo,1.d0)

!------------------------------------------------

	      gxx = fxx(:,cell)
	      gxy = fxy(:,cell)
	      gyx = fyx(:,cell)
	      gyy = fyy(:,cell)
!
!------------------------------------------------
! calcul de tenseur pour u
!
	   call cal_b(gx,gxx,agx,-by*detai)
	   call cal_b(gx,gxx,agy, ay*detai)

	   call cal_b(gx,gxy,agx, bx*detai)
	   call cal_b(gx,gxy,agy,-ax*detai)

	   call tri_f(cell,fbdxx,gx,detai)
	   call tri_f(cell,fbdxy,gx,detai)
!
!------------------------------------------------
! calcul de tenseur pour v
!
	   call cal_b(gy,gyx,agx,-by*detai)
	   call cal_b(gy,gyx,agy,ay*detai)

	   call cal_b(gy,gyy,agx,bx*detai)
	   call cal_b(gy,gyy,agy,-ax*detai)

	   call tri_f(cell,fbdyx,gy,detai)
	   call tri_f(cell,fbdyy,gy,detai)

	     call dtrsm4(nm,nm,amm,gx,dia1)
	     call dtrsm4(nm,nm,amm,gy,dia1)

	   endif

		u1(:,cell) = gx
		v1(:,cell) = gy

	enddo

!
!-------------------------------------------------
!--- h
!-------------------------------------------------
!
         uh = hh0 * u0
         vh = hh0 * v0

	call cal_flu(fbd,uh)
	call norm_flu(fbd,fbdx,1)
	call cal_flu(fbd,vh)
	call norm_flu(fbd,fbdy,2)

	call massbc(fbdx,fbdy)


	do cell=1,ne

         ax    = vect(1,cell)
         bx    = vect(2,cell)
         ay    = vect(3,cell)
         by    = vect(4,cell)
         detai = vect(5,cell)

         b =0.d0
         gx = uh(:,cell)
         gy = vh(:,cell)

         tk=pbdl_inv(cell)

!------------------------------------------------
! element sur le bord --- element courbe
!------------------------------------------------
       if (tk.gt.0) then

	 amm1=amm_bd(:,:,tk)
	 amm2=agx_bd(:,:,tk)
	 amm3=agy_bd(:,:,tk)

        do k=1,nm
         dia4(k) = 1.d0/amm1(k,k)
        enddo

         call cal_b(b,gx,amm2, by*detai)
         call cal_b(b,gx,amm3,-ay*detai)

         call cal_b(b,gy,amm2,-bx*detai)
         call cal_b(b,gy,amm3, ax*detai)

         call tri_fc(cell,fbdx,b,-detai)

         call dtrsm4(nm,nm,amm1,b,dia4)
	 
	 
!------------------------------------------------
! interior element
!------------------------------------------------
	 else

          call cal_b(b,gx,agx, by*detai)
          call cal_b(b,gx,agy,-ay*detai)

          call cal_b(b,gy,agx,-bx*detai)
          call cal_b(b,gy,agy, ax*detai)

          call tri_f(cell,fbdx,b,-detai)

	  call dtrsm4(nm,nm,amm,b,dia1)

	 endif
	  
	 h1(:,cell)=b

	enddo


	return
	end
!
!*****************************************************************
