!
!******************************************************************
! nonlinear version
! free-slip
! viscosity
! beta plane
! varying bottom topography
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
! --- some scalar
                 ax,by,ay,bx,deta,detai,dis,dnx,dny,mu,mv,&
                 mu0,mb0,fc(nm,ne),g0,time,pi,tday,hh0,&
! --- some other stuff
                 gxx(nm),gxy(nm),gyx(nm),gyy(nm),&
                 gnx(nm),gtx(nm),gny(nm),gty(nm),&
                 gnn(nm),gtn(nm),gnt(nm),gtt(nm),&
                 qxx(nm,ne),qxy(nm,ne),&
                 qyx(nm,ne),qyy(nm,ne),&
                 fxx(nm,ne),fxy(nm,ne),&
                 fyx(nm,ne),fyy(nm,ne),&
                 adu(nm,ne),adv(nm,ne),&
                 uh(nm,ne),vh(nm,ne),&
                 w1(nm,ne),z1(nm,ne),&
! --- some gauss tables
                 fbdxx(ngf_bd,nf),fbdxy(ngf_bd,nf),&
                 fbdyx(ngf_bd,nf),fbdyy(ngf_bd,nf),&
                 fbdx(ngf_bd,nf),fbdy(ngf_bd,nf),&
                 fbd(ngf_bd,nf),&
                 mff(ng_bd), mff1(ng_bd), &
		 mff2(ng_bd), mff3(ng_bd), mff4(ng_bd), &
		 jac0(ng_bd), jac1(ng_bd), jac2(ng_bd), &
		 jac11(ng_bd), jac12(ng_bd), &
		 jac21(ng_bd), jac22(ng_bd), &
		 ut(ng_bd), un(ng_bd), &
		 uxx(ng_bd), uxy(ng_bd), uyx(ng_bd), uyy(ng_bd), &
		 utx(ng_bd), uty(ng_bd), unx(ng_bd), uny(ng_bd), &
		 utt(ng_bd), utn(ng_bd), unt(ng_bd), unn(ng_bd)
      integer i,j,k,l,info,i1,i2,i3,ind,ns1,tk
      integer cell,face

      hh0=hb(1,1)

      
!----------- matrix inversion---------------

      do i=1,nm
         dia1(i) = 1.d0/amm(i,i)
      enddo
      do i=1,nm-nc-1
         dia2(i) = 1.d0/amml(i,i)
      enddo


!
!------------------------------------
! compute the gradient terms of u,v
!------------------------------------
!

	call cal_flu(fbd,u0)
	call norm_flu(fbd,fbdxx,1)
	call norm_flu(fbd,fbdxy,2)


	call cal_flu(fbd,v0)
	call norm_flu(fbd,fbdyx,1)
	call norm_flu(fbd,fbdyy,2)


	do cell=1,ne

         ax    = vect(1,cell)
         bx    = vect(2,cell)
         ay    = vect(3,cell)
         by    = vect(4,cell)
         detai = vect(5,cell)

	      uk = u0(:,cell)
	      vk = v0(:,cell)

!
!------------------------------------------------
! si la cellule est au bord, rotation
!------------------------------------------------
!
         tk=pbdl_inv(cell)

	 if (bdl(cell)) then

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

!------------------------------------------------
	      gx1=0.d0
	      gy1=0.d0

	   call cal_b(gx1,uk,amm2,-1.d0)
	   call cal_b(gy1,uk,amm3,-1.d0)

	   gxx = (  by * gx1 - ay * gy1) * detai
	   gxy = (- bx * gx1 + ax * gy1) * detai

           call tri_fc(cell,fbdxx,gxx,detai)
  	   call tri_fc(cell,fbdxy,gxy,detai)

!------------------------------------------------
	      gx1=0.d0
	      gy1=0.d0

	   call cal_b(gx1,vk,amm2,-1.d0)
	   call cal_b(gy1,vk,amm3,-1.d0)

	   gyx = (  by * gx1 - ay * gy1) * detai
	   gyy = (- bx * gx1 + ax * gy1) * detai

           call tri_fc(cell,fbdyx,gyx,detai)
  	   call tri_fc(cell,fbdyy,gyy,detai)
!
!------------------------------------------------
!       rotation des vitesse pour que la vitesse v soit 
!       normale a la face
!
         call dtrsm4(nm,nm,amm0,gxx,dia5)
         call dtrsm4(nm,nm,amm0,gxy,dia5)
         call dtrsm4(nm,nm,amm0,gyx,dia5)
         call dtrsm4(nm,nm,amm0,gyy,dia5)


!-------------------------------------------------
! cas normal sans rotation
!-------------------------------------------------

	   else

!------------------------------------------------
	      gx1=0.d0
	      gy1=0.d0

	   call cal_b(gx1,uk,agx,-1.d0)
	   call cal_b(gy1,uk,agy,-1.d0)

	   gxx = (  by * gx1 - ay * gy1) * detai
	   gxy = (- bx * gx1 + ax * gy1) * detai

           call tri_f(cell,fbdxx,gxx,detai)
  	   call tri_f(cell,fbdxy,gxy,detai)

!------------------------------------------------
	      gx1=0.d0
	      gy1=0.d0

	   call cal_b(gx1,vk,agx,-1.d0)
	   call cal_b(gy1,vk,agy,-1.d0)

	   gyx = (  by * gx1 - ay * gy1) * detai
	   gyy = (- bx * gx1 + ax * gy1) * detai

           call tri_f(cell,fbdyx,gyx,detai)
  	   call tri_f(cell,fbdyy,gyy,detai)

	   call dtrsm4(nm,nm,amm,gxx,dia1)
	   call dtrsm4(nm,nm,amm,gxy,dia1)
	   call dtrsm4(nm,nm,amm,gyx,dia1)
	   call dtrsm4(nm,nm,amm,gyy,dia1)

	   endif

	   do j=1,nm
		qxx(j,cell) = gxx(j)
		qxy(j,cell) = gxy(j)
		qyx(j,cell) = gyx(j)
		qyy(j,cell) = gyy(j)
	   enddo

	enddo


!------------------------------------------------
!--- compute the nonlinear terms
!------------------------------------------------

      do cell=1,ne

	   do j=1,nm
	      uk(j) = u0(j,cell)
	      vk(j) = v0(j,cell)
	      hk(j) = h0(j,cell)+hb(j,cell)
	      wk(j) = w0(j,cell)
	      zk(j) = z0(j,cell)
	      gxx(j) = qxx(j,cell)
	      gxy(j) = qxy(j,cell) - fc(j,cell)
	      gyx(j) = qyx(j,cell) + fc(j,cell)
	      gyy(j) = qyy(j,cell)
	   enddo


         tk=pbdl_inv(cell)

	 if (fflagt(cell).eq.1) then

           jac0 = jac(:,tk)

	   call spectoxy_c(mff1,uk)
	   call spectoxy_c(mff2,vk)

	   call spectoxy_c(mff3,gxx)
	   call spectoxy_c(mff4,gxy)

	   call multiply(ng_bd,mff,mff1,mff3)
	   call xytospec_c(uux,mff,jac0)
	   call multiply(ng_bd,mff,mff2,mff4)
	   call xytospec_c(vuy,mff,jac0)

	   call spectoxy_c(mff3,gyx)
	   call spectoxy_c(mff4,gyy)

	   call multiply(ng_bd,mff,mff1,mff3)
	   call xytospec_c(uvx,mff,jac0)
	   call multiply(ng_bd,mff,mff2,mff4)
	   call xytospec_c(vvy,mff,jac0)

!-------------------------------------------------
! calcul du flux de masse u*h,v*h
!

	   call spectoxy_c(mff3,hk)
	   call multiply(ng_bd,mff,mff1,mff3)
	   call xytospec_c(uhk,mff,jac0)
	   call multiply(ng_bd,mff,mff2,mff3)
	   call xytospec_c(vhk,mff,jac0)
!
! calcul du forcing en u (vent + friction)
!
           mff = sqrt( mff1 * mff1 + mff2 * mff2 ) + 1.d-1
	   call spectoxy_c(mff4,wk)
	   mff4 = mff4 - mb0 * mff * mff1
	   call divide(ng_bd,mff,mff4,mff3)
	   call xytospec_c(wk,mff,jac0)
!
! calcul du forcing en v (vent + friction)
!
	   call spectoxy_c(mff4,zk)
	   mff4 = mff4 - mb0 * mff * mff2
	   call divide(ng_bd,mff,mff4,mff3)
	   call xytospec_c(zk,mff,jac0)


!------------------------------------------------
!       rotation des vitesse pour que la vitesse v soit 
!       normale a la face

	 amm0=amm_bd(:,:,tk)
	 amm1=amm_bdl(:,:,tk)

         do k=1,nm-nc-1
          dia4(k) = 1.d0/amm1(k,k)
         enddo
         do k=1,nm
          dia5(k) = 1.d0/amm0(k,k)
         enddo

!------------------------------------------------
! reduction
!
         call dtrsm4(nm,nm,amm0,uhk,dia5)
         call dtrsm4(nm,nm,amm0,vhk,dia5)
	 
         call reduc_gn(uhk,bu)
         call reduc_gn(vhk,bv)
         call dtrsm4(nm,nm-nc-1,amm1,bu,dia4)
         call dtrsm4(nm,nm-nc-1,amm1,bv,dia4)
         call reduc_inv_gn(uhk,bu)
         call reduc_inv_gn(vhk,bv)


!********* curved element along an open boundary *************

	 else if (fflagt(cell).eq.5) then

           jac0 = jac(:,tk)

	   call spectoxy_c(mff1,uk)
	   call spectoxy_c(mff2,vk)

	   call spectoxy_c(mff3,gxx)
	   call spectoxy_c(mff4,gxy)

	   call multiply(ng_bd,mff,mff1,mff3)
	   call xytospec_c(uux,mff,jac0)
	   call multiply(ng_bd,mff,mff2,mff4)
	   call xytospec_c(vuy,mff,jac0)

	   call spectoxy_c(mff3,gyx)
	   call spectoxy_c(mff4,gyy)

	   call multiply(ng_bd,mff,mff1,mff3)
	   call xytospec_c(uvx,mff,jac0)
	   call multiply(ng_bd,mff,mff2,mff4)
	   call xytospec_c(vvy,mff,jac0)

!-------------------------------------------------
! calcul du flux de masse u*h,v*h
!

	   call spectoxy_c(mff3,hk)
	   call multiply(ng_bd,mff,mff1,mff3)
	   call xytospec_c(uhk,mff,jac0)
	   call multiply(ng_bd,mff,mff2,mff3)
	   call xytospec_c(vhk,mff,jac0)
!
! calcul du forcing en u (vent + friction)
!
           mff = sqrt( mff1 * mff1 + mff2 * mff2 ) + 1.d-1
	   call spectoxy_c(mff4,wk)
	   mff4 = mff4 - mb0 * mff * mff1
	   call divide(ng_bd,mff,mff4,mff3)
	   call xytospec_c(wk,mff,jac0)
!
! calcul du forcing en v (vent + friction)
!
	   call spectoxy_c(mff4,zk)
	   mff4 = mff4 - mb0 * mff * mff2
	   call divide(ng_bd,mff,mff4,mff3)
	   call xytospec_c(zk,mff,jac0)


!------------------------------------------------
!       rotation des vitesse pour que la vitesse v soit 
!       normale a la face

	 amm0=amm_bd(:,:,tk)

         do k=1,nm
          dia5(k) = 1.d0/amm0(k,k)
         enddo


         call dtrsm4(nm,nm,amm0,uhk,dia5)
         call dtrsm4(nm,nm,amm0,vhk,dia5)


! ******* interior element ********

	else

	   call spectoxy(mff1,uk)
	   call spectoxy(mff2,vk)

	   call spectoxy(mff3,gxx)
	   call spectoxy(mff4,gxy)

	   call multiply(ng,mff,mff1,mff3)
	   call xytospec(uux,mff)
	   call multiply(ng,mff,mff2,mff4)
	   call xytospec(vuy,mff)

	   call spectoxy(mff3,gyx)
	   call spectoxy(mff4,gyy)

	   call multiply(ng,mff,mff1,mff3)
	   call xytospec(uvx,mff)
	   call multiply(ng,mff,mff2,mff4)
	   call xytospec(vvy,mff)

!-------------------------------------------------
! calcul du flux de masse u*h,v*h
!

	   call spectoxy(mff3,hk)
	   call multiply(ng,mff,mff1,mff3)
	   call xytospec(uhk,mff)
	   call multiply(ng,mff,mff2,mff3)
	   call xytospec(vhk,mff)
!
! calcul du forcing en u (vent + friction)
!
           mff = sqrt( mff1 * mff1 + mff2 * mff2 ) + 1.d-1
	   call spectoxy(mff4,wk)
	   mff4 = mff4 - mb0 * mff * mff1
	   call divide(ng,mff,mff4,mff3)
	   call xytospec(wk,mff,jac0)
!
! calcul du forcing en v (vent + friction)
!
	   call spectoxy(mff4,zk)
	   mff4 = mff4 - mb0 * mff * mff2
	   call divide(ng,mff,mff4,mff3)
	   call xytospec(zk,mff,jac0)


	   call dtrsm4(nm,nm,amm,uhk,dia1)
	   call dtrsm4(nm,nm,amm,vhk,dia1)

	endif


	   do j=1,nm
	      adu(j,cell) = uux(j) + vuy(j)
	      adv(j,cell) = uvx(j) + vvy(j)
	      uh(j,cell)  = uhk(j)
	      vh(j,cell)  = vhk(j)
              w1(j,cell)  = wk (j)
              z1(j,cell)  = zk (j)
	   enddo
	enddo



!------------------------------------------------
!--- u,v
!------------------------------------------------

	fxx = mu0 * qxx
	fxy = mu0 * qxy
	fyx = mu0 * qyx
	fyy = mu0 * qyy

	call cal_flu(fbd,fxx)
	call norm_flu(fbd,fbdxx,1)
	call cal_flu(fbd,fxy)
	call norm_flu(fbd,fbdxy,2)
	call cal_flu(fbd,fyx)
	call norm_flu(fbd,fbdyx,1)
	call cal_flu(fbd,fyy)
	call norm_flu(fbd,fbdyy,2)

	fxx = fxx - g0 * h0
	fyy = fyy - g0 * h0

        call cal_flu(fbd,h0)
	call openbc1(fbd,time,1.d0)

	call norm_flu(fbd,fbdx,1)
	call norm_flu(fbd,fbdy,2)
	
	fbdxx = fbdxx - g0 * fbdx
	fbdyy = fbdyy - g0 * fbdy
	

	

	do cell=1,ne

         ax    = vect(1,cell)
         bx    = vect(2,cell)
         ay    = vect(3,cell)
         by    = vect(4,cell)
         detai = vect(5,cell)
	 
!------------------------------------------------
! calcul des termes sources (nonlin + vent + friction)
!
	 gx = - adu(:,cell) + w1(:,cell)
	 gy = - adv(:,cell) + z1(:,cell)

!------------------------------------------------

	      gxx = fxx(:,cell)
	      gxy = fxy(:,cell)
	      gyx = fyx(:,cell)
	      gyy = fyy(:,cell)

       tk=pbdl_inv(cell)

!
!------------------------------------------------
! si la cellule est au bord, noslip
!------------------------------------------------
!
	 if (fflagt(cell).eq.1) then
	 
	 
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
	 
         call reduc_gn(gx,bu)
         call reduc_gn(gy,bv)
         call dtrsm4(nm,nm-nc-1,amm1,bu,dia4)
         call dtrsm4(nm,nm-nc-1,amm1,bv,dia4)
         call reduc_inv_gn(gx,bu)
         call reduc_inv_gn(gy,bv)


!******* open boudanry condition *******

	 else if (fflagt(cell).eq.5) then


	 amm0=amm_bd(:,:,tk)
	 amm2=agx_bd(:,:,tk)
	 amm3=agy_bd(:,:,tk)

         do k=1,nm-nc-1
          dia4(k) = 1.d0/amm1(k,k)
         enddo
         do k=1,nm
          dia5(k) = 1.d0/amm0(k,k)
         enddo

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
