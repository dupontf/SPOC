!
!******************************************************************
! nonlinear version
! free-slip
! viscosity
! beta plane
! varying bottom topography
!******************************************************************
!
subroutine cal_var(u0,v0,h0,u1,v1,h1,mu0,mb0,fc,w0,z0,g0,hb,time)
     
    use mesh
    use gauss
    use utilities
    use boundary_condition
    
      implicit none
      double precision u0(nm,ne),v0(nm,ne),h0(nm,ne), &
                       u1(nm,ne),v1(nm,ne),h1(nm,ne),&
                 w0(nm,ne),z0(nm,ne),hb(nm,ne),&
                 dia1(nm),dia2(nm),&
                 ak(nm),b(nm),bu(nm),&
                 uk(nm),vk(nm),hk(nm),wk(nm),zk(nm),&
                 uhk(nm),vhk(nm),&
                 uux(nm),vuy(nm),&
                 uvx(nm),vvy(nm),&
                 gx(nm),gy(nm),gn(nm),gt(nm),&
                 gx1(nm),gy1(nm),&
                 ax,by,ay,bx,deta,detai,coeff,dis,dnx,dny,mu,mv,&
                 mu0,mb0,fc(nm,ne),g0,time,pi,tday,hh0,&
                 gxx(nm),gxy(nm),gyx(nm),gyy(nm),&
                 gnx(nm),gtx(nm),gny(nm),gty(nm),&
                 gnn(nm),gtn(nm),gnt(nm),gtt(nm),&
                 qxx(nm,ne),qxy(nm,ne),&
                 qyx(nm,ne),qyy(nm,ne),&
                 fxx(nm,ne),fxy(nm,ne),&
                 fyx(nm,ne),fyy(nm,ne),&
                 fbdxx(ngf,nf),fbdxy(ngf,nf),&
                 fbdyx(ngf,nf),fbdyy(ngf,nf),&
                 fbdx(ngf,nf),fbdy(ngf,nf),&
                 fbd(ngf,nf),&
                 adu(nm,ne),adv(nm,ne),&
                 uh(nm,ne),vh(nm,ne),&
                 w1(nm,ne),z1(nm,ne),&
                 mff(ng), mff1(ng), mff2(ng), mff3(ng), mff4(ng)
      integer i,j,k,l,info,i1,i2,i3,ind,ns1
      integer cell,face


      
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
	do face=1,nf
	   do j=1,ngf
	   fbdxx(j,face) = fbd(j,face)*norm(1,face)
	   fbdxy(j,face) = fbd(j,face)*norm(2,face)
	   enddo
	enddo


	call cal_flu(fbd,v0)
	do face=1,nf
	   do j=1,ngf
	   fbdyx(j,face) = fbd(j,face)*norm(1,face)
	   fbdyy(j,face) = fbd(j,face)*norm(2,face)
	   enddo
	enddo


	do cell=1,ne

         ax    = vect(1,cell)
         bx    = vect(2,cell)
         ay    = vect(3,cell)
         by    = vect(4,cell)
         detai = vect(5,cell)

	   do j=1,nm
	      uk(j) = u0(j,cell)
	      vk(j) = v0(j,cell)
	   enddo

!------------------------------------------------
	   do j=1,nm
	      gx1(j) =0.d0
	      gy1(j) =0.d0
	   enddo
	   call cal_b(gx1,uk,agx,-1.d0)
	   call cal_b(gy1,uk,agy,-1.d0)

	   do j=1,nm
		gxx(j) = (  by * gx1(j) - ay * gy1(j)) * detai
		gxy(j) = (- bx * gx1(j) + ax * gy1(j)) * detai
	   enddo

           call tri_f(cell,fbdxx,gxx,detai)
  	   call tri_f(cell,fbdxy,gxy,detai)

!------------------------------------------------
	   do j=1,nm
	      gx1(j) =0.d0
	      gy1(j) =0.d0
	   enddo
	   call cal_b(gx1,vk,agx,-1.d0)
	   call cal_b(gy1,vk,agy,-1.d0)

	   do j=1,nm
		gyx(j) = (  by * gx1(j) - ay * gy1(j)) * detai
		gyy(j) = (- bx * gx1(j) + ax * gy1(j)) * detai
	   enddo

           call tri_f(cell,fbdyx,gyx,detai)
  	   call tri_f(cell,fbdyy,gyy,detai)

!
!------------------------------------------------
! si la cellule est au bord, rotation
!------------------------------------------------
!
	 if (bdl(cell)) then
!
!------------------------------------------------
!       rotation des vitesse pour que la vitesse v soit 
!       normale a la face
!
	  dnx = normt(1,cell)
	  dny = normt(2,cell)
!
! rotation du gradient
!
	   call rot_vec(nm,gnx,gtx,gxx,gyx,dnx,dny)
	   call rot_vec(nm,gny,gty,gxy,gyy,dnx,dny)
	   call rot_vec(nm,gnn,gnt,gnx,gny,dnx,dny)
	   call rot_vec(nm,gtn,gtt,gtx,gty,dnx,dny)

!------------------------------------------------
! reduction de gn
!
	   call reduc_gn(gtn,bu)

!------------------------------------------------

	   call dtrsm4(nm,nm,     amm,gnn,dia1)
	   call dtrsm4(nm,nm,     amm,gnt,dia1)
	   call dtrsm4(nm,nm-nc-1,amml,bu,dia2)
	   call dtrsm4(nm,nm,     amm,gtt,dia1)

!------------------------------------------------
! reduction inverse de gn
!
	   call reduc_inv_gn(gtn,bu)

!------------------------------------------------
! rotation inverse de la vitesse
!
	   call rot_vec(nm,gnx,gny,gnn,gnt,dnx,dny)
	   call rot_vec(nm,gtx,gty,gtn,gtt,dnx,dny)
	   call rot_vec(nm,gxx,gyx,gnx,gtx,dnx,dny)
	   call rot_vec(nm,gxy,gyy,gny,gty,dnx,dny)

!-------------------------------------------------
! cas normal sans rotation
!-------------------------------------------------

	   else

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
!--- copmute the nonlinear terms
!------------------------------------------------

      do cell=1,ne

	   do j=1,nm
	      uk(j) = u0(j,cell)
	      vk(j) = v0(j,cell)
	      hk(j) = h0(j,cell)+hb(j,cell)
	      wk(j) = w0(j,cell) - mb0 * uk(j)
	      zk(j) = z0(j,cell) - mb0 * vk(j)
	      gxx(j) = qxx(j,cell)
	      gxy(j) = qxy(j,cell) - fc(j,cell)
	      gyx(j) = qyx(j,cell) + fc(j,cell)
	      gyy(j) = qyy(j,cell)
	   enddo

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
	   call spectoxy(mff4,wk)
	   call divide(ng,mff,mff4,mff3)
	   call xytospec(wk,mff)
!
! calcul du forcing en v (vent + friction)
!
	   call spectoxy(mff4,zk)
	   call divide(ng,mff,mff4,mff3)
	   call xytospec(zk,mff)


	 if (fflagt(cell).eq.1) then

!------------------------------------------------
!       rotation des vitesse pour que la vitesse v soit 
!       normale a la face

	  dnx = normt(1,cell)
	  dny = normt(2,cell)

	   call rot_vec(nm,gn,gt,uhk,vhk,dnx,dny)
	   call reduc_gn(gn,bu)

	   call dtrsm4(nm,nm-nc-1,amml,bu,dia2)
	   call dtrsm4(nm,nm,     amm, gt,dia1)

	   call reduc_inv_gn(gn,bu)
	   call rot_vec(nm,uhk,vhk,gn,gt,dnx,dny)

	else

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
	do face=1,nf
	   do j=1,ngf
	   fbdxx(j,face) = fbd(j,face)*norm(1,face)
	   enddo
	enddo
	call cal_flu(fbd,fxy)
	do face=1,nf
	   do j=1,ngf
	   fbdxy(j,face) = fbd(j,face)*norm(2,face)
	   enddo
	enddo
	call cal_flu(fbd,fyx)
	do face=1,nf
	   do j=1,ngf
	   fbdyx(j,face) = fbd(j,face)*norm(1,face)
	   enddo
	enddo
	call cal_flu(fbd,fyy)
	do face=1,nf
	   do j=1,ngf
	   fbdyy(j,face) = fbd(j,face)*norm(2,face)
	   enddo
	enddo

	fxx = fxx - g0 * h0
	fyy = fyy - g0 * h0

        call cal_flu(fbd,h0)
	call openbc1(fbd,time,1.d0)

	do face=1,nf
	   do j=1,ngf
	   fbdxx(j,face) = fbdxx(j,face) - g0 * fbd(j,face)*norm(1,face)
	   fbdyy(j,face) = fbdyy(j,face) - g0 * fbd(j,face)*norm(2,face)
	   enddo
	enddo
	

	

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

!
!------------------------------------------------
! si la cellule est au bord, rotation
!------------------------------------------------
!
	 if (fflagt(cell).eq.1) then
!
!------------------------------------------------
!       rotation des vitesse pour que la vitesse v soit 
!       normale a la face
!
	  dnx = normt(1,cell)
	  dny = normt(2,cell)
!
! rotation du gradient et reduction
!
	   call rot_vec(nm,gn,gt,gx,gy,dnx,dny)
	   call reduc_gn(gn,bu)
!
!------------------------------------------------
!
	   call dtrsm4(nm,nm-nc-1,amml,bu,dia2)
	   call dtrsm4(nm,nm,     amm, gt,dia1)
!
!------------------------------------------------
! reduction et rotation inverse
!
	   call reduc_inv_gn(gn,bu)
	   call rot_vec(nm,gx,gy,gn,gt,dnx,dny)
!
!-------------------------------------------------
! cas normal sans rotation
!-------------------------------------------------
!
	   else

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
	do face=1,nf
	   do j=1,ngf
	   fbdx(j,face) = fbd(j,face)*norm(1,face)
	   enddo
	enddo
	call cal_flu(fbd,vh)
	do face=1,nf
	   do j=1,ngf
	   fbdy(j,face) = fbd(j,face)*norm(2,face)
	   enddo
	enddo
	call massbc(fbdx,fbdy)


	do cell=1,ne

         ax    = vect(1,cell)
         bx    = vect(2,cell)
         ay    = vect(3,cell)
         by    = vect(4,cell)
         detai = vect(5,cell)

           do j=1,nm
              b(j) =0.d0
           enddo

	   do j=1,nm
	      ak(j) = uh(j,cell)
	   enddo
	   coeff = by*detai
	   call cal_b(b,ak,agx,coeff)
	   coeff =-ay*detai
	   call cal_b(b,ak,agy,coeff)

	   do j=1,nm
	      ak(j) = vh(j,cell)
	   enddo
	   coeff =-bx*detai
	   call cal_b(b,ak,agx,coeff)
	   coeff = ax*detai
	   call cal_b(b,ak,agy,coeff)


	   call tri_f(cell,fbdx,b,-detai)

	   call dtrsm4(nm,nm,amm,b,dia1)

	   do j=1,nm
		h1(j,cell) = b(j)
	   enddo

	enddo


	return
	end
!
!*****************************************************************
