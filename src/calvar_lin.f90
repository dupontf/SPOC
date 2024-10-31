!
!******************************************************************
! linear version
! no viscosity, no beta plane
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
	
	
	do face=1,nf
	   do j=1,ngf
	   fbdxx(j,face) = fbd(j,face)*norm(1,face)
	   fbdxy(j,face) = 0.d0
	   fbdyx(j,face) = 0.d0
	   fbdyy(j,face) = fbd(j,face)*norm(2,face)
	   enddo
	enddo


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
         uh = hh0 * u0
         vh = hh0 * v0

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
!	call massbc(fbdx,fbdy)


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
	   call tri_f(cell,fbdy,b,-detai)

	   call dtrsm4(nm,nm,amm,b,dia1)

	   do j=1,nm
		h1(j,cell) = b(j)
	   enddo

	enddo


	return
	end
!
!*****************************************************************
