c
c******************************************************************
c non linear dynamics
c viscosity: free-slip boundary condition
c wind forcing
c linear bottom friction (m/s)
c******************************************************************
c
      subroutine cal_var(nf,ne,u0,v0,h0,u1,v1,h1,
     &            tseg,pseg,ssit,itm,its,norm,nbdl,pbdl,nint,pint,
     &            bdl,cip,vect,normt,
     &            mu0,mb0,fc,beta,w0,z0,g0,hb,time)
      implicit none
      include '../include/maille.inc'
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision xi(ng_max),yi(ng_max),wi(ng_max)
      common /gauss/ wi,xi,yi
      double precision wif(ngf_max),tif(ngf_max)
      common /gaussf/ wif,tif
      integer ne,nf
      logical bdl(*)
      double precision u0(nnmod,*),v0(nnmod,*),h0(nnmod,*),
     &                 u1(nnmod,*),v1(nnmod,*),h1(nnmod,*),
     &            w0(nnmod,*),z0(nnmod,*),hb(nnmod,*),
     &            agx(nnmod,nnmod),agy(nnmod,nnmod),
     &            amm(nnmod,nnmod),amml(nnmod,nnmod),ammo(nnmod,nnmod),
     &            ammx(nnmod,nnmod),ammy(nnmod,nnmod),
     &            dia1(nnmod),dia2(nnmod),
     &            norm(2,*),normt(2,*),vect(5,*),
     &            ak(nnmod),b(nnmod),bu(nnmod),
     &            uk(nnmod),vk(nnmod),hk(nnmod),wk(nnmod),zk(nnmod),
     &            uhk(nnmod),vhk(nnmod),
     &            uux(nnmod),vuy(nnmod),
     &            uvx(nnmod),vvy(nnmod),
     &            gx(nnmod),gy(nnmod),gn(nnmod),gt(nnmod),
     &            gx1(nnmod),gy1(nnmod),
     &            ax,by,ay,bx,deta,detai,coeff,dis,dnx,dny,mu,mv,
     &            mu0,mb0,fc(nnmod,*),beta,g0,time,pi,tday,
     &            gxx(nnmod),gxy(nnmod),gyx(nnmod),gyy(nnmod),
     &            gnx(nnmod),gtx(nnmod),gny(nnmod),gty(nnmod),
     &            gnn(nnmod),gtn(nnmod),gnt(nnmod),gtt(nnmod),
     &            qxx(nnmod,nedim),qxy(nnmod,nedim),
     &            qyx(nnmod,nedim),qyy(nnmod,nedim),
     &            fxx(nnmod,nedim),fxy(nnmod,nedim),
     &            fyx(nnmod,nedim),fyy(nnmod,nedim),
     &            fbdxx(ngf_max,nfdim),fbdxy(ngf_max,nfdim),
     &            fbdyx(ngf_max,nfdim),fbdyy(ngf_max,nfdim),
     &            fbdx(ngf_max,nfdim),fbdy(ngf_max,nfdim),
     7            fbd(ngf_max,nfdim),
     &            adu(nnmod,nedim),adv(nnmod,nedim),
     &            uh(nnmod,nedim),vh(nnmod,nedim),
     &            w1(nnmod,nedim),z1(nnmod,nedim),
     &            mff(ng_max),mff1(ng_max),mff2(ng_max),
     &                       mff3(ng_max),mff4(ng_max),
     &            pico(ngf_max,nnmod,3),pi2(ng_max,nnmod)
      integer i,j,k,l,info,i1,i2,i3,ind,ns1
      integer tseg(2,*),pseg(2,*),ssit(2,*),itm(3,*),its(3,*)
      integer nbdl,pbdl(*),nint,pint(*)
      integer cip(2,0:nnmod)
      common /gauss2/ pi2
      common /gauss3/ pico
      common /matrices/ amm,amml,agx,agy,ammx,ammy,ammo

c
c fonction pour le vent
c
      pi=4.d0*atan(1.d0)
      tday=time/86400.d0
c
	do i=1,nm
         dia1(i) = 1.d0/amm(i,i)
      enddo
	do i=1,nm-nc-1
         dia2(i) = 1.d0/amml(i,i)
      enddo
c
c
c------------------------------------
c qxx-qxy--qyx-qyy
c------------------------------------
c
c
c calcul du la fonction le long de la frontiere
c------------------------------------
c
	call cal_flu(nf,tseg,fbd,u0,ssit)
	do i=1,nf
	   do j=1,ngf
	   fbdxx(j,i) = fbd(j,i)*norm(1,i)
	   fbdxy(j,i) = fbd(j,i)*norm(2,i)
	   enddo
	enddo
	call cal_flu(nf,tseg,fbd,v0,ssit)
	do i=1,nf
	   do j=1,ngf
	   fbdyx(j,i) = fbd(j,i)*norm(1,i)
	   fbdyy(j,i) = fbd(j,i)*norm(2,i)
	   enddo
	enddo
c
c------------------------------------
c solution pour chaque triangle
c------------------------------------
c
c
c$omp parallel share(u0,v0,qxx,qxy,qyx,qyy,fbdxx,fbdxy,fbdyx,fbdyy,
c$&                  nc,nm,normt,itm,its,bdl,vect,dia1,dia2,amm,amml),
c$&            local(i,j,ax,bx,ay,by,detai,
c$&                  uk,vk,gx1,gy1,gxx,gxy,gyx,gyy,
c$&                  gnx,gny,gtx,gty,gnn,gnt,gtn,gtt,bu,
c$&                  coeff,dnx,dny)
c$omp do schedule(runtime)
	do i=1,ne
c
	   do j=1,nm
	      uk(j) = u0(j,i)
	      vk(j) = v0(j,i)
	   enddo
c
         ax = vect(1,i)
         bx = vect(2,i)
         ay = vect(3,i)
         by = vect(4,i)
         detai = vect(5,i)
c
c------------------------------------------------
	   do j=1,nm
	      gx1(j) =0.d0
	      gy1(j) =0.d0
	   enddo
	   call cal_b(gx1,uk,agx,-1.d0)
	   call cal_b(gy1,uk,agy,-1.d0)
c
	   do j=1,nm
		gxx(j) = (  by * gx1(j) - ay * gy1(j)) * detai
		gxy(j) = (- bx * gx1(j) + ax * gy1(j)) * detai
	   enddo
c
         call tri_f(i,itm,its,fbdxx,gxx,detai)
  	   call tri_f(i,itm,its,fbdxy,gxy,detai)
c
c------------------------------------------------
	   do j=1,nm
	      gx1(j) =0.d0
	      gy1(j) =0.d0
	   enddo
	   call cal_b(gx1,vk,agx,-1.d0)
	   call cal_b(gy1,vk,agy,-1.d0)
c
	   do j=1,nm
		gyx(j) = (  by * gx1(j) - ay * gy1(j)) * detai
		gyy(j) = (- bx * gx1(j) + ax * gy1(j)) * detai
	   enddo
c
         call tri_f(i,itm,its,fbdyx,gyx,detai)
  	   call tri_f(i,itm,its,fbdyy,gyy,detai)
c
c
c------------------------------------------------
c si la cellule est au bord, rotation
c------------------------------------------------
c
	 if (bdl(i)) then
c
c------------------------------------------------
c       rotation des vitesse pour que la vitesse v soit 
c       normale a la face
c
	  dnx = normt(1,i)
	  dny = normt(2,i)
c
c rotation du gradient
c
	   call rot_vec(nm,gnx,gtx,gxx,gyx,dnx,dny)
	   call rot_vec(nm,gny,gty,gxy,gyy,dnx,dny)
	   call rot_vec(nm,gnn,gnt,gnx,gny,dnx,dny)
	   call rot_vec(nm,gtn,gtt,gtx,gty,dnx,dny)
c
c------------------------------------------------
c reduction de gn
c
	   call reduc_gn(nc,gtn,bu,cip)
c
c------------------------------------------------
c
	   call dtrsm4(nm,     amm,gnn,dia1)
	   call dtrsm4(nm,     amm,gnt,dia1)
	   call dtrsm4(nm-nc-1,amml,bu,dia2)
	   call dtrsm4(nm,     amm,gtt,dia1)
c
c------------------------------------------------
c reduction inverse de gn
c
	   call reduc_inv_gn(nc,gtn,bu,cip)
c
c------------------------------------------------
c rotation inverse de la vitesse
c
	   call rot_vec(nm,gnx,gny,gnn,gnt,dnx,dny)
	   call rot_vec(nm,gtx,gty,gtn,gtt,dnx,dny)
	   call rot_vec(nm,gxx,gyx,gnx,gtx,dnx,dny)
	   call rot_vec(nm,gxy,gyy,gny,gty,dnx,dny)
c
c-------------------------------------------------
c cas normal sans rotation
c-------------------------------------------------
c
	   else
c
	   call dtrsm4(nm,amm,gxx,dia1)
	   call dtrsm4(nm,amm,gxy,dia1)
	   call dtrsm4(nm,amm,gyx,dia1)
	   call dtrsm4(nm,amm,gyy,dia1)
c
	   endif
c
	   do j=1,nm
		qxx(j,i) = gxx(j)
		qxy(j,i) = gxy(j)
		qyx(j,i) = gyx(j)
		qyy(j,i) = gyy(j)
	   enddo
c
	enddo
c$omp end do
c$omp end parallel
c
c
c
c-------------------------------------------------
c calcul des termes non-lineaires
c
c$omp parallel share(u0,v0,h0,w0,hb,qxx,qxy,qyx,qyy,w1,z1,uh,vh,
c$&                  nc,nm,normt,itm,mb0,bdl,dia1,dia2,amm,amml),
c$&            local(i,j,
c$&                  uk,vk,hk,wk,zk,gxx,gxy,gyx,gyy,uhk,vhk,gn,gt,bu,
c$&                  uux,uvx,vuy,vvy,
c$&                  mff,mff1,mff2,mff3,mff4,
c$&                  dnx,dny)
c$omp do schedule(runtime)
      do i=1,ne
c
	   do j=1,nm
	      uk(j) = u0(j,i)
	      vk(j) = v0(j,i)
	      hk(j) = h0(j,i)+hb(j,i)
	      wk(j) = w0(j,i) - mb0 * uk(j)
	      zk(j) = z0(j,i) - mb0 * vk(j)
	      gxx(j) = qxx(j,i)
	      gxy(j) = qxy(j,i) - fc(j,i)
	      gyx(j) = qyx(j,i) + fc(j,i)
	      gyy(j) = qyy(j,i)
	   enddo
c
	   call spectoxy(mff1,uk)
	   call spectoxy(mff2,vk)
c
	   call spectoxy(mff3,gxx)
	   call spectoxy(mff4,gxy)
c
	   call multiply(ng,mff,mff1,mff3)
	   call xytospec(uux,mff)
	   call multiply(ng,mff,mff2,mff4)
	   call xytospec(vuy,mff)
c
	   call spectoxy(mff3,gyx)
	   call spectoxy(mff4,gyy)
c
	   call multiply(ng,mff,mff1,mff3)
	   call xytospec(uvx,mff)
	   call multiply(ng,mff,mff2,mff4)
	   call xytospec(vvy,mff)
c
c-------------------------------------------------
c calcul du flux de masse u*h,v*h
c
c
	   call spectoxy(mff3,hk)
	   call multiply(ng,mff,mff1,mff3)
	   call xytospec(uhk,mff)
	   call multiply(ng,mff,mff2,mff3)
	   call xytospec(vhk,mff)
c
c calcul du forcing en u (vent + friction)
c
	   call spectoxy(mff4,wk)
	   call divide(ng,mff,mff4,mff3)
	   call xytospec(wk,mff)
c
c calcul du forcing en v (vent + friction)
c
	   call spectoxy(mff4,zk)
	   call divide(ng,mff,mff4,mff3)
	   call xytospec(zk,mff)
c
	 if (bdl(i)) then
c
c------------------------------------------------
c       rotation des vitesse pour que la vitesse v soit 
c       normale a la face
c
	  dnx = normt(1,i)
	  dny = normt(2,i)
c
	   call rot_vec(nm,gn,gt,uhk,vhk,dnx,dny)
	   call reduc_gn(nc,gn,bu,cip)
c
	   call dtrsm4(nm-nc-1,amml,bu,dia2)
	   call dtrsm4(nm,     amm, gt,dia1)
c
	   call reduc_inv_gn(nc,gn,bu,cip)
	   call rot_vec(nm,uhk,vhk,gn,gt,dnx,dny)
c
	else
c
	   call dtrsm4(nm,amm,uhk,dia1)
	   call dtrsm4(nm,amm,vhk,dia1)
c
	endif
c
c
	   do j=1,nm
	      adu(j,i) = uux(j) + vuy(j)
	      adv(j,i) = uvx(j) + vvy(j)
	      uh(j,i)  = uhk(j)
	      vh(j,i)  = vhk(j)
		w1(j,i)  = wk (j)
		z1(j,i)  = zk (j)
	   enddo
	enddo
c$omp end do
c$omp end parallel
c
c
c------------------------------------------------
c--- u,v
c------------------------------------------------
c
	do i=1,ne
	   do j=1,nm
	      fxx(j,i) = mu0 * qxx(j,i) - g0 * h0(j,i)
	      fxy(j,i) = mu0 * qxy(j,i)
	      fyx(j,i) = mu0 * qyx(j,i)
	      fyy(j,i) = mu0 * qyy(j,i) - g0 * h0(j,i)
	   enddo
	enddo
c
	call cal_flu(nf,tseg,fbd,fxx,ssit)
	do i=1,nf
	   do j=1,ngf
	   fbdxx(j,i) = fbd(j,i)*norm(1,i)
	   enddo
	enddo
	call cal_flu(nf,tseg,fbd,fxy,ssit)
	do i=1,nf
	   do j=1,ngf
	   fbdxy(j,i) = fbd(j,i)*norm(2,i)
	   enddo
	enddo
	call cal_flu(nf,tseg,fbd,fyx,ssit)
	do i=1,nf
	   do j=1,ngf
	   fbdyx(j,i) = fbd(j,i)*norm(1,i)
	   enddo
	enddo
	call cal_flu(nf,tseg,fbd,fyy,ssit)
	do i=1,nf
	   do j=1,ngf
	   fbdyy(j,i) = fbd(j,i)*norm(2,i)
	   enddo
	enddo
c
c------------------------------------------------
c$omp parallel share(u0,v0,fxx,fxy,fyx,fyy,adu,adv,u1,v1,w1,z1,fc,
c$&                  nc,nm,normt,itm,its,bdl,vect,dia1,dia2,amm,amml),
c$&            local(i,j,ax,bx,ay,by,detai,
c$&                  uk,vk,wk,zk,gx,gy,gxx,gxy,gyx,gyy,gn,gt,bu,
c$&                  coeff,dnx,dny)
c$omp do schedule(runtime)
	do i=1,ne
c
         ax = vect(1,i)
         bx = vect(2,i)
         ay = vect(3,i)
         by = vect(4,i)
         detai = vect(5,i)
c
c------------------------------------------------
c calcul des termes sources (nonlin + vent + friction)
c
	   do j=1,nm
	      gx(j) =  - adu(j,i) + w1(j,i)
	      gy(j) =  - adv(j,i) + z1(j,i)
	   enddo
c
c
c------------------------------------------------
	   do j=1,nm
	      gxx(j) = fxx(j,i)
	      gxy(j) = fxy(j,i)
	      gyx(j) = fyx(j,i)
	      gyy(j) = fyy(j,i)
	   enddo
c
c------------------------------------------------
c calcul de tenseur pour u
c
	   call cal_b(gx,gxx,agx,-by*detai)
	   call cal_b(gx,gxx,agy, ay*detai)
c
	   call cal_b(gx,gxy,agx, bx*detai)
	   call cal_b(gx,gxy,agy,-ax*detai)
c
	   call tri_f(i,itm,its,fbdxx,gx,detai)
	   call tri_f(i,itm,its,fbdxy,gx,detai)
c
c------------------------------------------------
c calcul de tenseur pour v
c
	   call cal_b(gy,gyx,agx,-by*detai)
	   call cal_b(gy,gyx,agy,ay*detai)
c
	   call cal_b(gy,gyy,agx,bx*detai)
	   call cal_b(gy,gyy,agy,-ax*detai)
c
	   call tri_f(i,itm,its,fbdyx,gy,detai)
	   call tri_f(i,itm,its,fbdyy,gy,detai)
c
c
c------------------------------------------------
c si la cellule est au bord, rotation
c------------------------------------------------
c
	 if (bdl(i)) then
c
c------------------------------------------------
c       rotation des vitesse pour que la vitesse v soit 
c       normale a la face
c
	  dnx = normt(1,i)
	  dny = normt(2,i)
c
c rotation du gradient et reduction
c
	   call rot_vec(nm,gn,gt,gx,gy,dnx,dny)
	   call reduc_gn(nc,gn,bu,cip)
c
c------------------------------------------------
c
	   call dtrsm4(nm-nc-1,amml,bu,dia2)
	   call dtrsm4(nm,     amm, gt,dia1)
c
c------------------------------------------------
c reduction et rotation inverse
c
	   call reduc_inv_gn(nc,gn,bu,cip)
	   call rot_vec(nm,gx,gy,gn,gt,dnx,dny)
c
c
c-------------------------------------------------
c cas normal sans rotation
c-------------------------------------------------
c
	   else
c
	     call dtrsm4(nm,amm,gx,dia1)
	     call dtrsm4(nm,amm,gy,dia1)
c
	   endif
c
	   do j=1,nm
		u1(j,i) = gx(j)
		v1(j,i) = gy(j)
	   enddo
c
	enddo
c$omp end do
c$omp end parallel
c
c
c-------------------------------------------------
c--- h
c-------------------------------------------------
c
	call cal_flu(nf,tseg,fbd,uh,ssit)
	do i=1,nf
	   do j=1,ngf
	   fbdx(j,i) = fbd(j,i)*norm(1,i)
	   enddo
	enddo
	call cal_flu(nf,tseg,fbd,vh,ssit)
	do i=1,nf
	   do j=1,ngf
	   fbdy(j,i) = fbd(j,i)*norm(2,i)
	   enddo
	enddo
c
c-------------------------------------------------
c
c$omp parallel share(h1,uh,vh,fbdx,fbdy,
c$&                  nc,nm,itm,its,vect,dia1,amm),
c$&            local(i,j,ax,bx,ay,by,detai,b,ak,coeff)
c$omp do schedule(runtime)
	do i=1,ne
c
c       write(*,*) 'tk',i
	   do j=1,nm
	      b(j) =0.d0
	   enddo
         ax = vect(1,i)
         bx = vect(2,i)
         ay = vect(3,i)
         by = vect(4,i)
         detai = vect(5,i)
c
	   do j=1,nm
	      ak(j) = uh(j,i)
	   enddo
	   coeff = by*detai
	   call cal_b(b,ak,agx,coeff)
	   coeff =-ay*detai
	   call cal_b(b,ak,agy,coeff)
c
	   do j=1,nm
	      ak(j) = vh(j,i)
	   enddo
	   coeff =-bx*detai
	   call cal_b(b,ak,agx,coeff)
	   coeff = ax*detai
	   call cal_b(b,ak,agy,coeff)
c
c
	   call tri_f(i,itm,its,fbdx,b,-detai)
	   call tri_f(i,itm,its,fbdy,b,-detai)
c
	   call dtrsm4(nm,amm,b,dia1)
c
	   do j=1,nm
		h1(j,i) = b(j)
	   enddo
c
	enddo
c$omp end do
c$omp end parallel
c
	return
	end
c
c*****************************************************************
