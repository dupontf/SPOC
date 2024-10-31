c
c******************************************************************
c linear version
c no viscosity, no beta plane
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
     &            mu0,mb0,fc(*),beta,g0,time,pi,tday,hh0,
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
c
      hh0=hb(1,1)
      do i=1,nm
         dia1(i) = 1.d0/amm(i,i)
         dia2(i) = 1.d0/amml(i,i)
      enddo
c------------------------------------------------
c--- u,v
c------------------------------------------------
c
	do i=1,ne
	   do j=1,nm
	      fxx(j,i) = - g0 * h0(j,i)
	      fxy(j,i) = 0.d0
	      fyx(j,i) = 0.d0
	      fyy(j,i) = - g0 * h0(j,i)
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
c$&                  nc,nm,normt,dia1,dia2,amm,amml),
c$&            local(i,j,ax,bx,ay,by,detai,
c$&                  gx,gy,gxx,gxy,gyx,gyy,gn,gt,bu,
c$&                  coeff,dnx,dny)
c$omp do schedule(runtime)
	do i=1,ne
c
c       write(*,*) 'tk',i
c
         ax = vect(1,i)
         bx = vect(2,i)
         ay = vect(3,i)
         by = vect(4,i)
         detai = vect(5,i)
c
c
	   do j=1,nm
	      gx(j) =  0.d0
	      gy(j) =  0.d0
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
	   coeff = by*detai
	   call cal_b(gy,gyx,agx,-coeff)
	   coeff =-ay*detai
	   call cal_b(gy,gyx,agy,-coeff)
c
	   coeff =-bx*detai
	   call cal_b(gy,gyy,agx,-coeff)
	   coeff = ax*detai
	   call cal_b(gy,gyy,agy,-coeff)
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
	call cal_flu(nf,tseg,fbd,u0,ssit)
	do i=1,nf
	   do j=1,ngf
	   fbdx(j,i) = fbd(j,i)*norm(1,i)*hh0
	   enddo
	enddo
	call cal_flu(nf,tseg,fbd,v0,ssit)
	do i=1,nf
	   do j=1,ngf
	   fbdy(j,i) = fbd(j,i)*norm(2,i)*hh0
	   enddo
	enddo
c
c-------------------------------------------------
c
c$omp parallel share(h1,uh,vh,fbdx,fbdy,
c$&                  nc,nm),
c$&            local(i,j,i1,i2,i3,ax,bx,ay,by,deta,detai,b,ak,
c$&                  coeff,info)
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
	      ak(j) = u0(j,i)*hh0
	   enddo
	   coeff = by*detai
	   call cal_b(b,ak,agx,coeff)
	   coeff =-ay*detai
	   call cal_b(b,ak,agy,coeff)
c
	   do j=1,nm
	      ak(j) = v0(j,i)*hh0
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
