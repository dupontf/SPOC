	program interpol
c******************************************************************
c interpolation entre spectral finite elements et FDM C-grid
c******************************************************************
	implicit none
c
C
C 2-D array declaration
C
      include '../include/maille.inc'
      include '../include/sp.inc'
      integer nn,ne,nf,nebd
      INTEGER IN(3,NEDIM),tseg(2,nfdim),pseg(2,nfdim),
     &itm(3,nedim),nmit(nndim),tmit(NFTRDIM*NNDIM),
     & imit(NFTRDIM*NNDIM),ssit(2,nfdim),sbd(nbedim),
     & its(3,nedim)
      double precision norm(2,nfdim),
     &                 scx,scy,xgr(nndim),ygr(nndim)
      logical nobd(nndim)
      logical bdl(nedim)
c
c modes spectraux
c
      double precision 
     &       amm(nnmod,nnmod),
     &       ammo(nnmod,nnmod),
     &       amml(nnmod,nnmod),
     &       ammx(nnmod,nnmod),
     &       ammy(nnmod,nnmod),
     &       agx(nnmod,nnmod),
     &       agy(nnmod,nnmod),
     &       u0(nnmod,nedim),v0(nnmod,nedim),h0(nnmod,nedim),
     &       u1(nnmod,nedim),v1(nnmod,nedim),h1(nnmod,nedim)
      double precision f2,dfx,dfy,fx,fy
      integer i,j,k,l,i1,tk,tk1,n1,n2,k1
      integer nc,nm,ip
      external f2,dfx,dfy,ip,fx,fy
      integer nite,ite,info,per,nper,per2
      double precision 
     &       dt,mu0,mb0,f0,beta,g0,hh0,time,ke,pe,mass,ke0,pe0,
     &       lx,ly,sum,dx,dy,c0,pi,dis
	character*20 vieux
c
c------------------------------------
c fdm variables
	integer nd,nx,ny
	double precision ud,vd,hd,x0,y0,x1,y1,uref,href,uc,hc
	double precision error1,error2,error3,error4,error5,error6,
     & error7
	double precision lab(3)
	double precision fsol
	external fsol
c
c------------------------------------
c construction de la base de polynome
c------------------------------------
c
c
	open(2,file='oc.inp',status='old')
      read(2,*) nc
      read(2,*) lx,ly
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,*) scx,scy
      read(2,*) g0
	close(2)
c
	nm = (nc+1)*(nc+2)/2
	write(*,*) 'rang matrice',nm
	write(*,*) scx,scy
c
c------------------------------------
c parametres definissant la
c solution analytique
c------------------------------------
c
      write(*,60) g0
      c0=sqrt(10.d0)
      write(*,60) c0
      pi=4.d0*atan(1.d0)
      uc=sin(0.2d0*pi)*g0/c0
      hc=cos(0.2d0*pi)*1.d0
      write(*,60) uc,hc
c
c------------------------------------
c lecture du maillage
c------------------------------------
c
      call maille(nn,ne,nf,nebd,in,itm,nmit,tmit,imit,
     1                  pseg,tseg,sbd,xgr,ygr)
      write(*,*) nn,ne,nf
      do i=1,nn
	  xgr(i) = xgr(i) * scx
	  ygr(i) = ygr(i) * scy
      enddo
c
c------------------------------------
c lecture du fichier
c------------------------------------
c
	write(*,*) 'entrer fichier'
	read(*,10) vieux
10	format(a20)
	open(1,file=vieux,form='unformatted',status='old')
	read(1) time
	read(1) ((u0(j,i),j=1,nm),i=1,ne)
	read(1) ((v0(j,i),j=1,nm),i=1,ne)
	read(1) ((h0(j,i),j=1,nm),i=1,ne)
	close(1)
c
C-----------------------------------------------------------------------
c interpolation sur maillage FDM A-grid (centre)
c
c pour u
c
      error1 = 0.d0
      error2 = 0.d0
      error3 = 0.d0
      error4 = 0.d0
      error5 = 0.d0
      error6 = 0.d0
      error7 = 0.d0

      write(*,*) 'nx,ny'
      read(*,*) nx,ny
      dx = lx /dfloat(nx-1)
      dy = ly /dfloat(ny-1)
      tk=1
      call interlin(in,itm,tseg,pseg,xgr,ygr,tk,1.d6,1.d6,lab)
      tk1 = tk
      write(*,*) tk
      sum = lab(1) + lab(2) + lab(3)
      write(*,*) sum

      k=0
      k1=0
	do j=1,ny-1
	   do i=1,nx-1
	      x0 = 0.5d0 * lx * dfloat(2*i-nx  )/dfloat(nx-1)
	      y0 = 0.5d0 * ly * dfloat(2*j-ny  )/dfloat(ny-1)
	      call interlin(in,itm,tseg,pseg,xgr,ygr,tk,x0,y0,lab)
            sum = lab(1) + lab(2) + lab(3)
            if (sum.lt.0.5d0) then
             k=k+1
             tk=tk1
	    else
             k1=k1+1
		 x1 = 2.d0 * lab(2) - 1.d0
		 y1 = 2.d0 * lab(3) - 1.d0
	       ud=fsol(nc,x1,y1,u0,tk)
	       vd=fsol(nc,x1,y1,v0,tk)
	       hd=fsol(nc,x1,y1,h0,tk)
               uref= uc * sin ( 2.d0 * pi * x0/lx )
               href= hc * cos ( 2.d0 * pi * x0/lx )
c	       if (j.eq.ny/2) write(*,*) i,sin ( 2.d0 * pi * x0/lx ),uref
	       error1 = error1 + abs(uref-ud)
	       error2 = error2 + abs(vd)
	       error3 = error3 + abs(href-hd)
	       error4 = error4 + (uref-ud)**2
	       error5 = error5 + vd**2
	       error6 = error6 + (href-hd)**2
	      endif
	   enddo
	enddo
      write(*,*) k,' points etaient en dehors de la grille'
      write(*,*) k1,' points etaient a l''interieur de la grille'
c
C-----------------------------------------------------------------------
c ecriture du fichier apres interpolation
c
c
c calcul de dx moyen
c
	sum=0.d0
	do i=1,nf
	n1=pseg(1,i)
	n2=pseg(2,i)
	dis = sqrt ((xgr(n1)-xgr(n2))**2+(ygr(n1)-ygr(n2))**2)
	sum=sum+dis
	enddo
	sum=sum/dfloat(nf)
	write(*,21) 'dx moy',sum
	write(*,23) 
	write(*,22) 'Abs_norm',error1/dfloat(k1),error2/dfloat(k1),
     &  error3/dfloat(k1)
	write(*,22) 'RMS_norm',sqrt(error4/dfloat(k1)),
     &  sqrt(error5/dfloat(k1)),
     &  sqrt(error6/dfloat(k1))
        error7 = error4 + error5
	write(*,22) 'RMS_norm_v',sqrt(error7/dfloat(k1))
c
c
20	format(1x,20(f10.6,1x))
21	format(1x,a10,1p,e13.6)
22	format(1x,a10,1p,3(e10.3,3x))
23	format(15x,'u',10x,'v',10x,'h')
60	format(3(e23.16,2x))
c
	end
c*************************************************************
c functions
c*************************************************************
c
	function fsol(nc,x,y,b,i)
	implicit none
	include '../include/sp.inc'
	double precision b(nnmod,*)
	double precision fsol,x,y,po
	integer kc1,kc2,nc,i,j
	external po
	j=0
	fsol=0.d0
	  do kc1=0,nc
	  do kc2=0,nc-kc1
	    j=j+1
	    fsol = fsol + b(j,i) * po(kc1,x) * po(kc2,y)
	  enddo
	  enddo
	return
	end
c
