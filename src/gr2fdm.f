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
      integer i,j,k,l,i1,tk,tk1,n1,n2
      integer nc,nm,ip
      external f2,dfx,dfy,ip,fx,fy
      integer nite,ite,info,per,nper,per2
      double precision 
     &       dt,mu0,mb0,f0,beta,g0,hh0,time,ke,pe,mass,ke0,pe0,
     &       lx,ly,sum,dx,dy,c0,uc,pi,dis
	character*20 vieux
c
c------------------------------------
c fdm variables
	integer nd,nx,ny
	parameter(nd=1000)
	double precision ud(nd,nd),vd(nd,nd),hd(nd,nd),x0,y0,x1,y1,uref(nd)
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
	do j=1,ny-1
	   do i=1,nx-1
	      x0 = 0.5d0 * lx * dfloat(2*i-nx  )/dfloat(nx-1)
	      y0 = 0.5d0 * ly * dfloat(2*j-ny  )/dfloat(ny-1)
	      call interlin(in,itm,tseg,pseg,xgr,ygr,tk,x0,y0,lab)
            sum = lab(1) + lab(2) + lab(3)
            if (sum.lt.0.5d0) then
             ud(i,j)=-999.9
             vd(i,j)=-999.9
             hd(i,j)=-999.9
             k=k+1
             tk=tk1
	      else
		 x1 = 2.d0 * lab(2) - 1.d0
		 y1 = 2.d0 * lab(3) - 1.d0
	       ud(i,j)=fsol(nc,x1,y1,u0,tk)
	       vd(i,j)=fsol(nc,x1,y1,v0,tk)
	       hd(i,j)=fsol(nc,x1,y1,h0,tk)
	      endif
	   enddo
	enddo
      write(*,*) k,' points etaient en dehors de la grille'
c
C-----------------------------------------------------------------------
c ecriture du fichier apres interpolation
c
	write(*,*) g0
      c0=sqrt(10.d0)
	write(*,*) c0
      pi=4.d0*atan(1.d0)
      uc=sin(0.2d0*pi)*g0/c0
      do i=1,nx-1
         x0 = 0.5d0 * dfloat(2*i-nx  )/dfloat(nx-1)
         uref(i)= uc * sin ( 2.d0 * pi * x0 )
      enddo
c
c integral de abs(u-uref)
c
      sum=0.d0
      do j=1,ny-1
         do i=1,nx-1
            sum=sum+abs(ud(i,j)-uref(i))
         enddo
      enddo
      sum=sum*dx*dy
      write(*,21) 'abs(u)',sum
c
c integral de abs(v)
c
      sum=0.d0
      do j=1,ny-1
         do i=1,nx-1
            sum=sum+abs(vd(i,j))
         enddo
      enddo
      sum=sum*dx*dy
      write(*,21) 'abs(v)',sum
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
c
	open(4,file='dcm2.bin',form='unformatted')
	write(4) ((float(ud(i,j)),float(uref(i)),float(hd(i,j))
     &           ,i=1,nx-1),j=1,ny-1)
	close(4)
      call ipr_jnl(nx,ny,lx,ly,dx,dy)
c
20	format(1x,20(f10.6,1x))
21	format(1x,a10,1p,e13.6)
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
c -----------------------------------------------------------
c   ecriture du fichier .jnl des resultats pour ferret
c -----------------------------------------------------------
c
	subroutine ipr_jnl(nx,ny,lx,ly,dx,dy)
	implicit none
      double precision lx,ly,dx,dy,ll1,ll2
      integer nx,ny,i,j
      character fjnl*8,fbin*8,anum*3
c
      fjnl = 'dcm2.jnl'
      fbin = 'dcm2.bin'
c
	open(20,file=fjnl)
      ll1 = - 0.5d-3*lx + 0.5d-3*dx 
      ll2 = - ll1
	write(20,20) 'x',ll1,ll2,1.d-3*dx,'xa'
      ll1 = - 0.5d-3*ly + 0.5d-3*dy 
      ll2 = - ll1
	write(20,20) 'y',ll1,ll2,1.d-3*dy,'ya'
	write(20,21) 
	write(20,22) -lx*0.5d-3,lx*0.5d-3,-ly*0.5d-3,ly*0.5d-3
	write(20,23) 'u,v,h',0,nx-1,ny-1,fbin
      write(20,*) 'set variable /bad=-999.9 u'
      write(20,*) 'set variable /bad=-999.9 v'
      write(20,*) 'set variable /bad=-999.9 h'
	write(20,25) 'U','u'
	write(20,25) 'V','v'
	write(20,25) 'H','h'
c
c      write(20,*) 'set window /aspect=1.:axis'
      write(20,*) 'set region /x=-500:500/y=-500:500'
c      write(20,*) 'let va = abs(v)'
c      write(20,*) 'shade  va'
c      write(20,*) 'list va[x=@din,y=@din]'
      write(20,*) 'let uref = 1.858740172301e-3',
     &            '*sin(6.2831853072e-3*x[g=ga]+0*y[g=ga])'
      write(20,*) 'let uerr = abs(u-uref)'
      write(20,*) 'list uerr[x=@din,y=@din]'
	close(20)
c
 19   format(i3)
20    format('define axis/',a1,'=',f7.1,':',f7.1,':',f6.1,
     &'/units=km ',a2)
21    format('define grid/x=xa/y=ya ga')
22    format('define region /x=',f7.1,':',f7.1,'/y=',f7.1,':',f7.1,
     &' basin')
23    format('file /var=',a5,'/grid=ga/FORMAT=unformatted/skip=',i1,
     &'/COLUMNS=`',i3,'*',i3,'*3','` ',a10)
24    format('let ',a2,'=',a1,'*mask[d=1]')
25    format('set variable /title="',a1,'" ',a2)
26    format('shade ',a2)
27    format('fill /over ',a2)

	return
	end	
