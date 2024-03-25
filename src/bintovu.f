c******************************************************************
c equations shallow waters
c spectral discontinuous elements
c rk4 integration in time
c conditions frontieres pour u 
c viscosite 
c friction
c plan beta 
c******************************************************************
	program shwater
	implicit none
c
C
C 2-D array declaration
C
      include '../include/maille.inc'
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      integer nn,ne,nf,nebd
      INTEGER IN(3,NEDIM),tseg(2,nfdim),pseg(2,nfdim),
     &itm(3,nedim),nmit(nndim),tmit(NFTRDIM*NNDIM),
     & imit(NFTRDIM*NNDIM),ssit(2,nfdim),sbd(nbedim),
     & its(3,nedim),nbdl,pbdl(nedim),nint,pint(nedim)
      double precision norm(2,nfdim),
     &                 scx,scy,xgr(nndim),ygr(nndim)
      logical nobd(nndim)
      logical bdl(nedim)
c
c modes spectraux
c
      double precision 
     &       u0(nnmod,nedim),v0(nnmod,nedim),h0(nnmod,nedim),
     &       w0(nnmod,nedim),z0(nnmod,nedim),hb(nnmod,nedim)
      integer i,j,k,l,i1,i2,i3,ns1
      integer nite,ite,info,per,nper,per2,yearmax,daysmax,nitemax,
     &npermax
      double precision 
     &       dt,mu0,mb0,f0,beta,g0,hh0,time,ke,pe,mass,ke0,pe0,
     &       lx,ly,time0,amplvent,ax,ay,bx,by,deta,detai,
     &       dnx,dny,dis
      character*10 filename
c
c graphique
c
      include '../include/graph.inc'
      integer ngraph,nne
      double precision pgraph(nnmod,ngradim,ngradim)
c
c
c------------------------------------
c lecture des donnees
c------------------------------------
c
	open(2,file='oc.inp',status='old')
	read(2,*) nc
	read(2,*) lx,ly
	read(2,*) yearmax,daysmax,nitemax
	read(2,*) dt
	read(2,*) mu0
	read(2,*) mb0
	read(2,*) f0
	read(2,*) beta
	read(2,*) scx,scy
	read(2,*) g0
	read(2,*) hh0
	read(2,*) npermax
	read(2,*) per2
	read(2,*) amplvent
	close(2)
c
c
c------------------------------------
c construction de la base de polynome
c------------------------------------
c
	nm = (nc+1)*(nc+2)/2
	write(*,*) 'rang matrice',nm
	ngraph=10
	call cal_pgraph(nc,ngraph,pgraph,ngradim)
c
c------------------------------------
c lecture du maillage
c------------------------------------
c
	write(*,*) 'filename for the mesh'
	read(*,23) filename
      open(1,file=filename,status='old')
      read(1,*) nn,ne
	do i=1,nn
	   read(1,*) j,xgr(i),ygr(i)
	enddo
	do i=1,ne
	   read(1,*) k,(in(j,i),j=1,3)
	enddo
      close(1)
c
      write(*,*) 'fin lecture du fichier cm.dat.'
	call maille2(nn,ne,nf,in,itm,nmit,tmit,imit,
     1                  pseg,tseg,xgr,ygr)
      write(*,*) nn,ne,nf
      do i=1,nn
	  xgr(i) = xgr(i) * scx
	  ygr(i) = ygr(i) * scy
	enddo
c
c
c------------------------------------
c
	write(*,*) 'filename'
	read(*,23) filename
	open(1,file=filename,form='unformatted',status='old')
	read(1) time0
	read(1) ((u0(j,i),j=1,nm),i=1,ne)
	read(1) ((v0(j,i),j=1,nm),i=1,ne)
	read(1) ((h0(j,i),j=1,nm),i=1,ne)
	read(1) ((w0(j,i),j=1,nm),i=1,ne)
	read(1) ((z0(j,i),j=1,nm),i=1,ne)
	read(1) ((hb(j,i),j=1,nm),i=1,ne)
	close(1)
c
	write(*,22) 'time',time0/86400.d0
c
	write(*,*)'entre numero fichier graphique en sortie'
	read(*,*) nper
	call connec(nne,ne,ngraph)
	call spectovu(nne,ne,nper,ngraph)
	call spectovubin(ngraph,ne,xgr,ygr,in,nc,u0,v0,h0,nper,pgraph)
c
c
20      format(20(e10.3,1x))
21      format(1x,i7,(e13.6,1x),i3)
22      format(a7,5(1p,e13.6,1x))
23	format(a10)
c
c
	end
c
