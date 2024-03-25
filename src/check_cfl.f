c******************************************************************
c check cfl
c******************************************************************
	program check
	implicit none
c
C
C 2-D array declaration
C
      include '../include/maille.inc'
      include '../include/sp.inc'
      integer nn,ne,nf,nebd
      INTEGER IN(3,NEDIM)
      double precision scx,scy,xgr(nndim),ygr(nndim)
c
c modes spectraux
c
      integer i,j
      double precision 
     &       u0(nnmod,nedim),v0(nnmod,nedim),h0(nnmod,nedim),
     &       w0(nnmod,nedim),z0(nnmod,nedim),hb(nnmod,nedim)
      double precision 
     &       dt,mu0,mb0,f0,beta,g0,hh0,time,
     &       yearmax,daysmax,nitemax,lx,ly,time0
c
c quadrature de Gauss a ng points
c
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision xi(ng_max),yi(ng_max),wi(ng_max)
      common /gauss/ wi,xi,yi
      double precision wif(ngf_max),tif(ngf_max)
      common /gaussf/ wif,tif
      save /const/
c
c graphique
c
      include '../include/graph.inc'
      integer ngraph
      double precision pgraph(nnmod,ngradim,ngradim)
c
c
c------------------------------------
c lecture des donnees
c------------------------------------
c
	open(2,file='oc.inp',status='old')
	read(2,*) nc,ng,ngf
	read(2,*) 
	read(2,*) 
	read(2,*) dt
	read(2,*) 
	read(2,*) 
	read(2,*) 
	read(2,*) 
	read(2,*) scx,scy
	read(2,*) g0
	close(2)
26	format(a10)
c
c
c------------------------------------
c build the basis functions
c------------------------------------
c
      call read_gauss(ng,xi,yi,wi,ngf,wif,tif,nc,nm)
c
c graphic output grid
c
      ngraph=10
      call cal_pgraph(nc,ngraph,pgraph,ngradim)
c
c
c------------------------------------
c read the mesh and calculate the pointers
c------------------------------------
c
      call read_mesh(nn,ne,in,xgr,ygr,scx,scy)
c
c
c------------------------------------
c read initial state and forcing
c
      open(1,file='init.bin',form='unformatted',status='old')
      read(1) time0
      read(1) ((u0(j,i),j=1,nm),i=1,ne)
      read(1) ((v0(j,i),j=1,nm),i=1,ne)
      read(1) ((h0(j,i),j=1,nm),i=1,ne)
      read(1) ((w0(j,i),j=1,nm),i=1,ne)
      read(1) ((z0(j,i),j=1,nm),i=1,ne)
      read(1) ((hb(j,i),j=1,nm),i=1,ne)
      close(1)
c
      write(*,*) 'STARTING TIME  :',time0/86400.d0
c
c      call checkcfl(nc,ne,nnmod,hb,in,dt,g0,xgr,ygr)
      call checkcfl2(ne,hb,u0,v0,in,dt,g0,xgr,ygr,
     &pgraph,ngradim,ngraph)

20      format(20(e10.3,1x))
30      format(1x,i5,3x,3(i5,1x))
21      format(1x,i7,(e13.6,1x),i3)
22      format(1x,5(1p,e13.6,1x))
23      format(1x,e13.6,1x,2(i5),1x,f9.3)
c
c
	end
c
