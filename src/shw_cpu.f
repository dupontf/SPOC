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
     &       amm(nnmod,nnmod),
     &       ammo(nnmod,nnmod),
     &       amml(nnmod,nnmod),
     &       ammx(nnmod,nnmod),
     &       ammy(nnmod,nnmod),
     &       agx(nnmod,nnmod),
     &       agy(nnmod,nnmod),
     &       u0(nnmod,nedim),v0(nnmod,nedim),h0(nnmod,nedim),
     &       u1(nnmod,nedim),v1(nnmod,nedim),h1(nnmod,nedim),
     &       du1(nnmod,nedim),dv1(nnmod,nedim),dh1(nnmod,nedim),
     &       du2(nnmod,nedim),dv2(nnmod,nedim),dh2(nnmod,nedim),
     &       du3(nnmod,nedim),dv3(nnmod,nedim),dh3(nnmod,nedim),
     &       du4(nnmod,nedim),dv4(nnmod,nedim),dh4(nnmod,nedim),
     &       w0(nnmod,nedim),z0(nnmod,nedim),hb(nnmod,nedim),
     &       vect(5,nedim),normt(2,nedim)
      double precision f2,dfx,dfy,fx,fy
      integer i,j,k,l,i1,i2,i3,ns1
      integer ip,nne
      external f2,dfx,dfy,ip,fx,fy
      integer nite,ite,info,per,nper,per2,yearmax,daysmax,nitemax,
     &npermax
      double precision fc(nnmod,nedim),
     &       dt,mu0,mb0,f0,beta,g0,hh0,time,ke,pe,mass,ke0,pe0,
     &       lx,ly,time0,amplvent,ax,ay,bx,by,deta,detai,
     &       dnx,dny,dis
c
c
c quadrature de Gauss a ng points
c
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision xi(ng_max),yi(ng_max),wi(ng_max)
      common /gauss/ wi,xi,yi
      double precision wif(ngf_max),tif(ngf_max)
      common /gaussf/ wif,tif
      double precision pi2(ng_max,nnmod),
     &                 pico(ngf_max,nnmod,3),cip(2,0:nnmod)
      common /gauss2/ pi2
      common /gauss3/ pico
      save /const/
      common /matrices/ amm,amml,agx,agy,ammx,ammy,ammo
c
c graphique
c
	include '../include/graph.inc'
	integer ngraph
	double precision pgraph(nnmod,ngradim,ngradim)
C
      real ts,dtime
      real tarray(2) 
      external dtime
      integer pass,npass
c
c
c------------------------------------
c read input parameter file
c------------------------------------
c
	open(2,file='oc.inp',status='old')
	read(2,*) nc,ng,ngf
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
c----------------------------------------------
c write message
c
      write(*,*)
      write(*,30)
      write(*,*) 'READ INPUT PARAMETERS'
      write(*,30)
      write(*,*)
      write(*,*) 'time step              (s)   :',dt
      write(*,*) 'horizontal viscosity   (m2/s):',mu0
      write(*,*) 'linear bottom friction (m/s) :',mb0
      write(*,*) 'Coriolis f0            (1/s) :',f0
      write(*,*) 'Coriolis beta          (1/ms):',beta
      write(*,*) 'reduced gravity        (m/s2):',g0
      write(*,*) 'depth at rest          (m)   :',hh0
      nite=nitemax+int(daysmax*86400.d0/dt)
     &            +int(yearmax*365.d0*86400.d0/dt)
      per=int(nite/npermax)
      write(*,*) 'number of iterations                           :',nite
      write(*,*) 'number of iterations between two graphic output:',per
c
c
c------------------------------------
c build the basis functions
c------------------------------------
c
      call read_gauss(ng,xi,yi,wi,ngf,wif,tif,nc,nm)
      call cal_p2(xi,yi,pi2)
      call cal_pico(pico,tif)
      call calip(nc,cip)
c
c graphic output grid
c
      ngraph=10
      call cal_pgraph(nc,ngraph,pgraph,ngradim)
c
c build the matrices
c
      call cal_mat(nnmod,ng,nc,amm ,f2,xi,yi,wi)
      call cal_mat(nnmod,ng,nc,ammx,fx,xi,yi,wi)
      call cal_mat(nnmod,ng,nc,ammy,fy,xi,yi,wi)
      call cal_mat(nnmod,ng,nc,agx,dfx,xi,yi,wi)
      call cal_mat(nnmod,ng,nc,agy,dfy,xi,yi,wi)
      call cal_matl(nnmod,nc,amm,amml)
      do j=1,nm
         do i=1,nm
	    ammo(i,j)=amm(i,j)
         enddo
      enddo
c
c LU Decomposition of the mass matrices
c
      call dpotrf('U',nm,amm,nnmod,info)
      if (info.ne.0) write(*,*) 'pb avec matrice'
      call dpotrf('U',nm-nc-1,amml,nnmod,info)
      if (info.ne.0) write(*,*) 'pb avec matrice'
c
c------------------------------------
c read the mesh and calculate the pointers
c------------------------------------
c
      call read_mesh(nn,ne,in,xgr,ygr,scx,scy)
      call pointer_mesh(nn,ne,nf,in,itm,nmit,tmit,imit,pseg,tseg)
      call check_mesh(ne,nf,in,itm,tseg)
      call cal_rotnf(nf,xgr,ygr,norm,pseg)
      call cal_its(nf,itm,its,tseg,ssit)
      call calpre(nn,ne,nf,in,xgr,ygr,vect,itm,norm,normt,bdl,
     &                  nint,pint,nbdl,pbdl,nobd,nebd,sbd,tseg,pseg)
c
c
c--------------------------------------------------------
c compute coriolis parameter
c
      call coriolis(nc,ne,in,xgr,ygr,fc,f0,beta)
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
      do i=1,ne
         do j=1,nm
            u1(j,i)=0.d0
            v1(j,i)=0.d0
            h1(j,i)=0.d0
         enddo
      enddo
c
c-------------------------------------------------------------
c defines time iterators
c
      time = time0
      ite=0
      nper=0
c
        write(*,*) 'enter how many pass?'
	read(*,*) npass
c------------------------------------------
c compute energy and output graphic
c
      call connec(nne,ne,ngraph)
      call write_gra(ne,u0,'u',nper,nm,ngraph,xgr,ygr,in,pgraph)
      call write_gra(ne,v0,'v',nper,nm,ngraph,xgr,ygr,in,pgraph)
      call write_gra(ne,h0,'h',nper,nm,ngraph,xgr,ygr,in,pgraph)
      call energy(ne,u1,v1,h1,g0,hb,ke0,pe0,mass,
     &                  in,xgr,ygr,ammo,pi2)
      call energy(ne,u0,v0,h0,g0,hb,ke,pe,mass,
     &                  in,xgr,ygr,ammo,pi2)
      open(24,file='oc.energ')
      write(24,*) 'days, kinetic, potential energies'
      write(24,22) time /86400.d0,ke,pe-pe0,mass
      close(24)
c
c-------------------------------------------------------------
c write message
c
      write(*,30)
      write(*,31)
      write(*,30)
      write(*,23) ite,time/86400.d0,nper,ke,pe-pe0
c
c------------------------------------
c start iterating in time
c------------------------------------
c
	
	ts = dtime(tarray)
	do pass=1,npass
        ite=0
	
	do while (ite.lt.nite)
	ite = ite + 1
c
c-------------------------------rk4------------------------------
c
	  call cal_var(nf,ne,u0,v0,h0,du1,dv1,dh1,
     &            tseg,pseg,ssit,itm,its,norm,nbdl,pbdl,nint,pint,
     &            bdl,cip,vect,normt,
     &            mu0,mb0,fc,beta,w0,z0,g0,hb,time)
	  do i=1,ne
	   do j=1,nm
	    u1(j,i) = u0(j,i) + dt * du1(j,i) * 0.5d0
	    v1(j,i) = v0(j,i) + dt * dv1(j,i) * 0.5d0
	    h1(j,i) = h0(j,i) + dt * dh1(j,i) * 0.5d0
	   enddo
	  enddo
	  call cal_var(nf,ne,u1,v1,h1,du2,dv2,dh2,
     &            tseg,pseg,ssit,itm,its,norm,nbdl,pbdl,nint,pint,
     &            bdl,cip,vect,normt,
     &            mu0,mb0,fc,beta,w0,z0,g0,hb,time)
	  do i=1,ne
	   do j=1,nm
	    u1(j,i) = u0(j,i) + dt * du2(j,i) * 0.5d0
	    v1(j,i) = v0(j,i) + dt * dv2(j,i) * 0.5d0
	    h1(j,i) = h0(j,i) + dt * dh2(j,i) * 0.5d0
	   enddo
	  enddo
	  call cal_var(nf,ne,u1,v1,h1,du3,dv3,dh3,
     &            tseg,pseg,ssit,itm,its,norm,nbdl,pbdl,nint,pint,
     &            bdl,cip,vect,normt,
     &            mu0,mb0,fc,beta,w0,z0,g0,hb,time)
	  do i=1,ne
	   do j=1,nm
	    u1(j,i) = u0(j,i) + dt * du3(j,i)
	    v1(j,i) = v0(j,i) + dt * dv3(j,i)
	    h1(j,i) = h0(j,i) + dt * dh3(j,i)
	   enddo
	  enddo
	  call cal_var(nf,ne,u1,v1,h1,du4,dv4,dh4,
     &            tseg,pseg,ssit,itm,its,norm,nbdl,pbdl,nint,pint,
     &            bdl,cip,vect,normt,
     &            mu0,mb0,fc,beta,w0,z0,g0,hb,time)
	  do i=1,ne
	   do j=1,nm
	    u0(j,i) = u0(j,i) + dt / 6.d0 * ( 
     &    du1(j,i) + 2.d0 * du2(j,i) + 2.d0 * du3(j,i) + du4(j,i) )
	    v0(j,i) = v0(j,i) + dt / 6.d0 * ( 
     &    dv1(j,i) + 2.d0 * dv2(j,i) + 2.d0 * dv3(j,i) + dv4(j,i) )
	    h0(j,i) = h0(j,i) + dt / 6.d0 * ( 
     &    dh1(j,i) + 2.d0 * dh2(j,i) + 2.d0 * dh3(j,i) + dh4(j,i) )
         enddo
	  enddo
c
	time = time0 + dt * dfloat(ite)
c
	enddo
	enddo         ! npass end loop
	ts = dtime(tarray)
c
c-------------------------------------------------------
c final output
c
      open(2,file='temp.bin',form='unformatted')
      write(2) time
      write(2) ((u0(j,i),j=1,nm),i=1,ne)
      write(2) ((v0(j,i),j=1,nm),i=1,ne)
      write(2) ((h0(j,i),j=1,nm),i=1,ne)
      write(2) ((w0(j,i),j=1,nm),i=1,ne)
      write(2) ((z0(j,i),j=1,nm),i=1,ne)
      write(2) ((hb(j,i),j=1,nm),i=1,ne)
      close(2)
c
c-------------------------------------------------------
c energy
c
         call energy(ne,u0,v0,h0,g0,hb,ke,pe,mass,
     &                  in,xgr,ygr,ammo,pi2)
         open(24,file='oc.energ',access='append')
         write(24,22) time /86400.d0,ke,pe-pe0,mass
         close(24)
c
c-------------------------------------------------------
c visualisation
c
         nper=1
         call write_gra(ne,u0,'u',nper,nm,ngraph,xgr,ygr,in,pgraph)
         call write_gra(ne,v0,'v',nper,nm,ngraph,xgr,ygr,in,pgraph)
         call write_gra(ne,h0,'h',nper,nm,ngraph,xgr,ygr,in,pgraph)
         write(*,23) ite,time/86400.d0,nper,ke,pe-pe0
	write(*,*) 'time elapsed:',ts,ts/float(npass)
c
22    format(1x,5(1p,e13.6,1x))
30    format(1x,60('-'))
31    format('iteration|    days   | output # |     TKE     |',
     &'     TPE')
23    format(i8,' |  ',f8.2,' |  ',i6,'  |   ',1P,E8.1,'  |   ',1P,E8.1)
c
c
      stop
	end
c
c
c******************************************************************
c defines other pointers relevant to the spectral element computation
c******************************************************************
c
      subroutine calpre(nn,ne,nf,in,xgr,ygr,vect,itm,norm,normt,bdl,
     &                  nint,pint,nbdl,pbdl,nobd,nebd,sbd,tseg,pseg)
      implicit none
      include '../include/maille.inc'
      integer nn,ne,nf,nebd
      integer in(3,*),itm(3,*),nint,pint(*),nbdl,pbdl(*),
     & tseg(2,*),pseg(2,*),sbd(*)
      logical bdl(*),nobd(*)
      double precision xgr(*),ygr(*),vect(5,*),norm(2,*),normt(2,*),
     &                 ax,bx,ay,by,deta,detai,dnx,dny,dis
      integer i,ns1,i1,i2,i3
      do i=1,ne
	  bdl(i)=.false.
      enddo
      do i=1,nn
	nobd(i)=.false.
      enddo
      nebd=0
      do i=1,nf
	   if (tseg(2,i).eq.0) then
	      nobd(pseg(1,i))=.true.
	      nobd(pseg(2,i))=.true.
		nebd=nebd+1
		if (nebd.gt.nbedim) then 
		 write(*,*) 'NBEDIM too small'
		 stop
		endif
		sbd(nebd)=i
	      bdl(tseg(1,i))=.true.
	   endif
      enddo
      write(*,*) 'nebd:',nebd

      do i=1,ne
         i1 = in(1,i)
         i2 = in(2,i)
         i3 = in(3,i)
         ax = (xgr(i2)-xgr(i1))*0.5d0
         bx = (xgr(i3)-xgr(i1))*0.5d0
         ay = (ygr(i2)-ygr(i1))*0.5d0
         by = (ygr(i3)-ygr(i1))*0.5d0
         deta = ax * by - ay * bx
         detai = 1.d0/deta
         vect(1,i)=ax
         vect(2,i)=bx
         vect(3,i)=ay
         vect(4,i)=by
         vect(5,i)=detai
         ns1 = itm(1,i)
         dnx = norm(1,ns1)
         dny = norm(2,ns1)
         dis = sqrt(dnx**2+dny**2)
         dnx = dnx / dis
         dny = dny / dis
         normt(1,i)=dnx
         normt(2,i)=dny
      enddo
      nbdl=0
      nint=0
      do i=1,ne
         if (bdl(i)) then
            nbdl=nbdl+1
            pbdl(nbdl)=i
         else
            nint=nint+1
            pint(nint)=i
         endif
      enddo
      write(*,*) 'nbdl:',nbdl
      write(*,*) 'nint:',nint
c
      return
      end
c
c******************************************************************
c
	subroutine testnan(ke,time)
	double precision ke,time
      character*20 name
	integer ideb,ifin
      write(name,20) ke
20    format(1p,e13.6)
	ideb=1
	do while (name(ideb:ideb).eq.' ')
	   ideb=ideb+1
	enddo
	ifin=20
	do while (name(ifin:ifin).eq.' ')
	   ifin=ifin-1
	enddo
	if (name(ideb:ifin).eq.'nan') then
	   write(*,*) 'j''ai un probleme!!',time /86400.d0
	   stop
	endif
	return
      end
c
