c******************************************************************
c initialization of the variables for the discontinuous spectral 
c element model (SPOC)
c******************************************************************
      program ch2
      implicit none
c
C
C 2-D array declaration
C
      INCLUDE '../include/maille.inc'
      include '../include/sp.inc'
c
      integer nn,ne,nf,nebd,nne
      INTEGER IN(3,NEDIM),tseg(2,nfdim),pseg(2,nfdim),
     &itm(3,nedim),nmit(nndim),tmit(NFTRDIM*NNDIM),
     & imit(NFTRDIM*NNDIM),ssit(2,nfdim),sbd(nbedim),
     & its(3,nedim)
      double precision norm(2,nfdim),dx(2,nedim),dy(2,nedim),ar(nedim),
     &                 scx,scy,xgr(nndim),ygr(nndim)
      logical nobd(nndim),pco(nndim)
c
c modes spectraux
c
      double precision 
     &       amm(nnmod,nnmod),
     &       u(nnmod,nedim),
     &       v(nnmod,nedim),
     &       h(nnmod,nedim),
     &       w(nnmod,nedim),
     &       z(nnmod,nedim),
     &       hb(nnmod,nedim)
      integer i,j,k,ind,k1,k2,ne1,ne2,cote,i1,i2
      integer kc1,kc2,kc3,kc4
      common /i4_cheb/ kc1,kc2,kc3,kc4
      integer info,toto
      double precision mu0,mb0,f0,beta,g0,hh0,lx,ly,lvent
c
c quadrature de Gauss a ng points
c
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision xi(ng_max),yi(ng_max),wi(ng_max)
      common /gauss/ wi,xi,yi
      double precision wif(ngf_max),tif(ngf_max)
      common /gaussf/ wif,tif
      double precision pi2(ng_max,nnmod)
      common /gauss2/ pi2
c
      double precision f2,fvaguex,fvaguey,ampl,lf,x0,y0,
     &                 ftourbu,ftourbv,ftourbh,fvent1,fvent2,fvent3,
     &                 vmax,b0,fond0,fond1,amplvent,fsinus,fkelvin
      external f2,fvaguex,fvaguey,ftourbu,ftourbv,ftourbh,fvent1,
     &           fvent2,fvent3,
     &           fond0,fond1,fsinus,fkelvin
      double precision a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
c
c graphique
c
      include '../include/graph.inc'
      integer ngraph
      double precision pgraph(nnmod,ngradim,ngradim)
      character*10 typewindx,typewindy,typefond
c
c------------------------------------
c read input parameter file
c------------------------------------
c
	open(2,file='oc.inp',status='old')
	read(2,*) nc,ng,ngf
        read(2,*) lx,ly
	read(2,*) 
	read(2,*) 
	read(2,*) mu0
	read(2,*) mb0
	read(2,*) f0
	read(2,*) beta
	read(2,*) scx,scy
	read(2,*) g0
	read(2,*) hh0
	read(2,*) 
	read(2,*) 
	read(2,*) amplvent,lvent
	read(2,26) typewindx
	read(2,26) typewindy
	read(2,26) typefond
	close(2)
26	format(a10)
c
c------------------------------------
c build the basis functions
c------------------------------------
c
      call read_gauss(ng,xi,yi,wi,ngf,wif,tif,nc,nm)
      call cal_mat(nnmod,ng,nc,amm,f2,xi,yi,wi)
      call cal_p2(xi,yi,pi2)
c
c LU Decomposition of the mass matrices
c
      call dpotrf('U',nm,amm,nnmod,info)
      if (info.ne.0) write(*,*) 'pb avec matrice'
c
c------------------------------------
c read the mesh and calculate the pointers
c------------------------------------
c
      call read_mesh(nn,ne,in,xgr,ygr,scx,scy)
      call pointer_mesh(nn,ne,nf,in,itm,nmit,tmit,imit,pseg,tseg)
      call check_mesh(ne,nf,in,itm,tseg)
c
c------------------------------------
c compute the initial state and the forcing
c------------------------------------
c
c
	do i=1,ne
	   do j=1,nm
	      u(j,i)=0.d0
	      v(j,i)=0.d0
	      h(j,i)=0.d0
	      w(j,i)=0.d0
	      z(j,i)=0.d0
	   enddo
	enddo
c
c-------------------------------------------------
c initial state
c
      write(*,*)
      write(*,30)
      write(*,*) 'CHOOSE INITIAL STATE'
      write(*,30)
      write(*,*)
	write(*,*) 'rest              =0'
	write(*,*) 'wave in x         =1'
	write(*,*) 'wave in y         =2'
	write(*,*) 'eddy exp QG       =3'
	write(*,*) 'eddy sinus        =4'
	write(*,*) 'Kelvin wave       =5'
	read(*,*) toto
	write(*,*) 'you chose:',toto
c
      if (toto.eq.1) then
	   write(*,*) 'amplitude (m)'
	   read(*,*) a1
	   write(*,*) 'wave length (m)'
	   read(*,*) a2
	   call calh(fvaguex,ne,h,in,xgr,ygr,amm)
      elseif (toto.eq.2) then
	   write(*,*) 'amplitude (m)'
	   read(*,*) a1
	   write(*,*) 'wave length (m)'
	   read(*,*) a2
	   call calh(fvaguey,ne,h,in,xgr,ygr,amm)
      elseif (toto.eq.3) then
	 write(*,*) 'vmax (m/s)'
	 read(*,*) vmax
         write(*,*) 'max speed in the eddy in (m/s):',vmax
	 write(*,*) 'radius (m)'
	 read(*,*) a2
         B0      = (1.d0/a2)**2
         a3 = 0.d0
         a4 = 0.d0
         a1 = EXP(0.5d0) * VMAX * SQRT( 2.d0 * B0)
         a2 = b0
	 call calh(ftourbu,ne,u,in,xgr,ygr,amm)
         call calh(ftourbv,ne,v,in,xgr,ygr,amm)
         a1 = EXP(0.5d0) * VMAX * f0 / ( G0 * SQRT( 2.d0 * B0 ) )
	 write(*,*) 'amplitude of the eddy in (m):', a1
	 call calh(ftourbh,ne,h,in,xgr,ygr,amm)
      elseif (toto.eq.4) then
	 write(*,*) 'ampl (m)'
	 read(*,*) a1
	 write(*,*) 'radius (m)'
	 read(*,*) a2
         a3 = 0.d0
         a4 = 0.d0
	 call calh(fsinus,ne,h,in,xgr,ygr,amm)
	 do i=1,ne
	     do j=1,nm
	        u(j,i)=h(j,i)
		  v(j,i)=h(j,i)
	     enddo
	 enddo
      elseif (toto.eq.5) then
	 write(*,*) 'amplitude (m)'
	 read(*,*) a1
	 write(*,*) 'amplitude of the Kelvin wave in (m):', a1
         a3 = 0.d0
c         a5 = 1.d-10
         a5 = 1.d0/300.d3**2
	 write(*,*) 'length (m)'
	 read(*,*) ampl
         write(*,*) 'length of the Kelvin wave in (m):',ampl
	 a5 = 1.d0 / ampl**2
         a4 = -ly*0.5d0
         a2 = f0 / sqrt ( g0 * hh0 )
	 call calh(fkelvin,ne,h,in,xgr,ygr,amm)
         a1 = a1 * g0 / sqrt ( g0 * hh0 )
         write(*,*) 'max speed in the Kelvin wave in (m/s):',a1
	 call calh(fkelvin,ne,u,in,xgr,ygr,amm)
      endif
c
      write(*,*)
      write(*,30)
      write(*,*) 'CHOOSE WIND FORCING'
      write(*,30)
      write(*,*) 'zero wind forcing              =0'
      write(*,*) 'sinus wind in x varying in y   =1'
      write(*,*) 'wind in x arctan varying in  y =2'
      write(*,*) 'cosinus wind in x varying in y =3'
      if (typewindx.eq.'zero') toto=0
      if (typewindx.eq.'fvent1') toto=1
      if (typewindx.eq.'fvent2') toto=2
      if (typewindx.eq.'fvent3') toto=3
      write(*,*) 'you chose:',toto
c
      a1 = amplvent
      a2 = lvent
      if (toto.eq.1) then
	  call calh(fvent1,ne,w,in,xgr,ygr,amm)
      elseif (toto.eq.2) then
	  call calh(fvent2,ne,w,in,xgr,ygr,amm)
      elseif (toto.eq.3) then
	  call calh(fvent3,ne,w,in,xgr,ygr,amm)
      endif
c
c
      write(*,*)
      write(*,30)
      write(*,*) 'CHOOSE BATHYMETRY'
      write(*,30)
      write(*,*)
      write(*,*) 'constant hh0      =0'
      write(*,*) 'haidvogel channel =1'
      write(*,*) typefond
      if (typefond.eq.'constant') toto=0
      if (typefond.eq.'fond1') toto=1
      write(*,*) 'you chose:',toto
c
	if (toto.eq.1) then
         a1 = lx
         a2 = ly
	   call calh(fond1,ne,hb,in,xgr,ygr,amm)
        else
           do i=1,ne
              hb(1,i)=hh0
              do j=2,nm
                 hb(j,i)=0.d0
              enddo
           enddo
	endif
c
c
c
	ngraph=10
	call cal_pgraph(nc,ngraph,pgraph,ngradim)
	call connec(nne,ne,ngraph)
	call write_gra(ne,u,'u',0,nm,ngraph,xgr,ygr,in,pgraph)
	call write_gra(ne,v,'v',0,nm,ngraph,xgr,ygr,in,pgraph)
	call write_gra(ne,h,'h',0,nm,ngraph,xgr,ygr,in,pgraph)
	write(*,*) 'u/v/h0.gra graphic output files have been created'
c
c
	open(1,file='init.bin',form='unformatted')
	write(1) 0.d0
	write(1) ((u(j,i),j=1,nm),i=1,ne)
	write(1) ((v(j,i),j=1,nm),i=1,ne)
	write(1) ((h(j,i),j=1,nm),i=1,ne)
	write(1) ((w(j,i),j=1,nm),i=1,ne)
	write(1) ((z(j,i),j=1,nm),i=1,ne)
	write(1) ((hb(j,i),j=1,nm),i=1,ne)
	close(1)
c
      write(*,*)
      write(*,30)
      write(*,*) 'INIT.BIN has been written'
      write(*,30)
      write(*,*)
c
20 	format(20(f13.10,x))
21 	format(x,3(i3,x),2(f13.10,x))
30    format(1x,60('-'))
c
      end
c
c******************************************************************
c
      subroutine calh(ff,ne,h,in,xgr,ygr,amm)
      implicit none
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      integer ne
      double precision xgr(*),ygr(*),h(nnmod,*),amm(nnmod,*),
     &                 b(nnmod),dia1(nnmod)
      double precision mff(ng_max),pi2(ng_max,nnmod)
      double precision agx,agy,bgx,bgy,xof,yof,x,y,x1,y1,deta,ans,ansx
      integer i,i1,i2,i3,j,k,l,info
      integer in(3,*)
      double precision ff,fsol,po
      external ff,po
      double precision a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
      integer kc1,kc2
      double precision wi(ng_max),xi(ng_max),yi(ng_max)
      common /gauss/ wi,xi,yi
      common /gauss2/ pi2
c
      do i=1,nm
         dia1(i) = 1.d0/amm(i,i)
      enddo
c
	do i=1,ne
	   i1 = in(1,i)
	   i2 = in(2,i)
	   i3 = in(3,i)
	   agx = (xgr(i2)-xgr(i1))*0.5d0
	   agy = (ygr(i2)-ygr(i1))*0.5d0
	   bgx = (xgr(i3)-xgr(i1))*0.5d0
	   bgy = (ygr(i3)-ygr(i1))*0.5d0
	   xof = xgr(i1)
	   yof = ygr(i1)
	   deta = agx * bgy - agy * bgx
c
	   do j=1,ng
	      x = xi(j)
	      y = yi(j)
	      x1 = agx * (x+1.d0) + bgx * (y+1.d0) + xof
	      y1 = agy * (x+1.d0) + bgy * (y+1.d0) + yof
	      mff(j) = ff(x1,y1)
	   enddo
c
           call xytospec(b,mff)
           call dtrsm4(nm,amm,b,dia1)

c	   call f07fef('U',nm,1,amm,nnmod,b,nm,info)
	   do j=1,nm
	      h(j,i) = b(j)
	   enddo
	enddo
	return
	end
     
