c******************************************************************
c equations shallow waters
c spectral discontinuous elements
c rk4 integration in time
c conditions frontieres pour u 
c viscosite 
c friction
c plan beta 
c raffinement de maillage en cas de probleme sur la vitesse
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
     &       w0(nnmod,nedim),z0(nnmod,nedim),hb(nnmod,nedim),
     &       du1(nnmod,nedim),dv1(nnmod,nedim),dh1(nnmod,nedim),
     &       du2(nnmod,nedim),dv2(nnmod,nedim),dh2(nnmod,nedim),
     &       du3(nnmod,nedim),dv3(nnmod,nedim),dh3(nnmod,nedim),
     &       du4(nnmod,nedim),dv4(nnmod,nedim),dh4(nnmod,nedim),
     &       vect(5,nedim),normt(2,nedim)
      double precision f2,dfx,dfy,fx,fy
      integer i,j,k,l,i1,i2,i3,ns1
      integer ip
      external f2,dfx,dfy,ip,fx,fy
      integer nite,ite,info,per,nper,per2,
     &npermax
      double precision fc(nnmod,nedim),
     &       dt,mu0,mb0,f0,beta,g0,hh0,time,ke,pe,mass,ke0,pe0,
     &       yearmax,daysmax,nitemax,
     &       lx,ly,time0,amplvent,timemax,timeper,
     &       ax,ay,bx,by,deta,detai,dnx,dny,dis
c
      double precision fvent1,fvent2,fvent3,fond1
      external fvent1,fvent2,fvent3,fond1
c
c variables pour le raffinement
c
      double precision
     &       xgrnew(nndim),ygrnew(nndim),
     &       u1(nnmod,nedim),v1(nnmod,nedim),h1(nnmod,nedim),
     &       w1(nnmod,nedim),z1(nnmod,nedim),hb1(nnmod,nedim),
     &       u2(nnmod,nedim),v2(nnmod,nedim),h2(nnmod,nedim),
     &       w2(nnmod,nedim),z2(nnmod,nedim),hb2(nnmod,nedim),
     &       unew(nnmod,nedim),vnew(nnmod,nedim),hnew(nnmod,nedim),
     &       wnew(nnmod,nedim),znew(nnmod,nedim),hbnew(nnmod,nedim)
c
      integer nn_ab,ne_ab,ndoub
      integer in_ab(3,2*nedim),tlevel(2*nedim),torigi(2*nedim),
     &tabsol(nedim),tmt(3,2*nedim),tenfan(4,2*nedim)
      double precision xgr_ab(2*nndim),ygr_ab(2*nndim)
      logical tdoubl(2*nedim),tutili(2*nedim)

      common /adap_const/ nn_ab,ne_ab,ndoub
      common /adap_int/ tlevel,torigi,tabsol,in_ab,tmt,tenfan
      common /adap_vec/ xgr_ab,ygr_ab
      common /adap_log/ tdoubl,tutili

      logical fderaff(2*nedim),
     &fconv(nedim),fconv_ab(2*nedim),fraff_ab(2*nedim),deraff,timedir
      integer nnnew,nenew,innew(3,nedim),nemaxnew,nemax
      integer tinter(nedim)
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
      integer nne
c
      logical conv
      character*10 typewindx,typewindy,typefond
      double precision fac1,fac2,fac3,lvent
      double precision a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
      integer toto,tk1,tk,perraff
c
c------------------------------------
c lecture des donnees
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
	read(2,*) amplvent,lvent
	read(2,26) typewindx
	read(2,26) typewindy
	read(2,26) typefond
	read(2,*) toto
	read(2,*) perraff
	read(2,*) fac1
	read(2,*) fac2
	read(2,*) fac3
	close(2)
26	format(a10)
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
      timemax = daysmax*86400.d0 + yearmax*365.d0*86400.d0
      timeper=timemax/dfloat(npermax)
      write(*,*) 'integration time         (days):',timemax/86400.d0
      write(*,*) 'time between two graphic output:',timeper/86400.d0
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
c--------------------------------------------------------
c---init pour adaptation
c
      if (toto.eq.1) then 
         write(*,*) 'maillage nouveau pou raffinement'
         call init_adap(nn,ne,xgr,ygr,in,tseg,itm)
      else
         call read_adap('raff.bin',ne)
      endif
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
      write(*,*) 'STARTING TIME  :',time0/86400.d0
c
c-------------------------------------------------------------
c defines time iterators
c
      time = time0
      ite=0
      nper=0
c
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
c prepare arrays for refinement
c------------------------------------
c
      call put1in2(ne,nm,u0,v0,h0,w0,z0,hb,u2,v2,h2,w2,z2,hb2)
c      call checkcfl(nc,ne,nnmod,hb,in,dt,g0,xgr,ygr)
      call checkcfl2(ne,hb,u0,v0,in,dt,g0,xgr,ygr,
     &pgraph,ngradim,ngraph)
      open(25,file='oc.maille')
      write(25,*) 'time           nn     ne        dt'
      write(25,23) time /86400.d0,nn,ne,dt
      close(25)
c
c------------------------------------
c iteration
c------------------------------------
c
	do while (time.lt.timemax.and.ite.lt.nitemax)
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
	time = time + dt
c
c	 write(*,*) ite,time
c-------------------------------------------------------
c energy
c
      if (mod(ite,per2).eq.0) then
         call energy(ne,u0,v0,h0,g0,hb,ke,pe,mass,
     &                  in,xgr,ygr,ammo,pi2)
         open(24,file='oc.energ',access='append')
         write(24,22) time /86400.d0,ke,pe-pe0,mass
         close(24)
         call testnan(ke,time)
      endif
c
c
c------------------------------------------------------------
c test on the convergence
c------------------------------------------------------------
c
      if (mod(ite,perraff).eq.0) then
         nper=nper+1
         write(*,23) ite,time/86400.d0,nper,ke,pe-pe0
	 call save_temp('temp1.bin',ne,nm,time,u0,v0,h0,w0,z0,hb)
	 call save_mesh('cm1.dat',nn,ne,in,xgr,ygr,scx,scy)
	 call save_adap('raff1.bin',ne)
c
c------------------------------------------------------------
c test the convergence criteria on the solution
c
         call cal_conv(nn,ne,nf,nc,pico,ssit,conv,fconv,fderaff,timedir,
     &           tseg,u0,v0,h0,w0,z0,hb,deraff,fac1,fac2,fac3)
         call preref(ne,fderaff,fconv,
     &           ne_ab,fconv_ab,fraff_ab,tutili,torigi,tabsol)
c
c
c------------------------------------------------------------
c if the convergence criteria are not met (or the derefine criterium are)
c
c
	 if (conv) then
c
            call deraffiner(nn,ne,fconv_ab,fderaff)
            call adapmaille(fconv_ab,nnnew,nenew,nemaxnew)
            call protridouble(nnnew,nenew,nemaxnew)
            call elimpptri(nnnew,nenew,xgrnew,ygrnew,innew,nemaxnew)
	    call inter_only_new(ne,nenew,in,innew,tinter)
c
            write(*,*) 'new mesh, nn,ne:',nnnew,nenew
c
c
c----------------------------------------------------------------------
c calcul des variables sur la nouvelle grille
c----------------------------------------------------------------------
c si l'erreur n'est pas trop grande, on continue l'integration en temps;
c autrement, on revient a la sauvegarde anterieur
c
	    if (.not.timedir) then
	       time=time0
	       call put1in2(ne,nm,u2,v2,h2,w2,z2,hb2,u0,v0,h0,w0,z0,hb)
	    endif
c
c-----------------------------------------------------------------
c interpolate the fields and replace the current mesh by the new one
c
            call newinter(nn,ne,nf,u0,v0,h0,w0,z0,hb,
     &           nnnew,nenew,innew,xgrnew,ygrnew,
     &           xgr,ygr,in,itm,tseg,pseg,
     &           typewindx,typewindy,typefond,
     &           amplvent,lvent,lx,ly,g0,hh0,f0,beta)
c
            call pointer_mesh(nn,ne,nf,in,itm,nmit,tmit,imit,pseg,tseg)
            call cal_rotnf(nf,xgr,ygr,norm,pseg)
            call cal_its(nf,itm,its,tseg,ssit)
            call calpre(nn,ne,nf,in,xgr,ygr,vect,itm,norm,normt,bdl,
     &                  nint,pint,nbdl,pbdl,nobd,nebd,sbd,tseg,pseg)
	    call coriolis(nc,ne,in,xgr,ygr,fc,f0,beta)
	    call verif_rot(ne,u0,v0,normt,ammo,nint,pint,nbdl,pbdl,cip)
c
            call save_temp('temp2.bin',ne,nm,time,u0,v0,h0,w0,z0,hb)
            call save_mesh('cm2.dat',nn,ne,in,xgr,ygr,scx,scy)
            call save_adap('raff2.bin',ne)
c
c------------------------------------------------------------
c end remeshing 
c
         endif
c
	 time0=time
         call put1in2(ne,nm,u0,v0,h0,w0,z0,hb,u2,v2,h2,w2,z2,hb2)

c        call checkcfl(nc,ne,nnmod,hbnew,in,dt,g0,xgr,ygr)
         call checkcfl2(ne,hb,u0,v0,in,dt,g0,xgr,ygr,
     &pgraph,ngradim,ngraph)
c
         open(25,file='oc.maille',access='append')
         write(25,23) time /86400.d0,nn,ne,dt
         close(25)
c
c------------------------------------------------------------
c end test on convergence
c------------------------------------------------------------
      endif
c
c------------------------------------------------------------
c end iteration
c------------------------------------------------------------
      enddo
c
      call save_temp('temp.bin',ne,nm,time,u0,v0,h0,w0,z0,hb)
c
21      format(1x,i7,(e13.6,1x),i3)
22      format(1x,5(1p,e13.6,1x))
23      format(1x,e13.6,1x,2(i5),1x,f9.3)
30    format(1x,60('-'))
31    format('iteration|    days   | output # |     TKE     |',
     &'     TPE')
c
c
      end
c
c
c*****************************************************************
c force the velocity field to be non normal to trhe boundary
c******************************************************************
c
      subroutine verif_rot(ne,u0,v0,normt,ammo,nint,pint,nbdl,pbdl,cip)
      implicit none
      include '../include/maille.inc'
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      integer ne
      double precision u0(nnmod,*),v0(nnmod,*),
     &            ammo(nnmod,nnmod),
     &            normt(2,*),
     &            b(nnmod),bu(nnmod),
     &            uk(nnmod),vk(nnmod),
     &            uh(nnmod,nedim),vh(nnmod,nedim),
     &            qx(nnmod,nedim),qy(nnmod,nedim),
     &            gx(nnmod),gy(nnmod),gn(nnmod),gt(nnmod),
     &            dnx,dny
      integer i,j,info,ik
      integer nint,pint(*),nbdl,pbdl(*),cip(2,0:nnmod)
      double precision amm(nnmod,nnmod),amml(nnmod,nnmod)
      common /matrices/ amm,amml
c
c
c rotation et projection pour eviter les erreurs d'interpolation
c sur les frontieres
c
	do i=1,ne
	   do j=1,nm
	      uk(j) = u0(j,i)
	      vk(j) = v0(j,i)
	      gx(j)=0.d0
	      gy(j)=0.d0
	   enddo
	   call cal_b(gx,uk,ammo,1.d0)
	   call cal_b(gy,vk,ammo,1.d0)
	   do j=1,nm
	      uh(j,i)  = gx(j)
	      vh(j,i)  = gy(j)
	   enddo
	enddo
c
c
c------------------------------------------------
c cas normal sans rotation
c------------------------------------------------
	do ik=1,nint
         i = pint(ik)
	   do j=1,nm
		qx(j,ik) = uh(j,i)
		qy(j,ik) = vh(j,i)
	   enddo
      enddo
c
c
      call dtrsm3(nm,nint,amm,qx)
      call dtrsm3(nm,nint,amm,qy)
c
	do ik=1,nint
         i = pint(ik)
	   do j=1,nm
		u0(j,i) = qx(j,ik)
		v0(j,i) = qy(j,ik)
	   enddo
	enddo
c
c
c------------------------------------------------
c si la cellule est au bord, rotation
c------------------------------------------------
	do ik=1,nbdl
         i=pbdl(ik)
c        write(*,*) 'tk',i,ik
	   dnx = normt(1,i)
	   dny = normt(2,i)
c
	   do j=1,nm
		gx(j) = uh(j,i)
		gy(j) = vh(j,i)
	   enddo
c
	   call rot_vec(nm,gn,gt,gx,gy,dnx,dny)
	   call reduc_gn(nc,gn,bu,cip)
c
	   do j=1,nm
		qx(j,ik) = bu(j)
		qy(j,ik) = gt(j)
	   enddo
	enddo
c------------------------------------------------
c
      call dtrsm3(nm-nc-1,nbdl,amml,qx)
      call dtrsm3(nm,     nbdl,amm ,qy)
c
c------------------------------------------------
c reduction et rotation inverse
c
	do ik=1,nbdl
         i=pbdl(ik)
	   do j=1,nm
		bu(j) = qx(j,ik)
		gt(j) = qy(j,ik)
	   enddo
c
	   dnx = normt(1,i)
	   dny = normt(2,i)
c
	   call reduc_inv_gn(nc,gn,bu,cip)
	   call rot_vec(nm,gx,gy,gn,gt,dnx,dny)
c
	   do j=1,nm
		u0(j,i) = gx(j)
		v0(j,i) = gy(j)
	   enddo
c
	enddo
c
	return
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
c******************************************************************
c
      subroutine save_temp(name,ne,nm,time,u0,v0,h0,w0,z0,hb)
      implicit none
c
      include '../include/sp.inc'
      integer nm,ne
      double precision u0(nnmod,*),v0(nnmod,*),h0(nnmod,*),
     &                 w0(nnmod,*),z0(nnmod,*),hb(nnmod,*)
      double precision time
      character*20 name
      integer i,j
      
      open(2,file=name,form='unformatted')
      write(2) time
      write(2) ((u0(j,i),j=1,nm),i=1,ne)
      write(2) ((v0(j,i),j=1,nm),i=1,ne)
      write(2) ((h0(j,i),j=1,nm),i=1,ne)
      write(2) ((w0(j,i),j=1,nm),i=1,ne)
      write(2) ((z0(j,i),j=1,nm),i=1,ne)
      write(2) ((hb(j,i),j=1,nm),i=1,ne)
      close(2)
      return
      end
c
c******************************************************************
c
      subroutine put1in2(ne,nm,u1,v1,h1,w1,z1,hb1,u2,v2,h2,w2,z2,hb2)
      implicit none
      include '../include/sp.inc'
      integer ne,nm
      integer i,j
      double precision 
     & u1(nnmod,*),v1(nnmod,*),h1(nnmod,*),
     & w1(nnmod,*),z1(nnmod,*),hb1(nnmod,*),
     & u2(nnmod,*),v2(nnmod,*),h2(nnmod,*),
     & w2(nnmod,*),z2(nnmod,*),hb2(nnmod,*)
     
      do i=1,ne
         do j=1,nm
            u2(j,i) =u1(j,i)
            v2(j,i) =v1(j,i)
            h2(j,i) =h1(j,i)
            w2(j,i) =w1(j,i)
            z2(j,i) =z1(j,i)
            hb2(j,i)=hb1(j,i)
         enddo
      enddo
      return
      end
c
