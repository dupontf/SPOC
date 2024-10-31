!******************************************************************
! equations shallow waters
! spectral discontinuous elements
! rk4 integration in time
! conditions frontieres pour u 
! viscosite 
! friction
! plan beta 
!******************************************************************
program shwater

 use mesh
 use gauss
 use graphic
 use boundary_condition
 use utilities
 
 implicit none

!
! modes spectraux
!
      double precision, allocatable :: &
       u0(:,:), v0(:,:), h0(:,:), &
       u1(:,:), v1(:,:), h1(:,:), &
       du1(:,:), dv1(:,:), dh1(:,:), &
       du2(:,:), dv2(:,:), dh2(:,:), &
       du3(:,:), dv3(:,:), dh3(:,:), &
       du4(:,:), dv4(:,:), dh4(:,:), &
       w0(:,:), z0(:,:), hb(:,:), fc(:,:)

      external f2,dfx,dfy,ip,fx,fy
      double precision f2,dfx,dfy,fx,fy

      integer i,j,k,l,i1,i2,i3,ns1
      integer nite,ite,info,per,nper,per2,yearmax,daysmax,nitemax,npermax
      double precision :: &
         dt,mu0,mb0,f0,beta,g0,hh0,time,ke,pe,mass,pe0, &
         lx,ly,time0,amplvent,ax,ay,bx,by,deta,detai,&
         dnx,dny,dis,dt6,dt3,dt2

      character*60 :: format1,format2,format3

!
!------------------------------------
! read input parameter file
!------------------------------------
!
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
	read(2,*) npermax,ngraph
	read(2,*) per2
	read(2,*) amplvent
	close(2)
!
!----------------------------------------------
! write message

format1='1x,60(''-'')'
format2='1x,5(1p,e13.6,1x)'
format3='i8,'' |  '',f8.2,'' |  '',i6,''  |   '',1P,E8.1,''  |   '',1P,E8.1'

      write(*,*)
      write(*,format1)
      write(*,*) 'READ INPUT PARAMETERS'
      write(*,format1)
      write(*,*)
      write(*,*) 'time step              (s)   :',dt
      write(*,*) 'horizontal viscosity   (m2/s):',mu0
      write(*,*) 'linear bottom friction (m/s) :',mb0
      write(*,*) 'Coriolis f0            (1/s) :',f0
      write(*,*) 'Coriolis beta          (1/ms):',beta
      write(*,*) 'reduced gravity        (m/s2):',g0
      write(*,*) 'depth at rest          (m)   :',hh0

    nite=nitemax+int(daysmax*86400.d0/dt)+int(yearmax*365.d0*86400.d0/dt)
      per=int(nite/npermax)
      dt6=dt/6.d0
      dt3=dt/3.d0
      dt2=dt/2.d0

      write(*,*) 'number of iterations                           :',nite
      write(*,*) 'number of iterations between two graphic output:',per

!
!------------------------------------
! build the basis functions
!------------------------------------
!
      call read_gauss
      call cal_p2
      call cal_pico
      call calip
     
!
! build the matrices
!
      allocate(amm(nm,nm), agx(nm,nm), agy(nm,nm), amml(nm,nm))
      call cal_mat(nm,ng,nc,amm,f2,xi,yi,wi)
      call cal_mat(nm,ng,nc,agx,dfx,xi,yi,wi)
      call cal_mat(nm,ng,nc,agy,dfy,xi,yi,wi)
      call cal_matl(nm,nc,amm,amml)
      allocate(ammo(nm,nm))
      do j=1,nm
         do i=1,nm
	    ammo(i,j)=amm(i,j)
         enddo
      enddo
!
! LU Decomposition of the mass matrices
!
      call dpotrf('U',nm,amm,nm,info)
      if (info.ne.0) write(*,*) 'pb avec matrice'
      call dpotrf('U',nm-nc-1,amml,nm,info)
      if (info.ne.0) write(*,*) 'pb avec matrice'
!
!------------------------------------
! read the mesh and calculate the pointers
!------------------------------------
!
      call read_mesh
      call pointer_mesh
      call check_mesh

      call readbel
      call periodic

      call cal_rotnf
      call cal_its
      call calpre
      call openbc1_read
      
      
!
!--------------------------------------------------------
! compute coriolis parameter
!
      allocate(fc(nm,ne))
      fc=0.d0
      call coriolis(fc,f0,beta)

!
!------------------------------------
! read initial state and forcing
!
      allocate(u0(nm,ne),v0(nm,ne),h0(nm,ne),w0(nm,ne),z0(nm,ne),hb(nm,ne))
      allocate(u1(nm,ne),v1(nm,ne),h1(nm,ne))
      allocate(du1(nm,ne),dv1(nm,ne),dh1(nm,ne))
      allocate(du2(nm,ne),dv2(nm,ne),dh2(nm,ne))
      allocate(du3(nm,ne),dv3(nm,ne),dh3(nm,ne))
      allocate(du4(nm,ne),dv4(nm,ne),dh4(nm,ne))

      open(1,file='init.bin',form='unformatted',status='old')
      read(1) time0
      read(1) u0
      read(1) v0
      read(1) h0
      read(1) w0
      read(1) z0
      read(1) hb
      close(1)

!
!-------------------------------------------------------------
! defines time iterators
!
      time = time0
      ite=0
      nper=0
!
!------------------------------------------
! compute energy and output graphic
!
      call cal_pgraph
      call connec
      call write_gra(u0,'u',nper)
      call write_gra(v0,'v',nper)
      call write_gra(h0,'h',nper)
      call energy(u0,v0,h0,g0,hb,ke,pe,mass)
      open(24,file='oc.energ')
      write(24,*) 'days, kinetic, potential energies'
      write(24,format2) time /86400.d0,ke,pe,mass
      close(24)
!
!-------------------------------------------------------------
! write message
!
      write(*,format1)
      write(*,*) 'iteration|    days   | output # |     TKE     |',&
      '     TPE'
      write(*,format1)
      write(*,format3) ite,time/86400.d0,nper,ke,pe-pe0
!
!------------------------------------
! start iterating in time
!------------------------------------
!
	do while (ite.lt.nite)
	ite = ite + 1
!
!-------------------------------rk4------------------------------
!
	  call cal_var(u0,v0,h0,du1,dv1,dh1,mu0,mb0,fc,w0,z0,g0,hb,time)
	    u1 = u0 + dt2 * du1
	    v1 = v0 + dt2 * dv1
	    h1 = h0 + dt2 * dh1

	  call cal_var(u1,v1,h1,du2,dv2,dh2,mu0,mb0,fc,w0,z0,g0,hb,&
	   time+dt*0.5d0)
	    u1 = u0 + dt2 * du2
	    v1 = v0 + dt2 * dv2
	    h1 = h0 + dt2 * dh2

	  call cal_var(u1,v1,h1,du3,dv3,dh3,mu0,mb0,fc,w0,z0,g0,hb,&
	   time+dt*0.5d0)
	    u1 = u0 + dt * du3
	    v1 = v0 + dt * dv3
	    h1 = h0 + dt * dh3

	  call cal_var(u1,v1,h1,du4,dv4,dh4,mu0,mb0,fc,w0,z0,g0,hb,&
	   time+dt)

            u0 = u0 + dt6 * (du1+du4) + dt3 * (du2+du3)  
            v0 = v0 + dt6 * (dv1+dv4) + dt3 * (dv2+dv3)  
            h0 = h0 + dt6 * (dh1+dh4) + dt3 * (dh2+dh3)  

	time = time0 + dt * dfloat(ite)
!
!-------------------------------------------------------
! energy
!
      if (mod(ite,per2).eq.0) then
         call energy(u0,v0,h0,g0,hb,ke,pe,mass)
         open(24,file='oc.energ',status='old',position='append')
         write(24,format2) time /86400.d0,ke,pe-pe0,mass
         close(24)
         call testnan(ke,time)
      endif
!
!-------------------------------------------------------
! visualisation
!
       if (mod(ite,per).eq.0) then
         nper=int(ite/per)
         call write_gra(u0,'u',nper)
         call write_gra(v0,'v',nper)
         call write_gra(h0,'h',nper)
         write(*,format3) ite,time/86400.d0,nper,ke,pe-pe0
	   open(2,file='temp.bin',form='unformatted')
	   write(2) time
	   write(2) ((u0(j,i),j=1,nm),i=1,ne)
	   write(2) ((v0(j,i),j=1,nm),i=1,ne)
	   write(2) ((h0(j,i),j=1,nm),i=1,ne)
	   write(2) ((w0(j,i),j=1,nm),i=1,ne)
	   write(2) ((z0(j,i),j=1,nm),i=1,ne)
	   write(2) ((hb(j,i),j=1,nm),i=1,ne)
	   close(2)
       endif

      enddo
!
!-------------------------------------------------------
! final output
!
      open(2,file='temp.bin',form='unformatted')
      write(2) time
      write(2) ((u0(j,i),j=1,nm),i=1,ne)
      write(2) ((v0(j,i),j=1,nm),i=1,ne)
      write(2) ((h0(j,i),j=1,nm),i=1,ne)
      write(2) ((w0(j,i),j=1,nm),i=1,ne)
      write(2) ((z0(j,i),j=1,nm),i=1,ne)
      write(2) ((hb(j,i),j=1,nm),i=1,ne)
      close(2)



      stop
end program shwater
!
!******************************************************************
!
	subroutine testnan(ke,time)
	double precision ke,time
      character*20 name
	integer ideb,ifin

      write(name,'1p,e13.6') ke

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

