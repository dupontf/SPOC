!******************************************************************
! initialization of the variables for the discontinuous spectral 
! element model (SPOC)
!******************************************************************
program init

 use mesh
 use gauss
 use graphic
 use boundary_condition
 use utilities
 use gauss_bd
 use boundary_bd
 
 implicit none

      double precision, allocatable ::&
       u(:,:), v(:,:), h(:,:), w(:,:), z(:,:), hb(:,:)
        
      integer i,j,k,ind,k1,k2,ne1,ne2,cote,i1,i2
      integer info,toto,permax
      double precision mu0,mb0,f0,beta,g0,hh0,lx,ly,lvent,angle,pi

      double precision fvaguex,fvaguey,ampl,lf,x0,y0, &
                      ftourbu,ftourbv,ftourbh,fvent1,fvent2,fvent3, &
                      vmax,b0,fond0,fond1,fond2,amplvent,fsinus,fkelvin, &
                      fvagueg
      external fvaguex,fvaguey,ftourbu,ftourbv,ftourbh,fvent1, &
                fvent2,fvent3, &
                fond0,fond1,fond2,fsinus,fkelvin,&
                      fvagueg

      double precision f2,dfx,dfy
      external f2,dfx,dfy
 
      double precision a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
      character*60 :: format1
      character :: filename*60
      character*10 typewindx,typewindy,typefond

      pi=4.d0*atan(1.d0)
!
!------------------------------------
! read input parameter file
!------------------------------------
!
	open(2,file='oc.inp',status='old')
	read(2,*) nc,ng,ngf,ng_bd,ngf_bd
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
	read(2,*) permax,ngraph
	read(2,*) 
	read(2,*) amplvent,lvent
	read(2,'(a10)') typewindx
	read(2,'(a10)') typewindy
	read(2,'(a10)') typefond
	close(2)
!
!------------------------------------
! build the basis functions
!------------------------------------
!
      call read_gauss
      allocate(amm(nm,nm))
      call cal_mat(nm,ng,nc,amm,f2,xi,yi,wi)
      call cal_p2

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

      call cal_rotnf
      call calpre

! ********* init of curved elements *************

      call cal_pbdl_inv
      call readcoeff
!      coeff=0.d0
      
      call read_gauss_bd
      call cal_p2_bd
      call cal_courb_jac
      call cal_mat_courb
      
!
!------------------------------------
! compute the initial state and the forcing
!------------------------------------
!
!
      allocate(u(nm,ne),v(nm,ne),h(nm,ne),w(nm,ne),z(nm,ne),hb(nm,ne))
	      u=0.d0
	      v=0.d0
	      h=0.d0
	      w=0.d0
	      z=0.d0

!-------------------------------------------------
! initial state
!
format1='1x,60(''-'')'

      write(*,*)
      write(*,format1)
      write(*,*) 'CHOOSE INITIAL STATE'
      write(*,format1)
      write(*,*)
	write(*,*) 'rest              =0'
	write(*,*) 'wave in x         =1'
	write(*,*) 'wave in y         =2'
	write(*,*) 'eddy exp QG       =3'
	write(*,*) 'eddy sinus        =4'
	write(*,*) 'Kelvin wave       =5'
	write(*,*) 'propagating wave  =6'
	write(*,*) 'read from file    =7'
	read(*,*) toto
	write(*,*) 'you chose:',toto

      if (toto.eq.1) then
	   write(*,*) 'amplitude (m)'
	   read(*,*) a1
	   write(*,*) 'wave length (m)'
	   read(*,*) a2
	   call calh_c(fvaguex,h)
      elseif (toto.eq.2) then
	   write(*,*) 'amplitude (m)'
	   read(*,*) a1
	   write(*,*) 'wave length (m)'
	   read(*,*) a2
	   call calh_c(fvaguey,h)
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
	 call calh_c(ftourbu,u)
         call calh_c(ftourbv,v)
         a1 = EXP(0.5d0) * VMAX * f0 / ( G0 * SQRT( 2.d0 * B0 ) )
	 write(*,*) 'amplitude of the eddy in (m):', a1
	 call calh_c(ftourbh,h)
      elseif (toto.eq.4) then
	 write(*,*) 'ampl (m)'
	 read(*,*) a1
	 write(*,*) 'radius (m)'
	 read(*,*) a2
         a3 = 0.d0
         a4 = 0.d0
	 call calh_c(fsinus,h)
!	 do i=1,ne
!	     do j=1,nm
!	        u(j,i)=h(j,i)
!		  v(j,i)=h(j,i)
!	     enddo
!	 enddo
      elseif (toto.eq.5) then
	 write(*,*) 'amplitude (m)'
	 read(*,*) a1
	 write(*,*) 'amplitude of the Kelvin wave in (m):', a1
         a3 = 0.d0
!         a5 = 1.d-10
         a5 = 1.d0/300.d3**2
	 write(*,*) 'length (m)'
	 read(*,*) ampl
         write(*,*) 'length of the Kelvin wave in (m):',ampl
	 a5 = 1.d0 / ampl**2
         a4 = -ly*0.5d0
         a2 = f0 / sqrt ( g0 * hh0 )
	 call calh_c(fkelvin,h)
         a1 = a1 * g0 / sqrt ( g0 * hh0 )
         write(*,*) 'max speed in the Kelvin wave in (m/s):',a1
	 call calh_c(fkelvin,u)
      elseif (toto.eq.6) then
	   write(*,*) 'amplitude (m)'
	   read(*,*) a1
	   write(*,*) 'wave length (m)'
	   read(*,*) a2
	   write(*,*) 'wave angle (deg)'
	   read(*,*) angle
	   a3 = cos(angle*pi/180.d0)
	   a4 = sin(angle*pi/180.d0)
	 write(*,*) 'amplitude (m)'
	 call calh_c(fvagueg,h)
         a1 = a1 * g0 / sqrt ( g0 * hh0 )
         write(*,*) 'max speed in the gravity wave in (m/s):',a1
	 call calh_c(fvagueg,v)
	 do i=1,ne
	    do j=1,nm
	       u(j,i) = v(j,i) * a3
	       v(j,i) = v(j,i) * a4
	    enddo
	 enddo
      elseif (toto.eq.7) then
           write(*,*) 'type filename for u'
           read(*,fmt='(a60)') filename
	   write(*,*) filename
	   call calh_c2(filename,u)
           write(*,*) 'type filename for v'
           read(*,fmt='(a60)') filename
	   write(*,*) filename
	   call calh_c2(filename,v)
           write(*,*) 'type filename for h'
           read(*,fmt='(a60)') filename
	   write(*,*) filename
	   call calh_c2(filename,h)
      endif


      write(*,*)
      write(*,format1)
      write(*,*) 'CHOOSE WIND FORCING'
      write(*,format1)
      write(*,*) 'zero wind forcing              =0'
      write(*,*) 'sinus wind in x varying in y   =1'
      write(*,*) 'wind in x arctan varying in  y =2'
      write(*,*) 'cosinus wind in x varying in y =3'
      if (typewindx.eq.'zero') toto=0
      if (typewindx.eq.'fvent1') toto=1
      if (typewindx.eq.'fvent2') toto=2
      if (typewindx.eq.'fvent3') toto=3
      write(*,*) 'you chose:',toto

      a1 = amplvent
      a2 = lvent
      if (toto.eq.1) then
	  call calh_c(fvent1,w)
      elseif (toto.eq.2) then
	  call calh_c(fvent2,w)
      elseif (toto.eq.3) then
	  call calh_c(fvent3,w)
      endif


      write(*,*)
      write(*,format1)
      write(*,*) 'CHOOSE BATHYMETRY'
      write(*,format1)
      write(*,*)
      write(*,*) 'constant hh0      =0'
      write(*,*) 'haidvogel channel =1'
      write(*,*) typefond
      toto=10
      if (typefond.eq.'constant'.or.typefond.eq.'zero') toto=0
      if (typefond.eq.'fond1') toto=1
      if (typefond.eq.'fond2') toto=2
      write(*,*) 'you chose:',toto

	if (toto.eq.1) then
         a1 = lx
         a2 = ly
	   call calh_c(fond1,hb)
	elseif (toto.eq.2) then
         a1 = lx
         a2 = ly
	   call calh_c(fond2,hb)
        elseif (toto.eq.0) then
           do i=1,ne
              hb(1,i)=hh0
              do j=2,nm
                 hb(j,i)=0.d0
              enddo
           enddo
	else 
	 call calh_c2(typefond,hb)
	endif

!
!-------------- graphic output --------------
!

      call cal_pgraph
      call connec
      call write_gra(u,'u',0)
      call write_gra(v,'v',0)
      call write_gra(h,'h',0)
      call write_gra(hb,'z',0)
      write(*,*) 'u/v/h0.gra graphic output files have been created'


      open(1,file='init.bin',form='unformatted')
      write(1) 0.d0
      write(1) ((u(j,i),j=1,nm),i=1,ne)
      write(1) ((v(j,i),j=1,nm),i=1,ne)
      write(1) ((h(j,i),j=1,nm),i=1,ne)
      write(1) ((w(j,i),j=1,nm),i=1,ne)
      write(1) ((z(j,i),j=1,nm),i=1,ne)
      write(1) ((hb(j,i),j=1,nm),i=1,ne)
      close(1)


      write(*,*)
      write(*,format1)
      write(*,*) 'INIT.BIN has been written'
      write(*,format1)
      write(*,*)

      end
