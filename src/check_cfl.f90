!******************************************************************
! check cfl
!******************************************************************
program check

 use mesh
 use gauss
 use graphic

	implicit none

      integer i,j,npermax
      double precision, allocatable :: u0(:,:),v0(:,:), h0(:,:), & 
      w0(:,:),z0(:,:),hb(:,:)

      double precision &
            dt,mu0,mb0,f0,beta,g0,hh0,time,& 
            yearmax,daysmax,nitemax,lx,ly,time0


!------------------------------------
! lecture des donnees
!------------------------------------
!
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
	read(2,*) hh0
	read(2,*) npermax,ngraph
	close(2)
!
!------------------------------------
! build the basis functions
!------------------------------------
!
      call read_gauss
!
!------------------------------------
! read the mesh and calculate the pointers
!------------------------------------
!
      call read_mesh
      call pointer_mesh
      call check_mesh
!
!------------------------------------
! read initial state and forcing
!
      allocate(u0(nm,ne),v0(nm,ne),h0(nm,ne),w0(nm,ne),z0(nm,ne),hb(nm,ne))
      open(1,file='init.bin',form='unformatted',status='old')
      read(1) time0
      read(1) u0
      read(1) v0
      read(1) h0
      read(1) w0
      read(1) z0
      read(1) hb
      close(1)

      write(*,*) 'STARTING TIME  :',time0/86400.d0

!      call checkcfl(hb,dt,g0)

      call cal_pgraph
      call checkcfl2(hb,u0,v0,dt,g0)
      
      ngraph=2
      deallocate(pgraph)
      call cal_pgraph
      call connec
      call write_gra(hb,'u',0)

end
