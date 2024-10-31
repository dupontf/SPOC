!******************************************************************
! equations shallow waters
! spectral discontinuous elements
! rk4 integration in time
! conditions frontieres pour u 
! viscosite 
! friction
! plan beta 
!******************************************************************
program extract_bin

 use mesh
 use gauss
 use mesh_curved
 use utilities
 use graphic
 
 implicit none

!
! modes spectraux
!
      double precision, allocatable :: h(:,:)

      integer i,j,k
      integer nper,npermax
      character*60 :: format1,format2,format3
      character filename*80,varname*1
      integer ivar,nvar
      double precision time

!
!------------------------------------
! read input parameter file
!------------------------------------
!
	open(2,file='oc.inp',status='old')
	read(2,*) nc,ng,ngf
	read(2,*) 
	read(2,*) 
	read(2,*) 
	read(2,*) 
	read(2,*) 
	read(2,*) 
	read(2,*) 
	read(2,*) scx,scy
	read(2,*) 
	read(2,*) 
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

      call cal_rotnf
      call cal_its
      call calpre

! ********* init of curved elements *************

      call cal_pbdl_inv
      call readcoeff

!
!------------------------------------
! read initial state and forcing
!
      allocate(h(nm,ne))

      write(*,*) 'enter bin filename'
      read(*,'(a80)') filename
      write(*,*) 'enter variable index in file'
      read(*,*) nvar

      open(1,file=filename,form='unformatted',status='old',err=110)
      read(1) time
      do i=1,nvar-1
       read(1,err=100)
      enddo
      read(1,err=100) h
      close(1)

      format1='(''time in file (days)'',f10.2)'
      write(*,format1) time/86400.d0

!
!------------------------------------------
! compute energy and output graphic
!
      call cal_pgraph
      call connec

      write(*,*) 'enter variable name for ascii output (1 char)'
      read(*,'(a1)') varname
      write(*,*) 'enter index for ascii output'
      read(*,*) nper

      call write_gra(h,varname,nper)

      stop

100   write(*,*) 'wrong variable index in file',filename
      stop
110   write(*,*) 'cannot open file',filename
end program extract_bin
