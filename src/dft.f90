module harmonic

 use mesh
 use gauss
 use boundary_bd
 
 implicit none
 
 double complex, allocatable :: f(:,:,:)
 double precision start_anal,anal_interval
 double precision, parameter :: localpi=3.14159265358979323844d0
 integer :: ncycle=0,itelocal=0

 contains


!***********************************************************************
!
! FUNDY4 FIXED SUBROUTINE WRITTEN BY PROFESSOR DANIEL R. LYNCH
!
 FUNCTION PHASELAGD(ZZ)
!
! COMPUTES THE NEGATIVE ARGUMENT OF A COMPLEX NUMBER Z, IN DEGREES.
! RESULT IS REPORTED IN THE RANGE 0. TO 360. 
! 
      double COMPLEX Z,ZZ
      double precision PHASELAGD,FACTOR

      Z = ZZ
      FACTOR=180.d0/localPI
      PHASELAGD = 0.0d0
      IF(ABS(Z).GT.0.d0)  PHASELAGD = -FACTOR*ATAN2(AIMAG(Z), dble(Z))
      IF(PHASELAGD.LT.0.d0) PHASELAGD = PHASELAGD + 360.d0  
      RETURN
 END FUNCTION PHASELAGD
!***********************************************************************


!***********************************************************************
 subroutine init_dft

 open(1,name='anal.inp',action='read',status='old')
  read(1,*) start_anal
  read(1,*) anal_interval
 close(1)
 write(*,*) 'Start anal time',start_anal
 write(*,*) 'Time between two outputs',anal_interval
 return
 end subroutine init_dft


!***********************************************************************
 SUBROUTINE DFTana(dt,time,u,v,h)
!-----------------------------------------------------------------------
! This subroutine calculates the first harmonic time-domain component of
!   nsamp time-domain samples of 3 variables using a DFT.
! 
!-----------------------------------------------------------------------
!

    double precision, intent(in) :: dt,time
    double precision, intent(in) ::  u(nm,ne),v(nm,ne),h(nm,ne)

    character ccycle*7,filename*80
    double complex carg
    integer cell,j
    double precision fac,phaselag
    double precision ur(nm,ne),ui(nm,ne)
    

! Check array sizes and initialize F
    if(time.le.start_anal.and.time+dt.gt.start_anal) then
       allocate(f(3,nm,ne))
       do cell=1,ne
        do j=1,nm
         f(1,j,cell) = cmplx(0.d0,0.d0)
         f(2,j,cell) = cmplx(0.d0,0.d0)
         f(3,j,cell) = cmplx(0.d0,0.d0)
        enddo
       enddo
    endif

! Perform summation to build F for j=0,nsamp-1

 if(time.ge.start_anal) then
      fac=-time*omega
      carg=CMPLX(0.0,fac)
      itelocal = itelocal + 1

       do cell=1,ne
        do j=1,nm
         f(1,j,cell) = f(1,j,cell)+u(j,cell)*EXP(carg)
         f(2,j,cell) = f(2,j,cell)+v(j,cell)*EXP(carg)
         f(3,j,cell) = f(3,j,cell)+h(j,cell)*EXP(carg)
        enddo
       enddo


!---------------------------------------------------------------------
! Write result

 if(mod(time,anal_interval).lt.dt.and.time.gt.start_anal+dt) then

  ncycle = ncycle + 1
  write(ccycle,'(i7)')ncycle+1000000
  fac=2.d0/dble(itelocal)
  write(*,*)'Harmonic output:'

  ur=fac*dble(f(1,:,:))
  ui=fac*aimag(f(1,:,:))

  filename='anal_u'//ccycle(5:7)//'.bin'
 write(*,*) filename
 open(1,file=filename,form='unformatted',action='write')
  write(1) time
  write(1) ur
  write(1) ui
 close(1)

  ur=fac*dble(f(2,:,:))
  ui=fac*aimag(f(2,:,:))

  filename='anal_v'//ccycle(5:7)//'.bin'
 write(*,*) filename
 open(1,file=filename,form='unformatted',action='write')
  write(1) time
  write(1) ur
  write(1) ui
 close(1)

  ur=fac*dble(f(3,:,:))
  ui=fac*aimag(f(3,:,:))

  filename='anal_h'//ccycle(5:7)//'.bin'
 write(*,*) filename
 open(1,file=filename,form='unformatted',action='write')
  write(1) time
  write(1) ur
  write(1) ui
 close(1)

 endif
 endif
    
    
 return
 end SUBROUTINE DFTana

end module harmonic
