!******************************************************************
program open_boundary_1

    use mesh
    use gauss
    use boundary_condition

      implicit none

      integer i,j,k
      integer i1,i2,face
      double precision ax,ay,xof,yof,x,x1,y1,pi,rad,c0,freq
      double precision g0,hh0,lx,ly
      character*60 :: format1,format2,format3
!
!------------------------------------
! read input parameter file
!------------------------------------
!
        open(2,file='oc.inp',status='old')
        read(2,*) nc,ng,ngf
        read(2,*) lx,ly
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) scx,scy
        read(2,*) g0
        read(2,*) hh0
        close(2)

      pi=4.d0*atan(1.d0)
      rad= lx
      freq=2.d0 * pi /rad
      c0 = sqrt(g0 * hh0)
      omega = 2.d0 * pi * c0 /rad
      
      call read_gauss

      call read_mesh
      call pointer_mesh
      call check_mesh

      call readbel
!      call periodic

format1='60(e20.13,1x)'
format2='i4,1x,2(i5,1x),3x,i1,2(1x,i5)'
format3='4(e20.13,1x)'

      allocate(cosng(ngf,nobc1),sinng(ngf,nobc1))

      open(1,file='obc1.dat')
      open(2,file='obc1.xy')
      write(1,*) omega,0.d0
      do i=1,nobc1
         face=pobc1(i)
         i1=pseg(1,face)
         i2=pseg(2,face)
         ax = (xgr(i2)-xgr(i1))*0.5d0
         ay = (ygr(i2)-ygr(i1))*0.5d0
         xof = xgr(i1)
         yof = ygr(i1)

           do j=1,ngf
              x = tif(j)
              x1 = ax * (x+1.d0) + xof
              y1 = ay * (x+1.d0) + yof
              cosng(j) = cos(freq * (x1 + y1) / sqrt(2.d0) )
              sinng(j) = sin(freq * (x1 + y1) / sqrt(2.d0)  )
!              cosng(j) = cos(freq * y1 )
!              sinng(j) = sin(freq * y1 )
	      write(1,format3) cosng(j),sinng(j)
	      write(2,format3) x1,y1
           enddo
!	 write(1,format1) (cosng(j),j=1,ngf)
!	 write(1,format1) (sinng(j),j=1,ngf)
      enddo
      close(1)
      close(2)

      end
