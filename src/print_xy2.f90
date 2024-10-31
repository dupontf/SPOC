!******************************************************************
program print_xy

 use mesh
 use gauss
 
      implicit none

      integer i,j,k
      integer i1,i2,i3,face,cell,tk
      double precision ax,ay,bx,by,xof,yof,x,y,x1,y1,x2,y2,af,bf

!------------------------------------
! read input parameter file
!------------------------------------

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
        close(2)
      
      call read_mesh
      call pointer_mesh
      call check_mesh

      ng=ng*(ng+1)/2
      allocate(xi(ng),yi(ng))
      ng=0
      do j=0,nc
         do i=0,nc-j
	    ng=ng+1
            yi(ng) = dfloat(2*j-nc)/dfloat(nc)
            xi(ng) = dfloat(2*i-nc)/dfloat(nc)
	    write(*,*) i,j,xi(ng),yi(ng)
	 enddo
      enddo
      
      open(1,file='cm.xy')
      tk=0
      do i=1,ne

       i1 = in(1,i)
       i2 = in(2,i)
       i3 = in(3,i)
       ax = (xgr(i2)-xgr(i1))*0.5d0
       bx = (xgr(i3)-xgr(i1))*0.5d0
       ay = (ygr(i2)-ygr(i1))*0.5d0
       by = (ygr(i3)-ygr(i1))*0.5d0
       xof = xgr(i1)
       yof = ygr(i1)

        do j =1,ng
           y1 = xi(j)
           x1 = yi(j)
           x = xof + ax * (x1+1.d0) + bx * (y1+1.d0)
           y = yof + ay * (x1+1.d0) + by * (y1+1.d0)
           write(1,*) x,y
        enddo

      enddo
      close(1)

      end
