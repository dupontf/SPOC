!******************************************************************
program print_xy

 use mesh
 
      implicit none

      integer cell,i,j,k
      integer nc,ng
      integer i1,i2,i3,face,cell,tk
      double precision ax,ay,bx,by,xof,yof,x,y,x1,y1,x2,y2,af,bf
      double precision, allocatable :: xi(:),yi(:), &
      coeff(:,:)
 integer, allocatable :: pbdl_inv(:),bdl(:)

!------------------------------------
! read input parameter file
!------------------------------------

        open(2,file='oc.inp',status='old')
        read(2,*) nc
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

      allocate(bdl(ne))
      do i=1,ne
	  bdl(i)=.false.
      enddo
      do i=1,nf
	   if (tseg(2,i).eq.0) then
	      bdl(tseg(1,i))=.true.
	   endif
      enddo
      

! ********* init of curved elements *************

      allocate(pbdl_inv(ne))
      call cal_pbdl_inv(bdl,pbdl_inv)
      allocate(coeff(2,nfbd))
      call readcoeff(coeff,pbdl_inv)
      


      ng=(nc+2)*(nc+1)/2
      write(*,*) nc,ng
      allocate(xi(ng),yi(ng))
      ng=0
      do j=0,nc
         do i=0,nc-j
	    ng=ng+1
            yi(ng) = dfloat(2*j-nc)/dfloat(nc)
            xi(ng) = dfloat(2*i-nc)/dfloat(nc)
!	    write(*,'3(i3,1x),3x,2(f10.7,1x)') ng,i,j,xi(ng),yi(ng)
	 enddo
      enddo
      
      open(1,file='cm.xy')
      tk=0
      do cell=1,ne

       i1 = in(1,cell)
       i2 = in(2,cell)
       i3 = in(3,cell)
       ax = (xgr(i2)-xgr(i1))*0.5d0
       bx = (xgr(i3)-xgr(i1))*0.5d0
       ay = (ygr(i2)-ygr(i1))*0.5d0
       by = (ygr(i3)-ygr(i1))*0.5d0
       xof = xgr(i1)
       yof = ygr(i1)

       tk = pbdl_inv(cell)

       if (bdl(cell)) then
       
	 af = coeff(1,tk)
	 bf = coeff(2,tk)

        do j =1,ng
           y1 = xi(j)
           x1 = yi(j)
           x2 = x1 -     bf  * 0.5d0 * (x1**2 + x1 + x1*y1 + y1)
           y2 = y1 + (af+bf) * 0.5d0 * (x1**2 + x1 + x1*y1 + y1)
           x = xof + ax * (x2+1.d0) + bx * (y2+1.d0)
           y = yof + ay * (x2+1.d0) + by * (y2+1.d0)
           write(1,*) x,y
        enddo

        else

        do j =1,ng
           y1 = xi(j)
           x1 = yi(j)
           x = xof + ax * (x1+1.d0) + bx * (y1+1.d0)
           y = yof + ay * (x1+1.d0) + by * (y1+1.d0)
           write(1,*) x,y
        enddo
	
       endif

      enddo
      close(1)

      end

!**************************************************************
      subroutine cal_pbdl_inv(bdl,pbdl_inv)
      use mesh
      implicit none
      integer cell,j,tk
      integer pbdl_inv(*),bdl(*)
            
      tk=0
      
      do cell=1,ne
         pbdl_inv(cell)=0
	 if (bdl(cell)) then 
	    tk=tk+1
	    pbdl_inv(cell)=tk
	 endif
      enddo

      return
      end subroutine cal_pbdl_inv


!----------------------------------------
! read the curvature file
!----------------------------------------

      subroutine readcoeff(coeff,pbdl_inv)
      use mesh
      implicit none

      integer i,j,k,n1,n2,ind,i1,i2,i3,i4,ele
      character*100 line
      integer pnodebdy(2,nn)
      double precision af,bf,coeff(2,*)
      integer pbdl_inv(*)
      

!----------------------------------------
! construct the sbd and pbnodebdy pointers

      do i=1,nn
         pnodebdy(1,i)=0
         pnodebdy(2,i)=0
      enddo
      do i=1,nf
         if (tseg(2,i).eq.0) then
	    n1=pseg(1,i)
	    n2=pseg(2,i)
            pnodebdy(1,n2)=i
            pnodebdy(2,n1)=i
	 endif
      enddo      


      open(2,file='coef.dat',status='old')
102   read(2,'a80',end=103) line
      read(line,*) j,i1,i2,af,bf
      n1=pnodebdy(2,i1)
      n2=pnodebdy(1,i2)
      if (n1.ne.n2) then
	  write(*,*) 'in read coeff:'
	  write(*,*) 'pb with face defined by nodes:',i1,i2
	  stop
      endif
      ele=tseg(1,n1)
      i=pbdl_inv(ele)
      coeff(1,i) = af
      coeff(2,i) = bf
      goto 102
103   close(2)

      return
      end subroutine readcoeff

