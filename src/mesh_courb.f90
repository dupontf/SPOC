module mesh_curved

use mesh
use utilities
use boundary_condition

implicit none

 double precision, allocatable :: coeff(:,:)
 integer, allocatable :: pbdl_inv(:)
  
  
 contains


!**************************************************************
      subroutine cal_pbdl_inv
      implicit none
      integer cell,j,tk
      
      allocate(pbdl_inv(ne))
      
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

      subroutine readcoeff
      implicit none

      integer i,j,k,n1,n2,ind,i1,i2,i3,i4,ele
      character*100 line
      integer pnodebdy(2,nn)
      double precision af,bf
      
      allocate(coeff(2,nfbd))

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


!******************************************************************

end module mesh_curved
