c******************************************************************
c
      subroutine cal_b(b,ak,mat,coeff)
      implicit none
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision b(nnmod),ak(nnmod),mat(nnmod,nnmod),ans
      double precision coeff
      integer i,j
c
c integration dans triangle
c
	do i=1,nm
	   ans = 0.d0
	   do j=1,nm
	      ans = ans + mat(j,i) * ak(j)
	   enddo
	   b(i) = b(i) + ans * coeff
	enddo
c
	return
	end
c******************************************************************

