c******************************************************************
c
	subroutine cal_b(nnmod,nc,b,ak,mat,coeff)
	implicit none
c
	integer nnmod
	double precision b(*),ak(*),mat(nnmod,*),ans
	double precision coeff
	integer nc,nm
	integer i,j
c
c integration dans triangle
c
	nm = (nc+1)*(nc+2)/2
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
c
c******************************************************************
