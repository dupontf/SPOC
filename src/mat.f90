!
!******************************************************************
! integration dans triangle
!******************************************************************
!
subroutine cal_mat(nnmod,ng,nc,mat,ff,xi,yi,wi)

	implicit none
	integer nc,nm
	integer nnmod,ng
	double precision mat(nnmod,nnmod),xi(*),yi(*),wi(*)
	double precision ans,ff
	integer kc1,kc2,kc3,kc4
	common /i4_cheb/ kc1,kc2,kc3,kc4
	integer i,j,i1,i2
	external ff
!
! integration dans triangle
!
	nm = (nc+1)*(nc+2)/2

	do j=1,nm
	do i=1,nm
	   mat(i,j)=0.d0
	enddo
	enddo

	i1=0
	  do kc1=0,nc
	  do kc2=0,nc-kc1
	    i1=i1+1
	    i2=0
	    do kc3=0,nc
	    do kc4=0,nc-kc3
	       i2=i2+1
	       call inttri(ng,ff,wi,xi,yi,ans)
	       mat(i1,i2)=ans
	    enddo
	    enddo
	  enddo
	  enddo

	return
	
end subroutine cal_mat


subroutine inttri(ng,ff,wi,xi,yi,ans)
	implicit none
	integer ng
	double precision ff,xi(*),yi(*),wi(*),ans,x,y,ansx
	external ff
	integer i,j
	ans=0.d0
	do j=1,ng
	  y = yi(j)
	  x = xi(j)
	  ans = ans + wi(j) * ff(x,y)
	enddo
	return
end subroutine inttri

!
!*************************************************************
! functions
!*************************************************************
!
function f2(x,y)
	implicit none
	integer kc1,kc2,kc3,kc4
	common /i4_cheb/ kc1,kc2,kc3,kc4
	double precision f2,x,y,po
	external po
	f2=po(kc1,x)*po(kc2,y)*po(kc3,x)*po(kc4,y)
	return
end function f2

function dfx(x,y)
	implicit none
	integer kc1,kc2,kc3,kc4
	common /i4_cheb/ kc1,kc2,kc3,kc4
	double precision dfx,x,y,po,dpo
	external po,dpo
	dfx= po(kc1,x)*po(kc2,y)*dpo(kc3,x)*po(kc4,y)
	return
end function dfx

function dfy(x,y)
	implicit none
	integer kc1,kc2,kc3,kc4
	common /i4_cheb/ kc1,kc2,kc3,kc4
	double precision dfy,x,y,po,dpo
	external po,dpo
	dfy= po(kc1,x)*po(kc2,y)*po(kc3,x)*dpo(kc4,y)
	return
end function dfy


!
!***********************************************************
! transforms the mass-matrix for boundary elements
!***********************************************************
!
	subroutine cal_matl(nn,nc,amm,amml)
	implicit none
	integer nn,nc,nm
	double precision amm(nn,*),amml(nn,*)
	integer i,j,k,l,i1,i2,j1,j2,i0,j0,ip
	external ip

	nm = (nc+1)*(nc+2)/2

	do i=0,nc-1
	   do j=0,nc-1-i,2
	      do k=0,nc-1
		   do l=0,nc-1-k,2
		 i1 = ip(nc-1,i)+j+1
		   i2 = ip(nc  ,i)+j+2
		 j1 = ip(nc-1,k)+l+1
		   j2 = ip(nc  ,k)+l+2
		   i0 = ip(nc  ,i)+1
		   j0 = ip(nc  ,k)+1
		 amml(i1,j1) = &
                  amm(i2,j2) + amm(i0,j2) + amm(i2,j0) + amm(i0,j0)
		 enddo
	      enddo
	   enddo
	enddo
	do i=0,nc-1
	   do j=1,nc-1-i,2
	      do k=0,nc-1
		   do l=0,nc-1-k,2
		 i1 = ip(nc-1,i)+j+1
		   i2 = ip(nc  ,i)+j+2
		 j1 = ip(nc-1,k)+l+1
		   j2 = ip(nc  ,k)+l+2
		   i0 = ip(nc  ,i)+1
		   j0 = ip(nc  ,k)+1
		 amml(i1,j1) = &
                  amm(i2,j2) - amm(i0,j2) + amm(i2,j0) - amm(i0,j0)
		 enddo
	      enddo
	   enddo
	enddo
	do i=0,nc-1
	   do j=0,nc-1-i,2
	      do k=0,nc-1
		   do l=1,nc-1-k,2
		 i1 = ip(nc-1,i)+j+1
		   i2 = ip(nc  ,i)+j+2
		 j1 = ip(nc-1,k)+l+1
		   j2 = ip(nc  ,k)+l+2
		   i0 = ip(nc  ,i)+1
		   j0 = ip(nc  ,k)+1
		 amml(i1,j1) = &
                  amm(i2,j2) + amm(i0,j2) - amm(i2,j0) - amm(i0,j0)
		 enddo
	      enddo
	   enddo
	enddo
	do i=0,nc-1
	   do j=1,nc-1-i,2
	      do k=0,nc-1
		   do l=1,nc-1-k,2
		 i1 = ip(nc-1,i)+j+1
		   i2 = ip(nc  ,i)+j+2
		 j1 = ip(nc-1,k)+l+1
		   j2 = ip(nc  ,k)+l+2
		   i0 = ip(nc  ,i)+1
		   j0 = ip(nc  ,k)+1
		 amml(i1,j1) = &
                  amm(i2,j2) - amm(i0,j2) - amm(i2,j0) + amm(i0,j0)
		 enddo
	      enddo
	   enddo
	enddo

	return
	end
!
!*****************************************************************
! function IP computes the index in matrix for transformation
!*****************************************************************
!
	function ip(nc,i)
	integer ip,nc,i,j
	ip = 0
	j=0
	do while(j.lt.i)
	ip = ip + nc+1-j
	j=j+1
	enddo
	return
	end

