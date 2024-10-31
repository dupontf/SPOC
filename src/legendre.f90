!
!*************************************************
! compute Legendre polynomial at x location
!*************************************************
!
function po(n,x)
	implicit none
	integer n,i
	double precision p0,p1,x,a0,a2,a3,po

	p0=1.d0
	p1=x

	if (n.eq.0) po=p0
	if (n.eq.1) then 
	   po=p1
	else
	   i=1
	   do while (i.lt.n)
	      i=i+1
	      a0 = dfloat(i)
		a2 = dfloat(2*i-1)
		a3 = dfloat(i-1)
		po = (a2*x*p1-a3*p0)/a0
		p0=p1
		p1=po
	   enddo
	endif

	return
end function po

!
!*************************************************
! compute first derivative of Legendre polynomial at x location
!*************************************************
!
function dpo(n,x)
	implicit none
	integer n,i
	double precision p0,p1,x,a0,a2,a3,p2,dpo,dp0,dp1

	p0=1.d0
	p1=x
	dp0=0.d0
	dp1=1.d0

	if (n.eq.0) dpo=dp0
	if (n.eq.1) then 
	   dpo=dp1
	else
	   i=1
	   do while (i.lt.n)
	      i = i + 1
	      a0 = dfloat(i)
		a2 = dfloat(2*i-1)
		a3 = dfloat(i-1)
		dpo = (a2*x*dp1+a2*p1-a3*dp0)/a0
		p2 = (a2*x*p1-a3*p0)/a0
		p0=p1
		p1=p2
		dp0=dp1
		dp1=dpo
	   enddo
	endif

	return
end function dpo
