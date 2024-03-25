c
c******************************************************************
c integration de flux aux faces
c******************************************************************
c
      subroutine cal_flu(nf,tseg,fbdy,a0,ssit)
      implicit none
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      integer nf
      double precision pico(ngf_max,nnmod,3),a0(nnmod,*)
      double precision fbdy(ngf_max,*)
      integer tseg(2,*),ssit(2,*)
      integer i,j,k1,k2,i1
      integer cote
      common /gauss3/ pico
c
c$omp parallel share(a0,fbdy,tseg,ssit,pico,nm),
c$&            local(i,j,k1,k2,i1,cote)
c$omp do schedule(runtime)
	do i=1,nf
c
	  do j=1,ngf
	     fbdy(j,i)=0.d0
	  enddo
c
	   k1 = tseg(1,i)
	   k2 = tseg(2,i)
c
	   if (k2.eq.0) then
c
c la face est sur la frontiere: seul triangle 1
c
	  do i1=1,nm
	    do j=1,ngf
	       fbdy(j,i) =  fbdy(j,i) + a0(i1,k1) * pico(j,i1,1)
	    enddo
	  enddo
c
	  else
c
c triangle 1
c
        cote = ssit(1,i)
	  do i1=1,nm
	    do j=1,ngf
	       fbdy(j,i) =  fbdy(j,i) + a0(i1,k1) * pico(j,i1,cote)
	    enddo
	  enddo
c
c triangle 2
c
	  cote = ssit(2,i)
	  do i1=1,nm
	    do j=1,ngf
	       fbdy(j,i) =  fbdy(j,i) + a0(i1,k2) * pico(ngf-j+1,i1,cote)
	    enddo
	  enddo
	  do j=1,ngf
	     fbdy(j,i) = 0.5d0 * fbdy(j,i)
	  enddo
c
	 endif
c
	enddo
c$omp end do
c$omp end parallel
c
	return
	end
c
