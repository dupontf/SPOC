module utilities

use mesh
use gauss

implicit none

      integer, allocatable :: ssit(:,:),its(:,:)
      double precision, allocatable :: norm(:,:)
      logical, allocatable :: bdl(:)
      double precision, allocatable :: vect(:,:),normt(:,:)
      integer, allocatable :: cip(:,:)

 contains

!
!***************************************************************
! calcul de la normal pour les faces des triangles
!***************************************************************
!
!
subroutine cal_rotnf
      implicit none
      integer i,ne1,ne2
      double precision nx,ny,d

      allocate(norm(2,nf))

      do i=1,nf
	 ne1=pseg(1,i)
	 ne2=pseg(2,i)
	 nx= (ygr(ne2)-ygr(ne1))
	 ny=-(xgr(ne2)-xgr(ne1))
	 norm(1,i)= nx
	 norm(2,i)= ny
      enddo

      return
end subroutine cal_rotnf

!
!
!******************************************************************
!
subroutine cal_its
	implicit none
	integer i,j,tk

       allocate(ssit(2,nf),its(3,ne))

        do i=1,nf
         ssit(1,i)=0
         ssit(2,i)=0
        enddo

	do i=1,nf
	   tk = tseg(1,i)
	   j = 1
	   do while (itm(j,tk).ne.i.and.j.lt.4)
	    j = j + 1
	   enddo
	   if (j.eq.4) then
	      write(*,*) 'problem in cal_its. bad connectivity table'
		write(*,*) 'face',i,' triangle',tk,' index in triangle',j
		stop
	   endif
	   its(j,tk) = 1
         ssit(1,i)=j
	   tk = tseg(2,i)
	   if (tk.gt.0) then
	   j = 1
	   do while (itm(j,tk).ne.i)
	    j = j + 1
	   enddo
	   its(j,tk) = -1
         ssit(2,i)=j
	   endif
	enddo

	return
end subroutine cal_its

!
!******************************************************************
! defines other pointers relevant to the spectral element computation
!******************************************************************
!
subroutine calpre
      implicit none

      integer nbdl
      double precision ax,bx,ay,by,deta,detai,dnx,dny,dis
      integer i,ns1,i1,i2,i3

      allocate(bdl(ne))
      allocate(vect(5,ne),normt(2,ne))

      do i=1,ne
	  bdl(i)=.false.
      enddo

      do i=1,nf
	   if (tseg(2,i).eq.0) then
	      bdl(tseg(1,i))=.true.
	   endif
      enddo

      do i=1,ne
         i1 = in(1,i)
         i2 = in(2,i)
         i3 = in(3,i)
         ax = (xgr(i2)-xgr(i1))*0.5d0
         bx = (xgr(i3)-xgr(i1))*0.5d0
         ay = (ygr(i2)-ygr(i1))*0.5d0
         by = (ygr(i3)-ygr(i1))*0.5d0
         deta = ax * by - ay * bx
         detai = 1.d0/deta
         vect(1,i)=ax
         vect(2,i)=bx
         vect(3,i)=ay
         vect(4,i)=by
         vect(5,i)=detai
         ns1 = itm(1,i)
         dnx = norm(1,ns1)
         dny = norm(2,ns1)
         dis = sqrt(dnx**2+dny**2)
         dnx = dnx / dis
         dny = dny / dis
         normt(1,i)=dnx
         normt(2,i)=dny
      enddo

      nbdl=0
      do i=1,ne
         if (bdl(i)) then
            nbdl=nbdl+1
         endif
      enddo
      write(*,*) 'nbdl:',nbdl
 
      return
end subroutine calpre

!
!
!******************************************************************
!
	subroutine coriolis(fc,f0,beta)
	implicit none
	double precision fc(nm,ne),f0,beta
	integer cell,n1,n2,n3
	
	fc=0.d0
	
	do cell=1,ne
	   n1 = in(1,cell)
	   n2 = in(2,cell)
	   n3 = in(3,cell)
	   fc(1,cell) = f0 + 0.5d0 * beta *( ygr(n3) + ygr(n2))
	   fc(2,cell) =      0.5d0 * beta * ( ygr(n3) - ygr(n1))
	   fc(nc+2,cell) =   0.5d0 * beta * ( ygr(n2) - ygr(n1))
	enddo
	return
	end

!
!*****************************************************************
! function CALIP computes the index in matrix for transformation
!*****************************************************************
!
subroutine calip

	implicit none
	integer ip,i,j,nc2

      allocate(cip(2,0:nc))

	do i=0,nc
	  ip = 0
	  j=0
	  do while(j.lt.i)
	    ip = ip + nc+1-j
	    j=j+1
	  enddo
	  cip(1,i)=ip
	enddo

	nc2=nc-1
	do i=0,nc2
	  ip = 0
	  j=0
	  do while(j.lt.i)
	    ip = ip + nc2+1-j
	    j=j+1
	  enddo
	  cip(2,i)=ip
	enddo

	return
end


end module utilities
