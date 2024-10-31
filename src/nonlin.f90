!*******************************************************************
! routines pour le calcul des termes nonlineaires
!*******************************************************************
!
!*******************************************************************
!
      subroutine xytospec(b,mff)
      use gauss
      implicit none
      double precision b(*)
      double precision mff(ng),ans
      integer tk,i,j,k,i1,kc1,kc2
 
	i1=0
	do kc1=0,nc
	  do kc2=0,nc-kc1
	    i1=i1+1
	    call inttri2(mff,ans,i1)
	    b(i1)=ans
	  enddo
	enddo

	return
	end

!
!*******************************************************************
!
      subroutine spectoxy(mff,ak)
      use gauss
      implicit none
      double precision mff(ng),ak(nm)
      integer j,k,kc1,kc2,i1

	     mff = 0.d0

	 i1=0
	 do kc1=0,nc
	  do kc2=0,nc-kc1
	   i1=i1+1
	   do j=1,ng
	     mff(j) = mff(j) + ak(i1) * pi2(j,i1)
	   enddo
	  enddo
	 enddo

	return
	end

!
!*******************************************************************
!
      subroutine inttri2(mff,ans,i1)
      use gauss
      double precision mff(ng),ans
      integer i,i1
      ans=0.d0
      do i=1,ng
         ans = ans + wi(i) * mff(i) * pi2(i,i1)
      enddo
      return
      end
!
!*******************************************************************
!
      subroutine multiply(ng,mff,mff1,mff2)
      implicit none
      integer ng
      double precision mff(*),mff1(*),mff2(*)
      integer j,k
      do j=1,ng
	   mff(j) = mff1(j) * mff2(j)
      enddo
      return
      end	

      subroutine divide(ng,mff,mff1,mff2)
      implicit none
      integer ng
      double precision mff(*),mff1(*),mff2(*)
      integer j,k
      do j=1,ng
	   mff(j) = mff1(j) / mff2(j)
      enddo
      return
      end	
 
