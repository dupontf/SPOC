c*******************************************************************
c routines pour le calcul des termes nonlineaires
c*******************************************************************
c
c
c**************************************************
c calcul des polynomes d'avances
c**************************************************
c
      subroutine cal_p2(xi,yi,pi2)
      implicit none
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision pi2(ng_max,nnmod),xi(*),yi(*),po,x,y
      integer i,j,kc1,kc2,i1
      external po
      i1 = 0
      do kc1=0,nc
	   do kc2=0,nc-kc1
	      i1 = i1 + 1
	      do i=1,ng
	         y = yi(i)
	         x = xi(i)
	         pi2(i,i1)=po(kc1,x)*po(kc2,y)
		enddo
	   enddo
      enddo
      return
      end
c
c
c*******************************************************************
c
      subroutine xytospec(b,mff)
      implicit none
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision xi(ng_max),yi(ng_max),wi(ng_max)
      common /gauss/ wi,xi,yi
      double precision b(nnmod),pi2(ng_max,nnmod)
      double precision mff(ng_max),ans
      integer tk,i,j,k,i1,kc1,kc2
      common /gauss2/ pi2
c
	i1=0
	do kc1=0,nc
	  do kc2=0,nc-kc1
	    i1=i1+1
	    call inttri2(mff,ans,i1)
	    b(i1)=ans
	  enddo
	enddo
c
	return
	end
c
c
c*******************************************************************
c
      subroutine spectoxy(mff,ak)
      implicit none
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision mff(ng_max),pi2(ng_max,nnmod),ak(nnmod)
      integer j,k,kc1,kc2,i1
      common /gauss2/ pi2
c
	 do j=1,ng
	     mff(j) = 0.d0
	 enddo
c
	 i1=0
	 do kc1=0,nc
	  do kc2=0,nc-kc1
	   i1=i1+1
	   do j=1,ng
	     mff(j) = mff(j) + ak(i1) * pi2(j,i1)
	   enddo
	  enddo
	 enddo
c
	return
	end
c
c
c*******************************************************************
c
      subroutine inttri2(mff,ans,i1)
      implicit none
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision wi(ng_max),xi(ng_max),yi(ng_max)
      common /gauss/ wi,xi,yi
      double precision mff(ng_max),ans
      double precision pi2(ng_max,nnmod)
      common /gauss2/ pi2
      integer i,i1
      ans=0.d0
      do i=1,ng
         ans = ans + wi(i) * mff(i) * pi2(i,i1)
      enddo
      return
      end
c
c*******************************************************************
c
      subroutine multiply(ng,mff,mff1,mff2)
      implicit none
      include '../include/sp.inc'
      integer ng
      double precision mff(ng_max),mff1(ng_max),mff2(ng_max)
      integer j,k
      do j=1,ng
	   mff(j) = mff1(j) * mff2(j)
      enddo
      return
      end	
c 
      subroutine divide(ng,mff,mff1,mff2)
      implicit none
      include '../include/sp.inc'
      integer ng
      double precision mff(ng_max),mff1(ng_max),mff2(ng_max)
      integer j,k
      do j=1,ng
	   mff(j) = mff1(j) / mff2(j)
      enddo
      return
      end	
c 
