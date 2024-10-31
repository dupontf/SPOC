!******************************************************************
program base2legendre

 use mesh
 use gauss
 use gauss_bd

      implicit none

      integer i,j,k
      integer i1,i2,i3,face,cell,tk
      double precision ax,ay,bx,by,xof,yof,x,y,x1,y1,x2,y2,af,bf,sol
      double precision, allocatable :: h0(:)
      double precision fsol,pob2
      external fsol,pob2
      double precision g0,hh0,lx,ly
 integer, allocatable :: bdl(:)

!------------------------------------
! read input parameter file
!------------------------------------

        open(2,file='oc.inp',status='old')
        read(2,*) nc,ng,ngf,ng_bd,ngf_bd
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) scx,scy
        close(2)
      
      call read_gauss
      call read_gauss_bd

      do j=1,ng
         xi(j)=(1.d0+xi(j))*0.5d0
         yi(j)=(1.d0+yi(j))*0.5d0
      enddo
      
      do j=1,ng_bd
         xi_bd(j)=(1.d0+xi_bd(j))*0.5d0
         yi_bd(j)=(1.d0+yi_bd(j))*0.5d0
      enddo
      
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


 
      allocate(h0(nm))
      open(1,file='cm.depth')
      open(2,file='cm.depth2')

      do cell=1,ne

         do j=1,nm
	     read(1,*) h0(j)
	 enddo

       if (bdl(cell)) then
       
       do j =1,ng_bd
           y1 = xi_bd(j)
           x1 = yi_bd(j)
	   sol=fsol(nc,x1,y1,h0)
           write(2,*) sol
        enddo

       else

        do j =1,ng
           y1 = xi(j)
           x1 = yi(j)
	   sol=fsol(nc,x1,y1,h0)
           write(2,*) sol
         enddo
	 
	endif
	 
      enddo
      close(1)
      close(2)

      end

	function fsol(nc,x,y,b)
	implicit none
	integer nc,nm
	integer j
	double precision fsol,x,y,b(*),pob2
	external pob2
	nm = (nc+1)*(nc+2)/2
	fsol=0.d0
	do j=1,nm
	   fsol = fsol + b(j) * pob2(nc,j,x,y)
	enddo
	return
	end



!****************************************************************
! construction d'une base cardinale
!****************************************************************
	function pob2(nc,i,x,y)
	implicit none
	double precision x,y,pob2,p1,p2,p3,d0,po0,x0,y0
	integer nc,i,k,m,l,i1
	external p1,p2,p3

	d0 = 1.d0/dfloat(nc)
!
! cal coordonnes dans base tronquee (m,k)
!
	i1 = nc
	k = 0
	m = 0
	do while (i-m.gt.0)
	m = m + nc + 1 - k
	k = k + 1
	enddo
	l = i - m + nc + 1 - k
	m = k - 1

	pob2 =1.d0
	po0=1.d0
	x0=dfloat(l)*d0
	y0=dfloat(m)*d0
	do i1 = 0,l-1
	   pob2 = pob2 * (p2(x)  -dfloat(i1)*d0)
	   po0= po0* (p2(x0)  -dfloat(i1)*d0)
	enddo
	do i1 = 0,nc-m-1-l
	   pob2 = pob2 * (p1(x,y)-dfloat(i1)*d0)
	   po0= po0* (p1(x0,y0)-dfloat(i1)*d0)
	enddo
	do i1 = 0,m-1
	   pob2 = pob2 * (p3(y )  -dfloat(i1)*d0)
	   po0= po0* (p3(y0)  -dfloat(i1)*d0)
	enddo
	pob2 = pob2/po0

	end



	function p1(x,y)
	implicit none
	double precision x,y,p1
	p1 = 1.d0 -x -y
	end
	function p2(x)
	implicit none
	double precision x,p2
	p2 = x
	end
	function p3(y)
	implicit none
	double precision y,p3
	p3 = y
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

