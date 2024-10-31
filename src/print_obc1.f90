!******************************************************************
program open_boundary_1

    use mesh
    use boundary_condition
    use gauss_bd

      implicit none

      integer i,j,k,tk
      integer cell,i1,i2,i3,face
      double precision ax,bx,ay,by,&
      xof,yof,x,y,x1,y1,x2,y2,pi,rad,c0,freq,af,bf
      double precision g0,hh0,lx,ly
      character*60 :: format1,format2,format3
 integer, allocatable :: pbdl_inv(:),bdl(:)
 double precision, allocatable :: cosng(:,:),sinng(:,:)
 double precision omega
      double precision, allocatable :: coeff(:,:)

!
!------------------------------------
! read input parameter file
!------------------------------------
!
        open(2,file='oc.inp',status='old')
        read(2,*) nc,ng,ngf,ng_bd,ngf_bd
        read(2,*) lx,ly
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) scx,scy
        read(2,*) g0
        read(2,*) hh0
        close(2)

      pi=4.d0*atan(1.d0)
      rad= lx
      freq=2.d0 * pi /rad
      c0 = sqrt(g0 * hh0)
      omega = 2.d0 * pi * c0 /rad
      
      call read_gauss_bd

      call read_mesh
      call pointer_mesh
      call check_mesh

      call readbel

      allocate(bdl(ne))
      do i=1,ne
	  bdl(i)=.false.
      enddo
      do i=1,nf
	   if (tseg(2,i).eq.0) then
	      bdl(tseg(1,i))=.true.
	   endif
      enddo
      

! ********* init of curved elements *************

      allocate(pbdl_inv(ne))
      call cal_pbdl_inv(bdl,pbdl_inv)
      allocate(coeff(2,nfbd))
      call readcoeff(coeff,pbdl_inv)


format1='60(e20.13,1x)'
format2='i4,1x,2(i5,1x),3x,i1,2(1x,i5)'
format3='4(e20.13,1x)'

      allocate(cosng(ngf,nobc1),sinng(ngf,nobc1))

      open(1,file='obc1.xy')
      write(1,*) omega,0.d0
      do i=1,nobc1
         face=pobc1(i)

	 cell=tseg(1,face)
         i1 = in(1,cell)
         i2 = in(2,cell)
         i3 = in(3,cell)
         ax = (xgr(i2)-xgr(i1))*0.5d0
         bx = (xgr(i3)-xgr(i1))*0.5d0
         ay = (ygr(i2)-ygr(i1))*0.5d0
         by = (ygr(i3)-ygr(i1))*0.5d0
         xof = xgr(i1)
         yof = ygr(i1)

	 tk = pbdl_inv(cell)
	 af = coeff(1,tk)
	 bf = coeff(2,tk)

           do j=1,ngf_bd
              x1 = tif_bd(j)
              x2 =  x1   - bf      * 0.5d0 * (x1**2 -1.d0)
              y2 = -1.d0 + (af+bf) * 0.5d0 * (x1**2 -1.d0)
              x = ax * (x2+1.d0) + bx * (y2+1.d0) + xof
              y = ay * (x2+1.d0) + by * (y2+1.d0) + yof
              cosng(j) = cos(freq * (x + y) / sqrt(2.d0) )
              sinng(j) = sin(freq * (x + y) / sqrt(2.d0)  )
	      write(1,*) x,y
           enddo
      enddo
      close(1)
      close(2)

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

