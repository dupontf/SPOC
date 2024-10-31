!******************************************************************
program get_bdy

 use mesh
 
      implicit none

      integer i,j,k
      logical, allocatable :: pnbd(:)
      integer nnbd
      integer, allocatable :: p2f(:,:)

!------------------------------------
! read input parameter file
!------------------------------------

      open(2,file='oc.inp',status='old')
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) scx,scy
      close(1)

      call read_mesh
      call pointer_mesh
      call check_mesh
      
      allocate(pnbd(nn))
      pnbd=.false.
      allocate(p2f(2,nn))

      do i=1,nf
         if (tseg(2,i).eq.0) then
	    k=pseg(1,i)
	    p2f(2,k)=i
	    pnbd(k)=.true.
	    k=pseg(2,i)
	    pnbd(k)=.true.
	    p2f(1,k)=i
         endif
      enddo

      nnbd=0
      do i=1,nn
         if (pnbd(i)) then
	    nnbd=nnbd+1
	 endif
      enddo
      write(*,*) nnbd
      
      open(1,file='dessin.xy')
	 write(1,'a1') '>'

! find the first point on the bdy
10    i=1
      do while (.not.pnbd(i))
	 i=i+1
      enddo

      do while (nnbd.gt.0)
         do while (pnbd(i))
	    pnbd(i)=.false.
	    nnbd=nnbd-1
	    write(1,*) xgr(i),ygr(i)
            k=p2f(2,i)
	    i=pseg(2,k)
	 enddo
	 write(1,'a1') '>'
	 goto 10
      enddo
      end

