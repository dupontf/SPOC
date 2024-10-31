!******************************************************************
! This programs changes the connectivity table so that 
! the boundary elements have their first and second vertices 
! on the boundary.
! For the elements having only one vertex on the boundary, 
! this vertex must be their first.
!******************************************************************
program chmaille

 use mesh
 implicit none

 integer :: i,j,k,face,n1,n2,face_start,istart,iend,flag,iroot
 character :: filename*20,root*20
 integer, allocatable :: pnodebdy(:,:),fflag(:),pdone(:)
 double precision :: xmin

      scx=1.d0
      scy=1.d0
      call read_mesh
      call pointer_mesh
      
      allocate (pnodebdy(2,nn))
      pnodebdy=0
      allocate (pdone(nf))
      pdone=0
      
      xmin=xgr(1)
      istart=1
      do i=1,nfbd
         face=sbd(i)
	 pdone(face)=1
	 n1 = pseg(1,face)
	 n2 = pseg(2,face)
	 pnodebdy(2,n1)=face
	 pnodebdy(1,n2)=face
	 if (xgr(n1).lt.xmin) then
	    xmin=xgr(n1)
	    istart=n1
	 endif
	 if (xgr(n2).lt.xmin) then
	    xmin=xgr(n2)
	    istart=n2
	 endif
      enddo
      
      allocate(fflag(nf))
	 n1=istart
	 write(*,*) 'starting from node',n1
	 write(*,*) 'enter next node and code'
	 read(*,*) iend,flag

         face=pnodebdy(2,n1)
	 pdone(face)=0
	 fflag(face)=flag
	 n1 = pseg(2,face)

      do while (istart.ne.iend)
        do while (n1.ne.iend)
         face=pnodebdy(2,n1)
	 pdone(face)=0
	 fflag(face)=flag
	 n1 = pseg(2,face)
        enddo
	 write(*,*) 'enter next node and code'
	 read(*,*) iend,flag
      enddo

        do while (n1.ne.iend)
         face=pnodebdy(2,n1)
	 pdone(face)=0
	 fflag(face)=flag
	 n1 = pseg(2,face)
        enddo

      

      open(1,file='cm.bel',status='unknown',action='write')
      i=0
	 n1=istart
         face=pnodebdy(2,n1)
	 n1 = pseg(1,face)
	 n2 = pseg(2,face)
	 i=i+1
	 write(1,*) i,n1,n2,fflag(face)
	 n1=n2
        do while (n1.ne.iend)
         face=pnodebdy(2,n1)
	 n1 = pseg(1,face)
	 n2 = pseg(2,face)
	 i=i+1
	 write(1,*) i,n1,n2,fflag(face)
	 n1=n2
        enddo
      
!
! detection of islands
!
      do i=1,nfbd
         face=sbd(i)
	 if (pdone(face).eq.1) then
	    n1 = pseg(1,face)
	    n2 = pseg(2,face)
	    write(1,*) i,n1,n2,1
	 endif
      enddo
      
      close(1)
      
      open(1,file='cm.ncd',status='unknown',action='write')
      do i=1,nn
         n1=pnodebdy(1,i)
         n2=pnodebdy(2,i)
	 if (n1+n2.gt.0) then
	    k = fflag(n1)+fflag(n2)
	    write(1,*) i,k
	 endif
      enddo
      close(1)
      
      
end program chmaille


