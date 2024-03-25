c******************************************************************
      program get_bdy
      implicit none
c
C
C 2-D array declaration
C
      include '../include/maille.inc'
      include '../include/sp.inc'
      integer nn,ne,nf,nebd
      INTEGER IN(3,NEDIM),tseg(2,nfdim),pseg(2,nfdim),
     &itm(3,nedim),nmit(nndim),tmit(NFTRDIM*NNDIM),
     & imit(NFTRDIM*NNDIM),sbd(NBEDIM)
      integer i,j,k
      double precision xgr(nndim),ygr(nndim),scx,scy
      logical pnbd(nndim)
      integer nnbd
      integer p2f(2,nndim)
c
c------------------------------------
c read input parameter file
c------------------------------------
c
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
      call read_mesh(nn,ne,in,xgr,ygr,scx,scy)
      call pointer_mesh(nn,ne,nf,in,itm,nmit,tmit,imit,pseg,tseg)

      nebd=0
      do i=1,nn
         pnbd(i)=.false.
      enddo

      do i=1,nf
         if (tseg(2,i).eq.0) then
            nebd=nebd+1
            sbd(nebd)=i
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
	 write(1,20)
c find the first point on the bdy
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
	 write(1,20)
	 goto 10
      enddo
20    format('>')
      end
c
c
c        
      subroutine cit(n,k,np,k1)
      integer n,k,np,k1
c
      k1=k+np
      if(k1.gt.n) k1=k1-n
      if(k1.lt.1) k1=k1+n
c
      return
      end
