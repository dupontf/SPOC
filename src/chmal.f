c******************************************************************
      program chmaille
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
      character*10 name
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
      close(2)
      call read_mesh(nn,ne,in,xgr,ygr,scx,scy)
      call pointer_mesh(nn,ne,nf,in,itm,nmit,tmit,imit,pseg,tseg)
      nebd=0
      do i=1,nf
         if (tseg(2,i).eq.0) then
            nebd=nebd+1
            sbd(nebd)=i
         endif
      enddo
      call changemaille(nn,ne,nf,nebd,sbd,in,itm,tseg,pseg)
      name='cm.dat'
      call save_mesh(name,nn,ne,in,xgr,ygr,scx,scy)
      end
c
c**************************************************************
c elimination de triangles qui ont deux faces sur la frontiere
c rotation des reperes locales dans le triangle pour que
c l'axe des x soit tangeant a la frontiere
c
      subroutine changemaille(nn,ne,nf,nfbd,sbd,in,itm,tseg,pseg)
      implicit none
      include '../include/maille.inc'
      integer nn,ne,nf
      integer sbd(*),in(3,*),itm(3,*),tseg(2,*),pseg(2,*)
      integer nfbd,pebd(nedim)
      integer i,j,k,l,ns1,j1,j2,i1,i2,i3,i4,i5,i6,k1,k2,k3,k4,tk1,tk2
      integer pp(3),pf(3)
      do i=1,ne
         pebd(i)=0
      enddo
      do i=1,nfbd
         ns1=sbd(i)
         tk1 = tseg(1,ns1)
         pebd(tk1)=pebd(tk1)+1
      enddo
c
c
c si pebd(i)>=2 alors triangle -i- a deux faces sur la frontiere
c
      do i=1,ne
         if (pebd(i).gt.1) then
            write(*,*) 'tk',i
c
c recherche de l'arete partagee
c
            j=1
            ns1=itm(j,i)
            do while(tseg(2,ns1).eq.0)
              j=j+1
              ns1=itm(j,i)
            enddo
            tk1 = tseg(1,ns1)
            tk2 = tseg(2,ns1)
            j1=1
            do while(itm(j1,tk1).ne.ns1)
              j1=j1+1
            enddo
            j2=1
            do while(itm(j2,tk2).ne.ns1)
              j2=j2+1
            enddo
            call cit(3,j1,1,k1)
            call cit(3,j2,-1,k2)
            call cit(3,j2,1,k3)
            call cit(3,j1,-1,k4)
            i1 = in(k1,tk1)
            i2 = in(k2,tk2)
            i3 = in(k3,tk2)
            i4 = in(k4,tk1)
            in(k1,tk1)=i2
            in(k3,tk2)=i4
            itm(j1,tk1)=itm(k3,tk2)
            itm(j2,tk2)=itm(k1,tk1)
            itm(k1,tk1)=ns1
            itm(k3,tk2)=ns1
            pseg(1,ns1)=i2
            pseg(2,ns1)=i4
            if (i.eq.tk1) then
            tseg(1,itm(j2,tk2))=tk2
            else
            tseg(1,itm(j1,tk1))=tk1
            endif
         endif
      enddo
c
c
c rotation des aretes jusqu'a ce que l'arete 1 soit sur la frontiere
c
      do i=1,nfbd
         ns1=sbd(i)
         tk1 = tseg(1,ns1)
            j=1
            do while(itm(j,tk1).ne.ns1)
               j=j+1
            enddo
              write(*,*) i,j
            do k=1,3
              call cit(3,k,j-1,l)
              pp(k) = in (l,tk1)
              pf(k) = itm(l,tk1)
            enddo
            do k=1,3
              in (k,tk1)= pp(k)
              itm(k,tk1)= pf(k)
            enddo
      enddo
c
c
      return
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
