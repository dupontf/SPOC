c******************************************************************
c initialization of the variables for the discontinuous spectral 
c element model (SPOC)
c******************************************************************
      program ch2
      implicit none
c
C
C 2-D array declaration
C
      INCLUDE '../include/maille.inc'
      integer nn,ne,nf,nebd
      INTEGER IN(3,NEDIM),tseg(2,nfdim),pseg(2,nfdim),
     &itm(3,nedim),nmit(nndim),tmit(NFTRDIM*NNDIM),
     & imit(NFTRDIM*NNDIM),ssit(2,nfdim),sbd(nbedim),
     & its(3,nedim)
      double precision scx,scy,xgr(nndim),ygr(nndim)
      double precision xgr2(nndim),ygr2(nndim)
      logical nobd(nndim),pco(nndim)
      integer i,j,k,l,ind,n1,n2
      integer nd
      parameter(nd=100)
      integer npw,npe
      integer pw(nd),pe(nd)
      integer pfw(nd),pfe(nd),ppf(nfdim)
      integer india
      external india
      double precision xmin,xmax,eps,rot,pi
      integer toto
      integer nfp,pfp(2,nd)
      integer pnodebdy(2,nbedim)
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
c
c------------------------------------
c read the mesh and calculate the pointers
c------------------------------------
c
      call read_mesh(nn,ne,in,xgr,ygr,scx,scy)
      call pointer_mesh(nn,ne,nf,in,itm,nmit,tmit,imit,pseg,tseg)
      call check_mesh(ne,nf,in,itm,tseg)
      do i=1,nn
         pnodebdy(1,i)=0
         pnodebdy(2,i)=0
      enddo
      do i=1,nf
         if (tseg(2,i).eq.0) then
	    n1 = pseg(1,i)
	    n2 = pseg(2,i)
	    pnodebdy(2,n1)=i
	    pnodebdy(1,n2)=i
	 endif
      enddo
c
c--------------------------------------------
c define eps,pi
c
      eps=1.e-6*scx
      pi=4.d0*atan(1.d0)
c
c------------------------------------------------
c
      write(*,*)
      write(*,30)
      write(*,*) 'Choose the periodicity: 1=West/East 2=North/South'
      read(*,*) toto
      if (toto.eq.1) then
          rot=0
	  call findnodeper(nn,xgr,npw,npe,pw,pe,eps)
          call tagface(npw,npe,pw,pe,nfp,pfp,pnodebdy,tseg,xgr,ygr)
      elseif (toto.eq.2) then
          rot = pi * 0.5d0
	  call rotate(nn,xgr,ygr,xgr2,ygr2,rot)
	  call findnodeper(nn,xgr2,npw,npe,pw,pe,eps)
          call tagface(npw,npe,pw,pe,nfp,pfp,pnodebdy,tseg,xgr2,ygr2)
	  rot = - rot
	  call rotate(nn,xgr2,ygr2,xgr,ygr,rot)
      else
        write(*,*) 'no periodic condition choosen'
	stop
      endif
c
c----------------------------------------
c print the periodic points
c
      write(*,*) 'western points (or fake west)'
      do i=1,npw
        write(*,*) pw(i),xgr(pw(i)),ygr(pw(i))
      enddo
      write(*,*) 'eastern point (or fake east)'
      do i=1,npe
        write(*,*) pe(i),xgr(pe(i)),ygr(pe(i))
      enddo
c
c----------------------------------------
c write the bel file and save the mesh
c----------------------------------------
c
      open(2,file='cm.bel')
      k=0
      do i=1,nfp
         k=k+1
         n1=pfp(1,i)
         n2=pfp(2,i)
         write(2,31) k,pseg(1,n1),pseg(2,n1),6,pseg(1,n2),pseg(2,n2)
         k=k+1
         n2=pfp(2,i)
         write(2,31) k,pseg(1,n2),pseg(2,n2),6,pseg(1,n1),pseg(2,n1)
      enddo
      close(2)
      name='cm.dat'
      write(*,*) name
      call save_mesh(name,nn,ne,in,xgr,ygr,scx,scy)

c      
      write(*,*)
      write(*,30)
      write(*,*) 'End of periodic condition program'
      write(*,30)
      write(*,*)
c
20 	format(20(f13.10,x))
21 	format(x,3(i3,x),2(f13.10,x))
30    format(1x,60('-'))
31    format(i4,1x,2(i5,1x),3x,i1,2(1x,i5))
c
      end
c
      subroutine tagface(npw,npe,pw,pe,nfp,pfp,pnodebdy,tseg,xgr,ygr)
      implicit none
      integer npw,npe,pw(*),pe(*),nfp,pfp(2,*),pnodebdy(2,*),tseg(2,*)
      integer i,j,k,i1,i2,n1,n2
      double precision xgr(*),ygr(*),buf
      
      if (npw.ne.npe) then
         write(*,*) 'probleme avec les points de frontiere periodique'
         stop
      endif

c
c---classer par ordre de y croissant les points sur les frontieres est et ouest
c
      do i=1,npw-1
         n1 = pw(i)
         do j=i+1,npw
            n2=pw(j)
            if (ygr(n1).gt.ygr(n2)) then
               pw(j) = n1
               pw(i) = n2
                  n1 = pw(i)
            endif
         enddo
      enddo
      do i=1,npe-1
         n1 = pe(i)
         do j=i+1,npe
            n2 = pe(j)
            if (ygr(n1).gt.ygr(n2)) then
               pe(j) = n1
               pe(i) = n2
               n1 = pe(i)
            endif
         enddo
      enddo
c
c------------------------------------------------
c make sure that the periodic nodes have the same y
c
      do i=1,npw
         n1 = pw(i)
         n2 = pe(i)
         buf=0.5d0*(ygr(n1)+ygr(n2))
	 ygr(n1)=buf
	 ygr(n2)=buf
      enddo
c
c
c---marquer les faces sur la frontiere ouest et est
c
      nfp=0
      do i=1,npw-1
         n1=pw(i)
	 n2=pe(i)
	 i1=pnodebdy(1,n1)
	 i2=pnodebdy(2,n2)
	 nfp=nfp+1
	 pfp(1,nfp)=i1
	 pfp(2,nfp)=i2
      enddo
      return
      end
c
c------------------------------------
c find west-east nodes of the mesh
c------------------------------------
c
      subroutine findnodeper(nn,xgr,npw,npe,pw,pe,eps)
      implicit none
      integer nn,npe,npw
      integer i
      integer pw(*),pe(*)
      double precision xgr(*),xmin,xmax,eps
c
      xmin=xgr(1)
      xmax=xgr(1)
      do i=2,nn
         xmin=min(xmin,xgr(i))
         xmax=max(xmax,xgr(i))
      enddo
c
      do i=1,nn
         if (abs(xgr(i)-xmin).lt.eps) then
	    npw=npw+1
	    pw(npw)=i
	 endif
         if (abs(xgr(i)-xmax).lt.eps) then
	    npe=npe+1
	    pe(npe)=i
	 endif
      enddo
      return
      end
c      
c
c------------------------------------
c rotate the mesh
c------------------------------------
c
      subroutine rotate(nn,xgr,ygr,xgr2,ygr2,rot)
      implicit none
      integer nn,i
      double precision xgr(*),ygr(*),xgr2(*),ygr2(*)
      double precision rot,cosa,sina
c
      cosa = cos(rot)
      sina = sin(rot)
      do i=1,nn
         xgr2(i) = cosa * xgr(i) + sina * ygr(i)
         ygr2(i) = cosa * ygr(i) - sina * xgr(i)
      enddo
      return
      end
      
