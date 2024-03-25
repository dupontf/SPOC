c*****************************************************************
c routines handling the ouput the variables for GMT use
c*****************************************************************
      subroutine cal_pgraph(nc,n,pgraph,ngradim)
      implicit none
      include '../include/sp.inc'
      integer kc1,kc2,nc,nm,n,ngradim
      double precision pgraph(nnmod,ngradim,ngradim),x,y,po
      external po
      integer k,l,i1
      do k =1,n
         do l =1,n-k+1
           y = dfloat(2*k-n-1)/dfloat(n-1)
           x = dfloat(2*l-n-1)/dfloat(n-1)
           i1=0
           do kc1=0,nc
            do kc2=0,nc-kc1
             i1=i1+1
             pgraph(i1,l,k) = po(kc1,x)*po(kc2,y)
            enddo
           enddo
         enddo
      enddo
      return
      end
c
c
c----------------------------------------------
c calcul de la super connectivite
c----------------------------------------------
c
      subroutine connec(nne,ne,ngr)
      implicit none
      include '../include/maille.inc'
      include '../include/graph.inc'
      integer ne,ngr,nm,nne,i,j,k,l
      integer in(3,2*nmg*nedim)
      nm = ngr*(ngr+1)/2
      nne=0
      do i=1,ne
         l = 0
         do k=1,ngr-1
            do j=1,ngr-k-1
               nne = nne + 1
               in(1,nne) = nm*(i-1) + l +j
               in(2,nne) = nm*(i-1) + l +j+1
               in(3,nne) = nm*(i-1) + l + ngr -k + 1 +j
               nne = nne + 1
               in(1,nne) = nm*(i-1) + l +j+1
               in(2,nne) = nm*(i-1) + l + ngr -k + 1 +j+1
               in(3,nne) = nm*(i-1) + l + ngr -k + 1 +j
            enddo
            j = ngr - k
               nne = nne + 1
               in(1,nne) = nm*(i-1) + l +j
               in(2,nne) = nm*(i-1) + l +j+1
               in(3,nne) = nm*(i-1) + l + ngr -k + 1 +j
               l = l + ngr -k + 1
         enddo
      enddo
c 
      open(21,file='connec.dat')
      do i=1,nne
        write(21,10) in(1,i)-1,in(2,i)-1,in(3,i)-1
      enddo
      close(21)
c
c------------------------------------
c write messages to standart output
c
      write(*,*)
      write(*,30)
      write(*,*) 'number of sub-elements for graphic output:',nne
      write(*,*) 'the graphic output files are named: *.gra'
c
10    format(4(i6,1x))
30    format(1x,60('-'))
      return
      end
c
c
c*****************************************************************
c visu pour les elements courbes
c*****************************************************************
c
      subroutine write_gra(ne,h,name,nper,nm,ngr,xgr,ygr,in,pgraph)
      implicit none
      include '../include/maille.inc'
      include '../include/sp.inc'
      include '../include/graph.inc'
      integer ngr,nm,ne,n,nbdl,nper
      integer kc1,kc2
      double precision h(nnmod,*),b(nnmod),xgr(*),ygr(*)
      double precision x,y,sol,pp,
     &                 pgraph(nnmod,ngradim,ngradim)
      integer i,j,k,l,i1,i2,i3,ind,tk
      integer in(3,*)
      double precision ax,bx,ay,by,x1,y1,xof,yof,af,bf,x2,y2
      integer ideb,ifin
      character name*1,anum*3,fbin*8
c
 22   format(i3)
      write(anum,22) nper
      ideb=1
      do while (anum(ideb:ideb).eq.' ') 
         ideb=ideb+1
      enddo
c         
      fbin = name//anum(ideb:3)//'.gra'
      ifin=9-ideb

c	write(*,*) fbin
      open(1,file=fbin)
      do i=1,ne
         do j=1,nm
	    b(j)=h(j,i)
	 enddo
         i1 = in(1,i)
         i2 = in(2,i)
         i3 = in(3,i)
         ax = (xgr(i2)-xgr(i1))*0.5d0
         bx = (xgr(i3)-xgr(i1))*0.5d0
         ay = (ygr(i2)-ygr(i1))*0.5d0
         by = (ygr(i3)-ygr(i1))*0.5d0
         xof = xgr(i1)
         yof = ygr(i1)
c
         do k =1,ngr
         do l =1,ngr-k+1
           y1 = dfloat(2*k-ngr-1)/dfloat(ngr-1)
           x1 = dfloat(2*l-ngr-1)/dfloat(ngr-1)
           x = xof + ax * (x1+1.d0) + bx * (y1+1.d0)
           y = yof + ay * (x1+1.d0) + by * (y1+1.d0)
           sol = 0.d0
           do i1=1,nm
             pp = pgraph(i1,l,k)
             sol = sol + b(i1) * pgraph(i1,l,k)
           enddo
           write(1,20) x,y,sol
         enddo
         enddo
c
      enddo
      close(1)

c
20    format(3(e10.3,1x))
      return
      end
c
