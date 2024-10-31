module graphic

      use gauss
      use mesh
      
      implicit none

      integer :: ngraph,nne
      double precision, allocatable :: pgraph(:,:,:),x,y,po

 contains

!*****************************************************************
! routines handling the ouput the variables for GMT use
!*****************************************************************
      subroutine cal_pgraph
      implicit none
      integer kc1,kc2
      double precision x,y,po
      external po
      integer k,l,i1

      allocate(pgraph(nm, ngraph, ngraph))

      do k =1,ngraph
         do l =1,ngraph-k+1
           y = dfloat(2*k-ngraph-1)/dfloat(ngraph-1)
           x = dfloat(2*l-ngraph-1)/dfloat(ngraph-1)
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

!
!----------------------------------------------
! calcul de la super connectivite
!----------------------------------------------
!
      subroutine connec
      implicit none
      integer ngrm,ngrmele,i,j,k,l
      integer, allocatable :: in3(:,:)
      character*60 :: format1,format2

      ngrmele = (ngraph-1)*(ngraph-1)
      ngrm = ngraph*(ngraph+1)/2
      allocate(in3(3,ngrmele*ne))

      nne=0
      do i=1,ne
         l = 0
         do k=1,ngraph-1
            do j=1,ngraph-k-1
               nne = nne + 1
               in3(1,nne) = ngrm*(i-1) + l +j
               in3(2,nne) = ngrm*(i-1) + l +j+1
               in3(3,nne) = ngrm*(i-1) + l + ngraph -k + 1 +j
               nne = nne + 1
               in3(1,nne) = ngrm*(i-1) + l +j+1
               in3(2,nne) = ngrm*(i-1) + l + ngraph -k + 1 +j+1
               in3(3,nne) = ngrm*(i-1) + l + ngraph -k + 1 +j
            enddo
            j = ngraph - k
               nne = nne + 1
               in3(1,nne) = ngrm*(i-1) + l +j
               in3(2,nne) = ngrm*(i-1) + l +j+1
               in3(3,nne) = ngrm*(i-1) + l + ngraph -k + 1 +j
               l = l + ngraph -k + 1
         enddo
      enddo
 
format1='1x,60(''-'')'
format2='(4(i6,1x))'

      open(21,file='connec.dat')
      do i=1,nne
        write(21,format2) in3(1,i)-1,in3(2,i)-1,in3(3,i)-1
      enddo
      close(21)

!------------------------------------
! write messages to standart output
!

      write(*,*)
      write(*,format1)
      write(*,*) 'number of points along one element edge:',ngraph
      write(*,*) 'number of sub-elements for graphic output:',nne
      write(*,*) 'the graphic output files are named: *.gra'

      return
      end

!
!*****************************************************************
! visu pour les elements courbes
!*****************************************************************
!
      subroutine write_gra(h,name,nper)
      implicit none
      integer nper
      integer kc1,kc2
      double precision h(nm,*),b(nm)
      double precision x,y,sol,pp
      integer i,j,k,l,i1,i2,i3,ind,tk
      double precision ax,bx,ay,by,x1,y1,xof,yof,af,bf,x2,y2
      integer ideb,ifin
      character name*1,anum*3,fbin*8
      character*60 :: format1

      format1='3(e10.3,1x)'
      
      write(anum,'i3') nper
      ideb=1
      do while (anum(ideb:ideb).eq.' ') 
         ideb=ideb+1
      enddo
         
      fbin = name//anum(ideb:3)//'.gra'
      ifin=9-ideb

!	write(*,*) fbin
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

         do k =1,ngraph
         do l =1,ngraph-k+1
           y1 = dfloat(2*k-ngraph-1)/dfloat(ngraph-1)
           x1 = dfloat(2*l-ngraph-1)/dfloat(ngraph-1)
           x = xof + ax * (x1+1.d0) + bx * (y1+1.d0)
           y = yof + ay * (x1+1.d0) + by * (y1+1.d0)
           sol = 0.d0
           do i1=1,nm
             pp = pgraph(i1,l,k)
             sol = sol + b(i1) * pgraph(i1,l,k)
           enddo
           write(1,format1) x,y,sol
         enddo
         enddo

      enddo
      close(1)


      return
      end

end module graphic
