      program rmal
      implicit none
c
      integer nn,nnt
      integer ne,nt
      double precision x,y,xmin,xmax,ymin,ymax
      integer nx,ny,i1,i2,i3
      integer i,j,k,ind,ne1
c
      write(*,*) 'nx,ny de la grille reguliere'
      read(*,*) nx,ny
      write(*,*) 'donner dimensions max, xmin,xmax,ymin,ymax'
      read(*,*) xmin,xmax,ymin,ymax
c
      open(1,file='cm.dat')
      write(1,*) nx*ny,2*(nx-1)*(ny-1)
      ind=0
      do j=1,ny
         do i=1,nx
            ind=ind+1
            x=xmin+dfloat(i-1)*(xmax-xmin)/dfloat(nx-1)
            y=ymin+dfloat(j-1)*(ymax-ymin)/dfloat(ny-1)
            write(1,20) i,x,y
         enddo
      enddo
20    format(i6,1x,2(e23.16,2x))
c
      ne=ind
c
      ind=0
      do j=1,ny-1
         do i=1,nx-1
            ind=ind+1
            ne1=i+nx*(j-1)
            i1=ne1
            i2=ne1+1
            i3=ne1+nx
            write(1,*) i,i1,i2,i3
            ind=ind+1
            i1=ne1+1
            i2=ne1+1+nx
            i3=ne1+nx
            write(1,*) i,i1,i2,i3
         enddo
      enddo
      nt=ind
      close(1)
      end
