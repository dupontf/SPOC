      program rmal
      implicit none
c
      integer nn,ne
      real x,y,xmin,xmax,ymin,ymax,dx,xoff
      integer nx,ny,i1,i2,i3
      integer i,j,k,ind,ne1,j0
c
      write(*,*) 'nx,ny de la grille reguliere'
      read(*,*) nx,ny
      if (mod(nx,2).eq.0.or.mod(ny,2).eq.0) then
        write(*,*) '!!!! nx and ny must be odd numbers !!!!'
	stop
      endif
      write(*,*) 'donner dimensions max, xmin,xmax,ymin,ymax'
      read(*,*) xmin,xmax,ymin,ymax
      dx = (xmax-xmin)/float(nx-1)
      xoff = dx*0.5
c
      nn=(ny+1)*nx/2+(ny-1)*(nx+1)/2
      ne=(ny-1)*(2*nx-1)
      open(1,file='cm.dat')
      write(1,*) nn,ne
      ind=0
      do j=1,ny-1,2
         do i=1,nx
            ind=ind+1
            x=xmin+float(i-1)*(xmax-xmin)/float(nx-1)
            y=ymin+float(j-1)*(ymax-ymin)/float(ny-1)
            write(1,19) ind,x,y
         enddo
         j0=j+1
         i=1
            ind=ind+1
            x=xmin+float(i-1)*(xmax-xmin)/float(nx-1)
            y=ymin+float(j0-1)*(ymax-ymin)/float(ny-1)
            write(1,19) ind,x,y
         do i=2,nx
            ind=ind+1
            x=xmin+float(i-1)*(xmax-xmin)/float(nx-1)-xoff
            y=ymin+float(j0-1)*(ymax-ymin)/float(ny-1)
            write(1,19) ind,x,y
         enddo
         i=nx
            ind=ind+1
            x=xmin+float(i-1)*(xmax-xmin)/float(nx-1)
            y=ymin+float(j0-1)*(ymax-ymin)/float(ny-1)
            write(1,19) ind,x,y
      enddo
      j=ny
         do i=1,nx
            ind=ind+1
            x=xmin+float(i-1)*(xmax-xmin)/float(nx-1)
            y=ymin+float(j-1)*(ymax-ymin)/float(ny-1)
            write(1,19) ind,x,y
         enddo
      ind=0
      do j=1,ny-1,2
         do i=1,nx-1
            ind=ind+1
            i1=i+(2*nx+1)*((j-1)/2)
            i2=i+(2*nx+1)*((j-1)/2)+nx+1
            i3=i+(2*nx+1)*((j-1)/2)+nx
            write(1,20) ind,i1,i2,i3
            ind=ind+1
            i1=i+(2*nx+1)*((j-1)/2)
            i2=i+(2*nx+1)*((j-1)/2)+1
            i3=i+(2*nx+1)*((j-1)/2)+nx+1
            write(1,20) ind,i1,i2,i3
         enddo
         i=nx
            ind=ind+1
            i1=i+(2*nx+1)*((j-1)/2)
            i2=i+(2*nx+1)*((j-1)/2)+nx+1
            i3=i+(2*nx+1)*((j-1)/2)+nx
            write(1,20) ind,i1,i2,i3
         do i=1,nx-1
            ind=ind+1
            i1=i+(2*nx+1)*((j-1)/2)+nx
            i2=i+(2*nx+1)*((j-1)/2)+nx+1
            i3=i+(2*nx+1)*((j-1)/2)+2*nx+1
            write(1,20) ind,i1,i2,i3
            ind=ind+1
            i1=i+(2*nx+1)*((j-1)/2)+nx+1
            i2=i+(2*nx+1)*((j-1)/2)+2*nx+2
            i3=i+(2*nx+1)*((j-1)/2)+2*nx+1
            write(1,20) ind,i1,i2,i3
         enddo
         i=nx
            ind=ind+1
            i1=i+(2*nx+1)*((j-1)/2)+nx
            i2=i+(2*nx+1)*((j-1)/2)+nx+1
            i3=i+(2*nx+1)*((j-1)/2)+2*nx+1
            write(1,20) ind,i1,i2,i3
      enddo
      close(1)
19    format(i6,1x,2(e23.16,2x))
20    format(1x,i5,3x,3(i4,1x))
      end
