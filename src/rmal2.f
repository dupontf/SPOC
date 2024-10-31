      program rmal
      implicit none
c
      integer nn,ne
      double precision x,y,xmin,xmax,ymin,ymax,dx,xoff,r
      integer nx,ny,i1,i2,i3
      integer i,j,k,ind,ne1,j0
      double precision rand
      integer s0
      external rand
c
      write(*,*) 'nx,ny de la grille reguliere'
      read(*,*) nx,ny
      write(*,*) 'donner dimensions max, xmin,xmax,ymin,ymax'
      read(*,*) xmin,xmax,ymin,ymax
      dx = (xmax-xmin)/float(nx-1)
      xoff = dx*0.5
      write(*,*) 'enter r random amplitude'
      read(*,*) r
      r=r*dx
c
        s0 = 14441023
        call srand(s0)

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
	    call add_rand(xmin,xmax,ymin,ymax,x,y,r)
            write(1,*) ind,x,y
         enddo
         j0=j+1
         i=1
            ind=ind+1
            x=xmin+float(i-1)*(xmax-xmin)/float(nx-1)
            y=ymin+float(j0-1)*(ymax-ymin)/float(ny-1)
	    call add_rand(xmin,xmax,ymin,ymax,x,y,r)
            write(1,*) ind,x,y
         do i=2,nx
            ind=ind+1
            x=xmin+float(i-1)*(xmax-xmin)/float(nx-1)-xoff
            y=ymin+float(j0-1)*(ymax-ymin)/float(ny-1)
	    call add_rand(xmin,xmax,ymin,ymax,x,y,r)
            write(1,*) ind,x,y
         enddo
         i=nx
            ind=ind+1
            x=xmin+float(i-1)*(xmax-xmin)/float(nx-1)
            y=ymin+float(j0-1)*(ymax-ymin)/float(ny-1)
	    call add_rand(xmin,xmax,ymin,ymax,x,y,r)
            write(1,*) ind,x,y
      enddo
      j=ny
         do i=1,nx
            ind=ind+1
            x=xmin+float(i-1)*(xmax-xmin)/float(nx-1)
            y=ymin+float(j-1)*(ymax-ymin)/float(ny-1)
	    call add_rand(xmin,xmax,ymin,ymax,x,y,r)
            write(1,*) ind,x,y
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
20    format(1x,i5,3x,3(i4,1x))
      end
      subroutine add_rand(xmin,xmax,ymin,ymax,x,y,r)
      double precision rand
      double precision xmin,xmax,ymin,ymax,x,y,x1,y1,eps,r
      external rand
      eps=1.e-3
      if ((abs(x-xmin).gt.eps.and.abs(x-xmax).gt.eps)
     & .and.(abs(y-ymin).gt.eps.and.abs(y-ymax).gt.eps) ) then
       x1 = r*(rand()-0.5d0)
       y1 = r*(rand()-0.5d0)
       x = x + x1
       y = y + y1
      endif
      return
      end
      
