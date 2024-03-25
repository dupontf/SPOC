c******************************************************************
      program get_range
      implicit none
      integer n
      real x,y,a
      real xmin,xmax,ymin,ymax,amin,amax,fact
      character*20 filename
c
c------------------------------------
c read input parameter file
c------------------------------------
c
      fact = 1.1
      read(*,*) filename
      open(1,file=filename,status='old')
      n=0
      read(1,*,end=110) x,y,a
      xmin = x
      ymin = y
      amin = a
      xmax = x
      ymax = y
      amax = a
      n=n+1
100   read(1,*,end=110) x,y,a
      n=n+1
      xmin = min(xmin,x)
      ymin = min(ymin,y)
      amin = min(amin,a)
      xmax = max(xmax,x)
      ymax = max(ymax,y)
      amax = max(amax,a)
      goto 100
110   write(*,*) 'read',n,' points'
      write(*,*) 'suggested region:'
      write(*,20) xmin*fact,xmax*fact,ymin*fact,ymax*fact
      write(*,*) 'suggested range for color:'
      write(*,21) amin,amax,(amax-amin)/10.
20    format(sp,1p,e9.2,3('/',e9.2))
21    format(sp,1p,e9.2,2('/',e9.2))
      end
c
