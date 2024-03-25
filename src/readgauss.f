      subroutine read_gauss(ng,xi,yi,wi,ngf,wif,tif,nc,nm)
      implicit none
      include '../include/sp.inc'
      integer ng,ngf,nc,nm
      double precision xi(ng_max),yi(ng_max),wi(ng_max)
      double precision wif(ngf_max),tif(ngf_max)
      character*30 fng,fngf
      integer ideb,ifin
      character name*1,anum*3
      integer i
c
c-----------------------------------
c define name files to read
c
 22   format(i3)
      write(anum,22) ng
      ideb=1
      do while (anum(ideb:ideb).eq.' ') 
         ideb=ideb+1
      enddo
      fng  = 'include/tri2_'//anum(ideb:3)//'.dat'
c
      write(anum,22) ngf
      ideb=1
      do while (anum(ideb:ideb).eq.' ') 
         ideb=ideb+1
      enddo
      fngf = 'include/legen'//anum(ideb:3)//'.dat'
c
c-----------------------------------
c read files
c
      open(1,file=fng,status='old')
      do i=1,ng
         read(1,*) wi(i),xi(i),yi(i)
      enddo
      close(1)
      open(1,file=fngf,status='old')
      do i=1,ngf
         read(1,*) wif(i),tif(i)
      enddo
      close(1)
c
c------------------------------------
c write messages to standart output
c
c
      write(*,*)
      write(*,30)
      write(*,*) 'READ GAUSSIAN FILES'
      write(*,30)
      write(*,*)
      write(*,*) 'load 2d Gauss points on triangle: ',fng
      write(*,*) 'load 1d Gauss points: ',fngf
c
      write(*,10) 'tri gauss points    :',ng,ng_max
      write(*,10) '1d gauss points     :',ngf,ngf_max
      write(*,10) 'max polynomial order:',nc,nspdim
      nm = (nc+1)*(nc+2)/2
      write(*,10) 'matrice rank        :',nm,nnmod
c
30    format(1x,60('-'))
10    format(1x,a21,i7,3x,'(max = ',i7,')')
      return
      end
