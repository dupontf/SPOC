module gauss

 implicit none
 integer :: ng,ngf,nc,nm
 double precision, allocatable :: xi(:),yi(:),wi(:)
 double precision, allocatable :: wif(:),tif(:)
 double precision, allocatable :: pi2(:,:),pico(:,:,:)
 double precision, allocatable :: &
       amm(:,:), ammo(:,:), amml(:,:), agx(:,:), agy(:,:)
 
 contains

!
!
! --------------------------------------- 
!    read the gauss files
! --------------------------------------- 
!
 
subroutine read_gauss
      implicit none
      character*30 fng,fngf
      integer ideb,ifin
      character name*1,anum*3
      integer i
      character*60 :: format1,format2
     
      allocate(xi(ng),yi(ng),wi(ng))
      allocate(wif(ngf),tif(ngf))
!
!-----------------------------------
! define name files to read
!
      write(anum,'i3') ng
      ideb=1
      do while (anum(ideb:ideb).eq.' ') 
         ideb=ideb+1
      enddo
      fng  = 'include/tri2_'//anum(ideb:3)//'.dat'

      write(anum,'i3') ngf
      ideb=1
      do while (anum(ideb:ideb).eq.' ') 
         ideb=ideb+1
      enddo
      fngf = 'include/legen'//anum(ideb:3)//'.dat'

!-----------------------------------
! read files
!
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

!------------------------------------
! write messages to standart output
!
!
format1='1x,60(''-'')'
format2='(1x,a21,i7,3x,'

      write(*,*)
      write(*,format1)
      write(*,*) 'READ GAUSSIAN FILES'
      write(*,format1)
      write(*,*)
      write(*,*) 'load 2d Gauss points on triangle: ',fng
      write(*,*) 'load 1d Gauss points: ',fngf

      write(*,format2) 'tri gauss points    :',ng
      write(*,format2) '1d gauss points     :',ngf
      write(*,format2) 'max polynomial order:',nc
      nm = (nc+1)*(nc+2)/2
      write(*,format2) 'matrice rank        :',nm

      return

end subroutine read_gauss

!
!**************************************************
! compute polynomial value at Gauss point locations
!**************************************************
!
subroutine cal_p2

 implicit none
 double precision po,x,y
 integer i,j,kc1,kc2,i1
 external po

      allocate(pi2(ng,nm))

      i1 = 0
      do kc1=0,nc
	   do kc2=0,nc-kc1
	      i1 = i1 + 1
	      do i=1,ng
	         y = yi(i)
	         x = xi(i)
	         pi2(i,i1)=po(kc1,x)*po(kc2,y)
		enddo
	   enddo
      enddo
      return
end subroutine cal_p2

!
!**************************************************************
! computes the value of the polynomes along the boundaries of 
! the element
!**************************************************************
!
subroutine cal_pico

      implicit none
      double precision po,x,y
      external po
      integer i1,kc1,kc2,i,j
      
      allocate(pico(ngf,nm,3))

	i1=0
	do kc1=0,nc
	   do kc2=0,nc-kc1
	    i1=i1+1
	    do j=1,ngf
	       x = tif(j)
	       y = -1.d0
	       pico(j,i1,1) = po(kc1,x)*po(kc2,y)
	    enddo
	   enddo
	enddo
	i1=0
	do kc1=0,nc
	   do kc2=0,nc-kc1
	    i1=i1+1
	    do j=1,ngf
	       x = -tif(j)
	       y =  tif(j)
	       pico(j,i1,2) = po(kc1,x)*po(kc2,y)
	    enddo
	   enddo
	enddo
	i1=0
	do kc1=0,nc
	   do kc2=0,nc-kc1
	    i1=i1+1
	    do j=1,ngf
	       x = -1.d0
	       y = -tif(j)
	       pico(j,i1,3) = po(kc1,x)*po(kc2,y)
	    enddo
	   enddo
	enddo

	return
end subroutine cal_pico

end module gauss
