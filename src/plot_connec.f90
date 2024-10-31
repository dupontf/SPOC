!******************************************************************
program get_bdy

 use mesh
 use boundary_condition
 
      implicit none

      integer i,i1,i2
      character*60 :: format1,format2,format3

!------------------------------------
! read input parameter file
!------------------------------------

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
      close(1)

      call read_mesh
      call pointer_mesh
      call check_mesh

      call readbel
      
      
format3='2(e13.6,1x)'
      open(1,file='connec.xy')
      
      do i=1,nf
      
         if (fflag(i).eq.0) then
 	  write(1,'a1') '>'
          i1 = pseg(1,i)
          i2 = pseg(2,i)
	  write(1,format3) xgr(i1),ygr(i1)
	  write(1,format3) xgr(i2),ygr(i2)
	 endif
      enddo
      close(1)

      end

