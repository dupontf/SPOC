!
!***********************************************************************
! condition periodique
!
! la periodicite rend caduque l'utilisation de pseg,nmit,imit,tmit
! reste correct: in,tseg,itm
!***********************************************************************
!

module boundary_condition

 use mesh
 use gauss
 
 implicit none
 
 integer, allocatable :: fflag(:),fflagt(:),pfp(:,:)
 integer, allocatable :: pobc1(:)
 integer nfp, nobc1

 double precision, allocatable :: cosng(:,:),sinng(:,:)
 double precision omega,ramp_period
  
 contains
 
!
!***********************************************************************
! read the *.BEL file
!***********************************************************************
!
subroutine readbel

 implicit none

      integer i,j,k,n1,n2,jfp,ind,i1,i2,i3,i4,ele
      integer flag
      character*80 line
      integer, allocatable :: pnodebdy(:,:)
      
      allocate(pfp(2,nfbd))
      allocate(fflag(nf))
      allocate(fflagt(ne))
      allocate(pnodebdy(2,nn))
      allocate(pobc1(nfbd))
!
!----------------------------------------
! construct the sbd and pbnodebdy pointers
!
      flag=0
      pnodebdy=0


      do i=1,nf
	 fflag(i)=0
         if (tseg(2,i).eq.0) then
	    fflag(i)=1
	    n1=pseg(1,i)
	    n2=pseg(2,i)
            pnodebdy(1,n2)=i
            pnodebdy(2,n1)=i
            ele=tseg(1,i)
            fflagt(ele)=flag
	 endif
      enddo      
!
!----------------------------------------
! read the bel file
!----------------------------------------
!
      nobc1=0
      nfp=0
      write(*,*)
      open(2,file='cm.bel',status='old')
101   read(2,'a80',end=100) line
      read(line,*) j,i1,i2,flag
      n1=pnodebdy(2,i1)
      n2=pnodebdy(1,i2)
      if (n1.ne.n2) then
	  write(*,*) 'in periodic:'
	  write(*,*) 'pb with face defined by nodes:',i1,i2
	  stop
      endif
      fflag(n1)=flag
      ele=tseg(1,n1)
      fflagt(ele)=flag
!
! -------------- definition of open boundary -------------
! 
      if (flag.eq.5) then
         nobc1=nobc1+1
	 pobc1(nobc1)=n1
!
! ------------ definition of periodic boundary -------------
! 
      else if (flag.eq.6) then 
          read(line,*) j,i1,i2,flag,i3,i4
	  i1=pnodebdy(2,i3)
	  n2=pnodebdy(1,i4)
	  if (i1.ne.n2) then
	     write(*,*) 'in periodic:'
	     write(*,*) 'pb with face defined by nodes:',i3,i4
	     stop
	  endif
	  if (fflag(n2).eq.1) then
             nfp=nfp+1
	     pfp(1,nfp)=n1
             pfp(2,nfp)=n2
	  endif
      endif
      goto 101
100   close(2)
!
! ------ print some information ---------
!
      if (nfp.gt.0) then
       write(*,*) 'PERIODIC FACES: '
       write(*,*) 'number =',nfp
      elseif (nobc1.gt.0) then
       write(*,*) 'OPEN BOUNDARY FACES (type 1): '
       write(*,*) 'number =',nobc1
      endif
      
      deallocate(pnodebdy)
      
      return
end subroutine readbel

!
!***********************************************************************
! condition periodique
!
! la periodicite rend caduque l'utilisation de pseg,nmit,imit,tmit
! reste correct: in,tseg,itm
!***********************************************************************
!
subroutine periodic

 implicit none

      integer i,j,k,n1,n2,jfp,ind,i1,i2,i3,i4
      integer ppf(nf)

!
!---trouver les points sur les bords ouest et est
!
      if (nfp.gt.0) then

      do i=1,nfp
	 n1 = pfp(1,i)
	 n2 = pfp(2,i)
         if (n1.gt.n2) then
             pfp(1,i) = n2
             pfp(2,i) = n1
         endif
         n1 = pfp(1,i)
         n2 = pfp(2,i)
	     write(*,*) nfp,n1,n2
         tseg(2,n1) = tseg(1,n2)
         tseg(2,n2) = tseg(1,n1)
      enddo
      pfp(1,nfp+1)=0
      pfp(2,nfp+1)=0
!
!---creation d'un pointeur vers la valeur final de chaque face
!---et nettoyage de tseg,pseg
!
      ind=0
      do i=1,nf

         jfp = 1
         do while (pfp(2,jfp).ne.i.and.jfp.le.nfp) 
            jfp = jfp + 1
         enddo
         if (pfp(2,jfp).eq.i) then
            ppf(i)=ppf(pfp(1,jfp))
	    write(*,*) i,ind
         else
            ind=ind+1
            ppf(i)=ind
	    write(*,*) i,ind
            tseg(1,ind) = tseg(1,i)
            tseg(2,ind) = tseg(2,i)
            pseg(1,ind) = pseg(1,i)
            pseg(2,ind) = pseg(2,i)
	    fflag(ind)=fflag(i)
         endif
      enddo

!
!---elimination des faces en trop dans itm
!
      do i=1,ne
         do j=1,3
            k = itm(j,i)
            itm(j,i)=ppf(k)
         enddo
      enddo
!
!------------------------------
      nf = nf - nfp
!------------------------------
! end if condition on nfp

      endif

      return
end subroutine periodic

!
!***********************************************************************
! open boundary condition
!
! the elevation is prescribed at the boundary
!***********************************************************************
!
subroutine openbc1_read

 implicit none

    integer i,j,k

    if (nobc1.gt.0) then

      allocate(cosng(ngf,nobc1),sinng(ngf,nobc1))

      open(1,file='obc1.dat',status='old')
      read(1,*) omega,ramp_period
      do i=1,nobc1
        do j=1,ngf
	 read(1,*) cosng(j,i),sinng(j,i)
	enddo
      enddo
      close(1)

    endif

    return

end subroutine openbc1_read

!
!***********************************************************************
! open boundary condition
!
! the elevation is prescribed at the boundary
!***********************************************************************
!
subroutine openbc1(fbd,time,fact)

 implicit none

      integer i,j,k,face
      double precision fbd(ngf,*),time,phase,cosa,sina,fact,ans
      double precision pi,per,ramp

      if (nobc1.gt.0) then
	 if (ramp_period.gt.1.d-12) then
	  ramp = min (time,ramp_period)
	  ramp = ramp / ramp_period
	 else
	  ramp=1.d0
	 endif
	 phase = time * omega
	 cosa = cos(phase) * fact * ramp
	 sina = sin(phase) * fact * ramp
         do i=1,nobc1
            face = pobc1(i)
 	    do j=1,ngf
	       fbd(j,face) = cosng(j,i) * cosa + sinng(j,i) * sina
	    enddo
         enddo
      endif

      return
end subroutine openbc1

!
!***********************************************************************
! open boundary condition
! condition for mass conservation
!
!***********************************************************************
!
subroutine massbc(fbdx,fbdy)

 implicit none

      integer i,j,k,face,flag
      double precision fbdx(ngf,*),fbdy(ngf,*)

      do face=1,nf
       flag = fflag(face)
       if (flag.eq.0) then
        do j=1,ngf
         fbdx(j,face) = fbdx(j,face) + fbdy(j,face)
        enddo
       else if (flag.eq.1) then
        do j=1,ngf
         fbdx(j,face) = 0.d0
        enddo
       else if (flag.eq.5) then
        do j=1,ngf
         fbdx(j,face) = fbdx(j,face) + fbdy(j,face)
        enddo
       endif
      enddo

      return
end subroutine massbc

end module boundary_condition
