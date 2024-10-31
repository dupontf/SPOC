module mesh

 implicit none

 integer :: nn,ne,nf,nfbd
 integer :: nftrdim=10
 INTEGER, allocatable :: in(:,:),tseg(:,:),pseg(:,:), &
       itm(:,:),nmit(:),tmit(:),imit(:),sbd(:)
 double precision, allocatable :: xgr(:),ygr(:)
 double precision :: scx,scy

 
 contains
!
!
! --------------------------------------- 
!    read the mesh from file CM.DAT 
! --------------------------------------- 
!
      subroutine read_mesh
      implicit none
      integer i,j,k
      character*20 :: format1
!
! -------------------------------------
! read the file
!
      open(1,file='cm.dat',status='old',action='read')
      read(1,*) nn,ne
      allocate(xgr(nn))
      allocate(ygr(nn))
      allocate(in(3,ne))

	do i=1,nn
	   read(1,*) j,xgr(i),ygr(i)
	enddo
	do i=1,ne
	   read(1,*) k,(in(j,i),j=1,3)
	enddo
      close(1)
!
      xgr = xgr * scx
      ygr = ygr * scy
!
!------------------------------------
! write messages to standart output
!
!
      format1='1x,60(''-'')'
      write(*,*)
      write(*,fmt=format1)
      write(*,*) 'READ MESH'
      write(*,fmt=format1)
      write(*,*)
      write(*,*) 'cm.dat have been successfully read'

      return
      end subroutine read_mesh
!
!
!
! --------------------------------------- 
!    save the mesh 
! --------------------------------------- 
!
      subroutine save_mesh(name)
!
      implicit none
      integer i
      character*10 name
!
      write(*,*) 'filename: ',name
	open(11,file=name,status='unknown',action='write')
	write(11,*) nn,ne
	do i=1,nn
	    write(11,*) i,xgr(i)/scx,ygr(i)/scy
	enddo
	do i=1,ne
	   write(11,*) i,in(1,i),in(2,i),in(3,i)
	enddo
	close(11)
!
      return
      end subroutine save_mesh
!
!
! --------------------------------------- 
!					  
!    calcul des pointeurs pour maillage   
!      quelconque                         
!					  
! --------------------------------------- 
!
      subroutine pointer_mesh
!
      implicit none
      integer pouet
      integer, allocatable :: ndeb(:),nfin(:)
!
! entiers de service
      integer i,j,ind,k,tk,ik,lk,pj,it,ne1,ne2,ne3,i1,i2,nc,nd,buf,j1,j2
      logical trouve
      character :: format1*20
!
      allocate(itm(3,ne))
      allocate(tseg(2,3*nn))
      allocate(pseg(2,3*nn))
      allocate(nmit(nn+1))
      allocate(imit(nn*nftrdim))
      allocate(tmit(nn*nftrdim))
      allocate(ndeb(nn+1))
      allocate(nfin(nn+1))

      write(*,*) 'compute mesh-related pointers'
!
      tseg=0.0
!
!
! triangles (elements) mitoyens d'un noeud
! on suppose que chaque noeud est entoure de nftrdim quadrilateres
! et on remplit les vecteurs tmit et imit
!
      ndeb(1)=1
      nfin(1)=1
        do i=1,nn-1
           ndeb(i+1)=ndeb(i)+nftrdim
           nfin(i+1)=ndeb(i+1)
        enddo
!
        do i=1,ne
           do k=1,3
              ne1=in(k,i)
              ind=nfin(ne1)
              tmit(ind)=i
              imit(ind)=k
              nfin(ne1)=ind+1
              do while (nfin(ne1).eq.ndeb(ne1+1))
                 ne1=ne1+1
                 do j1=nfin(ne1)-1,ndeb(ne1),-1
                    tmit(j1+1)=tmit(j1)
                    imit(j1+1)=imit(j1)
                 enddo
                 ndeb(ne1)=ndeb(ne1)+1
                 nfin(ne1)=nfin(ne1)+1
              enddo
           enddo
         enddo
!
! elimination des vides laisses dans tmit et imit
!
         ind=0
         do i=1,nn
            do j=ndeb(i),nfin(i)-1
               ind=ind+1
               tmit(ind)=tmit(j)
               imit(ind)=imit(j)
            enddo
        enddo
        nmit(1)=1
        do i=1,nn
           nmit(i+1)=nmit(i)+nfin(i)-ndeb(i)
        enddo
!
!
!
      nf=0
      do i=1,ne
         do j=1,3
            ne1=in(j,i)
            if (j.eq.3) then
               ne2=in(1,i)
            else
               ne2=in(j+1,i)
            endif
!
!   checker les anciennes faces des elements mitoyens
!
            k=nmit(ne2)
            tk=tmit(k)
            trouve=.true.
            ik=imit(k)            
            do while(tk.lt.i.and.trouve)
               if (ik+1.lt.4) then
                  ne3=in(ik+1,tk)
               else
                  ne3=in(1,tk)
               endif
!               write(*,*) tk,ne2,ne3
               if (ne1.eq.ne3) then
                  lk=itm(ik,tk)
                  tseg(2,lk)=i
                  itm(j,i)=lk
                  trouve=.false.
               endif
               k=k+1
               tk=tmit(k)
               ik=imit(k)            
            enddo
!
! la face est nouvelle!
!
            if (trouve) then
               nf=nf+1
               itm(j,i)=nf
               tseg(1,nf)=i
               pseg(1,nf)=ne1
               pseg(2,nf)=ne2
            endif
         enddo      
      enddo

      deallocate (ndeb)
      deallocate (nfin)


! --------- creation of the boundary pointer

      nfbd=0
      do i=1,nf
         if (tseg(2,i).eq.0) then
            nfbd=nfbd+1
         endif
      enddo
      allocate(sbd(nfbd))
      k=0
      do i=1,nf
         if (tseg(2,i).eq.0) then
            k=k+1
            sbd(k)=i
         endif
      enddo
      
 format1='1x,a10,i7'

      write(*,fmt=format1) 'nodes    :',nn
      write(*,fmt=format1) 'elements :',ne
      write(*,fmt=format1) 'faces    :',nf

      return
      end subroutine pointer_mesh

!
!******************************************************************
! check the mesh for obvious problems
!******************************************************************
!
      subroutine check_mesh
      implicit none
      integer pebd(:)
      integer i,j,ns,tk

      allocate(pebd(ne))
      pebd=0

      do i=1,nfbd
         ns=sbd(i)
         tk = tseg(1,ns)
         pebd(tk)=pebd(tk)+1
      enddo
!
!
! if pebd(i) > 1 then the i-th triangle has more than one face on the boundary
!
      do i=1,ne
         if (pebd(i).gt.1) then
            write(*,*) 'Triangle',i,'has more than one face on the',&
      ' boundary. Cannot continue'
            stop
         endif
      enddo

!
!
! if itm(1,i) must be a face on the boundary if the i-th element is on the bdy
!
      do i=1,nfbd
         ns=sbd(i)
         tk = tseg(1,ns)
         j=1
	 do while (itm(j,tk).ne.ns.and.j.le.3)
	    j=j+1
	 enddo
	 if (j.ne.1) then 
	    write(*,*) 'Pb with triangle',tk,': its 1st face is not on ',&
      'the boundary'
	 endif
      enddo

      deallocate(pebd)

      return
      end subroutine check_mesh

end module mesh
