!******************************************************************
! This programs changes the connectivity table so that 
! the boundary elements have their first and second vertices 
! on the boundary.
! For the elements having only one vertex on the boundary, 
! this vertex must be their first.
!******************************************************************
program chmaille

 use mesh
 implicit none

 integer :: i,j,k
 character :: filename*10

      scx=1.d0
      scy=1.d0
      call read_mesh
      call pointer_mesh
      call changemaille
      filename='cm.dat'
      call save_mesh(filename)

end program chmaille

!
!**************************************************************
! elimination de triangles qui ont deux faces sur la frontiere
! rotation des reperes locales dans le triangle pour que
! l'axe des x soit tangeant a la frontiere
!
subroutine changemaille
 use mesh
      implicit none
      integer, allocatable :: pebd(:)
      integer i,j,k,l,ns1,j1,j2,i1,i2,i3,i4,i5,i6,k1,k2,k3,k4,tk1,tk2
      integer pp(3),pf(3),cit
      
      allocate(pebd(ne))
      pebd=0

      do i=1,nfbd
         ns1=sbd(i)
         tk1 = tseg(1,ns1)
         pebd(tk1)=pebd(tk1)+1
      enddo
!
!
! si pebd(i)>=2 alors triangle -i- a deux faces sur la frontiere
!
      do i=1,ne
         if (pebd(i).eq.2) then
            write(*,*) 'tk',i
!
! recherche de l'arete partagee
!
            j1=1
            ns1=itm(j1,i)
            do while(tseg(2,ns1).eq.0)
              j1=j1+1
              ns1=itm(j1,i)
            enddo
            tk1 = tseg(1,ns1)
            tk2 = tseg(2,ns1)
            j1=1
            do while(itm(j1,tk1).ne.ns1)
              j1=j1+1
            enddo
            j2=1
            do while(itm(j2,tk2).ne.ns1)
              j2=j2+1
            enddo

!------- find the 4 nodes to be exchanged

!           i1 = in(cit(3,j1,1),tk1)
!           i2 = in(cit(3,j1,2),tk1)
!           i3 = in(j1,tk1)
!           i4 = in(j2,tk2)
!	    
!	    in(1,tk1) = i2
!	    in(2,tk1) = i3
!	    in(3,tk1) = i4
!
!	    in(1,tk2) = i1
!	    in(2,tk2) = i2
!	    in(3,tk2) = i4

!------- 
            k1 = cit(3,j1,1)
            k2 = cit(3,j2,-1)
            k3 = cit(3,j2,1)
            k4 = cit(3,j1,-1)

            i1 = in(k1,tk1)
            i2 = in(k2,tk2)
            i3 = in(k3,tk2)
            i4 = in(k4,tk1)

            in(k1,tk1)=i2
            in(k3,tk2)=i4
            itm(j1,tk1)=itm(k3,tk2)
            itm(j2,tk2)=itm(k1,tk1)
            itm(k1,tk1)=ns1
            itm(k3,tk2)=ns1
            pseg(1,ns1)=i2
            pseg(2,ns1)=i4
            if (i.eq.tk1) then
            tseg(1,itm(j2,tk2))=tk2
            else
            tseg(1,itm(j1,tk1))=tk1
            endif
!         else if (pebd.eq.1) then
!	 
         endif
      enddo
!
!
! rotation des aretes jusqu'a ce que l'arete 1 soit sur la frontiere
!
      do i=1,nfbd
         ns1=sbd(i)
         tk1 = tseg(1,ns1)
            j=1
            do while(itm(j,tk1).ne.ns1)
               j=j+1
            enddo
!              write(*,*) i,j
            do k=1,3
              l = cit(3,k,j-1)
              pp(k) = in (l,tk1)
              pf(k) = itm(l,tk1)
            enddo
            do k=1,3
              in (k,tk1)= pp(k)
              itm(k,tk1)= pf(k)
            enddo
      enddo
!
!
      return
end subroutine changemaille
        
function cit(n,k,np)
      integer cit,n,k,np,k1

      k1=k+np
      if(k1.gt.n) k1=k1-n
      if(k1.lt.1) k1=k1+n
      cit=k1
      
return
end function cit


