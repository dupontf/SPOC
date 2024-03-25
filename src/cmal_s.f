c
c
c --------------------------------------- 
c    read the mesh from file CM.DAT 
c --------------------------------------- 
c
      subroutine read_mesh(nn,ne,in,xgr,ygr,scx,scy)
c
      implicit none
      include '../include/maille.inc'
      integer nn,ne
      double precision xgr(*),ygr(*)
      double precision scx,scy
      integer in(3,*)
      integer i,j,k
c
c -------------------------------------
c read the file
c
      open(1,file='cm.dat',status='old')
      read(1,*) nn,ne
	do i=1,nn
	   read(1,*) j,xgr(i),ygr(i)
	enddo
	do i=1,ne
	   read(1,*) k,(in(j,i),j=1,3)
	enddo
      close(1)
c
      do i=1,nn
         xgr(i) = xgr(i) * scx
         ygr(i) = ygr(i) * scy
      enddo
c
c------------------------------------
c write messages to standart output
c
c
      write(*,*)
      write(*,30)
      write(*,*) 'READ MESH'
      write(*,30)
      write(*,*)
      write(*,*) 'cm.dat have been successfully read'
c
30    format(1x,60('-'))
      return
      end
c
c
c
c --------------------------------------- 
c    save the mesh 
c --------------------------------------- 
c
      subroutine save_mesh(name,nn,ne,in,xgr,ygr,scx,scy)
c
      implicit none
      integer nn,ne
      double precision xgr(*),ygr(*)
      double precision scx,scy
      integer in(3,*)
      integer i
      character*10 name
c
	open(11,file=name)
	write(*,*) name
	write(11,*) nn,ne
	do i=1,nn
	    write(11,*) i,xgr(i)/scx,ygr(i)/scy
	enddo
	do i=1,ne
	   write(11,*) i,in(1,i),in(2,i),in(3,i)
	enddo
	close(11)
c
      return
      end
c
c
c --------------------------------------- c
c					  c
c    calcul des pointeurs pour maillage   c
c      quelconque                         c
c					  c
c --------------------------------------- c
c
      subroutine pointer_mesh(nn,ne,nf,in,itm,nmit,tmit,imit,
     1                  pseg,tseg)
c
      implicit none
      include '../include/maille.inc'
      integer pouet
c
c parametres et dimensions des vecteurs
       integer nn,ne,nf
c
       integer in(3,*),itm(3,*)
       integer nmit(*),tmit(*),imit(*)
       integer pseg(2,*),tseg(2,*)
c
       integer ndeb(nndim),nfin(nndim)
c
c entiers de service
      integer i,j,ind,k,tk,ik,lk,pj,it,ne1,ne2,ne3,
     1   i1,i2,nc,nd,buf,j1,j2
	 logical trouve
c
      write(*,*) 'compute mesh-related pointers'
c
	do i=1,nfdim
	   tseg(1,i)=0
	   tseg(2,i)=0
	enddo
c
c
c triangles (elements) mitoyens d'un noeud
c on suppose que chaque noeud est entoure de nftrdim quadrilateres
c et on remplit les vecteurs tmit et imit
c
      ndeb(1)=1
      nfin(1)=1
        do i=1,nn-1
           ndeb(i+1)=ndeb(i)+nftrdim
           nfin(i+1)=ndeb(i+1)
        enddo
c
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
c
c elimination des vides laisses dans tmit et imit
c
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
c
c
c
      nf=0
      do i=1,ne
         do j=1,3
            ne1=in(j,i)
            if (j.eq.3) then
               ne2=in(1,i)
            else
               ne2=in(j+1,i)
            endif
c
c   checker les anciennes faces des elements mitoyens
c
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
c               write(*,*) tk,ne2,ne3
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
c
c la face est nouvelle!
c
            if (trouve) then
               nf=nf+1
               itm(j,i)=nf
               tseg(1,nf)=i
               pseg(1,nf)=ne1
               pseg(2,nf)=ne2
            endif
         enddo      
      enddo
c
      write(*,10) 'nodes    :',nn,nndim
      write(*,10) 'elements :',ne,nedim
      write(*,10) 'faces    :',nf,nfdim
c
10    format(1x,a10,i7,3x,'(max = ',i7,')')
      return
      end
c
c******************************************************************
c check the mesh for obvious problems
c******************************************************************
c
      subroutine check_mesh(ne,nf,in,itm,tseg)
      implicit none
      include '../include/maille.inc'
      integer ne,nf,nfbd
      integer sbd(NBEDIM),pebd(nedim)
      integer in(3,*),itm(3,*),tseg(2,*)
      integer i,j,ns,tk

      do i=1,ne
         pebd(i)=0
      enddo

      nfbd=0
      do i=1,nf
         if (tseg(2,i).eq.0) then
            nfbd=nfbd+1
            sbd(nfbd)=i
         endif
      enddo

      do i=1,nfbd
         ns=sbd(i)
         tk = tseg(1,ns)
         pebd(tk)=pebd(tk)+1
      enddo
c
c
c if pebd(i) > 1 then the i-th triangle has more than one face on the boundary
c
      do i=1,ne
         if (pebd(i).gt.1) then
            write(*,*) 'Triangle',i,'has more than one face on the',
     &' boundary. Cannot continue'
            stop
         endif
      enddo

c
c
c if itm(1,i) must be a face on the boundary if the i-th element is on the bdy
c
      do i=1,nfbd
         ns=sbd(i)
         tk = tseg(1,ns)
         j=1
	 do while (itm(j,tk).ne.ns.and.j.le.3)
	    j=j+1
	 enddo
	 if (j.ne.1) then 
	    write(*,*) 'Pb with triangle',tk,': its 1st face is not on ',
     &'the boundary'
	 endif
      enddo

      return
      end
 
