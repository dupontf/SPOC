c*****************************************************************
c routines handling the ouput the variables for VU (CERCA) graph interface
c*****************************************************************
	subroutine cal_pgraph(nc,n,pgraph,ngradim)
	implicit none
	include '../include/sp.inc'
	integer kc1,kc2,nc,nm,n,ngradim
	double precision pgraph(nnmod,ngradim,ngradim),x,y,po
	external po
	integer k,l,i1
	do k =1,n
	   do l =1,n-k+1
	     y = dfloat(2*k-n-1)/dfloat(n-1)
	     x = dfloat(2*l-n-1)/dfloat(n-1)
	     i1=0
	     do kc1=0,nc
	      do kc2=0,nc-kc1
		 i1=i1+1
	       pgraph(i1,l,k) = po(kc1,x)*po(kc2,y)
		enddo
	     enddo
	   enddo
	enddo
	return
	end
c
c*****************************************************************
c
	subroutine visu(n,ne,xgr,ygr,in,nc,u,v,h,nper,pgraph,ngradim)
	implicit none
	include '../include/sp.inc'
	integer n,nc,nm,ne
	integer kc1,kc2,ngradim
	double precision h(nnmod,*),u(nnmod,*),v(nnmod,*),xgr(*),ygr(*)
	double precision a(nnmod),b(nnmod),c(nnmod),x,y,
     &sol1,sol2,sol3,po,pp,pgraph(nnmod,ngradim,ngradim)
	integer i,j,k,l,i1,i2,i3
	integer in(3,*)
	double precision ax,bx,ay,by,x1,y1,xof,yof
	integer ideb,ifin,nper
	external po
      character fbin*8,anum*3
c
c
c -----------------------------------------------------------
c   ecriture du fichier binaire des resultats
c -----------------------------------------------------------
c
 22   format(i3)
      write(anum,22) nper
      ideb=1
      do while (anum(ideb:ideb).eq.' ') 
         ideb=ideb+1
      enddo
c         
      fbin = 'f'//anum(ideb:3)//'.gra'
      ifin=9-ideb
c
	nm = (nc+1)*(nc+2)/2
c
	open(1,file=fbin)
	write(1,*) ne,n
	do i=1,ne
	   i1 = in(1,i)
	   i2 = in(2,i)
	   i3 = in(3,i)
	   ax = (xgr(i2)-xgr(i1))*0.5d0
	   bx = (xgr(i3)-xgr(i1))*0.5d0
	   ay = (ygr(i2)-ygr(i1))*0.5d0
	   by = (ygr(i3)-ygr(i1))*0.5d0
	   xof = xgr(i1)
	   yof = ygr(i1)
c
	   do j=1,nm
	      a(j)= u(j,i)
	      b(j)= v(j,i)
	      c(j)= h(j,i)
	   enddo
c
	   do k =1,n
	   do l =1,n-k+1
	     y1 = dfloat(2*k-n-1)/dfloat(n-1)
	     x1 = dfloat(2*l-n-1)/dfloat(n-1)
	     x = xof + ax * (x1+1.d0) + bx * (y1+1.d0)
	     y = yof + ay * (x1+1.d0) + by * (y1+1.d0)
	     sol1 = 0.d0
	     sol2 = 0.d0
	     sol3 = 0.d0
	     do i1=1,nm
	       pp = pgraph(i1,l,k)
	       sol1 = sol1 + a(i1) * pp
	       sol2 = sol2 + b(i1) * pp
	       sol3 = sol3 + c(i1) * pp
	     enddo
	     write(1,21) x,y,sol1,sol2,sol3
	   enddo
	   enddo
	enddo
	close(1)
c
20      format(20(f13.10,1x))
21      format(20(e13.6,1x))
	return
	end
c
c
c----------------------------------------------
c calcul de la super connectivite
c----------------------------------------------
c
	subroutine connec(nne,ne,ngr)
	implicit none
	include '../include/maille.inc'
	include '../include/graph.inc'
	integer ne,ngr,nm,nne,i,j,k,l
	integer in(3,2*nmg*nedim)
	nm = ngr*(ngr+1)/2
	nne=0
	do i=1,ne
	   l = 0
	   do k=1,ngr-1
	      do j=1,ngr-k-1
		   nne = nne + 1
		   in(1,nne) = nm*(i-1) + l +j
		   in(2,nne) = nm*(i-1) + l +j+1
		   in(3,nne) = nm*(i-1) + l + ngr -k + 1 +j
		   nne = nne + 1
		   in(1,nne) = nm*(i-1) + l +j+1
		   in(2,nne) = nm*(i-1) + l + ngr -k + 1 +j+1
		   in(3,nne) = nm*(i-1) + l + ngr -k + 1 +j
		enddo
		j = ngr - k
		   nne = nne + 1
		   in(1,nne) = nm*(i-1) + l +j
		   in(2,nne) = nm*(i-1) + l +j+1
		   in(3,nne) = nm*(i-1) + l + ngr -k + 1 +j
		   l = l + ngr -k + 1
	   enddo
	enddo
	write(*,*) ne,nne
c 
      open(21,file='connec.bin',form='unformatted')
	write(21) (in(1,i),in(2,i),in(3,i),i=1,nne)
	close(21)
c
	return
	end
c
	subroutine spectovu(nne,ne,nper,ngr)
	implicit none
	include '../include/maille.inc'
	integer ne,nper,ngr,n
	character fpie*8,fbin*8,ligne*60,val*20,anum*3
	integer i,nm,ideb,ifin,nne,k,l,j,ind
c
	nm=ngr*(ngr+1)/2
	n=ne*nm
c
c -----------------------------------------------------------
c   ecriture du fichier binaire des resultats
c -----------------------------------------------------------
c
 22   format(i3)
      write(anum,22) nper
      ideb=1
      do while (anum(ideb:ideb).eq.' ') 
         ideb=ideb+1
      enddo
c         
      fpie = 'f'//anum(ideb:3)//'.pie'
      fbin = 'f'//anum(ideb:3)//'.bin'
      ifin=9-ideb
	write(*,*) fpie(1:ifin)
c
c
      open(4,file=fpie(1:ifin))
      write(4,*) 'MAILLAGE ma()='
      write(4,*) '{'
      write(4,*) '  ZONE z1(LagrTrian03,coord,ele);'
      write(4,*) '};'
      write(4,29) 2*n,4
      write(4,28) 3*nne,4
      write(4,*) 
 28     format(x,'CHAMP <int> ele("connec.bin",',i6,',',i6,');')
 29     format(x,'CHAMP coord("xy.bin",',i6,',',i6,');')
 37     format(' VARIABLE ',a2,'(LagrTrian03,',a2,',ele,z1);')
 39   format(',',i6,',',i7,');')
c
      write(val,39) n,4
      ligne='CHAMP u  ("'//fbin(1:ifin)//'"'//val
      write(4,*) ligne
      write(val,39) n,(8+8*n)+4
      ligne='CHAMP v  ("'//fbin(1:ifin)//'"'//val
      write(4,*) ligne
      write(val,39) n,2*(8+8*n)+4
      ligne='CHAMP h  ("'//fbin(1:ifin)//'"'//val
      write(4,*) ligne
      write(val,39) n,3*(8+8*n)+4
c
        write(4,*) 
        write(4,*) ' SOLUTION sol( )='
	write(4,*) ' {'
        write(4,37) 'u','u'
        write(4,37) 'v','v'
        write(4,37) 'h','h'
	write(4,*) '};'
	close(4)
c
	return
	end
c
c*****************************************************************
c
	subroutine spectovubin(ngr,ne,xgr,ygr,in,nc,u,v,h,nper,pgraph)
	implicit none
	include '../include/maille.inc'
	include '../include/sp.inc'
	include '../include/graph.inc'
	integer ngr,nc,nm,ne,n
	integer kc1,kc2
	double precision h(nnmod,*),u(nnmod,*),v(nnmod,*),xgr(*),ygr(*)
	double precision va(nmg*nedim),vb(nmg*nedim),vc(nmg*nedim),
     &                 vx(nmg*nedim),vy(nmg*nedim)
	double precision a(nnmod),b(nnmod),c(nnmod),x,y,
     &                 sol1,sol2,sol3,po,pp,
     &                 pgraph(nnmod,ngradim,ngradim)
	integer i,j,k,l,i1,i2,i3,ind
	integer in(3,*)
	double precision ax,bx,ay,by,x1,y1,xof,yof
	integer ideb,ifin,nper
	external po
      character fbin*8,anum*3
c
c
c -----------------------------------------------------------
c   ecriture du fichier binaire des resultats
c -----------------------------------------------------------
c
 22   format(i3)
      write(anum,22) nper
      ideb=1
      do while (anum(ideb:ideb).eq.' ') 
         ideb=ideb+1
      enddo
c         
      fbin = 'f'//anum(ideb:3)//'.bin'
      ifin=9-ideb
c
	nm = (nc+1)*(nc+2)/2
	ind=0
c
	do i=1,ne
	   i1 = in(1,i)
	   i2 = in(2,i)
	   i3 = in(3,i)
	   ax = (xgr(i2)-xgr(i1))*0.5d0
	   bx = (xgr(i3)-xgr(i1))*0.5d0
	   ay = (ygr(i2)-ygr(i1))*0.5d0
	   by = (ygr(i3)-ygr(i1))*0.5d0
	   xof = xgr(i1)
	   yof = ygr(i1)
c
	   do j=1,nm
	      a(j)= u(j,i)
	      b(j)= v(j,i)
	      c(j)= h(j,i)
	   enddo
c
	   do k =1,ngr
	   do l =1,ngr-k+1
	     y1 = dfloat(2*k-ngr-1)/dfloat(ngr-1)
	     x1 = dfloat(2*l-ngr-1)/dfloat(ngr-1)
	     x = xof + ax * (x1+1.d0) + bx * (y1+1.d0)
	     y = yof + ay * (x1+1.d0) + by * (y1+1.d0)
	     sol1 = 0.d0
	     sol2 = 0.d0
	     sol3 = 0.d0
	     do i1=1,nm
	       pp = pgraph(i1,l,k)
	       sol1 = sol1 + a(i1) * pp
	       sol2 = sol2 + b(i1) * pp
	       sol3 = sol3 + c(i1) * pp
	     enddo
	     ind=ind+1
	     vx(ind)=x
	     vy(ind)=y
	     va(ind)=sol1
	     vb(ind)=sol2
	     vc(ind)=sol3
	   enddo
	   enddo
	enddo
c
	n=ind
c
	open(1,file=fbin,form='unformatted')
	write(1) (va(i),i=1,n)
	write(1) (vb(i),i=1,n)
	write(1) (vc(i),i=1,n)
	close(1)
c
      open(2,file='xy.bin',form='unformatted')
	write(2) (vx(i),vy(i),i=1,n)
	close(2)
c
	return
	end
c
