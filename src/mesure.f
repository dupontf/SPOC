c******************************************************************
c mesure de la longueur moyenne d'un maillage
c******************************************************************
	program ch2
	implicit none
c
C
C 2-D array declaration
C
      INCLUDE '../include/maille.inc'
c
	integer nn,ne,nf,nebd,nne
      INTEGER IN(3,NEDIM),tseg(2,nfdim),pseg(2,nfdim),
     &itm(3,nedim),nmit(nndim),tmit(NFTRDIM*NNDIM),
     & imit(NFTRDIM*NNDIM),ssit(2,nfdim),sbd(nbedim),
     & its(3,nedim)
	double precision xgr(nndim),ygr(nndim),dd,sum
c
	integer i,j,k,ind,i1,i2
c
c
c------------------------------------
c lecture du maillage
c------------------------------------
c
	call maille(nn,ne,nf,nebd,in,itm,nmit,tmit,imit,
     1                  pseg,tseg,sbd,xgr,ygr)
      write(*,*) nn,ne,nf
c
	sum = 0.d0
	do i=1,nf
	   i1 = pseg(1,i)
	   i2 = pseg(2,i)
	   dd = sqrt ( (xgr(i1)-xgr(i2))**2 + (ygr(i1)-ygr(i2))**2 )
	   sum = sum + dd
	enddo
c
	sum = sum /dfloat(nf)
	write(*,*) sum,20.d0/sum
c
	end
c
c******************************************************************
