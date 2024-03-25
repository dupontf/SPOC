c
c ------------------------------------------------------
c  recherche des coordonnees locales dans le triangle tk
c  pour le point xdp,ydp
c ------------------------------------------------------
c
      subroutine interlin(itt,itm,tseg,pseg,x,y,tk,xdp,ydp,lab)
c
      implicit none
c
	integer ne
	integer itt(3,*),itm(3,*),tseg(2,*),pseg(2,*)
      double precision x(*),y(*),lab(3),xdp,ydp,
     1     x1,x2,x3,y1,y2,y3,aire,eps,pr,d
      integer i,j,k,ns1,tk,tk1,n1,n2
      logical ext
c
      ext=.true.
      eps=1.0d-9
c
c
c coord local
c
      do while (ext)
c
      x1=x(itt(1,tk))
      y1=y(itt(1,tk))
      x2=x(itt(2,tk))
      y2=y(itt(2,tk))
      x3=x(itt(3,tk))
      y3=y(itt(3,tk))
      aire  = (x2-x3) * (y2-y1) - (x2-x1) *( y2-y3)
      lab(1)= (xdp-x2)*(ydp-y3) - (xdp-x3)*(ydp-y2)
      lab(2)= (xdp-x3)*(ydp-y1) - (xdp-x1)*(ydp-y3)
      lab(3)= (xdp-x1)*(ydp-y2) - (xdp-x2)*(ydp-y1)
c
      do k=1,3
         lab(k)=lab(k)/aire
         if (abs(lab(k)).lt.eps) lab(k)=0.d0
      enddo
c
      if ((lab(1).ge.0).and.(lab(2).ge.0).and.(lab(3).ge.0)) then
         ext=.false.
         tk1=tk
      else
c
         if     (lab(1).lt.0.d0) then
            ns1=itm(2,tk)
         elseif (lab(2).lt.0.d0) then
            ns1=itm(3,tk)
         elseif (lab(3).lt.0.d0) then
            ns1=itm(1,tk)
         endif
c     
         tk1=tseg(1,ns1)
         if (tk1.eq.tk) then
            tk1=tseg(2,ns1)
         endif
c
      endif
c
c      write(*,103) tk,tk1,lab
c
      if (tk1.gt.0) then
         tk=tk1
      else
c         write(*,104) xdp,ydp
         lab(1)=0.d0
         lab(2)=0.d0
         lab(3)=0.d0
         ext=.false.
      endif
c
      enddo
c
 103  format(1x,2(i4,2x),3(e12.6,2x))
 104  format('point a l''exterieur du domaine ',2(1p,e13.6,1x))
      return
      end
