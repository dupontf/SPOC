c
c******************************************************************
c
	subroutine checkcfl(nc,ne,nnmod,hb,in,dt,g0,xgr,ygr)
	implicit none
	integer nc,ne,in(3,*),nnmod
	double precision hb(nnmod,*),xgr(*),ygr(*),
     &dt,g0,ax,ay,dd1,dd2,dd3,cc0,cflmin,ra
	logical ok
	integer i,i1,i2,i3
c
	cflmin=6000.d0
	ra=0.30d0
c
	do i=1,ne
	   i1 = in(1,i)
	   i2 = in(2,i)
	   i3 = in(3,i)
	   ax = (xgr(i2)-xgr(i1))
	   ay = (ygr(i2)-ygr(i1))
	   dd1 = sqrt(ax**2+ay**2) /dfloat(nc) * ra
	   ax = (xgr(i3)-xgr(i1))
	   ay = (ygr(i3)-ygr(i1))
	   dd2 = sqrt(ax**2+ay**2) /dfloat(nc) * ra
	   ax = (xgr(i3)-xgr(i2))
	   ay = (ygr(i3)-ygr(i2))
	   dd3 = sqrt(ax**2+ay**2) /dfloat(nc) * ra
	   cc0 = sqrt(g0*hb(1,i))
	   cflmin=min(cflmin,dd1/cc0)
	   cflmin=min(cflmin,dd2/cc0)
	   cflmin=min(cflmin,dd3/cc0)
	enddo
c
        if (cflmin.ge.100) then 
	     dt = dfloat(10*int(cflmin/10))
        elseif (cflmin.ge.10.and.cflmin.lt.100) then 
	     dt = dfloat(2*int(cflmin/2))
	  elseif (cflmin.lt.10.and.cflmin.ge.1) then
	     dt = dfloat(int(cflmin))
	  elseif (cflmin.lt.1) then
	     dt = dfloat(int(cflmin*10)/10)
	  endif
c        ok = .true.
c        do while (ok)
c         if (cflmin.lt.dt) dt = dfloat(2*int(cflmin/2))
c         if (dt.lt.cflmin) ok=.false.
c        enddo
c
	write(*,*) 'dtmax', cflmin, dt
c
	return
	end
c
c
c******************************************************************
c
      subroutine checkcfl2(ne,hb,u,v,in,dt,g0,xgr,ygr,pgraph,ngradim,n)
      implicit none
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      integer ne,in(3,*)
      integer ngradim,n
      double precision hb(nnmod,*),u(nnmod,*),v(nnmod,*),xgr(*),ygr(*),
     &      pp,pgraph(nnmod,ngradim,ngradim),
     &      mff(ngradim,ngradim),sol0,
     &      mff1(ngradim,ngradim),sol1,
     &      mff2(ngradim,ngradim),sol2,maxvel,
     &      dt,g0,ax,ay,dd1,dd2,dd3,cc0,cflmin,ra,maxdepth
      logical ok
      integer i,i1,i2,i3,j,l,k
c
      cflmin=6000.d0
      ra=0.35d0
c
      do i=1,ne
c
c projection de la profondeur sur la grille reguliere graphique
c
	   do k =1,n
	   do l =1,n-k+1
	     sol0 = 0.d0
	     sol1 = 0.d0
	     sol2 = 0.d0
	     do i1=1,nm
	       pp = pgraph(i1,l,k)
	       sol0 = sol0 + hb(i1,i) * pp
	       sol1 = sol1 +  u(i1,i) * pp
	       sol2 = sol2 +  v(i1,i) * pp
	     enddo
	     mff(l,k) = sol0
	     mff1(l,k) = sol1
	     mff2(l,k) = sol2
	   enddo
	   enddo
c
c recherche max profondeur
c
         maxdepth = 0.d0
         maxvel   = 0.d0
	   do k =1,n
	   do l =1,n-k+1
	      maxdepth = max(maxdepth,mff(l,k))
	      maxvel   = max(maxvel,sqrt(mff1(l,k)**2+mff2(l,k)**2))
	   enddo
	   enddo
c
	   i1 = in(1,i)
	   i2 = in(2,i)
	   i3 = in(3,i)
	   ax = (xgr(i2)-xgr(i1))
	   ay = (ygr(i2)-ygr(i1))
	   dd1 = sqrt(ax**2+ay**2) /dfloat(nc) * ra
	   ax = (xgr(i3)-xgr(i1))
	   ay = (ygr(i3)-ygr(i1))
	   dd2 = sqrt(ax**2+ay**2) /dfloat(nc) * ra
	   ax = (xgr(i3)-xgr(i2))
	   ay = (ygr(i3)-ygr(i2))
	   dd3 = sqrt(ax**2+ay**2) /dfloat(nc) * ra
	   cc0 = sqrt(g0*maxdepth)+maxvel
	   cflmin=min(cflmin,dd1/cc0)
	   cflmin=min(cflmin,dd2/cc0)
	   cflmin=min(cflmin,dd3/cc0)
	enddo
c
        if (cflmin.ge.100) then 
	     dt = dfloat(10*int(cflmin/10))
        elseif (cflmin.ge.10.and.cflmin.lt.100) then 
	     dt = dfloat(2*int(cflmin/2))
	  elseif (cflmin.lt.10.and.cflmin.ge.1) then
	     dt = dfloat(int(cflmin))
	  elseif (cflmin.lt.1) then
	     dt = dfloat(int(cflmin*10)/10)
	  endif
c
	write(*,*) 'dtmax', cflmin, dt
c
	return
	end
c
