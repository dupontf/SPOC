!
!******************************************************************
!
	subroutine checkcfl(hb,dt,g0)
	use mesh
	use gauss
	implicit none
	double precision hb(nm,ne)
	double precision dt,g0,ax,ay,dd1,dd2,dd3,cc0,cflmin,ra
	integer i,i1,i2,i3

	cflmin=6000.d0
	ra=0.30d0

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

        if (cflmin.ge.100) then 
	     dt = dfloat(10*int(cflmin/10))
        elseif (cflmin.ge.10.and.cflmin.lt.100) then 
	     dt = dfloat(2*int(cflmin/2))
	  elseif (cflmin.lt.10.and.cflmin.ge.1) then
	     dt = dfloat(int(cflmin))
	  elseif (cflmin.lt.1) then
	     dt = dfloat(int(cflmin*10)/10)
	  endif

	write(*,*) 'dtmax', cflmin, dt

	return
	end

!
!******************************************************************
!
      subroutine checkcfl2(hb,u0,v0,dt,g0)
	use mesh
	use gauss
	use graphic
      implicit none
      double precision u0(nm,ne),v0(nm,ne),h0(nm,ne),hb(nm,ne)
      double precision pp, &
           mff(ngraph,ngraph),sol0,&
           mff1(ngraph,ngraph),sol1,&
           mff2(ngraph,ngraph),sol2,maxvel,&
           dt,g0,ax,ay,dd1,dd2,dd3,cc0,cflmin,ra,maxdepth
      integer i,i1,i2,i3,j,l,k

      cflmin=6000.d0
      ra=0.35d0

      do i=1,ne
!
! projection de la profondeur sur la grille reguliere graphique
!
	   do k =1,ngraph
	   do l =1,ngraph-k+1
	     sol0 = 0.d0
	     sol1 = 0.d0
	     sol2 = 0.d0
	     do i1=1,nm
	       pp = pgraph(i1,l,k)
	       sol0 = sol0 + hb(i1,i) * pp
	       sol1 = sol1 + u0(i1,i) * pp
	       sol2 = sol2 + v0(i1,i) * pp
	     enddo
	     mff(l,k) = sol0
	     mff1(l,k) = sol1
	     mff2(l,k) = sol2
	   enddo
	   enddo
!
! recherche max profondeur
!
         maxdepth = 0.d0
         maxvel   = 0.d0
	   do k =1,ngraph
	   do l =1,ngraph-k+1
	      maxdepth = max(maxdepth,mff(l,k))
	      maxvel   = max(maxvel,sqrt(mff1(l,k)**2+mff2(l,k)**2))
	   enddo
	   enddo

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
	   hb(:,i) = 0.d0
	   hb(1,i) = min(dd1,dd2,dd3)/cc0
	   cflmin=min(cflmin,dd1/cc0)
	   cflmin=min(cflmin,dd2/cc0)
	   cflmin=min(cflmin,dd3/cc0)
	enddo

        if (cflmin.ge.100) then 
	     dt = dfloat(10*int(cflmin/10.d0))
        elseif (cflmin.ge.10.and.cflmin.lt.100) then 
	     dt = dfloat(2*int(cflmin*0.5d0))
	  elseif (cflmin.lt.10.and.cflmin.ge.1) then
	     dt = dfloat(int(cflmin))
	  elseif (cflmin.lt.1) then
	     dt = dfloat(int(cflmin*10)/10)
	  endif

	write(*,*) 'dtmax', cflmin, dt

	return
	end

