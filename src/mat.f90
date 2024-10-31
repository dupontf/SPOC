c
c***************************************************************
c calcul de la normal pour les faces des triangles
c***************************************************************
c
c
      subroutine cal_rotnf(nf,x,y,norm,pseg)
      implicit none
	integer nf
      integer i,ne1,ne2,pseg(2,*)
      double precision norm(2,*),x(*),y(*),nx,ny,d
c
      do i=1,nf
	 ne1=pseg(1,i)
	 ne2=pseg(2,i)
	 nx= (y(ne2)-y(ne1))
	 ny=-(x(ne2)-x(ne1))
	 norm(1,i)= nx
	 norm(2,i)= ny
      enddo
c
      return
      end
c
c
c******************************************************************
c
	subroutine cal_its(nf,itm,its,tseg,ssit)
	implicit none
	integer nf
	integer i,j,tk
	integer itm(3,*),tseg(2,*),its(3,*),ssit(2,*)
c
        do i=1,nf
         ssit(1,i)=0
         ssit(2,i)=0
        enddo
c
	do i=1,nf
	   tk = tseg(1,i)
	   j = 1
	   do while (itm(j,tk).ne.i.and.j.lt.4)
	    j = j + 1
	   enddo
	   if (j.eq.4) then
	      write(*,*) 'problem in cal_its. bad connectivity table'
		write(*,*) 'face',i,' triangle',tk,' index in triangle',j
		stop
	   endif
	   its(j,tk) = 1
         ssit(1,i)=j
	   tk = tseg(2,i)
	   if (tk.gt.0) then
	   j = 1
	   do while (itm(j,tk).ne.i)
	    j = j + 1
	   enddo
	   its(j,tk) = -1
         ssit(2,i)=j
	   endif
	enddo
c
	return
	end
c
c
c******************************************************************
c
	subroutine coriolis(nc,ne,in,xgr,ygr,fc,f0,beta)
	implicit none
	include '../include/sp.inc'
	double precision xgr(*),ygr(*),fc(nnmod,*),f0,beta
	integer nc,ne,i,n1,n2,n3,in(3,*)
	do i=1,ne
	   n1 = in(1,i)
	   n2 = in(2,i)
	   n3 = in(3,i)
	   fc(1,i) = f0 + 0.5d0 * beta *( ygr(n3) + ygr(n2))
	   fc(2,i) =      0.5d0 * beta * ( ygr(n3) - ygr(n1))
	   fc(nc+2,i) =   0.5d0 * beta * ( ygr(n2) - ygr(n1))
	enddo
	return
	end
c
c
c******************************************************************
c integration dans triangle
c******************************************************************
c
	subroutine cal_mat(nnmod,ng,nc,mat,ff,xi,yi,wi)
	implicit none
	integer nnmod,ng
	double precision mat(nnmod,nnmod),xi(*),yi(*),wi(*)
	double precision xa,xb,ans,x,ff
	integer nc,nm
	integer kc1,kc2,kc3,kc4
	common /i4_cheb/ kc1,kc2,kc3,kc4
	integer i,j,i1,i2
	external xa,xb,ff
c
c integration dans triangle
c
	nm = (nc+1)*(nc+2)/2
c
	do j=1,nm
	do i=1,nm
	   mat(i,j)=0.d0
	enddo
	enddo
c
	i1=0
	  do kc1=0,nc
	  do kc2=0,nc-kc1
	    i1=i1+1
	    i2=0
	    do kc3=0,nc
	    do kc4=0,nc-kc3
	       i2=i2+1
	       call inttri(ng,ff,wi,xi,yi,ans)
	       mat(i1,i2)=ans
	    enddo
	    enddo
	  enddo
	  enddo
c
c dectection des zeros
c
c	do j=1,nm
c	  do i=1,nm
c	   if (abs(mat(i,j)).lt.1.d-6) mat(i,j)=0.d0
c	  enddo
c	enddo
c
c       write(*,*) 'mat'
c       do i=1,nm
c	  write(*,20) (mat(i,j),j=1,nm)
c       enddo
c
20      format(1x,20(f8.4,1x))
21      format(1x,4(i3,1x),f13.10,1x,i4,1x,i2)
	return
	end
c
c
c******************************************************************
c
      subroutine cal_bt(b,ak,mat,coeff)
      implicit none
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision b(*),ak(*),mat(nnmod,*),ans
      double precision coeff
      integer i,j
c
c integration dans triangle
c
	do i=1,nm
	   ans = 0.d0
	   do j=1,nm
	      ans = ans + mat(i,j) * ak(j)
	   enddo
	   b(i) = b(i) + ans * coeff
	enddo
c
	return
	end
c
c
c******************************************************************
c
      subroutine tri_f(tk,itm,its,fbdy,b,coeff)
      implicit none
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision wif(ngf_max),tif(ngf_max)
      common /gaussf/ wif,tif
      double precision pico(ngf_max,nnmod,3),
     &    b(*),fbdy(ngf_max,*),
     &    ans,var,coeff
      integer itm(3,*),its(3,*)
      integer tk,i,j,i1,k,kc1,kc2,j1
	integer isg,ns1,isg2
      double precision sig
      integer cote
      common /gauss3/ pico
c
	do k=1,3
c
	 ns1= itm(k,tk)
	 isg = its(k,tk)
	 isg2= (ngf+1)*(1-isg)/2
	 sig = dfloat(isg) * coeff *0.5d0
c
	  do i1=1,nm
	    ans = 0.d0
	    do j=1,ngf
	       j1 = isg2+isg*j
	       var = pico(j1,i1,k) * fbdy(j,ns1)
	       ans = ans + wif(j) * var
	    enddo
	    b(i1) = b(i1) + ans * sig
	  enddo
c
	enddo
c
	return
	end
c
c
c********************************************************************
c
      subroutine cal_pico(pico,ti)
      implicit none
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision ti(*),pico(ngf_max,nnmod,3),po,x,y
      external po
      integer i1,kc1,kc2,i,j
c
	nm=(nc+1)*(nc+2)/2
	i1=0
	do kc1=0,nc
	   do kc2=0,nc-kc1
	    i1=i1+1
	    do j=1,ngf
	       x = ti(j)
	       y = -1.d0
	       pico(j,i1,1) = po(kc1,x)*po(kc2,y)
	    enddo
	   enddo
	enddo
	i1=0
	do kc1=0,nc
	   do kc2=0,nc-kc1
	    i1=i1+1
	    do j=1,ngf
	       x = -ti(j)
	       y =  ti(j)
	       pico(j,i1,2) = po(kc1,x)*po(kc2,y)
	    enddo
	   enddo
	enddo
	i1=0
	do kc1=0,nc
	   do kc2=0,nc-kc1
	    i1=i1+1
	    do j=1,ngf
	       x = -1.d0
	       y = -ti(j)
	       pico(j,i1,3) = po(kc1,x)*po(kc2,y)
	    enddo
	   enddo
	enddo
c
	return
	end	
c
c
c**************************************************
c integration des polynomes sur le triangle
c**************************************************
c
	subroutine inttri(ng,ff,wi,xi,yi,ans)
	implicit none
	integer ng
	double precision ff,xi(*),yi(*),wi(*),ans,x,y,ansx
	external ff
	integer i,j
	ans=0.d0
	do j=1,ng
	  y = yi(j)
	  x = xi(j)
	  ans = ans + wi(j) * ff(x,y)
	enddo
	return
	end
c
c
c***********************************************************
c
c
	subroutine cal_matl(nn,nc,amm,amml)
	implicit none
	integer nn,nc,nm
	double precision amm(nn,*),amml(nn,*)
	integer i,j,k,l,i1,i2,j1,j2,i0,j0,ip
	external ip
c
	nm = (nc+1)*(nc+2)/2
c
	do i=0,nc-1
	   do j=0,nc-1-i,2
	      do k=0,nc-1
		   do l=0,nc-1-k,2
		 i1 = ip(nc-1,i)+j+1
		   i2 = ip(nc  ,i)+j+2
		 j1 = ip(nc-1,k)+l+1
		   j2 = ip(nc  ,k)+l+2
		   i0 = ip(nc  ,i)+1
		   j0 = ip(nc  ,k)+1
		 amml(i1,j1) = 
     &     amm(i2,j2) + amm(i0,j2) + amm(i2,j0) + amm(i0,j0)
		 enddo
	      enddo
	   enddo
	enddo
	do i=0,nc-1
	   do j=1,nc-1-i,2
	      do k=0,nc-1
		   do l=0,nc-1-k,2
		 i1 = ip(nc-1,i)+j+1
		   i2 = ip(nc  ,i)+j+2
		 j1 = ip(nc-1,k)+l+1
		   j2 = ip(nc  ,k)+l+2
		   i0 = ip(nc  ,i)+1
		   j0 = ip(nc  ,k)+1
		 amml(i1,j1) = 
     &     amm(i2,j2) - amm(i0,j2) + amm(i2,j0) - amm(i0,j0)
		 enddo
	      enddo
	   enddo
	enddo
	do i=0,nc-1
	   do j=0,nc-1-i,2
	      do k=0,nc-1
		   do l=1,nc-1-k,2
		 i1 = ip(nc-1,i)+j+1
		   i2 = ip(nc  ,i)+j+2
		 j1 = ip(nc-1,k)+l+1
		   j2 = ip(nc  ,k)+l+2
		   i0 = ip(nc  ,i)+1
		   j0 = ip(nc  ,k)+1
		 amml(i1,j1) = 
     &     amm(i2,j2) + amm(i0,j2) - amm(i2,j0) - amm(i0,j0)
		 enddo
	      enddo
	   enddo
	enddo
	do i=0,nc-1
	   do j=1,nc-1-i,2
	      do k=0,nc-1
		   do l=1,nc-1-k,2
		 i1 = ip(nc-1,i)+j+1
		   i2 = ip(nc  ,i)+j+2
		 j1 = ip(nc-1,k)+l+1
		   j2 = ip(nc  ,k)+l+2
		   i0 = ip(nc  ,i)+1
		   j0 = ip(nc  ,k)+1
		 amml(i1,j1) = 
     &     amm(i2,j2) - amm(i0,j2) - amm(i2,j0) + amm(i0,j0)
		 enddo
	      enddo
	   enddo
	enddo
c
c       write(*,*) 'mat'
c       do i=1,nm-nc-1
c          write(*,20) (amml(i,j),j=1,nm-nc-1)
c       enddo
20      format(20(f13.10,1x))
24      format(4(i2),20(f13.10,1x))
26      format(6(i2),20(f13.10,1x))
c
	return
	end
c
c*****************************************************************
c
	function ip(nc,i)
	integer ip,nc,i,j
	ip = 0
	j=0
	do while(j.lt.i)
	ip = ip + nc+1-j
	j=j+1
	enddo
	return
	end
c
	subroutine calip(nc,cip)
	implicit none
	include '../include/sp.inc'
	integer ip,nc,i,j,cip(2,0:nnmod),nc2
c
	do i=0,nc
	  ip = 0
	  j=0
	  do while(j.lt.i)
	    ip = ip + nc+1-j
	    j=j+1
	  enddo
	  cip(1,i)=ip
	enddo
c
	nc2=nc-1
	do i=0,nc2
	  ip = 0
	  j=0
	  do while(j.lt.i)
	    ip = ip + nc2+1-j
	    j=j+1
	  enddo
	  cip(2,i)=ip
	enddo
c
	return
	end
c
c
c*****************************************************************
c routines simples
c------------------------------------------------
c reduction de gn
c
	subroutine reduc_gn(nc,gn,gnr,cip)
	implicit none
	include '../include/sp.inc'
	integer nc,i,j,k,l,ind,ind0,cip(2,0:nnmod)
	double precision gn(*),gnr(*)
c
	do j=0,nc-1
	   ind0= cip(1,j)+1
	   do l=0,nc-1-j,2
	      k = cip(1,j)+l+2
	      ind = cip(2,j)+l+1
	      gnr(ind) = gn(k) + gn(ind0)
	   enddo
	   do l=1,nc-1-j,2
	      k = cip(1,j)+l+2
	      ind = cip(2,j)+l+1
	      gnr(ind) = gn(k) - gn(ind0) 
	   enddo
	enddo
c
	return
	end
c
c
c------------------------------------------------
c reduction inverse de gn
c
	subroutine reduc_inv_gn(nc,gn,gnr,cip)
	implicit none
	include '../include/sp.inc'
	integer nc,i,j,k,l,ind,nm,cip(2,0:nnmod)
	double precision gn(*),gnr(*)
c
	nm = (nc+1)*(nc+2)/2
c
	do j=0,nc-1
	   k = cip(1,j)+1
	   gn(k)=0.d0
	   do l=0,nc-j-1,2
	      ind = cip(2,j)+l+1
		gn(k) = gn(k) + gnr(ind)
	   enddo
	   do l=1,nc-j-1,2
	      ind = cip(2,j)+l+1
		gn(k) = gn(k) - gnr(ind)
	   enddo
	   do l=0,nc-j-1
	      k = cip(1,j)+l+2
	      ind = cip(2,j)+l+1
		gn(k) = gnr(ind)
	   enddo
	enddo
	gn(nm)=0.d0
c
	return
	end
c
c------------------------------------------------
c rotation du gradient
c
	subroutine rot_vec(nm,gn,gt,gx,gy,dnx,dny)
	integer nm
	double precision gn(*),gt(*),gx(*),gy(*),dnx,dny
c
	do j=1,nm
	   gn(j) = gx(j) * dnx + gy(j) * dny
	   gt(j) = gx(j) * dny - gy(j) * dnx
	enddo
c
	return
	end
c
c
c*************************************************************
c functions
c*************************************************************
c
	function fsol(nc,x,y,b,i)
	implicit none
	include '../include/sp.inc'
	double precision b(nnmod,*)
	double precision fsol,x,y,po
	integer kc1,kc2,nc,i,j
	external po
	j=0
	fsol=0.d0
	  do kc1=0,nc
	  do kc2=0,nc-kc1
	    j=j+1
	    fsol = fsol + b(j,i) * po(kc1,x) * po(kc2,y)
	  enddo
	  enddo
	return
	end
c
	function f2(x,y)
	implicit none
	integer kc1,kc2,kc3,kc4
	common /i4_cheb/ kc1,kc2,kc3,kc4
	double precision f2,x,y,po
	external po
	f2=po(kc1,x)
     &  *po(kc2,y)
     &  *po(kc3,x)
     &  *po(kc4,y)
	return
	end
c
	function dfx(x,y)
	implicit none
	integer kc1,kc2,kc3,kc4
	common /i4_cheb/ kc1,kc2,kc3,kc4
	double precision dfx,x,y,po,dpo
	external po,dpo
	dfx= po(kc1,x)
     &    *po(kc2,y)
     &    *dpo(kc3,x)
     &    *po(kc4,y)
	return
	end
c
	function dfy(x,y)
	implicit none
	integer kc1,kc2,kc3,kc4
	common /i4_cheb/ kc1,kc2,kc3,kc4
	double precision dfy,x,y,po,dpo
	external po,dpo
	dfy= po(kc1,x)
     &    *po(kc2,y)
     &    *po(kc3,x)
     &    *dpo(kc4,y)
	return
	end
c
	function fx(x,y)
	implicit none
	integer kc1,kc2,kc3,kc4
	common /i4_cheb/ kc1,kc2,kc3,kc4
	double precision fx,x,y,po,dpo
	external po,dpo
	fx= po(kc1,x)
     &    *po(kc2,y)
     &    *po(kc3,x)
     &    *po(kc4,y)*x
	return
	end
c
	function fy(x,y)
	implicit none
	integer kc1,kc2,kc3,kc4
	common /i4_cheb/ kc1,kc2,kc3,kc4
	double precision fy,x,y,po,dpo
	external po,dpo
	fy= po(kc1,x)
     &    *po(kc2,y)
     &    *po(kc3,x)
     &    *po(kc4,y)*y
	return
	end
c
