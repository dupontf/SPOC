!
!******************************************************************
! calcul des flux aux faces
!******************************************************************
!
subroutine cal_flu(fbdy,a0)

      use mesh
      use gauss
      use utilities

      implicit none
      double precision a0(nm,*)
      double precision fbdy(ngf,*)
      integer i,j,k1,k2,i1
      integer cote

	do i=1,nf

	  do j=1,ngf
	     fbdy(j,i)=0.d0
	  enddo

	   k1 = tseg(1,i)
	   k2 = tseg(2,i)

	   if (k2.eq.0) then
!
! la face est sur la frontiere: seul triangle 1
!
	  do i1=1,nm
	    do j=1,ngf
	       fbdy(j,i) =  fbdy(j,i) + a0(i1,k1) * pico(j,i1,1)
	    enddo
	  enddo

	  else
!
! triangle 1
!
        cote = ssit(1,i)
	  do i1=1,nm
	    do j=1,ngf
	       fbdy(j,i) =  fbdy(j,i) + a0(i1,k1) * pico(j,i1,cote)
	    enddo
	  enddo
!
! triangle 2
!
	  cote = ssit(2,i)
	  do i1=1,nm
	    do j=1,ngf
	       fbdy(j,i) =  fbdy(j,i) + a0(i1,k2) * pico(ngf-j+1,i1,cote)
	    enddo
	  enddo
	  do j=1,ngf
	     fbdy(j,i) = 0.5d0 * fbdy(j,i)
	  enddo

	 endif

	enddo

	return
end subroutine cal_flu

!
!******************************************************************
! Add the boundary element terms to the RHS
!******************************************************************
!
subroutine tri_f(tk,fbdy,b,coeff)

      use mesh
      use gauss
      use utilities

      implicit none
      double precision b(*),fbdy(ngf,*),ans,var,coeff
      integer tk,i,j,i1,k,kc1,kc2,j1
      integer isg,ns1,isg2
      double precision sig
      integer cote

	do k=1,3

	 ns1= itm(k,tk)
	 isg = its(k,tk)
	 isg2= (ngf+1)*(1-isg)/2
	 sig = dfloat(isg) * coeff *0.5d0

	  do i1=1,nm
	    ans = 0.d0
	    do j=1,ngf
	       j1 = isg2+isg*j
	       var = pico(j1,i1,k) * fbdy(j,ns1)
	       ans = ans + wif(j) * var
	    enddo
	    b(i1) = b(i1) + ans * sig
	  enddo

	enddo

	return

end subroutine tri_f


!
!*****************************************************************
! transformation for boundary elements for vectors
!------------------------------------------------
! reduction de gn
!
	subroutine reduc_gn(gn,gnr)
	use gauss
	use utilities
	implicit none
	integer i,j,k,l,ind,ind0
	double precision gn(*),gnr(*)

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

	return
	end

!
!------------------------------------------------
! reduction inverse de gn
!
	subroutine reduc_inv_gn(gn,gnr)
	use gauss
	use utilities
	implicit none
	integer i,j,k,l,ind
	double precision gn(*),gnr(*)

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

	return
	end
!
!------------------------------------------------
! rotation du gradient
!
	subroutine rot_vec(nm,gn,gt,gx,gy,dnx,dny)
	integer nm
	double precision gn(*),gt(*),gx(*),gy(*),dnx,dny

	do j=1,nm
	   gn(j) = gx(j) * dnx + gy(j) * dny
	   gt(j) = gx(j) * dny - gy(j) * dnx
	enddo

	return
	end
!******************************************************************
!
      subroutine cal_b(b,ak,mat,coeff)
      use gauss
      implicit none
      double precision b(nm),ak(nm),mat(nm,nm),ans
      double precision coeff
      integer i,j
!
! integration dans triangle
!
	do i=1,nm
	   ans = 0.d0
	   do j=1,nm
	      ans = ans + mat(j,i) * ak(j)
	   enddo
	   b(i) = b(i) + ans * coeff
	enddo

	return
	end
!******************************************************************

