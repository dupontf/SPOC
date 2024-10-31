!******************************************************************
! equations shallow waters
! spectral discontinuous elements
! rk4 integration in time
! conditions frontieres pour u 
! viscosite 
! friction
! plan beta 
!******************************************************************
program shwater

 use mesh
 use gauss
 use graphic
 use boundary_condition
 use utilities
 use gauss_bd
 use boundary_bd
 
 implicit none

!
! modes spectraux
!
      double precision, allocatable :: &
       u0(:,:), v0(:,:), h0(:,:), div(:,:), vor(:,:)

      external f2,dfx,dfy,ip,fx,fy
      double precision f2,dfx,dfy,fx,fy

      integer i,j,k,l,i1,i2,i3,ns1
      integer nite,ite,info,per,nper,per2,yearmax,daysmax,nitemax,npermax
      double precision :: &
         dt,mu0,mb0,f0,beta,g0,hh0,time,ke,pe,mass,pe0, &
         lx,ly,time0,amplvent,ax,ay,bx,by,deta,detai,&
         dnx,dny,dis,dt6,dt3,dt2

      character*60 :: format1,format2,format3,filename

!
!------------------------------------
! read input parameter file
!------------------------------------
!
	open(2,file='oc.inp',status='old')
	read(2,*) nc,ng,ngf,ng_bd,ngf_bd
	read(2,*) lx,ly
	read(2,*) yearmax,daysmax,nitemax
	read(2,*) dt
	read(2,*) mu0
	read(2,*) mb0
	read(2,*) f0
	read(2,*) beta
	read(2,*) scx,scy
	read(2,*) g0
	read(2,*) hh0
	read(2,*) npermax,ngraph
	read(2,*) per2
	read(2,*) amplvent
	close(2)
!
!----------------------------------------------
! write message

format1='1x,60(''-'')'
format2='1x,5(1p,e13.6,1x)'
format3='i8,'' |  '',f8.2,'' |  '',i6,''  |   '',1P,E8.1,''  |   '',1P,E8.1'

!
!------------------------------------
! build the basis functions
!------------------------------------
!
      call read_gauss
      call cal_p2
      call cal_pico
      call calip
     
!
! build the matrices
!
      allocate(amm(nm,nm), agx(nm,nm), agy(nm,nm), amml(nm,nm))
      call cal_mat(nm,ng,nc,amm,f2,xi,yi,wi)
      call cal_mat(nm,ng,nc,agx,dfx,xi,yi,wi)
      call cal_mat(nm,ng,nc,agy,dfy,xi,yi,wi)
      call cal_matl(nm,nc,amm,amml)
      allocate(ammo(nm,nm))
      do j=1,nm
         do i=1,nm
	    ammo(i,j)=amm(i,j)
         enddo
      enddo
!
! LU Decomposition of the mass matrices
!
      call dpotrf('U',nm,amm,nm,info)
      if (info.ne.0) write(*,*) 'pb avec matrice'
      call dpotrf('U',nm-nc-1,amml,nm,info)
      if (info.ne.0) write(*,*) 'pb avec matrice'
!
!------------------------------------
! read the mesh and calculate the pointers
!------------------------------------
!
      call read_mesh
      call pointer_mesh
      call check_mesh

      call readbel
      call periodic

      call cal_rotnf
      call cal_its
      call calpre
      

! ********* init of curved elements *************

      call cal_pbdl_inv
      call readcoeff
      
      call read_gauss_bd
      call cal_p2_bd
      call cal_pico_bd
      call cal_normc
      
      call cal_courb_jac
      call cal_mat_courb
      call cal_courb_nx
      call cal_courb_nx2
      
      call openbc1_read
      
!
!------------------------------------
! read initial state and forcing
!
      allocate(u0(nm,ne),v0(nm,ne),h0(nm,ne),div(nm,ne),vor(nm,ne))

      write(*,*) 'enter binary file name'
      read(*,fmt='a60') filename
      open(1,file=filename,form='unformatted',status='old',action='read')
      read(1) time0
      read(1) u0
      read(1) v0
      read(1) h0
      close(1)
      
      write(*,*) 'time',time0

!
!------------------------------------
! compute divergence and vorticity
!
	  call cal_var(u0,v0,h0,div,vor,time)

      write(*,*) 'enter binary file name for writing div/vort'
      read(*,fmt='a60') filename
	   open(2,file=filename,form='unformatted',action='write')
	   write(2) time
	   write(2) ((div(j,i),j=1,nm),i=1,ne)
	   write(2) ((vor(j,i),j=1,nm),i=1,ne)
	   write(2) ((h0(j,i),j=1,nm),i=1,ne)
	   close(2)

      stop
end program shwater
!
!******************************************************************
!
!
!******************************************************************
! nonlinear version
! free-slip
! viscosity
! beta plane
! varying bottom topography
! curved elements
! open boudary
!******************************************************************
!
subroutine cal_var(u0,v0,h0,div,vor,time)
     
    use mesh
    use gauss
    use utilities
    use boundary_condition
    use gauss_bd
    use boundary_bd
    
      implicit none
      double precision u0(nm,ne),v0(nm,ne),h0(nm,ne), &
                       div(nm,ne),vor(nm,ne),&
                 dia1(nm),dia2(nm),&
                 dia3(nm),dia4(nm), dia5(nm), &
                 ak(nm),b(nm),bu(nm),&
                 uk(nm),vk(nm),&
                 gx(nm),gy(nm),gn(nm),gt(nm),&
                 gx1(nm),gy1(nm),&
		 amm0(nm,nm),amm1(nm,nm),amm2(nm,nm), &
		 amm3(nm,nm), amm4(nm,nm), &
! --- some scalar
                 ax,by,ay,bx,deta,detai,dis,dnx,dny,mu,mv,&
                 time,pi,tday, &
! --- some other stuff
                 gxx(nm),gxy(nm),gyx(nm),gyy(nm),&
                 qxx(nm,ne),qxy(nm,ne),&
                 qyx(nm,ne),qyy(nm,ne),&
! --- some gauss tables
                 fbdxx(ngf_bd,nf),fbdxy(ngf_bd,nf),&
                 fbdyx(ngf_bd,nf),fbdyy(ngf_bd,nf),&
                 fbd(ngf_bd,nf)
      integer i,j,k,l,info,i1,i2,i3,ind,ns1,tk
      integer cell,face

      
!----------- matrix inversion---------------

      do i=1,nm
         dia1(i) = 1.d0/amm(i,i)
      enddo
      do i=1,nm-nc-1
         dia2(i) = 1.d0/amml(i,i)
      enddo


!
!------------------------------------
! compute the gradient terms of u,v
!------------------------------------
!

	call cal_flu(fbd,u0)
	call norm_flu(fbd,fbdxx,1)
	call norm_flu(fbd,fbdxy,2)


	call cal_flu(fbd,v0)
	call norm_flu(fbd,fbdyx,1)
	call norm_flu(fbd,fbdyy,2)


	do cell=1,ne

         ax    = vect(1,cell)
         bx    = vect(2,cell)
         ay    = vect(3,cell)
         by    = vect(4,cell)
         detai = vect(5,cell)

	      uk = u0(:,cell)
	      vk = v0(:,cell)

!
!------------------------------------------------
! boundary cell
!------------------------------------------------
!
         tk=pbdl_inv(cell)

	 if (bdl(cell)) then

	 amm0=amm_bd(:,:,tk)
	 amm1=amm_bdl(:,:,tk)
	 amm2=agx_bd(:,:,tk)
	 amm3=agy_bd(:,:,tk)

         do k=1,nm-nc-1
          dia4(k) = 1.d0/amm1(k,k)
         enddo
         do k=1,nm
          dia5(k) = 1.d0/amm0(k,k)
         enddo

!------------------------------------------------
	      gx1=0.d0
	      gy1=0.d0

	   call cal_b(gx1,uk,amm2,-1.d0)
	   call cal_b(gy1,uk,amm3,-1.d0)

	   gxx = (  by * gx1 - ay * gy1) * detai
	   gxy = (- bx * gx1 + ax * gy1) * detai

           call tri_fc(cell,fbdxx,gxx,detai)
  	   call tri_fc(cell,fbdxy,gxy,detai)

!------------------------------------------------
	      gx1=0.d0
	      gy1=0.d0

	   call cal_b(gx1,vk,amm2,-1.d0)
	   call cal_b(gy1,vk,amm3,-1.d0)

	   gyx = (  by * gx1 - ay * gy1) * detai
	   gyy = (- bx * gx1 + ax * gy1) * detai

           call tri_fc(cell,fbdyx,gyx,detai)
  	   call tri_fc(cell,fbdyy,gyy,detai)
!
!------------------------------------------------
!
         call dtrsm4(nm,nm,amm0,gxx,dia5)
         call dtrsm4(nm,nm,amm0,gxy,dia5)
         call dtrsm4(nm,nm,amm0,gyx,dia5)
         call dtrsm4(nm,nm,amm0,gyy,dia5)


!-------------------------------------------------
! cas normal sans rotation
!-------------------------------------------------

	   else

!------------------------------------------------
	      gx1=0.d0
	      gy1=0.d0

	   call cal_b(gx1,uk,agx,-1.d0)
	   call cal_b(gy1,uk,agy,-1.d0)

	   gxx = (  by * gx1 - ay * gy1) * detai
	   gxy = (- bx * gx1 + ax * gy1) * detai

           call tri_f(cell,fbdxx,gxx,detai)
  	   call tri_f(cell,fbdxy,gxy,detai)

!------------------------------------------------
	      gx1=0.d0
	      gy1=0.d0

	   call cal_b(gx1,vk,agx,-1.d0)
	   call cal_b(gy1,vk,agy,-1.d0)

	   gyx = (  by * gx1 - ay * gy1) * detai
	   gyy = (- bx * gx1 + ax * gy1) * detai

           call tri_f(cell,fbdyx,gyx,detai)
  	   call tri_f(cell,fbdyy,gyy,detai)

	   call dtrsm4(nm,nm,amm,gxx,dia1)
	   call dtrsm4(nm,nm,amm,gxy,dia1)
	   call dtrsm4(nm,nm,amm,gyx,dia1)
	   call dtrsm4(nm,nm,amm,gyy,dia1)

	   endif

	   do j=1,nm
		qxx(j,cell) = gxx(j)
		qxy(j,cell) = gxy(j)
		qyx(j,cell) = gyx(j)
		qyy(j,cell) = gyy(j)
	   enddo
	   

	enddo

!------------------------------------------------
! divergence and vorticity
!------------------------------------------------
        div = qxx + qyy
        vor = qyx - qxy
	

	return
	end
!
!*****************************************************************
