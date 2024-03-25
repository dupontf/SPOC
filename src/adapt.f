c**************************************************************
c global subroutine for mesh refinment
c**************************************************************
c
      subroutine global_refinment(nn,ne,nf,
     &           u0,v0,h0,w0,z0,hb,fc,
     &           u2,v2,h2,w2,z2,hb2,
     &           xgr,ygr,in,itm,tseg,pseg,nmit,tmit,imit,
     &           norm,vect,its,ssit,normt,bdl,nint,pint,nbdl,pbdl,
     &           nobd,nebd,sbd,cip,scx,scy,typewindx,typewindy,typefond,
     &           fac1,fac2,fac3,amplvent,lvent,lx,ly,time,time0,
     7           dt,g0,hh0,ngraph,pgraph,f0,beta)
c
      implicit none
      include '../include/sp.inc'
      include '../include/maille.inc'
c
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
c
      integer nn,ne,nf,nebd
      integer nnnew,nenew,innew(3,nedim),nemaxnew,nemax
      integer tinter(nedim)
      integer in(3,*),itm(3,*),tseg(2,*),pseg(2,*),nmit,tmit(*),imit(*),
     & ssit(2,*),sbd(*),
     & its(3,*),nbdl,pbdl(*),nint,pint(*)
      logical nobd(*),bdl(*)
      double precision 
     7       u0(nnmod,*),v0(nnmod,*),h0(nnmod,*),
     &       w0(nnmod,*),z0(nnmod,*),hb(nnmod,*),
     &       u2(nnmod,nedim),v2(nnmod,nedim),h2(nnmod,nedim),
     &       w2(nnmod,nedim),z2(nnmod,nedim),hb2(nnmod,nedim),
     &       xgr(*),ygr(*),norm(2,*),normt(2,*),vect(5,*),
     &       fc(nnmod,*)
      double precision scx,scy,time,time0,dt,g0,hh0

      double precision  xgrnew(nndim),ygrnew(nndim),
     &       unew(nnmod,nedim),vnew(nnmod,nedim),hnew(nnmod,nedim),
     &       wnew(nnmod,nedim),znew(nnmod,nedim),hbnew(nnmod,nedim)
      logical conv,deraff,timedir
      logical fderaff(2*nedim),
     &        fconv(nedim),fconv_ab(2*nedim),fraff_ab(2*nedim)

      double precision pi2(ng_max,nnmod),
     &                 pico(ngf_max,nnmod,3),cip(2,0:nnmod)
      common /gauss2/ pi2
      common /gauss3/ pico
      character*10 typewindx,typewindy,typefond
      double precision fac1,fac2,fac3,amplvent,lvent,lx,ly,
     & f0,beta
      double precision a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
      integer toto,tk1,tk,perraff
      double precision fvent1,fvent2,fvent3,fond1
      external fvent1,fvent2,fvent3,fond1
      include '../include/graph.inc'
      integer ngraph
      double precision pgraph(nnmod,ngradim,ngradim)
      double precision 
     &       amm(nnmod,nnmod),
     &       ammo(nnmod,nnmod),
     &       amml(nnmod,nnmod),
     &       ammx(nnmod,nnmod),
     &       ammy(nnmod,nnmod),
     &       agx(nnmod,nnmod),
     &       agy(nnmod,nnmod)
      common /matrices/ amm,amml,agx,agy,ammx,ammy,ammo
      integer i,j,k

      integer nn_ab,ne_ab,ndoub
      integer in_ab(3,2*nedim),tlevel(2*nedim),torigi(2*nedim),
     &tabsol(nedim),tmt(3,2*nedim),tenfan(4,2*nedim)
      double precision xgr_ab(2*nndim),ygr_ab(2*nndim)
      logical tdoubl(2*nedim),tutili(2*nedim)

      common /adap_const/ nn_ab,ne_ab,ndoub
      common /adap_int/ tlevel,torigi,tabsol,in_ab,tmt,tenfan
      common /adap_vec/ xgr_ab,ygr_ab
      common /adap_log/ tdoubl,tutili
c
c----------------------------------------------------------------
c save the previous fields
c
      call save_temp('temp1.bin',ne,nm,time,u0,v0,h0,w0,z0,hb)
      call save_mesh('cm1.dat',nn,ne,in,xgr,ygr,scx,scy)
      call save_adap('raff1.bin',ne)
c
c------------------------------------------------------------
c test the convergence criteria on the solution
c
      call cal_conv(nn,ne,nf,nc,pico,ssit,conv,fconv,fderaff,timedir,
     &           tseg,u0,v0,h0,w0,z0,hb,deraff,fac1,fac2,fac3)
c
	do i=1,ne_ab
	   fconv_ab(i)=.false.
	   fraff_ab(i)=.false.
	enddo
	do i=1,ne
	   tk=tabsol(i)
	   if (tk.gt.ne_ab) then
		  tk1 = torigi(tk)
		  tutili(tk1)=.true.
	   else
		  tk1 = tk
	   endif
	   if (fconv(i))   fconv_ab(tk1)=.true.
	   if (fderaff(i)) fraff_ab(tk )=.true.
	enddo
c
c
c
c test si il y a changement
c
c
	if (conv) then
c
c
         call deraffiner(nn,ne,fconv_ab,fderaff)
         call adapmaille(fconv_ab,nnnew,nenew,nemaxnew)
         call protridouble(nnnew,nenew,nemaxnew)
         call elimpptri(nnnew,nenew,xgrnew,ygrnew,innew,nemaxnew)
	 call inter_only_new(ne,nenew,in,innew,tinter)
c
         write(*,*) 'nouveau maillage',nnnew,nenew
c
c----------------------------------------------------------------------
c calcul des variables sur la nouvelle grille
c----------------------------------------------------------------------
c si l'erreur n'est pas trop grande, on continue l'integration en temps;
c autrement, on revient a la sauvegarde anterieur
c
	if (.not.timedir) then
	 time=time0
	 do i=1,nenew
	   do j=1,nm
	      u0(j,i)=u2(j,i)
	      v0(j,i)=v2(j,i)
	      h0(j,i)=h2(j,i)
	      w0(j,i)=w2(j,i)
	      z0(j,i)=z2(j,i)
	      hb(j,i)=hb2(j,i)
	   enddo
	 enddo
	endif
c
c-----------------------------------------------------------------
c interpolate the fields and replace the current mesh by the new one
c
        call newinter(nn,ne,nf,u0,v0,h0,w0,z0,hb,
     &           nnnew,nenew,innew,xgrnew,ygrnew,
     &           xgr,ygr,in,itm,tseg,pseg,
     &           typewindx,typewindy,typefond,
     &           amplvent,lvent,lx,ly,g0,hh0,f0,beta)
c
c      call checkcfl(nc,ne,nnmod,hbnew,in,dt,g0,xgr,ygr)
	call checkcfl2(ne,hbnew,unew,vnew,in,dt,g0,xgr,ygr,
     &pgraph,ngradim,ngraph)
c
c
        call pointer_mesh(nn,ne,nf,in,itm,nmit,tmit,imit,pseg,tseg)
        call cal_rotnf(nf,xgr,ygr,norm,pseg)
        call cal_its(nf,itm,its,tseg,ssit)
        call calpre(nn,ne,nf,in,xgr,ygr,vect,itm,norm,normt,bdl,
     &                  nint,pint,nbdl,pbdl,nobd,nebd,sbd,tseg,pseg)
	call coriolis(nc,ne,in,xgr,ygr,fc,f0,beta)
	call verif_rot(ne,u0,v0,normt,ammo,nint,pint,nbdl,pbdl,cip)
c
        call save_temp('temp2.bin',ne,nm,time,u0,v0,h0,w0,z0,hb)
        call save_mesh('cm2.dat',nn,ne,in,xgr,ygr,scx,scy)
        call save_adap('raff2.bin',ne)
c
c------------------------------------------------------------
       endif
c
	  time0=time
	  do i=1,ne
	   do j=1,nm
	      u2(j,i)=u0(j,i)
	      v2(j,i)=v0(j,i)
	      h2(j,i)=h0(j,i)
	      w2(j,i)=w0(j,i)
	      z2(j,i)=z0(j,i)
	      hb2(j,i)=hb(j,i)
	   enddo
	  enddo
      return
      end
c
c**************************************************************
c global subroutine for mesh refinment
c**************************************************************
c
      subroutine newinter(nn,ne,nf,u0,v0,h0,w0,z0,hb,
     &           nnnew,nenew,innew,xgrnew,ygrnew,
     &           xgr,ygr,in,itm,tseg,pseg,
     &           typewindx,typewindy,typefond,
     &           amplvent,lvent,lx,ly,g0,hh0,f0,beta)
c
      implicit none
      include '../include/sp.inc'
      include '../include/maille.inc'
c
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
c
      integer nn,ne,nf
      integer nnnew,nenew,innew(3,nedim),nemaxnew,nemax
      integer tinter(nedim)
      integer in(3,*),itm(3,*),tseg(2,*),pseg(2,*)
      double precision 
     7       u0(nnmod,*),v0(nnmod,*),h0(nnmod,*),
     &       w0(nnmod,*),z0(nnmod,*),hb(nnmod,*),
     &       xgr(*),ygr(*)

      double precision  xgrnew(nndim),ygrnew(nndim),
     &       unew(nnmod,nedim),vnew(nnmod,nedim),hnew(nnmod,nedim),
     &       wnew(nnmod,nedim),znew(nnmod,nedim),hbnew(nnmod,nedim)

      character*10 typewindx,typewindy,typefond
      double precision amplvent,lvent,lx,ly,g0,hh0,f0,beta

      double precision a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
      integer toto,tk1,tk,perraff
      double precision fvent1,fvent2,fvent3,fond1
      external fvent1,fvent2,fvent3,fond1
      integer i,j,k

	call interpol(ne,nn,xgr,ygr,in,tseg,pseg,itm,tinter,
     &                    u0,unew,nenew,nnnew,xgrnew,ygrnew,innew)
	call interpol(ne,nn,xgr,ygr,in,tseg,pseg,itm,tinter,
     &                    v0,vnew,nenew,nnnew,xgrnew,ygrnew,innew)
	call interpol(ne,nn,xgr,ygr,in,tseg,pseg,itm,tinter,
     &                    h0,hnew,nenew,nnnew,xgrnew,ygrnew,innew)
c
c recompute the wind forcing
c
        a1=amplvent
        a2 = lvent
	do i=1,nenew
	   do j=1,nm
	      wnew(j,i)=0.d0
	      znew(j,i)=0.d0
	   enddo
	enddo
	if (typewindx.eq.'fvent1') then
	  call calh2(fvent1,nenew,w0,wnew,innew,xgrnew,ygrnew,tinter)
	elseif (typewindx.eq.'fvent2') then
	  call calh2(fvent2,nenew,w0,wnew,innew,xgrnew,ygrnew,tinter)
	elseif (typewindx.eq.'fvent3') then
	  call calh2(fvent3,nenew,w0,wnew,innew,xgrnew,ygrnew,tinter)
	endif
c
c recompute the bathyemtry
c
        a1 = lx
        a2 = ly
	if (typefond.eq.'fond1') then
	 call calh2(fond1,nenew,hb,hbnew,innew,xgrnew,ygrnew,tinter)
        else
	 do i=1,nenew
	   do j=1,nm
	      hbnew(j,i)=0.d0
	   enddo
	   hbnew(1,i)=hh0
	 enddo
	endif
c
	do i=1,nenew
	   do j=1,nm
	      u0(j,i)=unew(j,i)
	      v0(j,i)=vnew(j,i)
	      h0(j,i)=hnew(j,i)
	      w0(j,i)=wnew(j,i)
	      z0(j,i)=znew(j,i)
	      hb(j,i)=hbnew(j,i)
	   enddo
	enddo
c
c---------------------------------------------------------------------
c l'ancienne grille est ecrasee par la nouvelle;
c le pas de temps est modifie si besoin est;
c creations des pointeurs dans le nouveau maillage;
c application de la condition periodique
c
	nn=nnnew
	ne=nenew
	do i=1,nn
	   xgr(i)=xgrnew(i)
	   ygr(i)=ygrnew(i)
	enddo
	do i=1,ne
	   do j=1,3
	      in(j,i)=innew(j,i)
	   enddo
	enddo

      return
      end
c
c******************************************************************
c
      subroutine calh2(ff,ne,h0,h1,in,xgr,ygr,tinter)
      implicit none
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision xi(ng_max),yi(ng_max),wi(ng_max)
      common /gauss/ wi,xi,yi
      double precision wif(ngf_max),tif(ngf_max)
      common /gaussf/ wif,tif
      integer ne
      double precision xgr(*),ygr(*),h0(nnmod,*),h1(nnmod,*),b(nnmod)
      double precision agx,agy,bgx,bgy,xof,yof,x,y,x1,y1,deta,ans,ansx
      integer i,i1,i2,i3,j,k,l,info,tk
      integer in(3,*),tinter(*)
      double precision ff,fsol,po
      external ff,fsol,po
      double precision a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
      integer kc1,kc2
      double precision pi2(ng_max,nnmod),mff(ng_max)
      double precision amm(nnmod,nnmod),amml(nnmod,nnmod)
      common /gauss2/ pi2
      common /matrices/ amm,amml
c
	do i=1,ne
c
	 if (tinter(i).eq.0) then
c
	   i1 = in(1,i)
	   i2 = in(2,i)
	   i3 = in(3,i)
	   agx = (xgr(i2)-xgr(i1))*0.5d0
	   agy = (ygr(i2)-ygr(i1))*0.5d0
	   bgx = (xgr(i3)-xgr(i1))*0.5d0
	   bgy = (ygr(i3)-ygr(i1))*0.5d0
	   xof = xgr(i1)
	   yof = ygr(i1)
	   deta = agx * bgy - agy * bgx
c
	   do j=1,ng
	      x = xi(j)
	      y = yi(j)
	      x1 = agx * (x+1.d0) + bgx * (y+1.d0) + xof
	      y1 = agy * (x+1.d0) + bgy * (y+1.d0) + yof
	      mff(j) = ff(x1,y1)
	   enddo
c
	  i1=0
	  do kc1=0,nc
	    do kc2=0,nc-kc1
	       i1=i1+1
	       call inttri2(mff,ans,i1)
	       b(i1)=ans
	    enddo
	   enddo
c
         call dtrsm3(nm,1,amm,b)
	   do j=1,nm
	      h1(j,i) = b(j)
	   enddo
c
	 else
	   tk = tinter(i)
	   do j=1,nm
	      h1(j,i)=h0(j,tk)
	   enddo
	 endif
c
	enddo
	return
	end
c
c
c******************************************************************
c routines pour le raffinement de maillage
c******************************************************************
c
      subroutine init_adap(nn,ne,xgr,ygr,in,tseg,itm)
      implicit none
      include '../include/maille.inc'
      double precision xgr(*),ygr(*)
      integer in(3,*),tseg(2,*),itm(3,*)
      integer nn,ne
      integer i,j,ne1,ne2
c
      integer nn_ab,ne_ab,ndoub
      integer in_ab(3,2*nedim),tlevel(2*nedim),torigi(2*nedim),
     &tabsol(nedim),tmt(3,2*nedim),tenfan(4,2*nedim)
      double precision xgr_ab(2*nndim),ygr_ab(2*nndim)
      logical tdoubl(2*nedim),tutili(2*nedim)

      common /adap_const/ nn_ab,ne_ab,ndoub
      common /adap_int/ tlevel,torigi,tabsol,in_ab,tmt,tenfan
      common /adap_vec/ xgr_ab,ygr_ab
      common /adap_log/ tdoubl,tutili
c
	  do i=1,ne
	   tlevel(i)=0
	   torigi(i)=0
	   tutili(i)=.true.
	   tabsol(i)=i
	   tdoubl(i)=.false.
	  enddo
	  ne_ab=ne
	  nn_ab=nn
c
	  do i=1,nn
	     xgr_ab(i)=xgr(i)
	     ygr_ab(i)=ygr(i)
	  enddo
c
	  do i=1,ne
	    do j=1,3
	      in_ab(j,i)=in(j,i)
	    enddo
	  enddo
c
	  ndoub=0
c
c construction tmt
c
	  do i=1,ne
	     do j=1,3
	        ne1=itm(j,i)
		  if (tseg(1,ne1).eq.i) tmt(j,i)=tseg(2,ne1)
		  if (tseg(2,ne1).eq.i) tmt(j,i)=tseg(1,ne1)
	     enddo
c	     write(*,*) (tmt(j,i),j=1,3)
	  enddo
c
c
c construction tenfan
c
	  do i=1,ne
	     do j=1,4
	        tenfan(j,i)=0
	     enddo
c	     write(*,*) (tenfan(j,i),j=1,4)
	  enddo
c
      return
      end
c
c
c******************************************************************
c read the adaptive arrays
c******************************************************************
c
      subroutine read_adap(name,ne)
      implicit none
      include '../include/maille.inc'
      integer ne,ne2
      integer i,j
      character*20 name
c
      integer nn_ab,ne_ab,ndoub
      integer in_ab(3,2*nedim),tlevel(2*nedim),torigi(2*nedim),
     &tabsol(nedim),tmt(3,2*nedim),tenfan(4,2*nedim)
      double precision xgr_ab(2*nndim),ygr_ab(2*nndim)
      logical tdoubl(2*nedim),tutili(2*nedim)

      common /adap_const/ nn_ab,ne_ab,ndoub
      common /adap_int/ tlevel,torigi,tabsol,in_ab,tmt,tenfan
      common /adap_vec/ xgr_ab,ygr_ab
      common /adap_log/ tdoubl,tutili
c
      open(10,file=name,status='old',form='unformatted')
      read(10) nn_ab,ne_ab,ndoub,ne2
      if (ne2.ne.ne) then
         write(*,*) 'pb avec fichier raff.bin'
         stop
      endif
	write(*,*) nn_ab,ne_ab,ndoub
	read(10) (tlevel(i),i=1,ne_ab)
	read(10) (torigi(i),i=1,ne_ab+ndoub)
	read(10) (tabsol(i),i=1,ne)
	read(10) (tutili(i),i=1,ne_ab)
	read(10) (tdoubl(i),i=1,ne_ab+ndoub)
	read(10) ((tmt(j,i),j=1,3),i=1,ne_ab)
	read(10) ((tenfan(j,i),j=1,4),i=1,ne_ab)
	read(10) ((in_ab(j,i),j=1,3),i=1,ne_ab)
	read(10) (xgr_ab(i),ygr_ab(i),i=1,nn_ab)
      close(10)
c	
c	
      return
      end
c
c
c******************************************************************
c save the adaptive arrays
c******************************************************************
c
      subroutine save_adap(name,ne)
      implicit none
      include '../include/maille.inc'
      integer ne
      integer i,j
      character*20 name
c
      integer nn_ab,ne_ab,ndoub
      integer in_ab(3,2*nedim),tlevel(2*nedim),torigi(2*nedim),
     &tabsol(nedim),tmt(3,2*nedim),tenfan(4,2*nedim)
      double precision xgr_ab(2*nndim),ygr_ab(2*nndim)
      logical tdoubl(2*nedim),tutili(2*nedim)

      common /adap_const/ nn_ab,ne_ab,ndoub
      common /adap_int/ tlevel,torigi,tabsol,in_ab,tmt,tenfan
      common /adap_vec/ xgr_ab,ygr_ab
      common /adap_log/ tdoubl,tutili
c
      open(10,file=name,form='unformatted')
	write(10) nn_ab,ne_ab,ndoub,ne
	write(10) (tlevel(i),i=1,ne_ab)
	write(10) (torigi(i),i=1,ne_ab+ndoub)
	write(10) (tabsol(i),i=1,ne)
	write(10) (tutili(i),i=1,ne_ab)
	write(10) (tdoubl(i),i=1,ne_ab+ndoub)
	write(10) ((tmt(j,i),j=1,3),i=1,ne_ab)
	write(10) ((tenfan(j,i),j=1,4),i=1,ne_ab)
	write(10) ((in_ab(j,i),j=1,3),i=1,ne_ab)
	write(10) (xgr_ab(i),ygr_ab(i),i=1,nn_ab)
      close(10)
c	
      return
      end
c
c*********************************************************************
c
      subroutine adapmaille(fconv_ab,nnnew,nenew,nemaxnew)
c
      implicit none
      include '../include/maille.inc'
      INTEGER ppremp(2*nndim),nemaxnew,
     &        peraff(nedim),neraff
      integer nnnew,nenew,nfnew,nlevel,llevel
c
      integer nn_ab,ne_ab,ndoub
      integer in_ab(3,2*nedim),tlevel(2*nedim),torigi(2*nedim),
     &tabsol(nedim),tmt(3,2*nedim),tenfan(4,2*nedim)
      double precision xgr_ab(2*nndim),ygr_ab(2*nndim)
      logical tdoubl(2*nedim),tutili(2*nedim)

      common /adap_const/ nn_ab,ne_ab,ndoub
      common /adap_int/ tlevel,torigi,tabsol,in_ab,tmt,tenfan
      common /adap_vec/ xgr_ab,ygr_ab
      common /adap_log/ tdoubl,tutili
c
      integer i,j,k,l,tk,i1,i2,i3,i4,i5,i6,j1,j2,tk1,tk2,k1,k2,ind
	logical conv,fconv_ab(*),ppelim(2*nndim)
      character*1 trans
	external trans
c
c---------------------------------------------------------------------
	nnnew=nn_ab
	nemaxnew = ne_ab
c
c---------------------------------------------------------------------
c rechercher triangles mitoyens et test sur raffinement
c entre deux trangles, il ne doit pas y avoir de decalage de plus 1
c dans le degre de raffinement
c---------------------------------------------------------------------
c
	call check_level(nemaxnew,tmt,tlevel,fconv_ab)
c
c-----------------------------------------------------------------
c test sur le nombre de triangles a raffiner
c-----------------------------------------------------------------
	neraff=0
	do i=1,ne_ab
	   if (fconv_ab(i)) then 
	      neraff=neraff+1
		peraff(neraff)=i
	   endif
	enddo
	write(*,*) 'nb triangles a raffiner',neraff
c
100	continue
c     write (*,21) (trans(fconv_ab(j)),j=1,ne_ab)
21	format(1x,50(a1),1x)
c
c	
c-----------------------------------------------------------------
c raffinement
c-----------------------------------------------------------------
c
c classement de peraff par level:
	do i = 1,neraff-1
	   do j = i+1,neraff
	      if (tlevel(peraff(i)).gt.tlevel(peraff(j))) then
		   ind = peraff(i)
		   peraff(i)=peraff(j)
		   peraff(j)=ind
		endif
	   enddo
	enddo
c        write(*,*) (peraff(ind),ind=1,neraff)
c
	do ind = 1,neraff
	   i = peraff(ind)
c--------------------------------------------
c-------stop if nnnew grow over nndim
c
	      if (nnnew+3.gt.2*nndim) then 
                write(*,*) 'depassement de 2*nndim',nnnew
		stop
		endif
c--------------------------------------------
c
		tutili(i)=.false.
		i1 = in_ab(1,i)
		i2 = in_ab(2,i)
		i3 = in_ab(3,i)
		nnnew=nnnew+1
		xgr_ab(nnnew)=0.5d0*(xgr_ab(i1)+xgr_ab(i2))
		ygr_ab(nnnew)=0.5d0*(ygr_ab(i1)+ygr_ab(i2))
		i4 = nnnew
		nnnew=nnnew+1
		xgr_ab(nnnew)=0.5d0*(xgr_ab(i2)+xgr_ab(i3))
		ygr_ab(nnnew)=0.5d0*(ygr_ab(i2)+ygr_ab(i3))
		i5 = nnnew
		nnnew=nnnew+1
		xgr_ab(nnnew)=0.5d0*(xgr_ab(i3)+xgr_ab(i1))
		ygr_ab(nnnew)=0.5d0*(ygr_ab(i3)+ygr_ab(i1))
		i6 = nnnew

c--------------------------------------------
c-------stop if nemaxnew grow over nedim
c
	      if (nemaxnew+4.gt.2*nedim) then 
		write(*,*) 'depassement de 2*nedim'
		stop
		endif
c--------------------------------------------
c
	      do j=1,4
		   tenfan(j,i)=nemaxnew+j
		enddo
	      do k=1,4
	       do j=1,4
		   tenfan(j,nemaxnew+k)=0
		 enddo
		enddo
		nlevel = max( nlevel , tlevel(i)+1 )
c
	      nemaxnew=nemaxnew+1
		tlevel(nemaxnew)=tlevel(i)+1
		torigi(nemaxnew)=i
		tutili(nemaxnew)=.true.
	      in_ab(1,nemaxnew) = i1
	      in_ab(2,nemaxnew) = i4
	      in_ab(3,nemaxnew) = i6
		tmt(1,nemaxnew)=tmt(1,i)
		tmt(2,nemaxnew)=nemaxnew+3
		tmt(3,nemaxnew)=tmt(3,i)
c		
	      nemaxnew=nemaxnew+1
		tlevel(nemaxnew)=tlevel(i)+1
		torigi(nemaxnew)=i
		tutili(nemaxnew)=.true.
	      in_ab(1,nemaxnew) = i4
	      in_ab(2,nemaxnew) = i2
	      in_ab(3,nemaxnew) = i5
		tmt(1,nemaxnew)=tmt(1,i)
		tmt(2,nemaxnew)=tmt(2,i)
		tmt(3,nemaxnew)=nemaxnew+2
c		
	      nemaxnew=nemaxnew+1
		tlevel(nemaxnew)=tlevel(i)+1
		torigi(nemaxnew)=i
		tutili(nemaxnew)=.true.
	      in_ab(1,nemaxnew) = i6
	      in_ab(2,nemaxnew) = i5
	      in_ab(3,nemaxnew) = i3
		tmt(1,nemaxnew)=nemaxnew+1
		tmt(2,nemaxnew)=tmt(2,i)
		tmt(3,nemaxnew)=tmt(3,i)
c		
	      nemaxnew=nemaxnew+1
		tlevel(nemaxnew)=tlevel(i)+1
		torigi(nemaxnew)=i
		tutili(nemaxnew)=.true.
	      in_ab(1,nemaxnew) = i5
	      in_ab(2,nemaxnew) = i6
	      in_ab(3,nemaxnew) = i4
		tmt(1,nemaxnew)=nemaxnew-1
		tmt(2,nemaxnew)=nemaxnew-3
		tmt(3,nemaxnew)=nemaxnew-2
c		
c---------------------------------------------------------------------
c verification de tmt. cas de nouveaux triangles mitoyens
c regle: les triangles nouveaux ont un tmt correct
c        les vieux triangles ont une connectivite ancienne
c---------------------------------------------------------------------
c
	     tk2=i
c
           do j=1,3
   		 tk1 = tmt(j,i)
		 j2=j
c
		 if (tk1.gt.0.and.(.not.tutili(tk1))) then
c
c trouver les anciens triangles mitoyens
c
                j1=1
                do while (tmt(j1,tk1).ne.tk2)
   		       j1 = j1 + 1
		    end do
c
                call cit(3,j2,1,k2)
                call cit(3,j1,1,k1)
c
c trouver le triangle remplacant
c
	          tmt(j1,tenfan(j1,tk1))=tenfan(k2,tk2)
	          tmt(j1,tenfan(k1,tk1))=tenfan(j2,tk2)
	          tmt(j2,tenfan(j2,tk2))=tenfan(k1,tk1)
	          tmt(j2,tenfan(k2,tk2))=tenfan(j1,tk1)
c
		  end if
	      end do
	 end do
c
c
c      write(*,*) 'tmt'
c	do i=1,nemaxnew
c	   write(*,*) i,(tmt(j,i),j=1,3)
c	enddo
c
c---------------------------------------------------------------------
c	
	do i=1,nemaxnew
	   fconv_ab(i)=.false.
	enddo
c
c
c---------------------------------------------------------------------
c rechercher triangles mitoyens et test sur raffinement
c un triangle est raffine si deux triangles mitoyens sont raffines
c
	do i=1,nemaxnew
	   if (tutili(i)) then
	   k=0
	   do j=1,3
	      tk=tmt(j,i)
		if (tk.gt.0) then 
		 if (.not.tutili(tk)) then 
		    k=k+1
		 endif
		endif
	   enddo
c	   if (k.ge.2) write(*,*) 'mit',i,k
	   if (k.ge.2) fconv_ab(i)=.true.
	   endif
	enddo
c
c---------------------------------------------------------------------
c rechercher triangles mitoyens et test sur raffinement
c entre deux trangles, il ne doit pas y avoir de decalage de plus 1
c dans le degre de raffinement
c
	call check_level(nemaxnew,tmt,tlevel,fconv_ab)
c
c-----------------------------------------------------------------
c test sur le nombre de triangles a raffiner
c-----------------------------------------------------------------
	neraff=0
	do i=1,nemaxnew
	   if (fconv_ab(i)) then 
	      neraff=neraff+1
		peraff(neraff)=i
	   endif
	enddo
	write(*,*) 'nb triangles a raffiner',neraff
	if (neraff.gt.0) goto 100
c
c
	return
	end
c
c
c*********************************************************************
c creation des triangles coupes en deux
c*********************************************************************
c
      subroutine protridouble(nnnew,nenew,nemaxnew)
c
      implicit none
      include '../include/maille.inc'
      INTEGER ppremp(2*nndim),nemaxnew
      integer nnnew,nenew,nfnew,nlevel,llevel
c
      integer nn_ab,ne_ab,ndoub
      integer in_ab(3,2*nedim),tlevel(2*nedim),torigi(2*nedim),
     &tabsol(nedim),tmt(3,2*nedim),tenfan(4,2*nedim)
      double precision xgr_ab(2*nndim),ygr_ab(2*nndim)
      logical tdoubl(2*nedim),tutili(2*nedim)

      common /adap_const/ nn_ab,ne_ab,ndoub
      common /adap_int/ tlevel,torigi,tabsol,in_ab,tmt,tenfan
      common /adap_vec/ xgr_ab,ygr_ab
      common /adap_log/ tdoubl,tutili
c
      integer i,j,k,l,tk,i1,i2,i3,i4,i5,i6,j1,j2,tk1,tk2,k1,k2
      logical conv,ppelim(2*nndim)
      character*1 trans
      external trans
c
c---------------------------------------------------------------------
c creation des triangles doubles
c---------------------------------------------------------------------
c
	do i=1,nemaxnew
	   tdoubl(i)=.false.
	enddo
c
	ndoub = 0
c
	do i=1,nemaxnew
	   if (tutili(i)) then
	      do j=1,3
		   tk1=tmt(j,i)
c
c--- si i a tlevel plus grand que celui de tk1, alors tk1 coupe en deux
c
		   if (tk1.gt.0.and.tlevel(i).gt.tlevel(tk1)) then
c	write(*,*) i,tk1
c
	          i1 = in_ab(1,tk1)
	          i2 = in_ab(2,tk1)
	          i3 = in_ab(3,tk1)
c
c--------------------------------------------
c-------stop if nnnew grow over nndim
c
	      if (nnnew+2.gt.2*nndim) then 
		write(*,*) 'depassement de 2*nndim'
		stop
		endif
c--------------------------------------------
c--------------------------------------------
c-------stop if nemaxnew + ndoub grow over nedim
c
	      if (nemaxnew+ndoub+3.gt.2*nedim) then 
		write(*,*) 'depassement de 2*nedim'
		stop
		endif
c--------------------------------------------
c tk2: triangle origine de i
c
	          tk2 = torigi(i)
c
c j1 est la face de tk1 donnant sur tk2
c j2 est la face de tk2 donnant sur tk1
c
                j1=1
                do while (tmt(j1,tk1).ne.tk2)
   		       j1 = j1 + 1
		    end do
c
                j2=1
                do while (tmt(j2,tk2).ne.tk1)
   		       j2 = j2 + 1
		    end do
c
                call cit(3,j2,1,k2)
c
                if (tenfan(j2,tk2).eq.i) then
                    i4 = in_ab(k2,i)
c
                 if (j1.eq.1) then
c
		      nnnew=nnnew+1
		      xgr_ab(nnnew)=0.5d0*(xgr_ab(i1)+xgr_ab(i2))
		      ygr_ab(nnnew)=0.5d0*(ygr_ab(i1)+ygr_ab(i2))
		      i4 = nnnew
c
	             ndoub = ndoub + 1
			 k = nemaxnew + ndoub
			 tdoubl(k) = .true.
			 torigi(k) = tk1
	             in_ab(1,k) = i4
	             in_ab(2,k) = i2
	             in_ab(3,k) = i3
			 tutili(tk1)=.false.
			 tutili(k)=.true.
c
                 elseif (j1.eq.2) then
c
		      nnnew=nnnew+1
		      xgr_ab(nnnew)=0.5d0*(xgr_ab(i3)+xgr_ab(i2))
		      ygr_ab(nnnew)=0.5d0*(ygr_ab(i3)+ygr_ab(i2))
		      i4 = nnnew
c
	             ndoub = ndoub + 1
			 k = nemaxnew + ndoub
			 tdoubl(k) = .true.
			 torigi(k) = tk1
	             in_ab(1,k) = i4
	             in_ab(2,k) = i3
	             in_ab(3,k) = i1
			 tutili(tk1)=.false.
			 tutili(k)=.true.
c
                 elseif (j1.eq.3) then
c
		      nnnew=nnnew+1
		      xgr_ab(nnnew)=0.5d0*(xgr_ab(i1)+xgr_ab(i3))
		      ygr_ab(nnnew)=0.5d0*(ygr_ab(i1)+ygr_ab(i3))
		      i4 = nnnew
c
	             ndoub = ndoub + 1
			 k = nemaxnew + ndoub
			 tdoubl(k) = .true.
			 torigi(k) = tk1
	             in_ab(1,k) = i1
	             in_ab(2,k) = i2
	             in_ab(3,k) = i4
			 tutili(tk1)=.false.
			 tutili(k)=.true.
c
		     endif
c
                else if (tenfan(k2,tk2).eq.i) then
                    i4 = in_ab(j2,i)
c
                 if (j1.eq.1) then
c
		      nnnew=nnnew+1
		      xgr_ab(nnnew)=0.5d0*(xgr_ab(i1)+xgr_ab(i2))
		      ygr_ab(nnnew)=0.5d0*(ygr_ab(i1)+ygr_ab(i2))
		      i4 = nnnew
c
	             ndoub = ndoub + 1
			 k = nemaxnew + ndoub
			 tdoubl(k) = .true.
			 torigi(k) = tk1
	             in_ab(1,k) = i1
	             in_ab(2,k) = i4
	             in_ab(3,k) = i3
			 tutili(tk1)=.false.
			 tutili(k)=.true.
c
                 elseif (j1.eq.2) then
c
		      nnnew=nnnew+1
		      xgr_ab(nnnew)=0.5d0*(xgr_ab(i3)+xgr_ab(i2))
		      ygr_ab(nnnew)=0.5d0*(ygr_ab(i3)+ygr_ab(i2))
		      i4 = nnnew
c
	             ndoub = ndoub + 1
			 k = nemaxnew + ndoub
			 tdoubl(k) = .true.
			 torigi(k) = tk1
	             in_ab(1,k) = i1
	             in_ab(2,k) = i2
	             in_ab(3,k) = i4
			 tutili(tk1)=.false.
			 tutili(k)=.true.
c
                 elseif (j1.eq.3) then
c
		      nnnew=nnnew+1
		      xgr_ab(nnnew)=0.5d0*(xgr_ab(i1)+xgr_ab(i3))
		      ygr_ab(nnnew)=0.5d0*(ygr_ab(i1)+ygr_ab(i3))
		      i4 = nnnew
c
	             ndoub = ndoub + 1
			 k = nemaxnew + ndoub
			 tdoubl(k) = .true.
			 torigi(k) = tk1
	             in_ab(1,k) = i4
	             in_ab(2,k) = i2
	             in_ab(3,k) = i3
			 tutili(tk1)=.false.
			 tutili(k)=.true.
c
		     endif
c
		    end if
c
		   endif
		enddo
	   endif
	enddo
c
c
	return
	end
c
c******************************************************************
c simplification de la table de connectivite
c elimination des points en double
c*********************************************************************
c
      subroutine elimpptri(nnnew,nenew,xgrnew,ygrnew,innew,nemaxnew)
c
      implicit none
      include '../include/maille.inc'
      double precision dd
      integer nnnew,nenew,innew(3,*)
      INTEGER ppremp(2*nndim),nemaxnew
c
      integer nn_ab,ne_ab,ndoub
      integer in_ab(3,2*nedim),tlevel(2*nedim),torigi(2*nedim),
     &tabsol(nedim),tmt(3,2*nedim),tenfan(4,2*nedim)
      double precision xgr_ab(2*nndim),ygr_ab(2*nndim)
      logical tdoubl(2*nedim),tutili(2*nedim)

      common /adap_const/ nn_ab,ne_ab,ndoub
      common /adap_int/ tlevel,torigi,tabsol,in_ab,tmt,tenfan
      common /adap_vec/ xgr_ab,ygr_ab
      common /adap_log/ tdoubl,tutili
c
      double precision xgrnew(*),ygrnew(*)
      integer i,j,k,i2
      logical ppelim(2*nndim)
      character*1 trans
      external trans
c
c
c---------------------------------------------------------------------
c verifier que les nouveaux points ne sont pas repetes
c
	do i=1,nnnew
	   ppelim(i)=.false.
	   ppremp(i)=0
	enddo
c
	write(*,*) nn_ab,nnnew
	do i=1,nnnew-1
	   if (.not.ppelim(i)) then
	   i2 = max(nn_ab,i)
	   do j=i2+1,nnnew
	      dd = (xgr_ab(i)-xgr_ab(j))**2+(ygr_ab(i)-ygr_ab(j))**2
		if (dd.lt.1.d-6) then
		   ppelim(j)=.true.
		   ppremp(j)=i
		endif
	   enddo
	   endif
c	   write(*,*) i,ppremp(i),xgr_ab(i),ygr_ab(i)
	enddo
	write(*,*) 'enlev pp'
c	write (*,21) (trans(ppelim(j)),j=1,nnnew)
21	format(1x,50(a1),1x)
c
c-------------------------------------------------------------
c elimination des noeuds en trop
c
	j=0
	do i=1,nnnew
	   if (.not.ppelim(i)) then
	      j=j+1
	      ppremp(i)=j
		xgr_ab(j)=xgr_ab(i)
		ygr_ab(j)=ygr_ab(i)
	   else
	      ppremp(i)=ppremp(ppremp(i))
	   endif
	enddo
	nnnew=j
c--------------------------------------------
c-------stop if nnnew grow over nndim
c
	      if (nnnew.gt.nndim) then 
		write(*,*) 'depassement de nndim'
		stop
		endif
c--------------------------------------------
c
	nn_ab=nnnew
	do i=1,nnnew
	  xgrnew(i)=xgr_ab(i)
	  ygrnew(i)=ygr_ab(i)
	enddo
c
	do i=ne_ab+1,nemaxnew+ndoub
	   do j=1,3
	      in_ab(j,i)=ppremp(in_ab(j,i))
	   enddo
	enddo
	ne_ab=nemaxnew
c
c
c
c----------------------------------------------------------
c elimination des triangles en trop
c
	nenew=0
	do i=1,ne_ab+ndoub
	   if (tutili(i)) then
	      nenew=nenew+1
c--------------------------------------------
c-------stop if nenew grow over nedim
c
	      if (nenew.gt.nedim) then 
		write(*,*) 'depassement de nedim'
		stop
		endif
c--------------------------------------------
c
		do j=1,3
		   innew(j,nenew)=in_ab(j,i)
		enddo
		tabsol(nenew)=i
	   endif
	enddo
c
c
	return
	end
c
c******************************************************************
c
c
c******************************************************************
c
      subroutine cal_conv(nn,ne,nf,nc,pico,ssit,conv,fconv,fderaff,timedir,
     &           tseg,u0,v0,h0,w0,z0,hb,deraff,fac1,fac2,fac3)
c
      implicit none
      include '../include/maille.inc'
      include '../include/sp.inc'
      integer nn,ne,nf
      INTEGER tseg(2,*),ssit(2,*)
      double precision 
     &       u0(nnmod,*),v0(nnmod,*),h0(nnmod,*),
     &       w0(nnmod,*),z0(nnmod,*),hb(nnmod,*)
      integer i,j,k,l
      integer nc,nm
      double precision lx,ly,ampl,lf,x0,y0
	double precision pico(ngf_max,nnmod,3)
	logical conv,fderaff(*),fconv(*),conv1,conv2,conv3,timedir,deraff
	double precision umax(nedim),vmax(nedim),hmax(nedim),supmax
	double precision fac1,fac2,fac3
c
c
	nm = (nc+1)*(nc+2)/2
c
c
	do i=1,ne
	   fconv(i)=.false.
	   fderaff(i)=.false.
	enddo
	conv = .false.
	deraff=.false.
c
c recherche max
c
	call checkconv(nf,ne,u0,pico,tseg,ssit,umax)
	call checkconv(nf,ne,v0,pico,tseg,ssit,vmax)
	call checkconv(nf,ne,h0,pico,tseg,ssit,hmax)
c
	supmax=umax(1)
	do i=1,ne
	   conv1 = (umax(i).gt.fac1)
	   conv2 = (vmax(i).gt.fac1)
	   conv3 = (hmax(i).gt.fac1)
	   if (conv1.or.conv2.or.conv3) then 
	      fconv(i)=.true.
         endif
	   conv1 = (umax(i).lt.fac2)
	   conv2 = (vmax(i).lt.fac2)
	   conv3 = (hmax(i).lt.fac2)
	   if (conv1.and.conv2.and.conv3) then 
	      fderaff(i)=.true.
		deraff=.true.
         endif
	   supmax = max( supmax , umax(i))
	   supmax = max( supmax , vmax(i))
	   supmax = max( supmax , hmax(i))
	enddo
c
c critere pour revenir en arriere dans le temps
c
	if (supmax.gt.fac3) then 
	   timedir=.false.
	else
	   timedir=.true.
	endif
c
	if (supmax.gt.fac1) conv=.true.
	if (deraff) conv=.true.
c
	write(*,*) conv,deraff,timedir
c
	return
	end
c
c
c
c******************************************************************
c prepare the refinment
c******************************************************************
      subroutine preref(ne,fderaff,fconv,
     &           ne_ab,fconv_ab,fraff_ab,tutili,torigi,tabsol)
      implicit none

      integer ne
      integer ne_ab
      integer i,tk,tk1
      integer torigi(*),tabsol(*)
      logical fderaff(*),fconv(*)
      logical fconv_ab(*),fraff_ab(*),tutili(*)
c
	do i=1,ne_ab
	   fconv_ab(i)=.false.
	   fraff_ab(i)=.false.
	enddo
	do i=1,ne
	   tk=tabsol(i)
	   if (tk.gt.ne_ab) then
		  tk1 = torigi(tk)
		  tutili(tk1)=.true.
	   else
		  tk1 = tk
	   endif
	   if (fconv(i))   fconv_ab(tk1)=.true.
	   if (fderaff(i)) fraff_ab(tk )=.true.
	enddo
      return
      end
c
c
c
c******************************************************************
c test sur les elements pour limiter l'interpolation aux elements nouveaux
c******************************************************************
      subroutine inter_only_new(ne,nenew,in,innew,tinter)
      implicit none

      integer ne
      integer nenew
      integer i,j,k,tk
      integer in(3,*),innew(3,*),tinter(*)
c
	do i=1,nenew
	   tinter(i)=0
	enddo
	tk=1
	do i=1,ne
	   k=0
	   do j=1,3
	      if (in(j,i).eq.innew(j,tk)) k=k+1
	   enddo
	   if (k.eq.3) then 
		tinter(tk)=i
	      tk=tk+1
	   endif
	enddo
      return
      end
c
c******************************************************************
c
      subroutine checkconv(nf,ne,a0,pico,tseg,ssit,amax)
      implicit none
      include '../include/sp.inc'
      integer nf,ne
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision a0(nnmod,*),pico(ngf_max,nnmod,3),
     & dif,sum,fbdy1(ngf_max),fbdy2(ngf_max),amax(*),vmax
      integer i,j,k1,k2,cote,i1
      integer tseg(2,*),ssit(2,*)
c
	nm = (nc+1)*(nc+2)/2
c
	do i=1,ne
	 amax(i)=0.d0
	enddo

	do i=1,nf
c
	  do j=1,ngf
	     fbdy1(j)=0.d0
	     fbdy2(j)=0.d0
	  enddo
c
	   k1 = tseg(1,i)
	   k2 = tseg(2,i)
c
	   if (k2.gt.0) then
c
c triangle 1
c
        cote = ssit(1,i)
	  do i1=1,nm
	    do j=1,ngf
	       fbdy1(j) =  fbdy1(j) + a0(i1,k1) * pico(j,i1,cote)
	    enddo
	  enddo
c
c triangle 2
c
	  cote = ssit(2,i)
	  do i1=1,nm
	    do j=1,ngf
	       fbdy2(j) =  fbdy2(j) + a0(i1,k2) * pico(ngf-j+1,i1,cote)
	    enddo
	  enddo
c
	  do j=1,ngf
	    dif = abs( fbdy2(j) - fbdy1(j) )
	    if (amax(k1).lt.dif) amax(k1) = dif
	    if (amax(k2).lt.dif) amax(k2) = dif
	  enddo
c
	 endif
c
	enddo
c
	vmax=abs(a0(1,1))
	do i=2,ne
           vmax=max(abs(a0(1,i)),vmax)
	enddo
	if (vmax.gt.1.d-6) then
	  do i=1,ne
	    amax(i)=amax(i)/vmax
	  enddo
	else
	  do i=1,ne
	    amax(i)=0.d0
	  enddo
	endif
c
	return
	end
c
c**************************************************************
c        
      subroutine cit(n,k,np,k1)
      integer n,k,np,k1
c
      k1=k+np
      if(k1.gt.n) k1=k1-n
      if(k1.lt.1) k1=k1+n
c
      return
      end
c
c
c******************************************************************
c
c
      subroutine check_level(ne_ab,tmt,tlevel,fconv_ab)
      implicit none
      include '../include/maille.inc'
      integer ne_ab,tmt(3,*),tlevel(*)
      logical fconv_ab(*)
      integer pe(nedim),ne
      integer i,j,k,l,ind,tk
c
	ne=0
c
	do i=1,ne_ab
	   if (fconv_ab(i)) then
		ne=ne+1
		pe(ne)=i
	   endif
	enddo
c
	ind=0
c
	do while(ind.lt.ne)
c
	   ind = ind + 1
	   i=pe(ind)
	   do j=1,3
	      tk=tmt(j,i)
		if (tk.gt.0) then
		 if (.not.fconv_ab(tk)) then
		  if (tlevel(tk).lt.tlevel(i)) then 
		     fconv_ab(tk)=.true.
c	           write(*,*) 'lvel',i,tk
		     ne=ne+1
		     pe(ne)=tk
		  endif
	       endif
	      endif
	   enddo
	enddo
c
      return
      end
c
c*********************************************************************
c
      subroutine deraffiner(nn,ne,fconv_ab,fderaff)
c
      implicit none
      include '../include/maille.inc'
      integer nn,ne
      integer nemax,
     &        ppremp(2*nndim),
     &        peremp(0:2*nedim)
      integer nlevel,llevel
c
      integer nn_ab,ne_ab,ndoub
      integer in_ab(3,2*nedim),tlevel(2*nedim),torigi(2*nedim),
     &tabsol(nedim),tmt(3,2*nedim),tenfan(4,2*nedim)
      double precision xgr_ab(2*nndim),ygr_ab(2*nndim)
      logical tdoubl(2*nedim),tutili(2*nedim)

      common /adap_const/ nn_ab,ne_ab,ndoub
      common /adap_int/ tlevel,torigi,tabsol,in_ab,tmt,tenfan
      common /adap_vec/ xgr_ab,ygr_ab
      common /adap_log/ tdoubl,tutili
c
      integer i,j,k,l,tk,i1,i2,i3,i4,i5,i6,j1,j2,tk1,tk2,k1,k2
      logical conv,fconv_ab(*),fderaff(*),ppelim(2*nndim),
     &        fenlev(2*nedim)
c
      character trans*1
      external trans
c
c---------------------------------------------------------------------
c
c
c	nlevel=0
c	do i=1,ne_ab
c	   nlevel = max( nlevel , tlevel(i) )
c	enddo
c      write(*,*) nlevel
c
c---------------------------------------------------------------------
c
c
c	write(*,*) 'test croise avec raff'
c	write (*,21) (trans(fderaff(j)),j=1,ne_ab)
	do i=1,ne_ab
	   if (fderaff(i).and.fconv_ab(i)) fderaff(i)=.false.
	enddo
c	write (*,21) (trans(fderaff(j)),j=1,ne_ab)
c
c
c--- test sur triangle origine----
c
	do i=1,ne_ab
	   fenlev(i)=.false.
	enddo
c
	do i=1,ne_ab
	  k=0
	  if (tenfan(1,i).gt.0) then
	   do j=1,4
	      tk = tenfan(j,i)
	      if (fderaff(tk)) k=k+1
	   enddo
	   if (k.eq.4) then
	    do j=1,4
	      tk = tenfan(j,i)
	      fenlev(tk)=.true.
	    enddo
	       tutili(i)=.true.
		 do j=1,4
		    tenfan(j,i)=0
		 enddo
	   endif
	  endif
	enddo
c
	write(*,*) 'enlev'
c	write (*,21) (trans(fenlev(j)),j=1,ne_ab)
c
c---------------------------------------------------------------------
c verification de tmt. cas de triangles mitoyens de meme ordre
c dont un disparait
c regle: les triangles nouveaux ont un tmt correct
c        les vieux triangles ont une connectivite ancienne
c---------------------------------------------------------------------
c
	do i=1,ne_ab
	  if (.not.fenlev(i)) then
           do j=1,3
   		 tk = tmt(j,i)
		 if (tlevel(tk).eq.tlevel(i).and.fenlev(tk)) then
c
c trouver le triangle remplacant
c
                 tmt(j,i) = torigi(tk)
		 end if
	     end do
        end if
	end do
c
c
c      write(*,*) 'tmt'
c	do i=1,ne_ab
c	   write(*,*) i,(tmt(j,i),j=1,3)
c	enddo
c
c
c--- enlever triangles marques----
c
	tk = 0
	peremp(0)=0
	do i=1,ne_ab
	   if (fenlev(i)) then
		peremp(i)=0
	   else
	      tk = tk + 1
		peremp(i)=tk
	   endif
	enddo
c
	nemax = tk
	write(*,*) 'nemax',nemax
c
	tk=0
	do i=1,ne_ab
	   if (.not.fenlev(i)) then
	      tk = tk+1
c		write(*,*) i,tk
	      do j=1,3
		   in_ab(j,tk)=in_ab(j,i)
		enddo
	      do j=1,3
		   tmt(j,tk)=peremp(tmt(j,i))
		enddo
	      do j=1,4
		   tenfan(j,tk)=peremp(tenfan(j,i))
		enddo
	      torigi(tk)=peremp(torigi(i))
	      tutili(tk)=tutili(i)
	      tlevel(tk)=tlevel(i)
	      fconv_ab(tk)=fconv_ab(i)
	   endif
	enddo
	do i=1,ne
	   tabsol(i)=peremp(tabsol(i))
	enddo
c
	ne_ab = nemax
	write(*,*) 'ne_ab',ne_ab
c
c
c-------- enlever les points en trop---------------
c
	do i=1,nn_ab
	   ppelim(i)=.true.
	enddo
	do i=1,ne_ab
	   do j=1,3
	      tk = in_ab(j,i)
		ppelim(tk)=.false.
	   enddo
	enddo
c	write(*,*) 'enlev point'
c	write (*,21) (trans(ppelim(j)),j=1,nn_ab)
c
	k=0
	do i=1,nn_ab
	   if (.not.ppelim(i)) then
	      k=k+1
	      ppremp(i)=k
		xgr_ab(k) = xgr_ab(i)
		ygr_ab(k) = ygr_ab(i)
	   endif
	enddo
	nn_ab=k
c
	do i=1,ne_ab
	   do j=1,3
	      in_ab(j,i) = ppremp(in_ab(j,i))
	   enddo
	enddo
	write(*,*) 'fin deraff',nn_ab,ne_ab
c
c
c
20	format(1x,200(i4,1x))
21	format(1x,50(a1),1x)
	return
	end
c
c
c*********************************************************************
c
      subroutine testderaff(ne,fderaff,ne_ab,tabsol,
     &           u,v,h,w,z,hb,nc,deraff)
c
      implicit none
      include '../include/maille.inc'
      include '../include/sp.inc'
      integer nn,ne,nc,nm
      integer tabsol(*)
      integer ne_ab
      integer i,j,k,l,tk
      logical fderaff(*),uconv,vconv,deraff
c
      double precision u(nnmod,*),v(nnmod,*),h(nnmod,*),
     &                 w(nnmod,*),z(nnmod,*),hb(nnmod,*),
     &                 a(nnmod),mat,fac1,fac2
c
c---------------------------------------------------------------------
c definite constantes
c fac1 = rapport entre plus grand coefficient polynome et constant
c fac2 = valeur limite de constante
c
	nm = (nc+1)*(nc+2)/2
	fac1= 1.d-6
	fac2= 1.d-10
c
	do i=1,ne_ab
	   fderaff(i)=.false.
	enddo
c
c---------------------------------------------------------------------
c
	do i=1,ne
	   tk = tabsol(i)
	   uconv=.false.
	   vconv=.false.
	   do j=1,nm
	      a(j)=(u(j,i))**2
	   enddo
	   mat = a(1) *fac1
	   if (a(nm).lt.mat.and.a(nc+1).lt.mat) then 
	      uconv=.true.
	   endif
	   if (abs(a(1)).lt.fac2) uconv=.true.
	   do j=1,nm
	      a(j)=(v(j,i))**2
	   enddo
	   mat = a(1) *fac1
	   if (a(nm).lt.mat.and.a(nc+1).lt.mat) then 
	      vconv=.true.
	   endif
	   if (abs(a(1)).lt.fac2) vconv=.true.
	   if (uconv.and.vconv) fderaff(tk)=.true.
	   if (fderaff(tk)) deraff=.true.
	enddo
21	format(10(e13.6,1x))
c
	return
	end
c
c  ----------------------------------------- c
c cree un charactere pour definir la frontiere
c  ----------------------------------------- c
c
      function trans(ival)
      logical ival
      character trans*1
c
      if (.not.ival) trans='.'
      if (ival)      trans='+'
c
      return
      end
c
c
c******************************************************************
c interpolation entre 2 spectral finite elements
c******************************************************************
c
      subroutine interpol(ne,nn,xgr,ygr,in,tseg,pseg,itm,tinter,
     &                    u,unew,nenew,nnnew,xgrnew,ygrnew,innew)
      implicit none
      include '../include/maille.inc'
      include '../include/sp.inc'
      integer nc,nm,ng,ngf
      common /const/ nc,nm,ng,ngf
      double precision xi(ng_max),yi(ng_max),wi(ng_max)
      common /gauss/ wi,xi,yi
c
      integer nn,ne,nf
      INTEGER IN(3,*),tseg(2,*),pseg(2,*),itm(3,*),tinter(*)
      double precision xgr(*),ygr(*)
      double precision amm(nnmod,nnmod),amml(nnmod,nnmod),
     &       u(nnmod,nedim),b(nnmod),mff(ng_max),
     &       pi2(ng_max,nnmod),mat
      integer i,j,k,l,i1,tk,i2,i3,info,tk1,kc1,kc2
      common /gauss2/ pi2
      common /matrices/ amm,amml
c
c------------------------------------
      double precision lab(3),xgrnew(*),ygrnew(*),unew(nnmod,*)
      double precision fsol
      external fsol
      integer innew(3,*),nnnew,nenew
      double precision agx,agy,bgx,bgy,deta,xof,yof,x,y,x0,y0,x1,y1
c
c
C-----------------------------------------------------------------------
c	write(*,*) 'gr'
c
	tk = 1
	do i=1,nenew
c
	 if (tinter(i).eq.0) then
c
	   i1 = innew(1,i)
	   i2 = innew(2,i)
	   i3 = innew(3,i)
	   agx = (xgrnew(i2)-xgrnew(i1))*0.5d0
	   bgx = (xgrnew(i3)-xgrnew(i1))*0.5d0
	   agy = (ygrnew(i2)-ygrnew(i1))*0.5d0
	   bgy = (ygrnew(i3)-ygrnew(i1))*0.5d0
	   xof = xgrnew(i1)
	   yof = ygrnew(i1)
c
	   do j=1,ng
	      y = yi(j)
	      x = xi(j)
	      x0 = agx * (x+1.d0) + bgx * (y+1.d0) + xof
	      y0 = agy * (x+1.d0) + bgy * (y+1.d0) + yof
	      call interlin(in,itm,tseg,pseg,xgr,ygr,tk,x0,y0,lab)
	      x1 = 2.d0 * lab(2) - 1.d0
	      y1 = 2.d0 * lab(3) - 1.d0
	      mff(j) = fsol(nc,x1,y1,u,tk)
	   enddo
c
	  call xytospec(b,mff)
c	  call dpotrs('U',nm,1,amm,nnmod,b,nnmod,info)
        call dtrsm3(nm,1,amm,b)
c
	  do j=1,nm
	     unew(j,i) = b(j)
	  enddo
c
	 else
	  tk1=tinter(i)
	  do j=1,nm
	     unew(j,i) = u(j,tk1)
	  enddo
	 endif
c
      enddo
c
21	format(1x,i3,10(e12.5,1x))
22	format(1x,2(i3,1x),10(e12.5,1x))
	return
	end
c***********************************************************************
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
         n1=pseg(1,ns1)
         n2=pseg(2,ns1)
         x1=x(n1)
         y1=y(n1)
         x2=x(n2)
         y2=y(n2)
         d=(x1-x2)**2+(y1-y2)**2
         pr=(x1-xdp)*(x1-x2)+(y1-ydp)*(y1-y2)
         pr=pr/d
         xdp=x1*(1.d0-pr)+x2*pr
         ydp=y1*(1.d0-pr)+y2*pr
      endif
c
      enddo
c
 103  format(x,2(i4,2x),3(e12.6,2x))
      return
      end
