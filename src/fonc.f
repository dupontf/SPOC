c
c**************************************************
c functions
c**************************************************
c
      function fvaguex(x,y)
      implicit none
      double precision fvaguex,x,y,pi,a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
      pi = 4.d0 * atan(1.d0)
      fvaguex = a1 * cos(x * pi / a2 * 2.d0)
      return
      end
c
      function fvaguey(x,y)
      implicit none
      double precision fvaguey,x,y,pi,a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
      pi = 4.d0 * atan(1.d0)
      fvaguey = a1 * cos(y * pi / a2 * 2.d0)
      return
      end
c
      function ftourbu(x,y)
      implicit none
      double precision ftourbu,x,y,a1,a2,a3,a4,a5,b0,rad
      common /fvague/ a1,a2,a3,a4,a5
      rad = ( x-a3 )**2 + ( y-a4 )**2
      ftourbu =  a1 * ( y - a4 ) *  EXP( -a2 * RAD )
      return
      end
c
      function ftourbv(x,y)
      implicit none
      double precision ftourbv,x,y,a1,a2,a3,a4,a5,b0,rad
      common /fvague/ a1,a2,a3,a4,a5
      rad = ( x-a3 )**2 + ( y-a4 )**2
      ftourbv = -a1 * ( x - a3 ) *  EXP( -a2 * RAD )
      return
      end
c
      function ftourbh(x,y)
      implicit none
      double precision ftourbh,x,y,a1,a2,a3,a4,a5,b0,rad
      common /fvague/ a1,a2,a3,a4,a5
      rad = ( x-a3 )**2 + ( y-a4 )**2
      ftourbh = a1 * EXP( -a2 * RAD )
      return
      end
c
      function fsinus(x,y)
      implicit none
      double precision fsinus,x,y,a1,a2,a3,a4,a5,b0,rad,pi
      common /fvague/ a1,a2,a3,a4,a5
      pi = 4.d0 * atan(1.d0)
      rad = sqrt ( ( x-a3 )**2 + ( y-a4 )**2 )
      fsinus = a1 * ( 1.d0 + cos( 2.d0 * pi * rad /a2 ) )
c      if (rad.gt.a2*0.5d0) fsinus=0.d0
      return
      end
c
      function fkelvin(x,y)
      implicit none
      double precision fkelvin,x,y,a1,a2,a3,a4,a5,b0,rad0,rad1,pi,b1
      common /fvague/ a1,a2,a3,a4,a5
      pi = 4.d0 * atan(1.d0)
      rad0 = ( x-a3 )**2
      rad1 = y-a4
c      a5 = 1.d-10
      fkelvin = a1 * exp( - rad0 * a5 ) * exp( - rad1 * a2 )
      return
      end
c
      function fvent1(x,y)
      implicit none
      double precision fvent1,x,y,a1,a2,a3,a4,a5,pi
      common /fvague/ a1,a2,a3,a4,a5
      pi = 4.d0 * atan(1.d0)
      fvent1 = a1 * sin ( pi * y / a2 )
      return
      end
c
      function fvent2(x,y)
      implicit none
      double precision x,y,fvent2,a1,a2,a3,a4,a5
        common /fvague/ a1,a2,a3,a4,a5
      fvent2 = a1 * 0.5d0 * (1.d0-tanh(y/10.d3))
      return
      end
c
      function fvent3(x,y)
      implicit none
      double precision fvent3,x,y,a1,a2,a3,a4,a5,pi
      common /fvague/ a1,a2,a3,a4,a5
      pi = 4.d0 * atan(1.d0)
      fvent3 = a1 * cos ( pi * y / a2 )
      return
      end
c
      function fvent4(x,y)
      implicit none
      double precision fvent4,x,y,a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
      fvent4 = a1 * y / a2 
      return
      end
c
      function fond1(x,y)
      implicit none
      double precision x,y,fond1,hmin,hmax,ls,lc,y0,yc,h0,lx,ly,pi,
     &                 a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
      pi = 4.d0 *atan(1.d0)
      hmin = 20.d0
      hmax = 4000.d0
      lx   = a1
      ly   = a2
      ls   = 10.d3
      lc   = 16.d3
      y0   = 32.d3 - ly * 0.5d0
      yc = y0 - lc * (cos(x/lx*pi))**24
      fond1 = hmin + 0.5d0*(hmax-hmin)*(1.d0+tanh((y-yc)/ls))
      return
      end
c
      function fond2(x,y)
      implicit none
      double precision x,y,fond2
      double precision a1,a2,a3,a4,a5
      common /fvague/ a1,a2,a3,a4,a5
      fond2 = a4 + 0.5d0 * a1 * ( 1.d0 + tanh( ( y - a2 ) * a3) )
      return
      end
c
