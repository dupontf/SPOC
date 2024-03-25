	program stat
	real x(100),y(100),a,b,c,d,e,x0,y0
	integer i,j,n
	character*20 namefile
	write(*,*) 'entrer file'
	read(*,20) namefile
20	format(a20)
	open(1,file=namefile,status='old')
        n=0
100     read(1,*,end=110) a,b
           n=n+1
	   x(n) = log(a)
	   y(n) = log(b)
           goto 100
c
c calcul des moyennes
c
110	write(*,*) n,' points lues'
        x0 = 0.
	y0 = 0.
	do i=1,n
         x0 = x0 + x(i)
         y0 = y0 + y(i)
	enddo
	x0 = x0 /float(n)
	y0 = y0 /float(n)
c	write(*,*) 'moy:',x0,y0
c
c on retranche la moyenne des series
c
	do i=1,n
	   x(i) = x(i) - x0
	   y(i) = y(i) - y0
	enddo
c
c calcul du coefficicent
c
	b = 0.
	c = 0.
	do i=1,n
	   b = b + x(i)**2
	   c = c + x(i)*y(i)
	enddo
	write(*,*) 'coefficient de la courbe en log/log:',-c/b
c
	end
