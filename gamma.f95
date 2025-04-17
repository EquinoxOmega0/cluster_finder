


FUNCTION gammp(a,x)
REAL a,gammp,x
REAL gammcf,gamser,gln
if(x.lt.0..or.a.le.0.) WRITE(*,*) 'bad arguments in gammp'
if(x.lt.a+1.)then
call gser(gamser,a,x,gln)
gammp=gamser
else
call gcf(gammcf,a,x,gln)
gammp=1.-gammcf
endif
return
END




FUNCTION gammq(a,x)

REAL a,gammq,x

REAL gammcf,gamser,gln

! IF (x.lt.0..or.a.le.0) WRITE(*,*) 'bad arguments in gammq'
IF (x.lt.a+1.) THEN
call gser(gamser,a,x,gln)
! WRITE(*,*) gamser,a,x,gln,'gammq'
gammq=1.-gamser
! WRITE(*,*) gammq,'gammq'
ELSE
call gcf(gammcf,a,x,gln)
gammq=gammcf
END IF

RETURN
END



SUBROUTINE gser(gamser,a,x,gln)
INTEGER ITMAX
REAL a,gamser,gln,x,EPS
PARAMETER (ITMAX=100,EPS=3.e-7)
LOGICAL looping
INTEGER n
REAL ap,del,summe,gammln
gln=gammln(a)
! WRITE(*,*) gln,a
if (x.le.0) then
if (x.lt.0) WRITE(*,*) 'x < 0 in gser'
gamser=0.
return
endif
ap=a
summe=1./a
del=summe
 looping=.true.
do n=1,ITMAX
 IF (looping) THEN
ap=ap+1
del=del*x/ap
summe=summe+del
if (abs(del).lt.abs(summe)*EPS) looping=.false.
END IF
END DO
 IF (looping) THEN
WRITE(*,*) 'a too large, ITMAX tooo small in gser'
END IF
gamser=summe*exp(-x+a*log(x)-gln)

RETURN
END





SUBROUTINE gcf(gammcf,a,x,gln)
INTEGER ITMAX
REAL a,gammcf,gln,x,EPS,FPMIN
PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)

LOGICAL looping
INTEGER i
REAL an,b,c,del,h,gammln
gln=gammln(a)
b=x+1.-a
 c=1./FPMIN
 d=1./b
 h=d
 looping=.true.
 DO i=1,ITMAX
 IF (looping) THEN
 an=-i*(i-a)
 b=b+2
 d=an*d+b
 if (abs(d).lt.FPMIN) d=FPMIN
 c=b+an/c
 if (ABS(c).lt.FPMIN) c=FPMIN
 d=1./d
 del=d*c
 h=h*del 
 if (abs(del-1.).lt.EPS) looping=.false.
 
 END IF
 END DO
 IF (looping) THEN
 WRITE(*,*) 'a too large, ITMAX to small in gcf'
 END IF
gammcf=exp(-x+a*log(x)-gln)*h
 RETURN
 END
 
 
 FUNCTION gammln(xx)
 REAL gammln,xx
 INTEGER j
 DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
 SAVE cof,stp
 DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,24.01409824083091d0,&
 -1.231739572450155d0,.1208650973866179d-2,-.5395239384953d-5,2.5066282746310005d0/
 x=xx
 y=x
!  WRITE(*,*) x,xx,y,'gammln'
 tmp=x+5.5D0
 tmp=(x+0.5D0)*log(tmp)-tmp
!  WRITE(*,*) tmp,'tmp gammln'
 ser=1.000000000190015D0
 DO j=1,6
 y=y+1.D0
 ser=ser+cof(j)/y
!  WRITE(*,*) y,ser,'y ser gammln'
 END DO
 gammln=tmp+log(stp*ser/x) 
!  WRITE(*,*) gammln,'2 gammln' 
 RETURN
 END
 
 
 