SUBROUTINE gauleg(x1,x2,x,w,n)
implicit  none
real(16) pi
parameter(pi=3.141592653589793238462643383279Q+00)
INTEGER n
real(16) x1,x2,x(n),w(n)
real(16) EPS
PARAMETER (EPS=3.Q-14)
INTEGER i,j,m
real(16) p1,p2,p3,pp,xl,xm,z,z1

m=(n+1)/2
xm=0.5Q+00*(x2+x1)
xl=0.5Q+00*(x2-x1)
do i=1,m
  z=cos(pi*(i-0.25Q+00)/(n+0.5Q+00))
  do
    p1=1.Q+00
    p2=0.Q+00
    do j=1,n
      p3=p2
      p2=p1
      p1=((2.Q+00*j-1.Q+00)*z*p2-(j-1.Q+00)*p3)/j
    end do
    pp=n*(z*p1-p2)/(z*z-1.Q+00)
    z1=z
    z=z1-p1/pp
    if(.not.abs(z-z1)>EPS) exit
  end do
  x(i)=xm-xl*z
  x(n+1-i)=xm+xl*z
  w(i)=2.Q+00*xl/((1.Q+00-z*z)*pp*pp)
  w(n+1-i)=w(i)
end do
END
