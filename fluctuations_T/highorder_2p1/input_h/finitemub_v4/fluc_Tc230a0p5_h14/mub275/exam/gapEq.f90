subroutine gapEqi(n, x, fvec)
!This subroutine is to calculate the partial derivative of Vall at initial
!wenrui wrote 2017.10.16

  implicit none

  integer n
  !n=2
  real(16) x(n), fvec(n)                     
  real(16) Sl,Ss
  real(16) lamd10,lamd20,lamd01
  real(16) ck,jl,js

  common /gapPara/ lamd10,lamd20,lamd01,ck,jl,js
  
  Sl=x(1)
  Ss=x(2)

  fvec(1)=-jl + lamd10*Sl - (ck*Sl*Ss)/Sqrt(2.Q+00) +                           &
    (lamd01*Sl*(Sl**2 - 2*Ss**2))/6.Q+00 + (lamd20*Sl*(Sl**2 + Ss**2))/2.Q+00

  fvec(2)=-js - (ck*Sl**2)/(2.Q+00*Sqrt(2.Q+00)) + lamd10*Ss -                  &
    (lamd01*Ss*(Sl**2 - 2*Ss**2))/3.Q+00 + (lamd20*Ss*(Sl**2 + Ss**2))/2.Q+00  
end
