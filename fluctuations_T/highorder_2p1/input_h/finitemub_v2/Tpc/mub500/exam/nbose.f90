subroutine nbose(k,mu,T,mboson2,dnbodm0,dnbodm1,dnbodm2,dnbodm3,dnbodm4,dnbodm5)

  implicit none

  real(16) k 
  real(16) T,mu(18)
  real(16) mboson2(18)
  real(16) dnbodm0(18),dnbodm1(18),dnbodm2(18),dnbodm3(18),dnbodm4(18),dnbodm5(18)
  real(16) x(18)
  real(16) nb(18)
  real(16) nbd1x(18),nbd2x(18),nbd3x(18),nbd4x(18),nbd5x(18)

  x=Sqrt(k**2 + mboson2)+mu
  call Fnb(x,T,nb)
  call nbdx(nb,T,nbd1x,nbd2x,nbd3x,nbd4x,nbd5x)
  x=Sqrt(k**2 + mboson2)
  dnbodm0=nb
  dnbodm1=nbd1x/(2.Q+00*x)
  dnbodm2=-(nbd1x - nbd2x*x)/(4.Q+00*x**3)
  dnbodm3=(3*nbd1x - 3*nbd2x*x + nbd3x*x**2)/(8.Q+00*x**5)
  dnbodm4=(-15*nbd1x + x*(15*nbd2x + x*(-6*nbd3x + nbd4x*x)))/(16.Q+00*x**7)
  dnbodm5=(105*nbd1x-105*nbd2x*x+45*nbd3x*x**2-10*nbd4x*x**3+nbd5x*x**4)        &
           /(32.Q+00*x**9)
end
