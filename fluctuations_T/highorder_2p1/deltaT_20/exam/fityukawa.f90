subroutine fityukawa(x,h,hs)
!This subroutine calculate the kurtosis analytically

  implicit none

  real(16) x,h,hs
  real(16) pdata(10,250)
  real(16) T,mu,muq,mus,muK
  integer iT
  REAL(16) :: p(10)

  common /Tmu/ T,mu,muq,mus,muK,iT
  common /parainput/ pdata

  P=Pdata(:,iT)

  h=1.Q+00!p(10)+p(9)*x+p(8)*x**2+p(7)*x**3+p(6)*x**4+p(5)*x**5+p(4)*x**6+p(3)*x**7+p(2)*x**8+p(1)*x**9

  P=Pdata(:,1)

  hs=1.Q+00!p(10)+p(9)*x+p(8)*x**2+p(7)*x**3+p(6)*x**4+p(5)*x**5+p(4)*x**6+p(3)*x**7+p(2)*x**8+p(1)*x**9
end
