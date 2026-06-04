subroutine nfplmimu(k,mu,T,l,lb,mfermion2,dnfplmudm,dnfmimudm)

  implicit none

  real(16) k 
  real(16) T,mu(3)
  real(16) mfermion2(3)
  real(16) l,lb 
  real(16) dnfplmudm(3,6),dnfmimudm(3,6)
  real(16) x(3)
  real(16) nf0(3),nf1(3),nf2(3)
  real(16) nfd0x(3),nfd1x(3),nfd2x(3),nfd3x(3),nfd4x(3),nfd5x(3)

  x=Sqrt(k**2 + mfermion2)+mu
  call Fnf0(x,T,lb,l,nf0,nf1,nf2)
  call nfdx(T,lb,l,nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x)
  x=Sqrt(k**2 + mfermion2)
  dnfplmudm(:,1)=nfd0x
  dnfplmudm(:,2)=nfd1x/(2.Q+00*x)
  dnfplmudm(:,3)=-(nfd1x - nfd2x*x)/(4.Q+00*x**3)
  dnfplmudm(:,4)=(3*nfd1x - 3*nfd2x*x + nfd3x*x**2)/(8.Q+00*x**5)
  dnfplmudm(:,5)=(-15*nfd1x + x*(15*nfd2x + x*(-6*nfd3x + nfd4x*x)))/(16.Q+00*x**7)
  dnfplmudm(:,6)=(105*nfd1x - 105*nfd2x*x + 45*nfd3x*x**2 -                     &
                 10*nfd4x*x**3 + nfd5x*x**4)/(32.Q+00*x**9)

  x=Sqrt(k**2 + mfermion2)-mu
  call Fnf0(x,T,l,lb,nf0,nf1,nf2)
  call nfdx(T,l,lb,nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x)
  x=Sqrt(k**2 + mfermion2)
  dnfmimudm(:,1)=nfd0x
  dnfmimudm(:,2)=nfd1x/(2.Q+00*x)
  dnfmimudm(:,3)=-(nfd1x - nfd2x*x)/(4.Q+00*x**3)
  dnfmimudm(:,4)=(3*nfd1x - 3*nfd2x*x + nfd3x*x**2)/(8.Q+00*x**5)
  dnfmimudm(:,5)=(-15*nfd1x + x*(15*nfd2x + x*(-6*nfd3x + nfd4x*x)))/(16.Q+00*x**7)
  dnfmimudm(:,6)=(105*nfd1x - 105*nfd2x*x + 45*nfd3x*x**2 -                     &
                 10*nfd4x*x**3 + nfd5x*x**4)/(32.Q+00*x**9)
  end
