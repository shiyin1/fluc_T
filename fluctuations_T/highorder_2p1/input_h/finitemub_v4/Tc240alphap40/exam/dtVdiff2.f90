subroutine dtVdiff2(kk,l,lb,dtVd1l,dtVd1lb)
!This subroutine calculate various differential kernal corresponding k>k_UV

  implicit none

  real(16) kk
  real(16) l,lb
  real(16) k,etaphi,etapsi
  real(16) mf2(3),mf(3)
  real(16) T,mu,muq,mus
  real(16) muF(3)
  real(16) xff(3),xfa(3)
  real(16) nf0(3),nf1(3),nf2(3)
  real(16) Nc
  parameter(Nc=3.Q+00)
  real(16) pi
  parameter(pi=3.141592653589793238462643383279Q+00)
  real(16) nfd1l(3),nfd1lb(3)
  real(16) nfd1lf(3),nfd1lbf(3),nfd1la(3),nfd1lba(3)
  real(16) dtVd1l,dtVd1lb

  common /Tmu/ T,mu,muq,mus
  common /mfUK/ mf

  k=kk

  etaphi=0.Q+00
  etapsi=0.Q+00

  muF=(/mu+2.Q+00/3.Q+00*muq,mu-1.Q+00/3.Q+00*muq,mu-1.Q+00/3.Q+00*muq-mus/)

  mf2=mf**2.Q+00

  xff=-muF + Sqrt(k**2 + mf2)
  call Fnf0(xff,T,l,lb,nf0,nf1,nf2)

  nfd1l=-(nf2*(-1 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
  nfd1lb=-(nf1*(-2 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))/2.Q+00

  nfd1lf=nfd1l
  nfd1lbf=nfd1lb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  xfa=muF + Sqrt(k**2 + mf2)
  call Fnf0(xfa,T,lb,l,nf0,nf1,nf2)

  nfd1l=-(nf1*(-2 + 3*nf0 + 3*l*nf1 + 3*lb*nf2))/2.Q+00
  nfd1lb=-(nf2*(-1 + 3*nf0 + 3*l*nf1 + 3*lb*nf2))

  nfd1la=nfd1l
  nfd1lba=nfd1lb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dtVd1l =(1.Q+00-etapsi/4.Q+00)*k**5*Nc*Sum((nfd1la +nfd1lf)/Sqrt(k**2+mf2))/(3.Q+00*Pi**2)
  dtVd1lb=(1.Q+00-etapsi/4.Q+00)*k**5*Nc*Sum((nfd1lba +nfd1lbf)/Sqrt(k**2+mf2))/(3.Q+00*Pi**2)
end
