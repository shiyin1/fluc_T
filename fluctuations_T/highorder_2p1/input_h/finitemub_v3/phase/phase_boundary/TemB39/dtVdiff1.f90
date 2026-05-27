subroutine dtVdiff1(kk,yflowk,rhok,l,lb,dtVd1l,dtVd1lb)
!This subroutine calculate various differential kernal corresponding 0<k<k_UV
!wenrui changed 2017.11.16

  implicit none

  real(16) pi
  parameter(pi=3.141592653589793238462643383279Q+00)
  real(16) kk,yflowk(50),rhok(2)
  real(16) l,lb,l_con,lb_con
  real(16) k,hl,hs,etaphi,etapsi,rho(2)
  real(16) sigmal,sigmas
  real(16) mf2(3)
  real(16) T,mu,muq,mus
  real(16) muF(3)
  real(16) xff(3),xfa(3)
  real(16) nf0(3),nf1(3),nf2(3)
  real(16) Nc
  parameter(Nc=3.Q+00)
  real(16) nfd1l(3),nfd1lb(3)
  real(16) nfd1lf(3),nfd1lbf(3),nfd1la(3),nfd1lba(3)
  real(16) dtVd1l,dtVd1lb
  real(16) h_mod,hs_mod
  REAL*16 :: h_const,lamd10_const,lamd20_const,lamd01_const,Tcglue_const,hs_const
  REAL*16 :: SLSS_IR(2)

  common /Tmu/ T,mu,muq,mus
  common /ini_const/ h_const,lamd10_const,lamd20_const,lamd01_const,Tcglue_const,hs_const
  COMMON /IR_RESULT/ SLSS_IR

  l_con=l
  lb_con=lb

  k=kk
!  hl=yflowk(21)
!  hs=yflowk(22)

  call fityukawa(k,h_mod,hs_mod)

  hl=h_mod*h_const
  hs=hs_mod*hs_const
!  hs=h_mod*h_const

  etaphi=0.Q+00
  etapsi=0.Q+00

  rho=rhok
!  sigmal=Sqrt(2/3.Q+00)*Sqrt(2*rho(1) - Sqrt(6.Q+00)*Sqrt(rho(2)))
!  sigmas=Sqrt(2/3.Q+00)*Sqrt(rho(1) + Sqrt(6.Q+00)*Sqrt(rho(2)))

  sigmal=SLSS_IR(1)
  sigmas=SLSS_IR(2)
  mf2(1)=(hl*sigmal/2.Q+00)**2
  mf2(2)=(hl*sigmal/2.Q+00)**2
  mf2(3)=(hs*sigmas)**2/2.Q+00

  muF=(/mu+2.Q+00/3.Q+00*muq,mu-1.Q+00/3.Q+00*muq,mu-1.Q+00/3.Q+00*muq-mus/)

  xff=-muF+ Sqrt(k**2 + mf2)
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
  dtVd1l =((1-etapsi/4.Q+00)*k**5*Nc*Sum((nfd1la +nfd1lf)/Sqrt(k**2+mf2)))/(3.Q+00*Pi**2)
  dtVd1lb=((1-etapsi/4.Q+00)*k**5*Nc*Sum((nfd1lba +nfd1lbf)/Sqrt(k**2+mf2)))/(3.Q+00*Pi**2)
end
