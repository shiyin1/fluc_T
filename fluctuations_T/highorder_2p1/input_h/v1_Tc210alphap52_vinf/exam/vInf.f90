subroutine vInf(k_UV,Vinfi)
!calculate the contribution to thermaldynamical potential arising from
!momenta above UV cutoff

  implicit none
  real(16) k_UV,Vinfi
  REAL(16) Vinfi2
  real(16) pi
  parameter(pi=3.141592653589793238462643383279Q+00)
  real(16) T,mu,muq,mus
  real(16) muF(3)
  real(16) l_com,lb_com
  real(16) l,lb,l_con,lb_con
  integer  npoint, iin
  parameter(npoint=256)
  real(16) w(npoint),y(npoint) 
  !Guass integral
  real(16) k,k_infi
  real(16) Nc
  parameter(Nc=3.Q+00)
  real(16) nff(3),nfa(3)
  real(16) nf0f(3),nf1f(3),nf2f(3)
  real(16) nf0a(3),nf1a(3),nf2a(3)
  real(16) xff(3),xfa(3)
  real(16) mf2(3),mf(3)
  real(16) mboson2_UV(18)
  real(16) muB(18)
  real(16) x(18)
  real(16) nbb(18),nba(18)

  common /Tmu/ T,mu,muq,mus
  common /polyakov_com/ l_com,lb_com
  common /mfUK/ mf,mboson2_UV


  l=l_com
  lb=lb_com

  muF=(/mu+2.Q+00/3.Q+00*muq,mu-1.Q+00/3.Q+00*muq,mu-1.Q+00/3.Q+00*muq-mus/)

  mf2=mf**2

  l_con=l
  lb_con=lb
  k_infi=k_UV*7.Q+00

  call gauleg(k_UV, k_infi, y, w, npoint)

  Vinfi=0.Q+00
  Vinfi2=0.Q+00
  do iin=1, npoint
  
  k=y(iin)

  xff=-muF + Sqrt(k**2+mf2)
  xfa=muF + Sqrt(k**2+mf2)

  l=l_con
  lb=lb_con

  call Fnf0(xff,T,l,lb,nf0f,nf1f,nf2f)

  nff=nf0f+lb*nf1f+l*nf2f

    l=l_con
    lb=lb_con

  call Fnf0(xfa,T,lb,l,nf0a,nf1a,nf2a)

    l=l_con
    lb=lb_con

    nfa=nf0a+l*nf1a+lb*nf2a

    muB=(/0.Q+00,0.Q+00,muq,muq,mus,mus,muq+mus,muq+mus,0.Q+00,           &
          0.Q+00,0.Q+00,muq,muq,mus,mus,muq+mus,muq+mus,0.Q+00 /)
    x=Sqrt(k**2 + mboson2_UV)+muB
    call Fnb(x,T,nbb)
    x=Sqrt(k**2 + mboson2_UV)-muB
    call Fnb(x,T,nba)
    Vinfi=Vinfi-w(iin)*k**4*Sum((nfa+nff)/Sqrt(k**2 + mf2))
    Vinfi2=Vinfi2-w(iin)*k**4*Sum((nba+nbb)/Sqrt(k**2 + mboson2_UV))
  end do
  Vinfi=(Vinfi2+4*Vinfi*Nc)/(12.Q+0*pi**2)
end
