subroutine PolyakovEq(n, x, fvec)

  implicit none

  integer n
  real(16) x(n), fvec(n)                     !x(1)
  real(16) l,lb
  integer i,j
  integer npoint1, ii
  parameter(npoint1=512)
  real(16) w1(npoint1),y1(npoint1) !Guass integral
  integer npoint2
  parameter(npoint2=512)
  real(16) w2(npoint2),y2(npoint2) !Guass integral

  real(16) kA(npoint1),lam1A(npoint1),lam2A(npoint1),lam3A(npoint1),lam4A(npoint1),lam5A(npoint1),lam0A(npoint1),hA(npoint1)
  real(16) ZphiA(npoint1),ZpsiA(npoint1),cA(npoint1),kappaA(npoint1),etaphiA(npoint1),etapsiA(npoint1)
  real(16) kk,lam1k,lam2k,lam3k,lam4k,lam5k,lam0k,hk,Zphik,Zpsik,ck,kappak,etaphik,etapsik,rhok
  real(16) rho0k0,Zphik0
  real(16) rho0
  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) dtVd1l,dtVd1lb
  real(16) Vd1l_1,Vd1lb_1,Vd1l_2,Vd1lb_2,Vd1l_4,Vd1lb_4,Vd1l,Vd1lb
  real(16) T,mu
  real(16) a1,a2,a3,a4,a5,b1,b2,b3,b4,c1,c2,c3,c4,c5,d1,d2,d3,d4,d5,T0,T_r,aPolya,cPolya,dPolya,bPolya


  common /y1w1/ y1, w1
  common /y2w2/ y2, w2
  common /flowconfi/ kA,lam1A,lam2A,lam3A,lam4A,lam5A,lam0A,hA,ZphiA,ZpsiA,cA,kappaA,etaphiA,etapsiA
  common /rho0Zphi/ rho0k0,Zphik0
  common /Tmu/ T,mu


  l=x(1)
  lb=x(2)


  Vd1l_1=0.Q+0
  Vd1lb_1=0.Q+0
  do ii=1, npoint1
    kk=kA(ii)
    lam1k=lam1A(ii)
    lam2k=lam2A(ii)
    lam3k=lam3A(ii)
    lam4k=lam4A(ii)
    lam5k=lam5A(ii)
    lam0k=lam0A(ii)
    hk=hA(ii)
    Zphik=ZphiA(ii)
    Zpsik=ZpsiA(ii)
    ck=cA(ii)
    kappak=kappaA(ii)
    etaphik=etaphiA(ii)
    etapsik=etapsiA(ii)

    rhok=kappak!rho0k0*Zphik/Zphik0

    call dtVdiff1(kk,lam1k,lam2k,lam3k,lam4k,lam5k,lam0k,hk,Zphik,Zpsik,ck,kappak,etaphik,etapsik,rhok,l,lb,dtVd1l,dtVd1lb)

    Vd1l_1=Vd1l_1+w1(ii)/kk*dtVd1l
    Vd1lb_1=Vd1lb_1+w1(ii)/kk*dtVd1lb

  end do

  Vd1l_2=0.Q+0
  Vd1lb_2=0.Q+0
  do ii=1, npoint2
    kk=y2(ii)

    call dtVdiff2(kk,l,lb,dtVd1l,dtVd1lb)

    Vd1l_2=Vd1l_2-w2(ii)/kk*dtVd1l
    Vd1lb_2=Vd1lb_2-w2(ii)/kk*dtVd1lb

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!glue potential
  a1=-44.14Q+0
  a2=151.4Q+0
  a3=-90.0677Q+0
  a4=2.77173Q+0
  a5=3.56403Q+0

  b1=-0.32665Q+0
  b2=-82.9823Q+0
  b3=3.Q+0
  b4=5.85559Q+0

  c1=-50.7961Q+0
  c2=114.038Q+0
  c3=-89.4596Q+0
  c4=3.08718Q+0
  c5=6.72812Q+0

  d1=27.0885Q+0
  d2=-56.0859Q+0
  d3=71.2225Q+0
  d4=2.9715Q+0
  d5=6.61433Q+0

!  T0=208./hc
  T0=270.Q+0/hc

!  T_r=T/T0
  T_r=0.52Q+0*(T-T0)/T0+1.Q+0

  aPolya=(a1+a2/T_r+a3/T_r**2)/(1.Q+0+a4/T_r+a5/T_r**2)
  cPolya=(c1+c2/T_r+c3/T_r**2)/(1.Q+0+c4/T_r+c5/T_r**2)
  dPolya=(d1+d2/T_r+d3/T_r**2)/(1.Q+0+d4/T_r+d5/T_r**2)

  bPolya=b1*T_r**(-b4)*(1.Q+0-exp(b2/T_r**b3))

  Vd1l_4=((3*cPolya*l**2)/2.Q+0 - (aPolya*lb)/2.Q+0 + 2*dPolya*l*lb**2 +             &
    (bPolya*(12*l**2 - 6*lb - 6*l*lb**2))/                                     &
     (1 - 6*l*lb - 3*l**2*lb**2 + 4*(l**3 + lb**3)))*T**4

  Vd1lb_4=(-(aPolya*l)/2.Q+0 + 2*dPolya*l**2*lb + (3*cPolya*lb**2)/2.Q+0 +           &
    (bPolya*(-6*l - 6*l**2*lb + 12*lb**2))/                                    &
     (1 - 6*l*lb - 3*l**2*lb**2 + 4*(l**3 + lb**3)))*T**4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Vd1l=Vd1l_1+Vd1l_2+Vd1l_4
  Vd1lb=Vd1lb_1+Vd1lb_2+Vd1lb_4

  fvec(1)=Vd1l
  fvec(2)=Vd1lb

end





