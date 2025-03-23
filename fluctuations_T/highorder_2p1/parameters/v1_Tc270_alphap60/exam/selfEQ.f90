subroutine selfEQ(l_i,lb_i,l,lb,rho0,mboson,mfermion,Vtotal,fpifK,h,Z_wave,ck)
!This subroutine solve the FRG equation and the equation for the Polyakov loop self-consistently,
!Inputing guess l_i,lb_i, outputing l,lb and other physical variables
!wenrui try to change 2017.11.10

  implicit none
  real(16) l_i,lb_i
  real(16) l,lb
  real(16) rho0(2)
  real(16) mboson(8),mfermion(2)
  real(16) Vtotal,fpifK(2),h(2),Z_wave(4),ck
  real(16) l_new,lb_new,l_old,lb_old
  real(16) l_com,lb_com
  real(16) T,mu,muq,mus
  !temperature and chemical potential
  integer i,j
  integer jmax  
  !maximal number of loops
  parameter(jmax=600)
  INTEGER kmax,kount 
  !variables in common block of subroutine odeint
  INTEGER KMAXX,NMAX
  PARAMETER (NMAX=50,KMAXX=200)
  real(16) dxsav,xp(KMAXX),yp(NMAX,KMAXX) 
  !variables in common block of subroutine odeint
  real(16) x_bef(KMAXX),y_bef(KMAXX) 
  !used for interpolation
  real(16) k_UV,k_IR,t_UV,t_IR
  real(16) epsi_l,epsi_lb,epsi_err
  logical stopp
  real(16) Vall,Vinfi,Vglue
  external PolyakovEq
  integer nn
  parameter(nn=2)
  real(16) x(nn)
  integer npoint1, ii
  parameter(npoint1=256)
  real(16) w1(npoint1),y1(npoint1) 
  !Guass integral
  integer npoint2
  parameter(npoint2=256)
  real(16) w2(npoint2),y2(npoint2) 
  !Guass integral
  real(16) x_aft(npoint1),y_aft(npoint1) 
  !used for interpolation
  real(16) kA(npoint1),yflowA(npoint1,50)
  real(16) rho0k0(2),Zphik0
  real(16) a1,a2,a3,a4,a5,b1,b2,b3,b4,c1,c2,c3,c4,c5,d1,d2,d3,d4,d5,T0,T_r,aPolya,cPolya,dPolya,bPolya
  logical check
  real(16) pi
  parameter(pi=3.141592653589793238462643383279Q+00)
  real(16) k_infi

  COMMON /path/ kmax,kount,dxsav,xp,yp
  common /polyakov_com/ l_com,lb_com
  common /kRange/k_UV,k_IR,t_UV,t_IR
  common /y1w1/ y1, w1
  common /y2w2/ y2, w2
  common /flowconfi/ kA,yflowA
  common /rho0Zphi/ rho0k0,Zphik0
  common /Tmu/ T,mu,muq,mus

  epsi_err=1.Q-10

  l_new=l_i
  lb_new=lb_i

  j=0
  !start of loops
  stopp=.false.
  do while((.not. stopp).and.(j < jmax))
    j=j+1

    l_old=l_new
    lb_old=lb_new

    l_com=l_old
    lb_com=lb_old

    call FRG(rho0,mboson,mfermion,Vall,fpifK,h,Z_wave,ck)

    rho0k0=rho0
    Zphik0=Z_wave(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!interpolation
    if(j==1)then
      call gauleg(k_UV, k_IR, y1, w1, npoint1)
      k_infi=k_UV*10.Q+00
      call gauleg(k_UV, k_infi, y2, w2, npoint2)
    end if

    x_aft=y1
    kA=x_aft

    x_bef=k_UV*exp(xp)

    do i=1,50
    y_bef=yp(i,:)
    call intLin(KMAXX,kount,npoint1,x_bef,y_bef,x_aft,y_aft)
    yflowA(:,i)=y_aft
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    x(1)=l_com
    x(2)=lb_com
    call newt(x, nn, check, PolyakovEq)
    l_new=x(1)
    lb_new=x(2)

    epsi_l=abs(l_new-l_old)/abs(l_new)
    epsi_lb=abs(lb_new-lb_old)/abs(lb_new)

    if(epsi_l<epsi_err.and.epsi_lb<epsi_err)then
!    if(j==1)then
      stopp=.true. 
    end if
  end do

  l=l_old
  lb=lb_old

  l_com=l
  lb_com=lb

  call vInf(k_UV,Vinfi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!glue potential

  a1=-44.14Q+00
  a2=151.4Q+00
  a3=-90.0677Q+00
  a4=2.77173Q+00
  a5=3.56403Q+00

  b1=-0.32665Q+00
  b2=-82.9823Q+00
  b3=3.Q+00
  b4=5.85559Q+00

  c1=-50.7961Q+00
  c2=114.038Q+00
  c3=-89.4596Q+00
  c4=3.08718Q+00
  c5=6.72812Q+00

  d1=27.0885Q+00
  d2=-56.0859Q+00
  d3=71.2225Q+00
  d4=2.9715Q+00
  d5=6.61433Q+00

  T0=270.Q+00

  T_r=0.60Q+00*(T-T0)/T0+1.Q+00

  aPolya=(a1+a2/T_r+a3/T_r**2)/(1.Q+00+a4/T_r+a5/T_r**2)
  cPolya=(c1+c2/T_r+c3/T_r**2)/(1.Q+00+c4/T_r+c5/T_r**2)
  dPolya=(d1+d2/T_r+d3/T_r**2)/(1.Q+00+d4/T_r+d5/T_r**2)

  bPolya=b1*T_r**(-b4)*(1.Q+00-exp(b2/T_r**b3))

  Vglue=T**4*(-(aPolya*l*lb)/2.Q+00 + dPolya*l**2*lb**2 + (cPolya*(l**3 + lb**3))/2.Q+00 +  &
    bPolya*Log(1 - 6*l*lb - 3*l**2*lb**2 + 4*(l**3 + lb**3)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Vtotal=Vall+Vinfi+Vglue
end
