subroutine derivs(x,y,dydx)
!Calculating the right hand side of differential equations
!wenrui wrote 2017.09.22
!wenrui checked 2017.10.16

  implicit none

  real(16) pi
  parameter(pi=3.141592653589793238462643383279Q+00)
  !constant
  integer NMAX 
  !maximal number of differential equations
  parameter(NMAX=50)
  real(16) x,y(NMAX),dydx(NMAX)
  integer N_str(4) 
  !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck
  !number of lam h Z ck
  real(16) k 
  ! IR cutoff in flow equations
  real(16) lam00,lam10,lam20,lam30,lam40,lam50,lam60,lam70
  real(16) lam01,lam11,lam21,lam31,lam41,lam51
  real(16) lam02,lam12,lam22,lam32
  real(16) lam03,lam13
  !constant of Taylor expansion of potential V
  real(16) dlam00dt,dlam10dt,dlam20dt,dlam30dt,dlam40dt,dlam50dt,dlam60dt
  real(16) dlam70dt,dlam01dt,dlam11dt,dlam21dt,dlam31dt,dlam41dt,dlam51dt
  real(16) dlam02dt,dlam12dt,dlam22dt,dlam32dt
  real(16) dlam03dt,dlam13dt
  real(16) hlk,hsk,dhldt,dhsdt
  !Yukawa coupling constants
  real(16) Z_pi,Z_K,Z_l,Z_s
  real(16) kappa1,kappa2
  !expansion points
  real(16) eta_pi,eta_K,eta_l,eta_s
  !meson and quark anomanous dimension
  real(16) Nc,Nf
  parameter(Nc=3.Q+00,Nf=3.Q+00)
  !number of color and flavor
  real(16) rho,rho2
  !1/2(sigma_l^2+sigma_s^2)
  !1/24(sigma_l^2-28sigma_s^2)^2
  real(16) T,mu,muq,mus,muK
  !temperature and chemical potential
  real(16) muB(18),muF(3)
  real(16) k_UV,k_IR,t_UV,t_IR
  real(16) ck,dckdt
  !anomaly breaking of axial U(1)A symmetry constant (KMT term)
  real(16) mboson2(18)
  !mf02,ma02*3,mkap2*4,msig2
  !meta2,mpi2*3,mK2*4,meta12
  real(16) dmbo2drho(18,12)
  !mf02,ma02*3,mkap2*4,msig2
  !meta2,mpi2*3,mK2*4,meta12
  !mbo2d00 is mass^2 of mboson
  !00,10,20,30,40,50,01,11,21,31,02,12
  real(16) sin2phi(2),cos2phi(2)
  !sin^2 and cos^2 of mixing angles
  !scalar,psendoscalar
  real(16) mfermion2(3)
  !ml2,ms2
  real(16) mbo2d10rho(18),mbo2d20rho(18),mbo2d30rho(18),mbo2d40rho(18),mbo2d50rho(18)
  real(16) mbo2d01rho(18),mbo2d11rho(18),mbo2d21rho(18),mbo2d31rho(18)
  real(16) mbo2d02rho(18),mbo2d12rho(18)  
  real(16) mfe2d10rho(3),mfe2d01rho(3),mfe2d02rho(3)
  !derivatives of mass
  real(16) dnbopldm0(18),dnbopldm1(18),dnbopldm2(18),dnbopldm3(18),dnbopldm4(18),dnbopldm5(18)
  real(16) dnbomidm0(18),dnbomidm1(18),dnbomidm2(18),dnbomidm3(18),dnbomidm4(18),dnbomidm5(18)
  real(16) nfminumu(3),nfplusmu(3)
  real(16) dnfmimudm1(3),dnfmimudm2(3),dnfmimudm3(3),dnfmimudm4(3),dnfmimudm5(3)
  real(16) dnfplmudm1(3),dnfplmudm2(3),dnfplmudm3(3),dnfplmudm4(3),dnfplmudm5(3)
  real(16) dnfplmudm(3,6),dnfmimudm(3,6)
  real(16) l_com,lb_com
  real(16) l,lb 
  !polyakov loop
  real(16) p0
  !p0,ext 
  complex(16) FB21lcomplex(18),FB12lcomplex(18),FB21scomplex(18),FB12scomplex(18)
  real(16) FB21l(18),FB12l(18),FB21s(18),FB12s(18)
  !terms of function B10 
  real(16) capL11l(18),capL11s(18)
  !running Yukawa coupling
  !f0,a0,kappa,sigma,eta,pi,K,eta1
  real(16) lb0(18),lbm(18),lbm2(18),lbm3(18),lbm4(18),lbm5(18)
  real(16) lf0(3),lfm(3),lfm2(3),lfm3(3),lfm4(3),lfm5(3)
  !bosonic/fermionic loop
  real(16) lbt10(18),lbt20(18),lbt30(18),lbt40(18),lbt50(18)
  real(16) lbt01(18),lbt11(18),lbt21(18),lbt31(18)
  real(16) lbt02(18),lbt12(18) 
  real(16) lft10(3),lft20(3),lft30(3),lft40(3),lft50(3)
  real(16) lft01(3),lft11(3),lft21(3),lft31(3)
  real(16) lft02(3),lft12(3)
  !terms of drdtV ,terms of function 26
  real(16) k44pi2
  !calculate accelerate,no sence
!  real(16) dvbconst(8)
!  parameter(dvbconst=(/1.Q+00,3.Q+00,4.Q+00,1.Q+00,1.Q+00,3.Q+00,4.Q+00,1.Q+00/))
!  real(16) dvfconst(2)
!  parameter(dvfconst=(/2.Q+00,1.Q+00/))
  !function 26, constants of flow equation
  real(16) dhlconst(18),dhsconst(18)
  real(16) dr00dtv,dr10dtV,dr20dtV,dr30dtV,dr40dtV,dr50dtV
  real(16) dr01dtV,dr11dtV,dr21dtV,dr31dtV
  real(16) dr02dtV,dr12dtV
  real(16) dkappa1dt,dkappa2dt
  real(16) sl,ss,dsldt,dssdt
  real(16) jl,js,djldt,djsdt

  real*16 :: h_mod,hs_mod
  REAL*16 :: h_const,lamd10_const,lamd20_const,lamd01_const,Tcglue_const,hs_const
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  common /strucFun/ N_str
  common /Tmu/ T,mu,muq,mus,muK
  common /kRange/k_UV,k_IR,t_UV,t_IR
  common /polyakov_com/ l_com,lb_com
  common /ini_const/ h_const,lamd10_const,lamd20_const,lamd01_const,Tcglue_const,hs_const

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)

  k=k_UV*exp(x)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  lam10=y(1)
  lam20=y(2)
  lam30=y(3)
  lam40=y(4)
  lam50=y(5)
  lam60=y(6)
  lam70=y(7)
  lam01=y(8)
  lam11=y(9)
  lam21=y(10)
  lam31=y(11)
  lam41=y(12)
  lam51=y(13)
  lam02=y(14)
  lam12=y(15)
  lam22=y(16)
  lam32=y(17)
  lam03=y(18)
  lam13=y(19)
  lam00=y(Nv+1)
  hlk=y((Nv+1)+1)
  hsk=y((Nv+1)+2)
  Z_pi=y((Nv+1)+(Nh+2)+1)
  Z_K=y((Nv+1)+(Nh+2)+2)
  Z_l=y((Nv+1)+(Nh+2)+3)
  Z_s=y((Nv+1)+(Nh+2)+4)
  ck=y((Nv+1)+(Nh+2)+Nz+1)
  Sl=y((Nv+1)+(Nh+2)+Nz+2)
  Ss=y((Nv+1)+(Nh+2)+Nz+3)
  jl=y((Nv+1)+(Nh+2)+Nz+4)
  js=y((Nv+1)+(Nh+2)+Nz+5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  kappa1=(Sl**2 + Ss**2)/2.Q+00
  kappa2=(Sl**2 - 2*Ss**2)**2/24.Q+00
  rho=kappa1
  rho2=kappa2
  !calculations are performed at expansion point kappa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call fityukawa(k,h_mod,hs_mod)

!  write(*,*)h_mod
  hlk=h_mod*h_const
  hsk=hs_mod*hs_const
!  hsk=h_mod*h_const
!  stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !square of boson mass and their derivatives

  call massbose2(y,dmbo2drho,cos2phi,sin2phi)

  mboson2  =dmbo2drho(:,1)
  mbo2d10rho=dmbo2drho(:,2)
  mbo2d20rho=dmbo2drho(:,3)
  mbo2d30rho=dmbo2drho(:,4)
  mbo2d40rho=dmbo2drho(:,5)
  mbo2d50rho=dmbo2drho(:,6)
  mbo2d01rho=dmbo2drho(:,7)
  mbo2d11rho=dmbo2drho(:,8)
  mbo2d21rho=dmbo2drho(:,9)
  mbo2d31rho=dmbo2drho(:,10)
  mbo2d02rho=dmbo2drho(:,11)
  mbo2d12rho=dmbo2drho(:,12)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !scalar,psendoscalar
  !f0,a0,kap,sig,eta,pi,K,eta1
  dhlconst=(/cos2phi(1),1.Q+00,1.Q+00,1.Q+00,0.Q+00,0.Q+00,0.Q+00,0.Q+00,&
          sin2phi(1),-cos2phi(2),-1.Q+00,-1.Q+00,-1.Q+00,0.Q+00,0.Q+00,0.Q+00,0.Q+00,-sin2phi(2)/)
  dhsconst=(/sin2phi(1),0.Q+00,0.Q+00,0.Q+00,0.Q+00,0.Q+00,0.Q+00,0.Q+00,&
          cos2phi(1),-sin2phi(2),0.Q+00,0.Q+00,0.Q+00,0.Q+00,0.Q+00,0.Q+00,0.Q+00,-cos2phi(2)/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !square of fermion mass and their derivatives
  !ml2,ms2
  mfermion2(1)=(hlk**2*Sl**2)/4.Q+00
  mfe2d10rho(1)=hlk**2/3.Q+00
  mfe2d01rho(1)=hlk**2/(Sl**2 - 2*Ss**2)
  mfe2d02rho(1)=(-12.Q+00*hlk**2)/(Sl**2 - 2.Q+00*Ss**2)**3
  mfermion2(2)=(hlk**2*Sl**2)/4.Q+00
  mfe2d10rho(2)=hlk**2/3.Q+00
  mfe2d01rho(2)=hlk**2/(Sl**2 - 2.Q+00*Ss**2)
  mfe2d02rho(2)=(-12.Q+00*hlk**2)/(Sl**2 - 2.Q+00*Ss**2)**3
  mfermion2(3)=(hsk**2*Ss**2)/2.Q+00
  mfe2d10rho(3)=hsk**2/3.Q+00
  mfe2d01rho(3)=(-2.Q+00*hsk**2)/(Sl**2 - 2.Q+00*Ss**2)
  mfe2d02rho(3)=(24.Q+00*hsk**2)/(Sl**2 - 2.Q+00*Ss**2)**3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  l=l_com
  lb=lb_com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  eta_pi=0.Q+00
  eta_K=0.Q+00
  eta_l=0.Q+00
  eta_s=0.Q+00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calculate nboson and it's mass derivation 
  muB=(/0.Q+00,0.Q+00,muq,muq,mus,mus,muq+mus+muK,muq+mus+muK,0.Q+00,           &
        0.Q+00,0.Q+00,muq,muq,mus,mus,muq+mus+muK,muq+mus+muK,0.Q+00 /)
  call nbose(k,muB,T,mboson2,dnbopldm0,dnbopldm1,dnbopldm2,dnbopldm3,dnbopldm4,dnbopldm5)
  call nbose(k,-muB,T,mboson2,dnbomidm0,dnbomidm1,dnbomidm2,dnbomidm3,dnbomidm4,dnbomidm5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calculate nf(E-mu) and it's mass derivation 
!calculate nf(E+mu) and it's mass derivation 
!product by nfpl_nfmi.nb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  muF=(/mu+2.Q+00/3.Q+00*muq,mu-1.Q+00/3.Q+00*muq,mu-1.Q+00/3.Q+00*muq-mus/)
  call nfplmimu(k,muF,T,l,lb,mfermion2,dnfplmudm,dnfmimudm)
  nfminumu=dnfmimudm(:,1)
  dnfmimudm1=dnfmimudm(:,2)
  dnfmimudm2=dnfmimudm(:,3)
  dnfmimudm3=dnfmimudm(:,4)
  dnfmimudm4=dnfmimudm(:,5)
  dnfmimudm5=dnfmimudm(:,6)
  nfplusmu=dnfplmudm(:,1)
  dnfplmudm1=dnfplmudm(:,2)
  dnfplmudm2=dnfplmudm(:,3)
  dnfplmudm3=dnfplmudm(:,4)
  dnfplmudm4=dnfplmudm(:,5)
  dnfplmudm5=dnfplmudm(:,6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  p0=pi*T
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !running Yukawa coupling
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!the l_boson^0 and its mass derivation
!producted by lbm_product.nb
  lb0=-((1 + dnbomidm0 + dnbopldm0)*(-5 + eta_pi)*k)/                            &
    (15.Q+00*Sqrt(k**2 + mboson2))
  lbm=((1 - eta_pi/5.Q+00)*k*(-1 - dnbomidm0 - dnbopldm0 +                           &
        2*(dnbomidm1 + dnbopldm1)*(k**2 + mboson2)))/                           &
    (6.Q+00*(k**2 + mboson2)**1.5)
  lbm2=((1 - eta_pi/5.Q+00)*k*(-3*                                                   &
         (-1 - dnbomidm0 - dnbopldm0 +                                          &
           2*(dnbomidm1 + dnbopldm1)*(k**2 + mboson2)) +                        &
        2*(k**2 + mboson2)*(dnbomidm1 + dnbopldm1 +                             &
           2*(dnbomidm2 + dnbopldm2)*(k**2 + mboson2))))/                       &
    (12.Q+00*(k**2 + mboson2)**2.5)
  lbm3=-((-5 + eta_pi)*k*(-15 - 15*dnbomidm0 - 15*dnbopldm0 +                    &
         18*dnbomidm1*k**2 + 18*dnbopldm1*k**2 - 12*dnbomidm2*k**4 -            &
         12*dnbopldm2*k**4 + 8*dnbomidm3*k**6 + 8*dnbopldm3*k**6 +              &
         18*dnbomidm1*mboson2 + 18*dnbopldm1*mboson2 -                          &
         24*dnbomidm2*k**2*mboson2 - 24*dnbopldm2*k**2*mboson2 +                &
         24*dnbomidm3*k**4*mboson2 + 24*dnbopldm3*k**4*mboson2 -                &
         12*dnbomidm2*mboson2**2 - 12*dnbopldm2*mboson2**2 +                    &
         24*dnbomidm3*k**2*mboson2**2 + 24*dnbopldm3*k**2*mboson2**2 +          &
         8*dnbomidm3*mboson2**3 + 8*dnbopldm3*mboson2**3))/                     &
    (120.Q+00*(k**2 + mboson2)**3.5)
  lbm4=-((-5 + eta_pi)*k*(105 + 105*dnbomidm0 + 105*dnbopldm0 -                  &
         120*dnbomidm1*k**2 - 120*dnbopldm1*k**2 + 72*dnbomidm2*k**4 +          &
         72*dnbopldm2*k**4 - 32*dnbomidm3*k**6 - 32*dnbopldm3*k**6 +            &
         16*dnbomidm4*k**8 + 16*dnbopldm4*k**8 - 120*dnbomidm1*mboson2 -        &
         120*dnbopldm1*mboson2 + 144*dnbomidm2*k**2*mboson2 +                   &
         144*dnbopldm2*k**2*mboson2 - 96*dnbomidm3*k**4*mboson2 -               &
         96*dnbopldm3*k**4*mboson2 + 64*dnbomidm4*k**6*mboson2 +                &
         64*dnbopldm4*k**6*mboson2 + 72*dnbomidm2*mboson2**2 +                  &
         72*dnbopldm2*mboson2**2 - 96*dnbomidm3*k**2*mboson2**2 -               &
         96*dnbopldm3*k**2*mboson2**2 + 96*dnbomidm4*k**4*mboson2**2 +          &
         96*dnbopldm4*k**4*mboson2**2 - 32*dnbomidm3*mboson2**3 -               &
         32*dnbopldm3*mboson2**3 + 64*dnbomidm4*k**2*mboson2**3 +               &
         64*dnbopldm4*k**2*mboson2**3 + 16*dnbomidm4*mboson2**4 +               &
         16*dnbopldm4*mboson2**4))/(240.Q+00*(k**2 + mboson2)**4.5)
  lbm5=-((-5 + eta_pi)*k*(-945 - 945*dnbomidm0 - 945*dnbopldm0 +                 &
         1050*dnbomidm1*k**2 + 1050*dnbopldm1*k**2 - 600*dnbomidm2*k**4 -       &
         600*dnbopldm2*k**4 + 240*dnbomidm3*k**6 + 240*dnbopldm3*k**6 -         &
         80*dnbomidm4*k**8 - 80*dnbopldm4*k**8 + 32*dnbomidm5*k**10 +           &
         32*dnbopldm5*k**10 + 1050*dnbomidm1*mboson2 +                          &
         1050*dnbopldm1*mboson2 - 1200*dnbomidm2*k**2*mboson2 -                 &
         1200*dnbopldm2*k**2*mboson2 + 720*dnbomidm3*k**4*mboson2 +             &
         720*dnbopldm3*k**4*mboson2 - 320*dnbomidm4*k**6*mboson2 -              &
         320*dnbopldm4*k**6*mboson2 + 160*dnbomidm5*k**8*mboson2 +              &
         160*dnbopldm5*k**8*mboson2 - 600*dnbomidm2*mboson2**2 -                &
         600*dnbopldm2*mboson2**2 + 720*dnbomidm3*k**2*mboson2**2 +             &
         720*dnbopldm3*k**2*mboson2**2 - 480*dnbomidm4*k**4*mboson2**2 -        &
         480*dnbopldm4*k**4*mboson2**2 + 320*dnbomidm5*k**6*mboson2**2 +        &
         320*dnbopldm5*k**6*mboson2**2 + 240*dnbomidm3*mboson2**3 +             &
         240*dnbopldm3*mboson2**3 - 320*dnbomidm4*k**2*mboson2**3 -             &
         320*dnbopldm4*k**2*mboson2**3 + 320*dnbomidm5*k**4*mboson2**3 +        &
         320*dnbopldm5*k**4*mboson2**3 - 80*dnbomidm4*mboson2**4 -              &
         80*dnbopldm4*mboson2**4 + 160*dnbomidm5*k**2*mboson2**4 +              &
         160*dnbopldm5*k**2*mboson2**4 + 32*dnbomidm5*mboson2**5 +              &
         32*dnbopldm5*mboson2**5))/(480.Q+00*(k**2 + mboson2)**5.5)  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!the l_femi^0 and its mass derivation
!producted by lfm_product1.nb

  lf0=((1 - eta_l/4.Q+00)*k*(1 - nfminumu - nfplusmu))/                        &
    (3.Q+00*Sqrt(k**2 + mfermion2))

  lfm=((-dnfmimudm1 - dnfplmudm1)*(1 - eta_l/4.Q+00)*k)/                       &
     (3.Q+00*Sqrt(k**2 + mfermion2)) -                                          &
    ((1 - eta_l/4.Q+00)*k*(1 - nfminumu - nfplusmu))/                          &
     (6.Q+00*(k**2 + mfermion2)**1.5Q+00)

  lfm2=-((-dnfmimudm1 - dnfplmudm1)*(1 - eta_l/4.Q+00)*k)/                     &
     (3.Q+00*(k**2 + mfermion2)**1.5Q+00) +                                     &
    ((-dnfmimudm2 - dnfplmudm2)*(1 - eta_l/4.Q+00)*k)/                         &
     (3.Q+00*Sqrt(k**2 + mfermion2)) +                                          &
    ((1 - eta_l/4.Q+00)*k*(1 - nfminumu - nfplusmu))/                          &
     (4.Q+00*(k**2 + mfermion2)**2.5Q+00)

  lfm3=(3*(-dnfmimudm1 - dnfplmudm1)*(1 - eta_l/4.Q+00)*k)/                    &
     (4.Q+00*(k**2 + mfermion2)**2.5Q+00) -                                     &
    ((-dnfmimudm2 - dnfplmudm2)*(1 - eta_l/4.Q+00)*k)/                         &
     (2.Q+00*(k**2 + mfermion2)**1.5Q+00) +                                     &
    ((-dnfmimudm3 - dnfplmudm3)*(1 - eta_l/4.Q+00)*k)/                         &
     (3.Q+00*Sqrt(k**2 + mfermion2)) -                                          &
    (5*(1 - eta_l/4.Q+00)*k*(1 - nfminumu - nfplusmu))/                        &
     (8.Q+00*(k**2 + mfermion2)**3.5Q+00)

  lfm4=(-5*(-dnfmimudm1 - dnfplmudm1)*(1 - eta_l/4.Q+00)*k)/                   &
     (2.Q+00*(k**2 + mfermion2)**3.5Q+00) +                                     &
    (3*(-dnfmimudm2 - dnfplmudm2)*(1 - eta_l/4.Q+00)*k)/                       &
     (2.Q+00*(k**2 + mfermion2)**2.5Q+00) -                                     &
    (2*(-dnfmimudm3 - dnfplmudm3)*(1 - eta_l/4.Q+00)*k)/                        &
     (3.Q+00*(k**2 + mfermion2)**1.5Q+00) +                                     &
    ((-dnfmimudm4 - dnfplmudm4)*(1 - eta_l/4.Q+00)*k)/                         &
     (3.Q+00*Sqrt(k**2 + mfermion2)) +                                          &
    (35*(1 - eta_l/4.Q+00)*k*(1 - nfminumu - nfplusmu))/                       &
     (16.Q+00*(k**2 + mfermion2)**4.5Q+00)

  lfm5=(175*(-dnfmimudm1 - dnfplmudm1)*(1 - eta_l/4.Q+00)*k)/                  &
     (16.Q+00*(k**2 + mfermion2)**4.5Q+00) -                                    &
    (25*(-dnfmimudm2 - dnfplmudm2)*(1 - eta_l/4.Q+00)*k)/                       &
     (4.Q+00*(k**2 + mfermion2)**3.5Q+00) +                                     &
    (5*(-dnfmimudm3 - dnfplmudm3)*(1 - eta_l/4.Q+00)*k)/                       &
     (2.Q+00*(k**2 + mfermion2)**2.5Q+00) -                                     &
    (5*(-dnfmimudm4 - dnfplmudm4)*(1 - eta_l/4.Q+00)*k)/                       &
     (6.Q+00*(k**2 + mfermion2)**1.5Q+00) +                                     &
    ((-dnfmimudm5 - dnfplmudm5)*(1 - eta_l/4.Q+00)*k)/                         &
     (3.Q+00*Sqrt(k**2 + mfermion2)) -                                          &
    (315*(1 - eta_l/4.Q+00)*k*(1 - nfminumu - nfplusmu))/                      &
     (32.Q+00*(k**2 + mfermion2)**5.5Q+00)                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!the terms of flow function
!producted by lt_product.nb

  lbt10=lbm*mbo2d10rho

  lbt20=lbm2*mbo2d10rho**2 + lbm*mbo2d20rho

  lbt30=lbm3*mbo2d10rho**3 + 3*lbm2*mbo2d10rho*mbo2d20rho +                     &
    lbm*mbo2d30rho

  lbt40=lbm4*mbo2d10rho**4 + 6*lbm3*mbo2d10rho**2*mbo2d20rho +                  &
    3*lbm2*mbo2d20rho**2 + 4*lbm2*mbo2d10rho*mbo2d30rho + lbm*mbo2d40rho

  lbt50=lbm5*mbo2d10rho**5 + 10*lbm4*mbo2d10rho**3*mbo2d20rho +                 &
    15*lbm3*mbo2d10rho*mbo2d20rho**2 + 10*lbm3*mbo2d10rho**2*mbo2d30rho +       &
    10*lbm2*mbo2d20rho*mbo2d30rho + 5*lbm2*mbo2d10rho*mbo2d40rho +              &
    lbm*mbo2d50rho

  lbt01=lbm*mbo2d01rho

  lbt11=lbm2*mbo2d01rho*mbo2d10rho + lbm*mbo2d11rho

  lbt21=lbm3*mbo2d01rho*mbo2d10rho**2 + 2*lbm2*mbo2d10rho*mbo2d11rho +          &
    lbm2*mbo2d01rho*mbo2d20rho + lbm*mbo2d21rho

  lbt31=lbm4*mbo2d01rho*mbo2d10rho**3 +                                         &
    3*lbm3*mbo2d10rho**2*mbo2d11rho +                                           &
    3*lbm3*mbo2d01rho*mbo2d10rho*mbo2d20rho +                                   &
    3*lbm2*mbo2d11rho*mbo2d20rho + 3*lbm2*mbo2d10rho*mbo2d21rho +               &
    lbm2*mbo2d01rho*mbo2d30rho + lbm*mbo2d31rho

  lbt02=lbm2*mbo2d01rho**2 + lbm*mbo2d02rho

  lbt12=lbm3*mbo2d01rho**2*mbo2d10rho + lbm2*mbo2d02rho*mbo2d10rho +            &
    2*lbm2*mbo2d01rho*mbo2d11rho + lbm*mbo2d12rho 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  lft10=lfm*mfe2d10rho
  lft20=lfm2*mfe2d10rho**2 
  lft30=lfm3*mfe2d10rho**3 
  lft40=lfm4*mfe2d10rho**4 
  lft50=lfm5*mfe2d10rho**5
  lft01=lfm*mfe2d01rho
  lft11=lfm2*mfe2d01rho*mfe2d10rho 
  lft21=lfm3*mfe2d01rho*mfe2d10rho**2
  lft31=lfm4*mfe2d01rho*mfe2d10rho**3
  lft02=lfm2*mfe2d01rho**2 + lfm*mfe2d02rho
  lft12=lfm3*mfe2d01rho**2*mfe2d10rho + lfm2*mfe2d02rho*mfe2d10rho 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!the rho deviration of the dtV
!equation 26
  k44pi2=k**4/4.Q+00/pi/pi

  dr00dtV=k44pi2*(Sum(lb0)-4*Nc*(Sum(lf0)))
  dr10dtV=k44pi2*(Sum(lbt10)-4*Nc*(Sum(lft10)))
  dr20dtV=k44pi2*(Sum(lbt20)-4*Nc*(Sum(lft20)))
  dr30dtV=k44pi2*(Sum(lbt30)-4*Nc*(Sum(lft30)))
  dr40dtV=k44pi2*(Sum(lbt40)-4*Nc*(Sum(lft40)))
  dr50dtV=k44pi2*(Sum(lbt50)-4*Nc*(Sum(lft50)))
  dr01dtV=k44pi2*(Sum(lbt01)-4*Nc*(Sum(lft01)))
  dr11dtV=k44pi2*(Sum(lbt11)-4*Nc*(Sum(lft11)))
  dr21dtV=k44pi2*(Sum(lbt21)-4*Nc*(Sum(lft21)))
  dr31dtV=k44pi2*(Sum(lbt31)-4*Nc*(Sum(lft31)))
  dr02dtV=k44pi2*(Sum(lbt02)-4*Nc*(Sum(lft02)))
  dr12dtV=k44pi2*(Sum(lbt12)-4*Nc*(Sum(lft12)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dSldt=-(((dr01dtV + (eta_pi*js*Sl - eta_pi*jl*Ss)/                             &
            (Sl**3*Ss - 2*Sl*Ss**3))*                                           &
         (lam20*Ss - (lam11*Ss*(Sl**2 - 2*Ss**2))/3.Q+00 +                          &
           (4*js + Sqrt(2.Q+00)*ck*(Sl**2 - 4*Ss**2))/(12.Q+00*Ss**2)) -            &
        (dr10dtV - (eta_pi*(js*Sl + 2*jl*Ss))/(6.Q+00*Sl*Ss))*                      &
         (lam11*Ss - (lam02*Ss*(Sl**2 - 2*Ss**2))/3.Q+00 -                          &
           (16*jl*Ss**3 + Sqrt(2.Q+00)*ck*Sl*(Sl**2 - 2*Ss**2)**2 +             &
              4*js*(Sl**3 - 6*Sl*Ss**2))/(2.Q+00*Sl*Ss**2*(Sl**2 - 2*Ss**2)**2)     &
))/((lam11*Sl + (lam02*Sl*(Sl**2 - 2*Ss**2))/6.Q+00 +                               &
           (-4*js*Sl**3 + 6*jl*Sl**2*Ss - 4*jl*Ss**3)/                          &
            (Ss*(Sl**3 - 2*Sl*Ss**2)**2))*                                      &
         (lam20*Ss - (lam11*Ss*(Sl**2 - 2*Ss**2))/3.Q+00 +                          &
           (4*js + Sqrt(2.Q+00)*ck*(Sl**2 - 4*Ss**2))/(12.Q+00*Ss**2)) -            &
        ((2*jl)/(3.Q+00*Sl**2) + (Sl*                                               &
              (-(Sqrt(2.Q+00)*ck) + 6*lam20*Ss + lam11*Sl**2*Ss -               &
                2*lam11*Ss**3))/(6.Q+00*Ss))*                                       &
         (lam11*Ss - (lam02*Ss*(Sl**2 - 2*Ss**2))/3.Q+00 -                          &
           (16*jl*Ss**3 + Sqrt(2.Q+00)*ck*Sl*(Sl**2 - 2*Ss**2)**2 +             &
              4*js*(Sl**3 - 6*Sl*Ss**2))/(2.Q+00*Sl*Ss**2*(Sl**2 - 2*Ss**2)**2))))
  dSsdt=(-2*Sl*Ss*(-6*Ss*(dr10dtV*                                              &
            (24*js*Sl**3 + Ss*                                                  &
               (jl*(-36*Sl**2 + 24*Ss**2) -                                     &
                 Sl**3*(Sl**2 - 2*Ss**2)**2*                                    &
                  (6*lam11 + lam02*Sl**2 - 2*lam02*Ss**2))) +                   &
           dr01dtV*(Sl**2 - 2*Ss**2)**2*                                        &
            (-(Sqrt(2.Q+00)*ck*Sl**3) +                                         &
              Ss*(4*jl + Sl**3*(6*lam20 + lam11*Sl**2 - 2*lam11*Ss**2)))) +     &
        eta_pi*(24*js**2*Sl**3 -                                                &
           2*jl*Sl*Ss*(3*Sqrt(2.Q+00)*ck*Sl*(Sl**2 - 2*Ss**2) +                 &
              Ss*(24*jl + Sl*(Sl**2 - 2*Ss**2)*                                 &
                  (-18*lam20 +                                                  &
                    (Sl**2 - 2*Ss**2)*                                          &
                     (3*lam11 + lam02*Sl**2 - 2*lam02*Ss**2)))) +               &
           js*(6*Sqrt(2.Q+00)*ck*Sl**3*(Sl**2 - 2*Ss**2) +                      &
              Ss*(-12*jl*(Sl**2 - 6*Ss**2) -                                    &
                 Sl**3*(Sl**2 - 2*Ss**2)*                                       &
                  (36*lam20 + (Sl**2 - 2*Ss**2)*                                &
                     (12*lam11 + lam02*Sl**2 - 2*lam02*Ss**2)))))))/            &
    (-96*js**2*Sl**4 - 12*ck**2*Sl**4*(Sl**2 - 2*Ss**2)**2 +                    &
      Sqrt(2.Q+00)*ck*Sl*Ss*(12*jl*(5*Sl**4 - 30*Sl**2*Ss**2 + 16*Ss**4) +      &
         Sl**3*(Sl**2 - 2*Ss**2)**2*                                            &
          (36*lam20 + (Sl**2 - 2*Ss**2)*                                        &
             (12*lam11 + lam02*Sl**2 - 8*lam02*Ss**2))) +                       &
      4*Ss**4*(96*jl**2 - 9*(lam11**2 - lam02*lam20)*Sl**4*                     &
          (Sl**2 - 2*Ss**2)**3 +                                                &
         4*jl*Sl*(9*lam20*(7*Sl**2 - 2*Ss**2) +                                 &
            (Sl**2 - 2*Ss**2)**2*(-6*lam11 + lam02*(Sl**2 - 2*Ss**2)))) -       &
      4*js*Sl*(12*Sqrt(2.Q+00)*ck*Sl**3*(Sl**2 - 5*Ss**2) +                     &
         Ss*(jl*(-60*Sl**2 + 168*Ss**2) -                                       &
            Sl**3*(36*lam20*(Sl**2 - 8*Ss**2) +                                 &
               (Sl**2 - 2*Ss**2)**2*(12*lam11 + lam02*Sl**2 - 2*lam02*Ss**2)))))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dkappa1dt=dSldt*Sl+dSsdt*Ss
  dkappa2dt=((dSldt*Sl - 2*dSsdt*Ss)*(Sl**2 - 2*Ss**2))/6.Q+00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dlam00dt=dr00dtV+(dkappa1dt+eta_pi*kappa1)*lam10+(dkappa2dt+2*eta_pi*kappa2)*lam01+0*eta_pi*lam00
  dlam10dt=dr10dtV+(dkappa1dt+eta_pi*kappa1)*lam20+(dkappa2dt+2*eta_pi*kappa2)*lam11+1*eta_pi*lam10
  dlam20dt=dr20dtV+(dkappa1dt+eta_pi*kappa1)*lam30+(dkappa2dt+2*eta_pi*kappa2)*lam21+2*eta_pi*lam20
  dlam30dt=dr30dtV+(dkappa1dt+eta_pi*kappa1)*lam40+(dkappa2dt+2*eta_pi*kappa2)*lam31+3*eta_pi*lam30
  dlam40dt=dr40dtV+(dkappa1dt+eta_pi*kappa1)*lam50+(dkappa2dt+2*eta_pi*kappa2)*lam41+4*eta_pi*lam40
  dlam50dt=dr50dtV+(dkappa1dt+eta_pi*kappa1)*lam60+(dkappa2dt+2*eta_pi*kappa2)*lam51+5*eta_pi*lam50
  dlam01dt=dr01dtV+(dkappa1dt+eta_pi*kappa1)*lam11+(dkappa2dt+2*eta_pi*kappa2)*lam02+2*eta_pi*lam01
  dlam11dt=dr11dtV+(dkappa1dt+eta_pi*kappa1)*lam21+(dkappa2dt+2*eta_pi*kappa2)*lam12+3*eta_pi*lam11
  dlam21dt=dr21dtV+(dkappa1dt+eta_pi*kappa1)*lam31+(dkappa2dt+2*eta_pi*kappa2)*lam22+4*eta_pi*lam21
  dlam31dt=dr31dtV+(dkappa1dt+eta_pi*kappa1)*lam41+(dkappa2dt+2*eta_pi*kappa2)*lam32+5*eta_pi*lam31
  dlam02dt=dr02dtV+(dkappa1dt+eta_pi*kappa1)*lam12+(dkappa2dt+2*eta_pi*kappa2)*lam03+4*eta_pi*lam02
  dlam12dt=dr12dtV+(dkappa1dt+eta_pi*kappa1)*lam22+(dkappa2dt+2*eta_pi*kappa2)*lam13+5*eta_pi*lam12
  dlam60dt=0.Q+00
  dlam70dt=0.Q+00
  dlam41dt=0.Q+00
  dlam51dt=0.Q+00
  dlam22dt=0.Q+00
  dlam32dt=0.Q+00
  dlam03dt=0.Q+00
  dlam13dt=0.Q+00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dhldt=0.Q+0
  dhsdt=0.Q+0
  dckdt=0.Q+00
  djldt=jl*eta_pi/2.Q+00
  djsdt=js*eta_pi/2.Q+00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dydx(1)=dlam10dt
  dydx(2)=dlam20dt
  dydx(3)=dlam30dt
  dydx(4)=dlam40dt
  dydx(5)=dlam50dt
  dydx(6)=dlam60dt
  dydx(7)=dlam70dt
  dydx(8)=dlam01dt
  dydx(9)=dlam11dt
  dydx(10)=dlam21dt
  dydx(11)=dlam31dt
  dydx(12)=dlam41dt
  dydx(13)=dlam51dt
  dydx(14)=dlam02dt
  dydx(15)=dlam12dt
  dydx(16)=dlam22dt
  dydx(17)=dlam32dt
  dydx(18)=dlam03dt
  dydx(19)=dlam13dt
  dydx(Nv+1)=dlam00dt
!  dydx((Nv+1)+1)=dhldt
  dydx((Nv+1)+1)=0.Q+00
!  dydx((Nv+1)+2)=dhsdt
  dydx((Nv+1)+2)=0.Q+00
  dydx((Nv+1)+(Nh+2)+1)=0.Q+00
  dydx((Nv+1)+(Nh+2)+2)=0.Q+00
  dydx((Nv+1)+(Nh+2)+3)=0.Q+00
  dydx((Nv+1)+(Nh+2)+4)=0.Q+00
  dydx((Nv+1)+(Nh+2)+Nz+1)=dckdt
  dydx((Nv+1)+(Nh+2)+Nz+2)=dSldt
  dydx((Nv+1)+(Nh+2)+Nz+3)=dSsdt
  dydx((Nv+1)+(Nh+2)+Nz+4)=djldt
  dydx((Nv+1)+(Nh+2)+Nz+5)=djsdt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  write(*,*) k
end
