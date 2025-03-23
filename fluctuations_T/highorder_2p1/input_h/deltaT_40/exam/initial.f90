subroutine initial(Nflow,yflow)
!make the initialization
!wenrui checked 2017.10.16

  implicit none

  integer Nflow
  real(16) yflow(50) 
  integer N_str(4) 
  !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck
  real(16) lamd10,lamd20,lamd01
  real(16) lam00,lam10,lam20,lam30,lam40,lam50,lam60,lam70
  real(16) lam01,lam11,lam21,lam31,lam41,lam51
  real(16) lam02,lam12,lam22,lam32
  real(16) lam03,lam13
  real(16) hlk,hsk
  !Yukawa coupling constants
  real(16) Z_pi,Z_K,Z_l,Z_s
  real(16) ck,jl,js
  real(16) kappa(2)
  real(16) kappa1,kappa2
  real(16) sigmal,sigmas
  real(16) sigmal_old,sigmas_old
  real(16) sigmal_i,sigmas_i
  integer nn,i
  parameter(nn=2)
  real(16) x(nn)
  logical check,stopp
  external gapEqi
  real(16) delta,epsi
  parameter(epsi=1.Q-07)
  integer imax
  parameter(imax=10000)
  real(16) mf(3)
  REAL(16) :: MKaon,MKappa

  real*16 :: h_mod,hs_mod
  REAL*16 :: h_const,lamd10_const,lamd20_const,lamd01_const,Tcglue_const,hs_const

  real(16) mboson2_UV(18)
  real(16) dmbo2drho(18,12)
  real(16) sin2phi(2),cos2phi(2)


  common /strucFun/ N_str
  common /gapPara/ lamd10,lamd20,lamd01,ck,jl,js
  common /mfUK/ mf,mboson2_UV
  common /ini_const/ h_const,lamd10_const,lamd20_const,lamd01_const,Tcglue_const,hs_const

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)

  lamd10=lamd10_const**2
  !in unit of MeV^2
  lamd20=lamd20_const
  lamd01=lamd01_const

  ck=4808.Q+00

  !explicit chiral symmetry breaking term
  !jl=119.Q+00**3
  jl=121.Q+00**3
  !in unit of MeV^3
  js=337Q+00**3
  !in unit of MeV^3

  sigmal_i=40.Q+00
  sigmas_i=60.Q+00

  sigmal_old=sigmal_i
  sigmas_old=sigmas_i

  i=0
  !start of loops
  stopp=.false.
  do while((.not. stopp).and.(i < imax))
    i=i+1
    x(1)=sigmal_old
    x(2)=sigmas_old
    call newt(x, nn, check, gapEqi)
    sigmal=x(1)
    sigmas=x(2)
    delta=abs((sigmal-sigmal_old)/sigmal_old)+abs((sigmas-sigmas_old)/sigmas_old)
    if(abs(delta)<epsi)then
      stopp=.true.
    else
      sigmal_old=sigmal
      sigmas_old=sigmas
    end if
  end do

  kappa(1)=(sigmal**2 + sigmas**2)/2.Q+00
  !in unit of MeV^2
  kappa(2)=(sigmal**2 - 2*sigmas**2)**2/24.Q+00
  !in unit of MeV^4

  lam10=lamd10+lamd20*kappa(1)
  !in unit of MeV^2
  lam20=lamd20
  lam30=0.0Q+00
  lam40=0.0Q+00
  lam50=0.0Q+00
  lam60=0.0Q+00
  lam70=0.0Q+00
  lam01=lamd01
  lam11=0.0Q+00
  lam21=0.0Q+00
  lam31=0.0Q+00
  lam41=0.0Q+00
  lam51=0.0Q+00
  lam02=0.0Q+00
  lam12=0.0Q+00
  lam22=0.0Q+00
  lam32=0.0Q+00
  lam03=0.0Q+00
  lam13=0.0Q+00
  !expansion coefficients of effective potential V
  call fityukawa(1.Q+3,h_mod,hs_mod)
  hlk=h_mod*h_const
  hsk=hs_mod*hs_const
!  hsk=h_mod*h_const
  !expansion coefficients of Yukawa coupling
  Z_pi=1.0Q+00 
  Z_K=1.0Q+00 
  !meson wave function renormalization
  Z_l=1.0Q+00
  Z_s=1.0Q+00
  !quark wave function renormalization

  mf(1)=hlk*Sigmal/2.Q+00
  mf(2)=hlk*Sigmal/2.Q+00
  mf(3)=hsk*Sigmas/Sqrt(2.Q+00)

!  write(*,*) real(sigmal,kind=4),real(sigmas,kind=4)
!  write(*,*) 'ml_UV',real(mf(1),kind=4)
!  write(*,*) 'ms_UV',real(mf(3),kind=4)
  yflow(1)=lam10
  yflow(2)=lam20
  yflow(3)=lam30
  yflow(4)=lam40
  yflow(5)=lam50
  yflow(6)=lam60
  yflow(7)=lam70
  yflow(8)=lam01
  yflow(9)=lam11
  yflow(10)=lam21
  yflow(11)=lam31
  yflow(12)=lam41
  yflow(13)=lam51
  yflow(14)=lam02
  yflow(15)=lam12
  yflow(16)=lam22
  yflow(17)=lam32
  yflow(18)=lam03
  yflow(19)=lam13
  yflow(Nv+1)=lam00
  yflow((Nv+1)+1)=hlk
  yflow((Nv+1)+2)=hsk
  yflow((Nv+1)+(Nh+2)+1)=Z_pi
  yflow((Nv+1)+(Nh+2)+2)=Z_K
  yflow((Nv+1)+(Nh+2)+3)=Z_l
  yflow((Nv+1)+(Nh+2)+4)=Z_s
  yflow((Nv+1)+(Nh+2)+Nz+1)=ck
  yflow((Nv+1)+(Nh+2)+Nz+2)=Sigmal
  yflow((Nv+1)+(Nh+2)+Nz+3)=Sigmas
  yflow((Nv+1)+(Nh+2)+Nz+4)=jl
  yflow((Nv+1)+(Nh+2)+Nz+5)=js

  call massbose2(yflow,dmbo2drho,cos2phi,sin2phi)

  mboson2_UV=dmbo2drho(:,1)
  mkappa=sqrt(mboson2_UV(5))
  mkaon =sqrt(mboson2_UV(17))
!  WRITE(*,*) 'kappa=',mkappa
!  WRITE(*,*) 'kaon =',mkaon
!  WRITE(*,*) sqrt( mboson2)
!  STOP
end
