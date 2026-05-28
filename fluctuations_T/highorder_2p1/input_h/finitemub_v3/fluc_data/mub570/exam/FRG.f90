subroutine FRG(rho0,mboson,mfermion,Vall,fpifk,h,Z_wave,ck) 
!This program solve FRG flow equations with fixed expansion point.
!wenrui wrote 2017.10.17
!wenrui changed 2017.11.10
!wenrui changed 2018.06.03

  implicit none

  integer Nv,Nh
  !Nv:order of Tylor expansion for effective potential V
  !Nh: order of Yukawa coupling h
  parameter(Nv=19)
  parameter(Nh=0)
  integer Nz 
  !number of wave function renormalizations
  parameter(Nz=4)
  integer Nck
  !number of others:ck,sigmal,sigmas,jl,js
  parameter(Nck=5)
  integer Nflow
  !number of flow equations
  parameter(Nflow=(Nv+1)+(Nh+2)+Nz+Nck)
  real(16) yflow(Nflow)
  !dependent variables in flow equations
  integer N_str(4) 
  !store the structure of functions of ODE
  real(16) T,mu,muq,mus
  real(16) k_UV,k_IR,t_UV,t_IR
  external derivs,rkqs
  real(16) eps_ode,h1,hmin 
  !variables in subroutine odeint
  integer nok,nbad
  !variables in subroutine odeint
  INTEGER kmax,kount
  !variables in common block of subroutine odeint
  INTEGER KMAXX,NMAX
  PARAMETER (NMAX=50,KMAXX=200)
  real(16) dxsav,xp(KMAXX),yp(NMAX,KMAXX) 
  !variables in common block of subroutine odeint
  real(16) rho0(2)
  real(16) fpi,fK
  real(16) hlk,hsk
  real(16) Z_wave(4)
  real(16) ck,jl,js
  real(16) l_com,lb_com
  real(16) h(2)
  real(16) mboson(8),mfermion(2),Vall,fpifk(2)

  common /strucFun/ N_str
  common /Tmu/ T,mu,muq,mus
  common /kRange/k_UV,k_IR,t_UV,t_IR
  common /odeContr/ eps_ode,h1,hmin
  COMMON /path/ kmax,kount,dxsav,xp,yp
  common /polyakov_com/ l_com,lb_com

  N_str(1)=Nv
  N_str(2)=Nh
  N_str(3)=Nz
  N_str(4)=Nck

  k_UV=1000.Q+00
  !in unit of MeV
  k_IR=0.01Q+00
  !in unit of MeV
  t_UV=0.Q+00
  t_IR=log(k_IR/k_UV)

  eps_ode=1.Q-07
  h1=t_IR/200.Q+00
  hmin=0.Q+00
  kmax=KMAXX
  dxsav=t_IR/10000.Q+00

  call initial(Nflow,yflow)
  call odeint(yflow,Nflow,t_UV,t_IR,eps_ode,h1,hmin,nok,nbad,derivs,rkqs)
  call phypoint2(Nflow,yflow,rho0,mboson,mfermion,Vall,fpifk)
  Z_wave(1)=yflow((Nv+1)+(Nh+2)+1)
  Z_wave(2)=yflow((Nv+1)+(Nh+2)+2)
  Z_wave(3)=yflow((Nv+1)+(Nh+2)+3)
  Z_wave(4)=yflow((Nv+1)+(Nh+2)+4)
  ck=yflow((Nv+1)+(Nh+2)+Nz+1)
end 
