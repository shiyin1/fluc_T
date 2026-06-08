program QM

  implicit none

  real(16) pi
  parameter(pi=3.141592653589793238462643383279Q+00)
  integer i,j,Tnum,munum,js,Ti
  parameter(Tnum=21,munum=1,Ti=1.Q+00)
  real(16) Tmu(2,munum),dt
  real(16) T,mu,muq,mus,muK,muB,T_delta
  !temperature and chemical potential
  real(16) T_up,T_down
  real(16) l,lb 
  !polyakov loop
  real(16) l_i1,l_i2,lb_i1,lb_i2
  !polyakov loop
  real(16) rho0(2),mboson(8),mfermion(2),Vtotal,fpifK(2),h(2),Z_wave(4),ck
  real(16) Vtotal0
  real(16) k_UV_i(2,munum),k_UV(2,munum),k_UV_0(2)
  real(16) V_res(Tnum,munum),T_res(Tnum,munum),P_res(Tnum,munum),mu_res(Tnum,munum)
  real(16) pre_com(munum)
  integer  mmax
  parameter(mmax=10)
  !order of chebyshev polynomial
  real(16) dcd0mu(mmax),dcd1mu(mmax),dcd2mu(mmax),dcd3mu(mmax),dcd4mu(mmax),dcd5mu(mmax),dcd6mu(mmax)
  real(16) chi(mmax,Tnum)
  real(16) factorial
  real(16) chebev
  external chebev

  real(16) mubi

  real(16) pdata(10,250)
  integer  iT

  REAL*16 :: h_const,lamd10_const,lamd20_const,lamd01_const,Tcglue_const,hs_const

  common /Tmu/ T,mu,muq,mus,muK,iT
  common /parainput/ pdata
  common /prefit/ pre_com
  common /ini_const/ h_const,lamd10_const,lamd20_const,lamd01_const,Tcglue_const,hs_const

  open(unit=51,file='../iniconst.dat')
  READ(51,*)h_const
  READ(51,*)hs_const
  READ(51,*)lamd10_const
  READ(51,*)lamd20_const
  READ(51,*)lamd01_const
  READ(51,*)Tcglue_const
  close (51)

  open(unit=51,file='../fitres.dat')
  DO i=1,250
    read(51,*) pdata(:,i)
  END DO
  close (51)

  OPEN(UNIT=51,FILE='./M1.DAT')
    READ(51,*)J
    READ(51,*)JS
  CLOSE(51)
  WRITE(*,*)'LOAD OK',J

  muk=0.Q+00
  muBi=200.Q+0
  
  T_up=20.Q+0
  T_down=-T_up
  T_delta=cos(pi*(real(J,kind=16)-0.5Q+00)/real(Tnum,kind=16))*(0.5Q+00*(T_up-T_down))+0.5Q+00*(T_up+T_down)
  dT=1.Q+00
  muS=0.Q+0
  muB=muBi
  muq=0.Q+00
  mu=muB/3.Q+0

    l_i1=1.Q-10
    l_i2=l_i1
    lb_i1=1.Q-10
    lb_i2=lb_i1

    do i=1,250
      iT=i
      T=Ti+real(i-1,kind=16)*dT
      call selfEQ(l_i1,lb_i1,l,lb,rho0,mboson,mfermion,Vtotal,fpifK,h,Z_wave,ck)
      l_i1=2*l-l_i2
      l_i2=l
      lb_i1=2*lb-lb_i2
      lb_i2=lb
      !T_res(i,j)=T
      !mu_res(i,j)=mu
      !V_res(i,j)=Vtotal
      write(*,*)'i=',i,'j=',j,'js=',js
      write(*,"('Temperature=',F7.2)")T
      write(*,"('mu mus muq =',3F7.2)")mu,mus,muq
      write(*,"('muK =',F7.2)")muK
      write(*,"('loop       =',F20.16)")l
      write(*,"('loopbar    =',F20.16)")lb
      write(*,"('mbo=  ' 8F8.2)")mboson
      write(*,"('mfe=  ' 2F8.2)")mfermion
      write(*,"('fpifK=  ' 2F8.2)")fpifK
    
      OPEN(UNIT=51,FILE='./BUFFER/VTOTAL.DAT',POSITION='APPEND')
      WRITE(51,*)VTOTAL
      CLOSE(51)

      OPEN(UNIT=51,FILE='./BUFFER/MF.DAT',POSITION='APPEND')
      WRITE(51,*)REAL(mfermion,KIND=8)
      CLOSE(51)

      !OPEN(UNIT=51,FILE='./BUFFER/MBO.DAT',POSITION='APPEND')
      !WRITE(51,"(8F8.2)")REAL(mboson,KIND=4)
      !CLOSE(51)

      OPEN(UNIT=51,FILE='./BUFFER/FPIFK.DAT',POSITION='APPEND')
      WRITE(51,*)REAL(fpifK,KIND=8)
      CLOSE(51)

      OPEN(UNIT=51,FILE='./BUFFER/T.DAT',POSITION='APPEND')
      WRITE(51,*)REAL(T,KIND=8)
      CLOSE(51)

      !OPEN(UNIT=51,FILE='./BUFFER/MU.DAT',POSITION='APPEND')
      !WRITE(51,*)REAL(MU,KIND=8)
      !CLOSE(51)

      OPEN(UNIT=51,FILE='./BUFFER/l.DAT',POSITION='APPEND')
      WRITE(51,*)REAL(l,KIND=8)
      CLOSE(51)

      OPEN(UNIT=51,FILE='./BUFFER/lb.DAT',POSITION='APPEND')
      WRITE(51,*)REAL(lb,KIND=8)
      CLOSE(51)

    end do

end
