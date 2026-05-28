subroutine massbose2(y,dmbo2drho,cos2phi,sin2phi)
!this subroutine is just calculate the mass^2 and it's derivatives of bosons
!wenrui wrote 2017.10.20

  implicit none

  integer NMAX 
  !maximal number of differential equations
  parameter(NMAX=50)
  real(16) y(NMAX)
  integer N_str(4) 
  !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck
  !number of lam h Z ck
  real(16) lam00,lam10,lam20,lam30,lam40,lam50,lam60,lam70
  real(16) lam01,lam11,lam21,lam31,lam41,lam51
  real(16) lam02,lam12,lam22,lam32
  real(16) lam03,lam13
  !constant of Taylor expansion of potential V
  real(16) ck
  !anomaly breaking of axial U(1)A symmetry constant (KMT term)
  real(16) Sl,Ss
  !expansion points
  real(16) hd00s11,hd00s99,hd00s19,hd00s22,hd00s55
  real(16) hd00p11,hd00p99,hd00p19,hd00p22,hd00p55
  real(16) hd10s11,hd10s99,hd10s19,hd10s22,hd10s55
  real(16) hd10p11,hd10p99,hd10p19,hd10p22,hd10p55
  real(16) hd20s11,hd20s99,hd20s19,hd20s22,hd20s55
  real(16) hd20p11,hd20p99,hd20p19,hd20p22,hd20p55
  real(16) hd30s11,hd30s99,hd30s19,hd30s22,hd30s55
  real(16) hd30p11,hd30p99,hd30p19,hd30p22,hd30p55
  real(16) hd40s11,hd40s99,hd40s19,hd40s22,hd40s55
  real(16) hd40p11,hd40p99,hd40p19,hd40p22,hd40p55
  real(16) hd50s11,hd50s99,hd50s19,hd50s22,hd50s55
  real(16) hd50p11,hd50p99,hd50p19,hd50p22,hd50p55
  real(16) hd01s11,hd01s99,hd01s19,hd01s22,hd01s55
  real(16) hd01p11,hd01p99,hd01p19,hd01p22,hd01p55
  real(16) hd11s11,hd11s99,hd11s19,hd11s22,hd11s55
  real(16) hd11p11,hd11p99,hd11p19,hd11p22,hd11p55
  real(16) hd21s11,hd21s99,hd21s19,hd21s22,hd21s55
  real(16) hd21p11,hd21p99,hd21p19,hd21p22,hd21p55
  real(16) hd31s11,hd31s99,hd31s19,hd31s22,hd31s55
  real(16) hd31p11,hd31p99,hd31p19,hd31p22,hd31p55
  real(16) hd02s11,hd02s99,hd02s19,hd02s22,hd02s55
  real(16) hd02p11,hd02p99,hd02p19,hd02p22,hd02p55
  real(16) hd12s11,hd12s99,hd12s19,hd12s22,hd12s55
  real(16) hd12p11,hd12p99,hd12p19,hd12p22,hd12p55
  !hamilton and it's derivatives Appendix D 
  real(16) thetap,thetas
  !mixing angles of 08
  real(16) sin2phi(2),cos2phi(2)
  !sin^2 and cos^2 of mixing angles
  !scalar,psendoscalar
  real(16) a00,a10,a20,a30,a40,a50,a01,a11,a21,a31,a02,a12
  real(16) b00,b10,b20,b30,b40,b50,b01,b11,b21,b31,b02,b12
  real(16) mboson2(18)
  !mf02,ma02,mkap2,msig2,meta2,mpi2,mK2,meta12
  real(16) mbo2d10rho(18),mbo2d20rho(18),mbo2d30rho(18),mbo2d40rho(18),mbo2d50rho(18)
  real(16) mbo2d01rho(18),mbo2d11rho(18),mbo2d21rho(18),mbo2d31rho(18)
  real(16) mbo2d02rho(18),mbo2d12rho(18)  
  !derivatives of mass
  real(16) dmbo2drho(18,12)
  !mf02,ma02,mkap2,msig2,meta2,mpi2,mK2,meta12
  !mbo2d00 is mboson
  !00,10,20,30,40,50,01,11,21,31,02,12

  common /strucFun/ N_str

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)

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
  ck=y((Nv+1)+(Nh+2)+Nz+1)
  Sl=y((Nv+1)+(Nh+2)+Nz+2)
  Ss=y((Nv+1)+(Nh+2)+Nz+3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calculations are performed at expansion point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calculate H^(08)_ij Appendix D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calculate derivatives of H^(08)_ij
  !producted by hamilton_derivation.nb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  hd00p11=lam10 + (ck*(2*Sl + Sqrt(2.Q+00)*Ss))/3.Q+00
  hd00p22=lam10 - (ck*Ss)/Sqrt(2.Q+00) + (lam01*(Sl**2 - 2*Ss**2))/6.Q+00
  hd00p55=(6*lam10 + lam01*Sl**2 + 4*lam01*Ss**2 -                              &
      3*Sl*(ck + Sqrt(2.Q+00)*lam01*Ss))/6.Q+00
  hd00p99=(6*lam10 - 4*ck*Sl - lam01*Sl**2 + Sqrt(2.Q+00)*ck*Ss +                    &
      2*lam01*Ss**2)/6.Q+00
  hd00p19=(ck*(-(Sqrt(2.Q+00)*Sl) + 2*Ss) + Sqrt(2.Q+00)*lam01*(Sl**2 - 2*Ss**2))/6.Q+00
  hd10p11=lam20 + (ck*(4/Sl + Sqrt(2.Q+00)/Ss))/9.Q+00
  hd20p11=lam30 + (ck*(-8/Sl**3 - Sqrt(2.Q+00)/Ss**3))/27.Q+00
  hd30p11=lam40 + (ck*(16/Sl**5 + Sqrt(2.Q+00)/Ss**5))/27.Q+00
  hd40p11=lam50 + (5*ck*(-32/Sl**7 - Sqrt(2.Q+00)/Ss**7))/81.Q+00
  hd50p11=(35*ck*(64/Sl**9 + Sqrt(2.Q+00)/Ss**9))/243.Q+00
  hd01p11=-((2*Sqrt(2.Q+00)*ck*Sl - 4*ck*Ss - 3*lam11*Sl**3*Ss +                     &
        6*lam11*Sl*Ss**3)/(3*Sl**3*Ss - 6*Sl*Ss**3))
  hd11p11=(2*Sqrt(2.Q+00)*ck*Sl**3 - 8*ck*Ss**3 + 9*lam21*Sl**5*Ss**3 -              &
      18*lam21*Sl**3*Ss**5)/(9*Sl**5*Ss**3 - 18*Sl**3*Ss**5)
  hd21p11=(-2*Sqrt(2.Q+00)*ck*Sl**5 + (16*ck + 9*lam31*Sl**7)*Ss**5 -                &
      18*lam31*Sl**5*Ss**7)/(9.Q+00*Sl**5*Ss**5*(Sl**2 - 2*Ss**2))
  hd31p11=(10*ck*(Sqrt(2.Q+00)*Sl**7 - 16*Ss**7))/                                   &
    (27.Q+00*Sl**7*Ss**7*(Sl**2 - 2*Ss**2))
  hd02p11=(3*lam12*Ss**3*(Sl**3 - 2*Sl*Ss**2)**3 -                              &
      4*ck*(Sqrt(2.Q+00)*Sl**5 - 8*Sqrt(2.Q+00)*Sl**3*Ss**2 + 14*Sl**2*Ss**3 -            &
         4*Ss**5))/(3.Q+00*Sl**3*Ss**3*(Sl**2 - 2*Ss**2)**3)
  hd12p11=(4*ck*(Sqrt(2.Q+00)*Sl**7 - 4*Sqrt(2.Q+00)*Sl**5*Ss**2 +                        &
        12*Sl**2*Ss**5 - 8*Ss**7))/(3.Q+00*Sl**5*Ss**5*(Sl**2 - 2*Ss**2)**3)
  hd10p22=(6*lam20 - (Sqrt(2.Q+00)*ck)/Ss + lam11*(Sl**2 - 2*Ss**2))/6.Q+00
  hd20p22=lam30 + ((Sqrt(2.Q+00)*ck)/Ss**3 + 3*lam21*(Sl**2 - 2*Ss**2))/18.Q+00
  hd30p22=lam40 + (-((Sqrt(2.Q+00)*ck)/Ss**5) + 3*lam31*(Sl**2 - 2*Ss**2))/18.Q+00
  hd40p22=lam50 + (5*ck)/(27.Q+00*Sqrt(2.Q+00)*Ss**7)
  hd50p22=(-35*ck)/(81.Q+00*Sqrt(2.Q+00)*Ss**9)
  hd01p22=-(-6*Sqrt(2.Q+00)*ck - 12*lam01*Ss -                                       &
       Ss*(Sl**2 - 2*Ss**2)*(6*lam11 + lam02*Sl**2 - 2*lam02*Ss**2))/           &
    (6.Q+00*Ss*(Sl**2 - 2*Ss**2))
  hd11p22=  -(2*Sqrt(2.Q+00)*ck - Ss**3*(12*lam11 +                             &
          (Sl**2 - 2*Ss**2)*(6*lam21 + lam12*(Sl**2 - 2*Ss**2))))/              &
    (6.Q+00*Ss**3*(Sl**2 - 2*Ss**2))
  hd21p22=(Sqrt(2.Q+00)*ck + 6*lam21*Ss**5 + 3*lam31*Sl**2*Ss**5 -                   &
      6*lam31*Ss**7)/(3*Sl**2*Ss**5 - 6*Ss**7)
  hd31p22=(-5*Sqrt(2.Q+00)*ck + 18*lam31*Ss**7)/(9.Q+00*Ss**7*(Sl**2 - 2*Ss**2))
  hd02p22=(2*Sqrt(2.Q+00)*ck*(Sl**2 - 8*Ss**2) +                                     &
      Ss**3*(-24*lam01 + (Sl**2 - 2*Ss**2)**2*                                  &
          (4*lam02 + lam12*Sl**2 - 2*lam12*Ss**2)))/                            &
    (Ss**3*(Sl**2 - 2*Ss**2)**3)
  hd12p22=(2*(Sqrt(2.Q+00)*ck*(Sl**2 - 4*Ss**2) -                                    &
        2*Ss**5*(-6*lam11 + lam12*(Sl**2 - 2*Ss**2)**2)))/                      &
    (Ss**5*(-Sl**2 + 2*Ss**2)**3)
  hd10p55=(4*lam01 + 6*lam20 + lam11*Sl**2 + 4*lam11*Ss**2 -                    &
      (2*(ck + Sqrt(2.Q+00)*lam01*Ss))/Sl -                                          &
      (Sqrt(2.Q+00)*Sl*(lam01 + 3*lam11*Ss**2))/Ss)/6.Q+00
  hd20p55=(Sqrt(2.Q+00)*lam01*(Sl**2 - 2*Ss**2)**2 +                                 &
      Ss**2*(4*ck*Ss + 3*Sl**2*                                                 &
          (-2*Sqrt(2.Q+00)*lam11*Sl**2 +                                             &
            Sl*(8*lam11 + 6*lam30 + lam21*Sl**2)*Ss -                           &
            Sqrt(2.Q+00)*(4*lam11 + 3*lam21*Sl**2)*Ss**2 + 4*lam21*Sl*Ss**3)))/      &
    (18.Q+00*Sl**3*Ss**3)
  hd30p55=(-(Sqrt(2.Q+00)*lam01*(Sl**2 - 2*Ss**2)**2*(Sl**2 + 2*Ss**2)) +            &
      3*Sqrt(2.Q+00)*lam11*Ss**2*(Sl**3 - 2*Sl*Ss**2)**2 +                           &
      Ss**4*(-8*ck*Ss + 3*Sl**4*                                                &
          (lam31*Sl**3*Ss - 6*Sqrt(2.Q+00)*lam21*Ss**2 -                             &
            3*Sqrt(2.Q+00)*Sl**2*(lam21 + lam31*Ss**2) +                             &
            2*Sl*Ss*(6*lam21 + 3*lam40 + 2*lam31*Ss**2))))/(18.Q+00*Sl**5*Ss**5)
  hd40p55=(Sqrt(2.Q+00)*lam01*(Sl**2 - 2*Ss**2)**2*                                  &
       (5*Sl**4 + 12*Sl**2*Ss**2 + 20*Ss**4) +                                  &
      2*Ss**2*(-6*Sqrt(2.Q+00)*lam11*(Sl**2 + 2*Ss**2)*                              &
          (Sl**3 - 2*Sl*Ss**2)**2 +                                             &
         9*Sqrt(2.Q+00)*lam21*Sl**4*(Sl**2*Ss - 2*Ss**3)**2 +                        &
         Ss**4*(40*ck*Ss + 9*Sl**6*                                             &
             (3*lam50*Sl*Ss - 2*lam31*                                          &
                (Sqrt(2.Q+00)*Sl**2 - 4*Sl*Ss + 2*Sqrt(2.Q+00)*Ss**2)))))/                &
    (54.Q+00*Sl**7*Ss**7)
hd50p55=(-5*Sqrt(2.Q+00)*lam01*(Sl**2 - 2*Ss**2)**2*(Sl**2 + 2*Ss**2)*               &
       (7*Sl**4 + 4*Sl**2*Ss**2 + 28*Ss**4) +                                   &
      5*Ss**2*(-224*ck*Ss**7 +                                                  &
         18*Sqrt(2.Q+00)*lam31*Sl**6*Ss**4*(Sl**2 - 2*Ss**2)**2 -                    &
         18*Sqrt(2.Q+00)*lam21*Sl**4*(Sl**2 + 2*Ss**2)*                              &
          (Sl**2*Ss - 2*Ss**3)**2 +                                             &
         3*Sqrt(2.Q+00)*lam11*(Sl**3 - 2*Sl*Ss**2)**2*                               &
          (5*Sl**4 + 12*Sl**2*Ss**2 + 20*Ss**4)))/(162.Q+00*Sl**9*Ss**9)
  hd01p55=(6*lam11 + lam02*Sl**2 + (6*Sqrt(2.Q+00)*lam01)/(Sl*Ss) -                  &
      3*Sqrt(2.Q+00)*lam02*Sl*Ss + 4*lam02*Ss**2 -                                   &
      (6*(ck + 2*lam01*Sl - Sqrt(2.Q+00)*lam01*Ss))/(Sl**3 - 2*Sl*Ss**2))/6.Q+00
  hd11p55=(4*lam02 + 6*lam21 + lam12*Sl**2 -                                    &
      (2*Sqrt(2.Q+00)*lam01)/(Sl*Ss**3) -                                            &
      (Sqrt(2.Q+00)*(2*lam01 - 6*lam11*Sl**2 + lam02*Sl**4))/(Sl**3*Ss) -            &
      (Sqrt(2.Q+00)*(2*lam02 + 3*lam12*Sl**2)*Ss)/Sl + 4*lam12*Ss**2 +               &
   (4*ck + 6*lam11*Sl**2*(-2*Sl + Sqrt(2.Q+00)*Ss))/(Sl**5 - 2*Sl**3*Ss**2))/6.Q+00
  hd21p55=  (6*Sqrt(2.Q+00)*lam01*(Sl**6 - Sl**4*Ss**2 - 4*Ss**6) +            &
      Ss**2*(Sqrt(2.Q+00)*lam02*Sl**2*(Sl**2 - 2*Ss**2)**3 -                         &
         12*Sqrt(2.Q+00)*lam11*Sl**2*(Sl**2 - 2*Ss**2)*(Sl**2 + Ss**2) -             &
         6*Ss**2*(4*ck*Ss + Sl**4*                                              &
             (3*lam21*(-(Sqrt(2.Q+00)*Sl**2) + 2*Sl*Ss + Sqrt(2.Q+00)*Ss**2) +            &
               (Sl**2 - 2*Ss**2)*                                               &
                (Sqrt(2.Q+00)*lam12*Sl**2 - 4*lam12*Sl*Ss - 3*lam31*Sl*Ss +          &
                  2*Sqrt(2.Q+00)*lam12*Ss**2)))))/                                   &
    (18.Q+00*Sl**5*Ss**5*(Sl**2 - 2*Ss**2))
  hd31p55=(2*Sqrt(2.Q+00)*lam01*                                                     &
       (-5*Sl**8 + 5*Sl**6*Ss**2 + 2*Sl**4*Ss**4 - 4*Sl**2*Ss**6 +              &
         40*Ss**8) + Ss**2*(80*ck*Ss**5 - 36*lam31*Sl**7*Ss**5 +                &
         8*Sqrt(2.Q+00)*Sl**2*Ss**6*(-9*lam11 + 2*lam02*Ss**2) -                     &
         Sqrt(2.Q+00)*Sl**10*(lam02 - 3*lam12*Ss**2) -                               &
         4*Sqrt(2.Q+00)*Sl**4*Ss**6*(4*lam02 - 9*lam21 + 6*lam12*Ss**2) -            &
         18*Sqrt(2.Q+00)*Sl**6*Ss**2*                                                &
          (lam11 - lam21*Ss**2 + (-2*lam12 + lam31)*Ss**4) +                    &
         2*Sqrt(2.Q+00)*Sl**8*(9*lam11 +                                             &
            Ss**2*(2*lam02 - 9*lam21 + 9*(-lam12 + lam31)*Ss**2))))/            &
    (18.Q+00*Sl**7*Ss**7*(Sl**2 - 2*Ss**2))
  hd02p55=lam12 + (2*Sqrt(2.Q+00)*lam01)/(Sl**3*Ss**3) +                             &
    (2*Sqrt(2.Q+00)*lam02)/(Sl*Ss) +                                                 &
    (12*(ck + 2*lam01*Sl - Sqrt(2.Q+00)*lam01*Ss))/(Sl*(Sl**2 - 2*Ss**2)**3) +       &
    (2*(ck - 3*Sqrt(2.Q+00)*lam01*Ss))/(Sl**3*(Sl**2 - 2*Ss**2)**2) +                &
    (lam02*(-4*Sl + 2*Sqrt(2.Q+00)*Ss))/(Sl**3 - 2*Sl*Ss**2)
  hd12p55=(-2*(3*Sqrt(2.Q+00)*lam01*(Sl**2 - 2*Ss**2)**2*(Sl**4 - Ss**4) +           &
        Ss**2*(Sqrt(2.Q+00)*lam02*Sl**2*(Sl**2 - 2*Ss**2)**3*(Sl**2 + Ss**2) -       &
           3*lam11*(Sqrt(2.Q+00)*Sl**8 - 6*Sqrt(2.Q+00)*Sl**6*Ss**2 +                     &
              12*Sl**5*Ss**3 + 3*Sqrt(2.Q+00)*Sl**4*Ss**4 -                          &
              2*Sqrt(2.Q+00)*Sl**2*Ss**6) +                                          &
           3*Ss**2*(6*ck*Sl**2*Ss - 4*ck*Ss**3 -                                &
              lam12*Sl**4*(Sl**2 - 2*Ss**2)**2*                                 &
               (Sqrt(2.Q+00)*Sl**2 - 2*Sl*Ss - Sqrt(2.Q+00)*Ss**2)))))/                   &
    (3.Q+00*Sl**5*Ss**5*(Sl**2 - 2*Ss**2)**3)
  hd10p99=(18*lam20 - (8*ck)/Sl - 3*lam11*Sl**2 + (Sqrt(2.Q+00)*ck)/Ss +             &
      6*lam11*Ss**2)/18.Q+00
  hd20p99= (54*lam30 + (16*ck)/Sl**3 - 9*lam21*Sl**2 - (Sqrt(2.Q+00)*ck)/Ss**3 +     &
      18*lam21*Ss**2)/54.Q+00
  hd30p99= (54*lam40 - (32*ck)/Sl**5 - 9*lam31*Sl**2 + (Sqrt(2.Q+00)*ck)/Ss**5 +     &
      18*lam31*Ss**2)/54.Q+00
  hd40p99=lam50 + (5*ck*(64/Sl**7 - Sqrt(2.Q+00)/Ss**7))/162.Q+00
  hd50p99=(35*ck*(-128/Sl**9 + Sqrt(2.Q+00)/Ss**9))/486.Q+00
  hd01p99=-(2*ck*(Sqrt(2.Q+00)*Sl + 4*Ss) +                                          &
       Sl*Ss*(12*lam01 + (Sl**2 - 2*Ss**2)*                                     &
           (-6*lam11 + lam02*Sl**2 - 2*lam02*Ss**2)))/                          &
    (6.Q+00*Sl*Ss*(Sl**2 - 2*Ss**2))
  hd11p99=(2*ck*(Sqrt(2.Q+00)*Sl**3 + 8*Ss**3) -                                     &
      3*Sl**3*Ss**3*(12*lam11 +                                                 &
         (Sl**2 - 2*Ss**2)*(-6*lam21 + lam12*Sl**2 - 2*lam12*Ss**2)))/          &
    (18.Q+00*Sl**3*Ss**3*(Sl**2 - 2*Ss**2))
  hd21p99=(-(ck*(Sqrt(2.Q+00)*Sl**5 + 16*Ss**5)) +                                   &
      9*Sl**5*Ss**5*(-2*lam21 + lam31*(Sl**2 - 2*Ss**2)))/                      &
    (9.Q+00*Sl**5*Ss**5*(Sl**2 - 2*Ss**2))
  hd31p99=(5*Sqrt(2.Q+00)*ck*Sl**7 + 160*ck*Ss**7 - 54*lam31*Sl**7*Ss**7)/           &
    (27*Sl**9*Ss**7 - 54*Sl**7*Ss**9)
  hd02p99=(-2*ck*(Sqrt(2.Q+00)*Sl**5 - 8*Sqrt(2.Q+00)*Sl**3*Ss**2 -                       &
         28*Sl**2*Ss**3 + 8*Ss**5) +                                            &
      3*Sl**3*Ss**3*(24*lam01 +                                                 &
         (Sl**2 - 2*Ss**2)**2*(-4*lam02 + lam12*Sl**2 - 2*lam12*Ss**2)))/       &
    (3.Q+00*Sl**3*Ss**3*(Sl**2 - 2*Ss**2)**3)
  hd12p99=(2*(ck*(Sqrt(2.Q+00)*Sl**7 - 4*Sqrt(2.Q+00)*Sl**5*Ss**2 -                       &
           24*Sl**2*Ss**5 + 16*Ss**7) -                                         &
        6*Sl**5*Ss**5*(-6*lam11 + lam12*(Sl**2 - 2*Ss**2)**2)))/                &
    (3.Q+00*Sl**5*Ss**5*(Sl**2 - 2*Ss**2)**3)
  hd10p19=(2*ck*(Sl - Sqrt(2.Q+00)*Ss) +                                             &
      3*Sqrt(2.Q+00)*lam11*Sl*Ss*(Sl**2 - 2*Ss**2))/(18.Q+00*Sl*Ss)
  hd20p19=(ck*((4*Sqrt(2.Q+00))/Sl**3 - 2/Ss**3) +                                   &
      9*Sqrt(2.Q+00)*lam21*(Sl**2 - 2*Ss**2))/54.Q+00
  hd30p19=(2*ck*((-4*Sqrt(2.Q+00))/Sl**5 + Ss**(-5)) +                               &
      9*Sqrt(2.Q+00)*lam31*(Sl**2 - 2*Ss**2))/54.Q+00
  hd40p19=(5*ck*((8*Sqrt(2.Q+00))/Sl**7 - Ss**(-7)))/81.Q+00
  hd50p19=(35*ck*((-16*Sqrt(2.Q+00))/Sl**9 + Ss**(-9)))/243.Q+00
  hd01p19=(-2*ck*(2*Sl + Sqrt(2.Q+00)*Ss) +                                          &
      Sqrt(2.Q+00)*Sl*Ss*(12*lam01 + lam02*(Sl**2 - 2*Ss**2)**2))/                   &
    (6.Q+00*Sl*Ss*(Sl**2 - 2*Ss**2))
  hd11p19=(4*ck*(Sl**3 + Sqrt(2.Q+00)*Ss**3) +                                       &
      3*Sqrt(2.Q+00)*Sl**3*Ss**3*(12*lam11 + lam12*(Sl**2 - 2*Ss**2)**2))/           &
    (18.Q+00*Sl**3*Ss**3*(Sl**2 - 2*Ss**2))
  hd21p19=(-2*(-9*Sqrt(2.Q+00)*lam21*Sl**5*Ss**5 +                                   &
        ck*(Sl**5 + 2*Sqrt(2.Q+00)*Ss**5)))/(9.Q+00*Sl**5*Ss**5*(Sl**2 - 2*Ss**2))
  hd31p19=(10*ck*Sl**7 + 40*Sqrt(2.Q+00)*ck*Ss**7 +                                  &
      54*Sqrt(2.Q+00)*lam31*Sl**7*Ss**7)/(27*Sl**9*Ss**7 - 54*Sl**7*Ss**9)
  hd02p19=(2*(ck*(-2*Sl**5 + 16*Sl**3*Ss**2 + 7*Sqrt(2.Q+00)*Sl**2*Ss**3 -           &
           2*Sqrt(2.Q+00)*Ss**5) + 6*Sqrt(2.Q+00)*Sl**3*Ss**3*                            &
         (-6*lam01 + lam02*(Sl**2 - 2*Ss**2)**2)))/                             &
    (3.Q+00*Sl**3*Ss**3*(Sl**2 - 2*Ss**2)**3)
  hd12p19=(4*(ck*(Sl**7 - 4*Sl**5*Ss**2 - 3*Sqrt(2.Q+00)*Sl**2*Ss**5 +               &
           2*Sqrt(2.Q+00)*Ss**7) + 3*Sqrt(2.Q+00)*Sl**5*Ss**5*                            &
         (-6*lam11 + lam12*(Sl**2 - 2*Ss**2)**2)))/                             &
    (3.Q+00*Sl**5*Ss**5*(Sl**2 - 2*Ss**2)**3)
  hd00s11=(12*lam01*(Sl**2 - 2*Sqrt(2.Q+00)*Sl*Ss + 2*Ss**2) +                       &
      lam02*(Sl**2 - 2*Ss**2)**2*(Sl**2 - 2*Sqrt(2.Q+00)*Sl*Ss + 2*Ss**2) +          &
      6*(9*lam10 - 3*ck*(2*Sl + Sqrt(2.Q+00)*Ss) +                                   &
         lam11*(Sl**2 - 2*Ss**2)*(2*Sl**2 - Sqrt(2.Q+00)*Sl*Ss - 2*Ss**2) +          &
         3*lam20*(2*Sl**2 + 2*Sqrt(2.Q+00)*Sl*Ss + Ss**2)))/54.Q+00
  hd00s22=lam10 + (ck*Ss)/Sqrt(2.Q+00) + (lam01*(7*Sl**2 - 2*Ss**2))/6.Q+00
  hd00s55=(6*lam10 + 3*ck*Sl +                                                  &
      lam01*(Sl**2 + 3*Sqrt(2.Q+00)*Sl*Ss + 4*Ss**2))/6.Q+00
  hd00s99=(-6*lam01*(Sl**2 - 8*Sqrt(2.Q+00)*Sl*Ss - 22*Ss**2) +                      &
      lam02*(Sl**2 - 2*Ss**2)**2*(Sl**2 + 4*Sqrt(2.Q+00)*Sl*Ss + 8*Ss**2) +          &
      6*(18*lam10 + 3*ck*(4*Sl - Sqrt(2.Q+00)*Ss) +                                  &
         2*lam11*(Sl**2 + Sqrt(2.Q+00)*Sl*Ss - 4*Ss**2)*(Sl**2 - 2*Ss**2) +          &
         6*lam20*(Sl**2 - 2*Sqrt(2.Q+00)*Sl*Ss + 2*Ss**2)))/108.Q+00
  hd00s19=(18*ck*(Sqrt(2.Q+00)*Sl - 2*Ss) +                                          &
      6*lam01*(5*Sqrt(2.Q+00)*Sl**2 + 4*Sl*Ss - 14*Sqrt(2.Q+00)*Ss**2) +                  &
      36*lam20*(Sqrt(2.Q+00)*Sl**2 - Sl*Ss - Sqrt(2.Q+00)*Ss**2) +                        &
      (Sl**2 - 2*Ss**2)*(lam02*(Sl**2 - 2*Ss**2)*                               &
          (Sqrt(2.Q+00)*Sl**2 + 2*Sl*Ss - 4*Sqrt(2.Q+00)*Ss**2) +                         &
         6*lam11*(2*Sqrt(2.Q+00)*Sl**2 + Sl*Ss + 4*Sqrt(2.Q+00)*Ss**2)))/108.Q+00
  hd10s11=(96*lam01 + 342*lam20 + 3*lam12*Sl**6 - (18*Sqrt(2.Q+00)*ck)/Ss +          &
      6*(4*lam11 + 9*lam30)*Ss**2 + 8*(4*lam02 + 9*lam21)*Ss**4 +               &
      24*lam12*Ss**6 + 2*Sl**4*(4*lam02 + 18*lam21 - 3*lam12*Ss**2) -           &
      (2*Sqrt(2.Q+00)*Sl**5*(lam02 + 3*lam12*Ss**2))/Ss -                            &
      (4*Sqrt(2.Q+00)*Sl*(6*lam01 - 9*lam20 + 9*(2*lam11 - 3*lam30)*Ss**2 -          &
           (2*lam02 + 9*lam21)*Ss**4 + 6*lam12*Ss**6))/Ss +                     &
      4*Sl**2*(15*lam11 + 27*lam30 -                                            &
         Ss**2*(8*lam02 + 27*lam21 + 3*lam12*Ss**2)) +                          &
      (2*Sqrt(2.Q+00)*Sl**3*(-3*lam11 +                                              &
           Ss**2*(2*lam02 - 9*lam21 + 12*lam12*Ss**2)))/Ss +                    &
      (-72*ck - 8*Sqrt(2.Q+00)*Ss*(6*lam01 - 9*lam20 - 3*lam11*Ss**2 +               &
            2*lam02*Ss**4))/Sl)/162.Q+00
  hd20s11=(288*lam11 + 783*lam30 + 6*(4*lam12 + 9*lam31)*Sl**4 +                &
      (9*Sqrt(2.Q+00)*ck)/Ss**3 + 9*(-4*lam21 + 9*lam40)*Ss**2 +                     &
      12*(8*lam12 + 9*lam31)*Ss**4 +                                            &
      (Sqrt(2.Q+00)*Sl**5*(lam02 - 6*lam12*Ss**2))/Ss**3 +                           &
      6*Sl**2*(21*lam21 + 27*lam40 - (16*lam12 + 27*lam31)*Ss**2) +             &
      (Sqrt(2.Q+00)*Sl**3*(3*lam11 - 2*(4*lam02 + 9*lam21)*Ss**2 +                   &
           3*(4*lam12 - 9*lam31)*Ss**4))/Ss**3 -                                &
      (4*Sqrt(2.Q+00)*(12*lam01 - 18*lam20 + 27*(lam11 - 2*lam30)*Ss**2 +            &
           2*(4*lam02 - 9*lam21)*Ss**4 + 12*lam12*Ss**6))/(Sl*Ss) +             &
      (6*Sqrt(2.Q+00)*Sl*(2*lam01 - 3*lam20 + 3*(-5*lam11 + 6*lam30)*Ss**2 +         &
           (4*lam02 - 18*lam21 + 27*lam40)*Ss**4 +                              &
           (4*lam12 + 9*lam31)*Ss**6))/Ss**3 +                                  &
      (72*ck + 8*Sqrt(2.Q+00)*Ss*(6*lam01 - 9*lam20 - 3*lam11*Ss**2 +                &
            2*lam02*Ss**4))/Sl**3)/243.Q+00
  hd30s11=(-(Sqrt(2.Q+00)*Sl**5*                                                     &
         (9*ck + 12*lam01*Sl - 18*lam20*Sl + 3*lam11*Sl**3 + lam02*Sl**5))      &
+ 3*Sqrt(2.Q+00)*Sl**4*(8*lam01 - 12*lam20 + 2*(8*lam11 - 9*lam30)*Sl**2 +           &
         (2*lam02 + 3*lam21)*Sl**4 + lam12*Sl**6)*Ss**2 -                       &
      Sqrt(2.Q+00)*Sl**2*(-48*lam01 + 72*lam20 +                                     &
         72*(2*lam11 - 3*lam30)*Sl**2 +                                         &
         2*(4*lam02 + 81*lam21 - 81*lam40)*Sl**4 +                              &
         3*(8*lam12 + 9*lam31)*Sl**6)*Ss**4 +                                   &
      9*(-16*ck + 3*Sl**5*(16*lam21 + 39*lam40 +                                &
            6*(lam31 + lam50)*Sl**2))*Ss**5 -                                   &
      2*Sqrt(2.Q+00)*(48*lam01 - 72*lam20 + 12*(-4*lam11 + 9*lam30)*Sl**2 +          &
         2*(4*lam02 + 27*lam21 - 81*lam40)*Sl**4 -                              &
         9*(4*lam12 - 6*lam31 + 9*lam50)*Sl**6)*Ss**6 +                         &
      27*(-4*lam31 + 3*lam50)*Sl**5*Ss**7 +                                     &
      12*Sqrt(2.Q+00)*(4*lam11 + Sl**2*                                              &
          (4*lam02 - 6*lam21 + (-8*lam12 + 9*lam31)*Sl**2))*Ss**8 +             &
      16*Sqrt(2.Q+00)*(-2*lam02 + 3*lam12*Sl**2)*Ss**10)/(243.Q+00*Sl**5*Ss**5)
  hd40s11=(Sqrt(2.Q+00)*Sl**12*(5*lam02 - 12*lam12*Ss**2) +                          &
      Sqrt(2.Q+00)*Sl**10*(15*lam11 - 4*(7*lam02 + 9*lam21)*Ss**2 +                  &
         18*(4*lam12 + 3*lam31)*Ss**4) -                                        &
      16*Sqrt(2.Q+00)*Sl**2*Ss**6*(24*lam01 - 36*lam20 +                             &
         9*(5*lam11 - 12*lam30)*Ss**2 + 4*(7*lam02 - 9*lam21)*Ss**4 +           &
         24*lam12*Ss**6) + 2*Sqrt(2.Q+00)*Sl**8*                                     &
       (30*lam01 - 45*lam20 + 9*(-11*lam11 + 12*lam30)*Ss**2 +                  &
         2*(11*lam02 + 90*lam21 - 81*lam40)*Ss**4 -                             &
         6*(8*lam12 + 63*lam31 - 54*lam50)*Ss**6) +                             &
      9*Sl**7*(5*Sqrt(2.Q+00)*ck + 3*(64*lam31 + 147*lam50)*Ss**7) +                 &
      160*Ss**7*(9*ck + Sqrt(2.Q+00)*Ss*                                             &
          (6*lam01 - 9*lam20 - 3*lam11*Ss**2 + 2*lam02*Ss**4)) +                &
      16*Sqrt(2.Q+00)*Sl**4*Ss**4*(-6*lam01 + 9*lam20 +                              &
         Ss**2*(33*lam11 - 54*lam30 +                                           &
            (11*lam02 + 18*lam21 - 81*lam40)*Ss**2 +                            &
            9*(4*lam12 - 3*lam31)*Ss**4)) -                                     &
      8*Sqrt(2.Q+00)*Sl**6*Ss**2*(12*lam01 - 18*lam20 +                              &
         Ss**2*(-39*lam11 + 54*lam30 +                                          &
            2*(2*lam02 + 54*lam21 - 81*lam40)*Ss**2 +                           &
            3*(8*lam12 + 9*lam31 - 54*lam50)*Ss**4)))/(729.Q+00*Sl**7*Ss**7)
  hd50s11=(-5*(63*ck*(Sqrt(2.Q+00)*Sl**9 + 64*Ss**9) +                               &
        Sqrt(2.Q+00)*(Sl**2 - 2*Ss**2)**2*                                           &
         (7*Sl**6*(12*lam01 - 18*lam20 + 3*lam11*Sl**2 + lam02*Sl**4) -         &
           Sl**4*(-216*lam01 + 324*lam20 +                                      &
              6*(28*lam11 - 45*lam30)*Sl**2 +                                   &
              5*(2*lam02 + 9*lam21)*Sl**4 + 15*lam12*Sl**6)*Ss**2 +             &
           2*Sl**2*(216*lam01 - 324*lam20 +                                     &
              108*(-2*lam11 + 3*lam30)*Sl**2 +                                  &
              (-4*lam02 + 99*lam21 - 162*lam40)*Sl**4 +                         &
              3*(4*lam12 + 9*lam31)*Sl**6)*Ss**4 +                              &
           4*(168*lam01 - 252*lam20 + 6*(-32*lam11 + 45*lam30)*Sl**2 +          &
              (-4*lam02 + 117*lam21 - 162*lam40)*Sl**4 +                        &
              3*(2*lam12 - 18*lam31 + 27*lam50)*Sl**6)*Ss**6 -                  &
           8*(42*lam11 + 5*(2*lam02 - 9*lam21)*Sl**2 +                          &
              3*(-4*lam12 + 9*lam31)*Sl**4)*Ss**8 +                             &
           16*(14*lam02 - 15*lam12*Sl**2)*Ss**10)))/(2187.Q+00*Sl**9*Ss**9)
  hd01s11=(147*lam11 + 6*lam12*Sl**4 + 40*lam02*Ss**2 + 9*lam21*Ss**2 +         &
      12*lam12*Ss**4 + 2*Sl**2*(8*lam02 + 9*lam21 - 9*lam12*Ss**2) +            &
      (Sqrt(2.Q+00)*Sl**3*(2*lam02 - 3*lam12*Ss**2))/Ss +                            &
      (6*(3*ck*(Sl - Sqrt(2.Q+00)*Ss) +                                              &
           (Sqrt(2.Q+00)*Sl - 2*Ss)*Ss*(-2*lam01 + 3*lam20 + 6*lam11*Ss**2)))/       &
       (-(Sl**2*Ss**2) + 2*Ss**4) +                                             &
      (6*Sqrt(2.Q+00)*Sl*(lam11 + Ss**2*(-7*lam02 + 3*lam21 + lam12*Ss**2)))/        &
       Ss + (18*ck + 2*Sqrt(2.Q+00)*Ss*                                              &
     (6*lam01 - 9*lam20 - 3*lam11*Ss**2 + 2*lam02*Ss**4))/(Sl*Ss**2))/27.Q+00
  hd11s11= (6*(10*lam12 + 9*lam31)*Sl**2 +                                          &
      (Sqrt(2.Q+00)*Sl**3*(-2*lam02 + 3*lam12*Ss**2))/Ss**3 +                        &
      3*(48*lam02 + 177*lam21 + (32*lam12 + 9*lam31)*Ss**2) +                   &
      (18*(ck*(Sl - Sqrt(2.Q+00)*Ss) -                                               &
           (Sqrt(2.Q+00)*Sl - 2*Ss)*Ss**3*(2*lam11 + 3*lam30 + 6*lam21*Ss**2))       &
)/(Ss**4*(Sl**2 - 2*Ss**2)) - (6*Sqrt(2.Q+00)*Sl*                                    &
         (lam11 + Ss**2*(5*lam02 - 6*lam21 +                                    &
              3*(7*lam12 - 3*lam31)*Ss**2)))/Ss**3 +                            &
      (-36*ck + 4*Sqrt(2.Q+00)*Ss*(-6*lam01 + 9*lam20 + 3*lam11*Ss**2 -              &
            2*lam02*Ss**4))/(Sl**3*Ss**2) +                                     &
      (-18*ck + 6*Sqrt(2.Q+00)*Ss*(-4*lam01 + 6*lam20 +                              &
            (13*lam11 - 9*lam30)*Ss**2 + 3*(-4*lam02 + lam21)*Ss**4 +           &
            4*lam12*Ss**6))/(Sl*Ss**4))/81.Q+00
  hd21s11=(9*(32*lam12 + 69*lam31) +                                            &
      (Sqrt(2.Q+00)*Sl**3*(2*lam02 - 3*lam12*Ss**2))/Ss**5 +                         &
      (18*(ck*(Sl - Sqrt(2.Q+00)*Ss) +                                               &
           3*(Sqrt(2.Q+00)*Sl - 2*Ss)*Ss**5*(2*lam21 + lam40 + 2*lam31*Ss**2))       &
)/(-(Sl**2*Ss**6) + 2*Ss**8) +                                                  &
      (6*Sqrt(2.Q+00)*Sl*(lam11 + Ss**2*                                             &
            (lam02 - 3*lam21 + (-11*lam12 + 9*lam31)*Ss**2)))/Ss**5 +           &
      (72*ck + 8*Sqrt(2.Q+00)*Ss*(6*lam01 - 9*lam20 - 3*lam11*Ss**2 +                &
            2*lam02*Ss**4))/(Sl**5*Ss**2) +                                     &
      (36*ck + 4*Sqrt(2.Q+00)*Ss*(6*lam01 - 9*lam20 +                                &
            18*(-lam11 + lam30)*Ss**2 + 10*lam02*Ss**4 - 6*lam12*Ss**6))/       &
       (Sl**3*Ss**4) + (18*ck -                                                 &
         2*Sqrt(2.Q+00)*Ss*(-12*lam01 + 18*lam20 +                                   &
            3*(7*lam11 - 12*lam30)*Ss**2 +                                      &
            (22*lam02 - 72*lam21 + 27*lam40)*Ss**4 +                            &
            3*(22*lam12 - 9*lam31)*Ss**6))/(Sl*Ss**6))/81.Q+00
  hd31s11=(-10*Sqrt(2.Q+00)*Sl**7*                                                   &
       (9*ck + 12*lam01*Sl - 18*lam20*Sl + 3*lam11*Sl**3 + lam02*Sl**5)*Ss      &
+ Sqrt(2.Q+00)*Sl**6*(120*lam01 - 180*lam20 + 54*(5*lam11 - 6*lam30)*Sl**2 +         &
         2*(7*lam02 + 36*lam21)*Sl**4 + 15*lam12*Sl**6)*Ss**3 +                 &
      2*Sqrt(2.Q+00)*Sl**4*(24*lam01 - 36*lam20 +                                    &
         6*(-16*lam11 + 27*lam30)*Sl**2 +                                       &
         (34*lam02 - 171*lam21 + 162*lam40)*Sl**4 +                             &
         18*(lam12 - 3*lam31)*Sl**6)*Ss**5 -                                    &
      2*Sqrt(2.Q+00)*Sl**2*(48*lam01 - 72*lam20 - 48*lam11*Sl**2 +                   &
         Sl**4*(-8*lam02 + 54*(lam21 + 3*lam40) +                               &
            3*(88*lam12 - 63*lam31 + 54*lam50)*Sl**2))*Ss**7 +                  &
      36*(40*ck + 3*(10*lam31 + 3*lam50)*Sl**7)*Ss**8 +                         &
      4*Sqrt(2.Q+00)*(240*lam01 - 360*lam20 +                                        &
         108*(-2*lam11 + 3*lam30)*Sl**2 -                                       &
         2*(56*lam02 - 117*lam21 + 81*lam40)*Sl**4 +                            &
         3*(92*lam12 - 126*lam31 + 27*lam50)*Sl**6)*Ss**9 -                     &
      8*Sqrt(2.Q+00)*(60*lam11 - 2*(14*lam02 + 9*lam21)*Sl**2 +                      &
         27*(2*lam12 - lam31)*Sl**4)*Ss**11 +                                   &
      64*Sqrt(2.Q+00)*(5*lam02 - 6*lam12*Sl**2)*Ss**13)/                             &
    (243.Q+00*Ss**8*(Sl**9 - 2*Sl**7*Ss**2))
  hd02s11=(267*lam12 + (4*Sqrt(2.Q+00)*Sl*(lam02 + 3*lam12*Ss**2))/Ss**3 +           &
      (18*(Sqrt(2.Q+00)*(2*lam01 - 3*lam20)*Sl*Ss +                                  &
           2*Sqrt(2.Q+00)*lam11*Sl*Ss**3 - 16*lam11*Ss**4 +                          &
           ck*(-5*Sl + 2*Sqrt(2.Q+00)*Ss)))/(Ss**4*(Sl**2 - 2*Ss**2)**2) +           &
      (-216*ck*(Sl - Sqrt(2.Q+00)*Ss) -                                              &
         72*(Sqrt(2.Q+00)*Sl - 2*Ss)*Ss*(-2*lam01 + 3*lam20 + 6*lam11*Ss**2))/       &
       (Ss**2*(-Sl**2 + 2*Ss**2)**3) +                                          &
      (3*(12*ck*Sl + 15*Sqrt(2.Q+00)*lam11*Sl*Ss**3 -                                &
           4*(Sqrt(2.Q+00)*Sl - 2*Ss)*Ss**5*                                         &
            (-5*lam02 + 3*lam21 + 6*lam12*Ss**2)))/                             &
       (Ss**6*(Sl**2 - 2*Ss**2)) +                                              &
      (18*ck + 2*Sqrt(2.Q+00)*Ss*(6*lam01 - 9*lam20 - 3*lam11*Ss**2 +                &
            2*lam02*Ss**4))/(Sl**3*Ss**4) -                                     &
      (36*ck + Sqrt(2.Q+00)*Ss**3*(33*lam11 +                                        &
            4*Ss**2*(-17*lam02 + 9*lam21 + 3*lam12*Ss**2)))/(Sl*Ss**6))/27.Q+00
  hd12s11=(4*(Sqrt(2.Q+00)*lam02*Sl**12 + 6*(lam12 - 3*lam31)*Sl**9*Ss**5 +          &
        Sqrt(2.Q+00)*Sl**8*(12*lam01 - 18*lam20 + 18*(-lam11 + lam30)*Ss**2 +        &
           (-37*lam02 + 18*lam21)*Ss**4 + 9*(17*lam12 - 10*lam31)*Ss**6) +      &
        2*Sqrt(2.Q+00)*Sl**4*Ss**4*                                                  &
         (18*lam01 - 27*lam20 + 3*(-19*lam11 + 9*lam30)*Ss**2 +                 &
           4*(5*lam02 - 9*lam21)*Ss**4 + 6*(13*lam12 - 6*lam31)*Ss**6) +        &
        3*Sl**7*(3*Sqrt(2.Q+00)*ck + 24*lam21*Ss**5 -                                &
           8*(lam12 - 3*lam31)*Ss**7) -                                         &
        8*Ss**7*(9*ck + Sqrt(2.Q+00)*Ss*                                             &
            (6*lam01 - 9*lam20 - 3*lam11*Ss**2 + 2*lam02*Ss**4)) +              &
        12*Sl**5*Ss**2*(-3*Sqrt(2.Q+00)*ck +                                         &
           Ss**3*(6*lam11 + 9*lam30 + 6*lam21*Ss**2 +                           &
              2*(lam12 - 3*lam31)*Ss**4)) +                                     &
        4*Sl**2*Ss**5*(27*ck +                                                  &
           Sqrt(2.Q+00)*Ss*(12*lam01 - 18*lam20 + 9*(lam11 - lam30)*Ss**2 +          &
              (-16*lam02 + 9*lam21)*Ss**4 + 6*lam12*Ss**6)) +                   &
        Sqrt(2.Q+00)*Sl**10*(3*lam11 +                                               &
           Ss**2*(4*lam02 - 9*(lam21 + (3*lam12 - 2*lam31)*Ss**2))) +           &
        3*Sqrt(2.Q+00)*Sl**6*Ss**2*(-16*lam01 + 24*lam20 +                           &
           Ss**2*(11*lam11 - 36*lam30 +                                         &
              Ss**2*(16*lam02 - 27*lam21 - 94*lam12*Ss**2 + 48*lam31*Ss**2)     &
))))/(27.Q+00*Sl**5*Ss**5*(-Sl**2 + 2*Ss**2)**3)
  hd10s22=(Sqrt(2.Q+00)*ck + 8*lam01*Ss + 6*lam20*Ss + 7*lam11*Sl**2*Ss -            &
      2*lam11*Ss**3)/(6.Q+00*Ss)
  hd20s22=(8*lam11)/3.Q+00 + lam30 + (21*lam21*Sl**2 - (Sqrt(2.Q+00)*ck)/Ss**3 -      &
       6*lam21*Ss**2)/18.Q+00
  hd30s22= 4*lam21 + lam40 + (21*lam31*Sl**2 + (Sqrt(2.Q+00)*ck)/Ss**5 - 6*lam31*Ss**2)/     &
     18.Q+00
  hd40s22=(16*lam31)/3.Q+00 + lam50 - (5*ck)/(27.Q+00*Sqrt(2.Q+00)*Ss**7)
  hd50s22=(35*ck)/(81.Q+00*Sqrt(2.Q+00)*Ss**9)
  hd01s22=-(6*Sqrt(2.Q+00)*ck - 36*lam01*Ss -                                        &
       Ss*(Sl**2 - 2*Ss**2)*(6*lam11 + 7*lam02*Sl**2 - 2*lam02*Ss**2))/         &
    (6.Q+00*Ss*(Sl**2 - 2*Ss**2))
  hd11s22= -(-2*Sqrt(2.Q+00)*ck - Ss**3*(36*lam11 +                                  &
          (Sl**2 - 2*Ss**2)*(8*lam02 + 6*lam21 + 7*lam12*Sl**2 -                &
             2*lam12*Ss**2)))/(6.Q+00*Ss**3*(Sl**2 - 2*Ss**2))
  hd21s22=(Sqrt(2.Q+00)*ck + Ss**5*                                                  &
       (-18*lam21 - (8*lam12 + 3*lam31)*(Sl**2 - 2*Ss**2)))/                    &
    (-3*Sl**2*Ss**5 + 6*Ss**7)
  hd31s22=(5*Sqrt(2.Q+00)*ck + 54*lam31*Ss**7)/(9*Sl**2*Ss**7 - 18*Ss**9)
  hd02s22=(2*Sqrt(2.Q+00)*ck*(Sl**2 - 8*Ss**2) +                                     &
      Ss**3*(72*lam01 - (Sl**2 - 2*Ss**2)**2*                                   &
          (12*lam02 + lam12*Sl**2 - 2*lam12*Ss**2)))/                           &
    (-(Sl**2*Ss) + 2*Ss**3)**3
  hd12s22=(-2*Sqrt(2.Q+00)*ck*(Sl**2 - 4*Ss**2) -                                    &
      12*Ss**5*(-6*lam11 + lam12*(Sl**2 - 2*Ss**2)**2))/                        &
    (Ss**5*(-Sl**2 + 2*Ss**2)**3)
  hd10s55=(4*lam01 + 6*lam20 + lam11*Sl**2 + 4*lam11*Ss**2 +                    &
      (2*(ck + Sqrt(2.Q+00)*lam01*Ss))/Sl +                                          &
      (Sqrt(2.Q+00)*Sl*(lam01 + 3*lam11*Ss**2))/Ss)/6.Q+00
  hd20s55=(-(Sqrt(2.Q+00)*lam01*(Sl**2 - 2*Ss**2)**2) +                              &
      Ss**2*(-4*ck*Ss + 3*Sl**2*                                                &
          (2*Sqrt(2.Q+00)*lam11*Sl**2 +                                              &
            Sl*(8*lam11 + 6*lam30 + lam21*Sl**2)*Ss +                           &
            Sqrt(2.Q+00)*(4*lam11 + 3*lam21*Sl**2)*Ss**2 + 4*lam21*Sl*Ss**3)))/      &
    (18.Q+00*Sl**3*Ss**3)
  hd30s55=  (Sqrt(2.Q+00)*lam01*(Sl**2 - 2*Ss**2)**2*(Sl**2 + 2*Ss**2) -             &
      3*Sqrt(2.Q+00)*lam11*Ss**2*(Sl**3 - 2*Sl*Ss**2)**2 +                           &
      Ss**4*(8*ck*Ss + 3*Sl**4*                                                 &
          (lam31*Sl**3*Ss + 6*Sqrt(2.Q+00)*lam21*Ss**2 +                             &
            3*Sqrt(2.Q+00)*Sl**2*(lam21 + lam31*Ss**2) +                             &
            2*Sl*Ss*(6*lam21 + 3*lam40 + 2*lam31*Ss**2))))/(18.Q+00*Sl**5*Ss**5)
  hd40s55=(-(Sqrt(2.Q+00)*lam01*(Sl**2 - 2*Ss**2)**2*                                &
         (5*Sl**4 + 12*Sl**2*Ss**2 + 20*Ss**4)) +                               &
      2*Ss**2*(6*Sqrt(2.Q+00)*lam11*(Sl**2 + 2*Ss**2)*(Sl**3 - 2*Sl*Ss**2)**2 -      &
         9*Sqrt(2.Q+00)*lam21*Sl**4*(Sl**2*Ss - 2*Ss**3)**2 +                        &
         Ss**4*(-40*ck*Ss + 9*Sl**6*                                            &
             (2*Sqrt(2.Q+00)*lam31*Sl**2 + 8*lam31*Sl*Ss + 3*lam50*Sl*Ss +           &
               4*Sqrt(2.Q+00)*lam31*Ss**2))))/(54.Q+00*Sl**7*Ss**7)
  hd50s55=(5*Sqrt(2.Q+00)*lam01*(Sl**2 - 2*Ss**2)**2*(Sl**2 + 2*Ss**2)*              &
       (7*Sl**4 + 4*Sl**2*Ss**2 + 28*Ss**4) +                                   &
      5*Ss**2*(224*ck*Ss**7 - 18*Sqrt(2.Q+00)*lam31*Sl**6*Ss**4*                     &
          (Sl**2 - 2*Ss**2)**2 +                                                &
         18*Sqrt(2.Q+00)*lam21*Sl**4*(Sl**2 + 2*Ss**2)*                              &
          (Sl**2*Ss - 2*Ss**3)**2 -                                             &
         3*Sqrt(2.Q+00)*lam11*(Sl**3 - 2*Sl*Ss**2)**2*                               &
          (5*Sl**4 + 12*Sl**2*Ss**2 + 20*Ss**4)))/(162.Q+00*Sl**9*Ss**9)
  hd01s55=(6*lam11 + lam02*Sl**2 - (6*Sqrt(2.Q+00)*lam01)/(Sl*Ss) +                  &
      3*Sqrt(2.Q+00)*lam02*Sl*Ss + 4*lam02*Ss**2 +                                   &
      (6*(ck - lam01*(2*Sl + Sqrt(2.Q+00)*Ss)))/(Sl**3 - 2*Sl*Ss**2))/6.Q+00
  hd11s55=(4*lam02 + 6*lam21 + lam12*Sl**2 +                                    &
      (2*Sqrt(2.Q+00)*lam01)/(Sl*Ss**3) +                                            &
      (Sqrt(2.Q+00)*(2*lam01 - 6*lam11*Sl**2 + lam02*Sl**4))/(Sl**3*Ss) +            &
      (Sqrt(2.Q+00)*(2*lam02 + 3*lam12*Sl**2)*Ss)/Sl + 4*lam12*Ss**2 -               &
      (2*(2*ck + 3*lam11*Sl**2*(2*Sl + Sqrt(2.Q+00)*Ss)))/                           &
       (Sl**5 - 2*Sl**3*Ss**2))/6.Q+00
  hd21s55=(6*Sqrt(2.Q+00)*lam01*(-Sl**6 + Sl**4*Ss**2 + 4*Ss**6) +                   &
      Ss**2*(-(Sqrt(2.Q+00)*lam02*Sl**2*(Sl**2 - 2*Ss**2)**3) +                      &
         12*Sqrt(2.Q+00)*lam11*Sl**2*(Sl**2 - 2*Ss**2)*(Sl**2 + Ss**2) +             &
         6*Ss**2*(4*ck*Ss + Sl**4*                                              &
             (-3*lam21*(Sqrt(2.Q+00)*Sl**2 + 2*Sl*Ss - Sqrt(2.Q+00)*Ss**2) +              &
               (Sl**2 - 2*Ss**2)*                                               &
                (Sqrt(2.Q+00)*lam12*Sl**2 + 4*lam12*Sl*Ss + 3*lam31*Sl*Ss +          &
                  2*Sqrt(2.Q+00)*lam12*Ss**2)))))/                                   &
    (18.Q+00*Sl**5*Ss**5*(Sl**2 - 2*Ss**2))
  hd31s55=(2*Sqrt(2.Q+00)*lam01*(Sl**2 - 2*Ss**2)*                                   &
       (5*Sl**6 + 5*Sl**4*Ss**2 + 8*Sl**2*Ss**4 + 20*Ss**6) +                   &
      Ss**2*(-80*ck*Ss**5 - 36*lam31*Sl**7*Ss**5 +                              &
         8*Sqrt(2.Q+00)*Sl**2*Ss**6*(9*lam11 - 2*lam02*Ss**2) +                      &
         Sqrt(2.Q+00)*Sl**10*(lam02 - 3*lam12*Ss**2) +                               &
         4*Sqrt(2.Q+00)*Sl**4*Ss**6*(4*lam02 - 9*lam21 + 6*lam12*Ss**2) +            &
         18*Sqrt(2.Q+00)*Sl**6*Ss**2*                                                &
          (lam11 - lam21*Ss**2 + (-2*lam12 + lam31)*Ss**4) -                    &
         2*Sqrt(2.Q+00)*Sl**8*(9*lam11 +                                             &
            Ss**2*(2*lam02 - 9*lam21 + 9*(-lam12 + lam31)*Ss**2))))/            &
    (18.Q+00*Sl**7*Ss**7*(Sl**2 - 2*Ss**2))
  hd02s55=lam12 + (2*(-((Sqrt(2.Q+00)*lam01)/Ss**3) -                                &
         (Sqrt(2.Q+00)*lam02*Sl**2)/Ss +                                             &
         (6*Sl**2*(-ck + 2*lam01*Sl + Sqrt(2.Q+00)*lam01*Ss))/                       &
          (Sl**2 - 2*Ss**2)**3 +                                                &
         (-ck + 3*Sqrt(2.Q+00)*lam01*Ss)/(Sl**2 - 2*Ss**2)**2 -                 &
         (lam02*Sl**2*(2*Sl + Sqrt(2.Q+00)*Ss))/(Sl**2 - 2*Ss**2)))/Sl**3
  hd12s55=(2*(3*Sqrt(2.Q+00)*lam01*(Sl**2 - 2*Ss**2)**2*(Sl**4 - Ss**4) +       &
        Ss**2*(Sqrt(2.Q+00)*lam02*Sl**2*(Sl**2 - 2*Ss**2)**3*(Sl**2 + Ss**2) +  &
           lam11*(-3*Sqrt(2.Q+00)*Sl**8 + 18*Sqrt(2.Q+00)*Sl**6*Ss**2 +         &
              36*Sl**5*Ss**3 - 9*Sqrt(2.Q+00)*Sl**4*Ss**4 +                     &
              6*Sqrt(2.Q+00)*Sl**2*Ss**6) -                                     &
           3*Ss**2*(-6*ck*Sl**2*Ss + 4*ck*Ss**3 +                               &
              lam12*Sl**4*(Sl**2 - 2*Ss**2)**2*                                 &
               (Sqrt(2.Q+00)*Sl**2 + 2*Sl*Ss - Sqrt(2.Q+00)*Ss**2)))))/                   &
    (3.Q+00*Sl**5*Ss**5*(Sl**2 - 2*Ss**2)**3)
  hd10s99=(2*Sqrt(2.Q+00)*Sl*(-9*ck + 24*lam01*Sl - 36*lam20*Sl +                    &
         6*lam11*Sl**3 + 2*lam02*Sl**5) +                                       &
      (144*ck + 240*lam01*Sl + 612*lam20*Sl +                                   &
         6*(-11*lam11 + 18*lam30)*Sl**3 + 4*(5*lam02 + 9*lam21)*Sl**5 +         &
         3*lam12*Sl**7)*Ss + 4*Sqrt(2.Q+00)*                                         &
       (24*lam01 - 36*lam20 +                                                   &
         Sl**2*(36*lam11 - 54*lam30 +                                           &
            Sl**2*(-2*lam02 + 9*lam21 + 3*lam12*Sl**2)))*Ss**2 +                &
      4*Sl*(123*lam11 + 54*lam30 +                                              &
         Sl**2*(-20*lam02 - 54*lam21 + 3*lam12*Sl**2))*Ss**3 -                  &
      8*Sqrt(2.Q+00)*(6*lam11 + Sl**2*(2*lam02 + 9*lam21 + 6*lam12*Sl**2))*     &
       Ss**4 + 4*Sl*(20*lam02 + 72*lam21 - 21*lam12*Sl**2)*Ss**5 +              &
      16*Sqrt(2.Q+00)*(2*lam02 + 3*lam12*Sl**2)*Ss**6 + 96*lam12*Sl*Ss**7)/     &
    (324.Q+00*Sl*Ss)
  hd20s99=  (720*lam11 + 1350*lam30 + 6*(10*lam12 + 9*lam31)*Sl**4 +            &
      (9*Sqrt(2.Q+00)*ck)/Ss**3 + 18*(49*lam21 + 18*lam40)*Ss**2 +                   &
      48*(5*lam12 + 9*lam31)*Ss**4 -                                            &
      (2*Sqrt(2.Q+00)*Sl**5*(lam02 - 6*lam12*Ss**2))/Ss**3 -                         &
      3*Sl**2*(57*lam21 - 54*lam40 + 4*(20*lam12 + 27*lam31)*Ss**2) +           &
      (2*Sqrt(2.Q+00)*Sl**3*(-3*lam11 + 2*(4*lam02 + 9*lam21)*Ss**2 +                &
           3*(-4*lam12 + 9*lam31)*Ss**4))/Ss**3 +                               &
      (8*Sqrt(2.Q+00)*(12*lam01 - 18*lam20 + 27*(lam11 - 2*lam30)*Ss**2 +            &
           2*(4*lam02 - 9*lam21)*Ss**4 + 12*lam12*Ss**6))/(Sl*Ss) -             &
      (12*Sqrt(2.Q+00)*Sl*(2*lam01 - 3*lam20 + 3*(-5*lam11 + 6*lam30)*Ss**2 +        &
           (4*lam02 - 18*lam21 + 27*lam40)*Ss**4 +                              &
           (4*lam12 + 9*lam31)*Ss**6))/Ss**3 +                                  &
      (-144*ck - 16*Sqrt(2.Q+00)*Ss*                                                 &
          (6*lam01 - 9*lam20 - 3*lam11*Ss**2 + 2*lam02*Ss**4))/Sl**3)/486.Q+00
  hd30s99=(Sqrt(2.Q+00)*Sl**5*(-9*ck + 24*lam01*Sl - 36*lam20*Sl +                   &
         6*lam11*Sl**3 + 2*lam02*Sl**5) -                                       &
      6*Sqrt(2.Q+00)*Sl**4*(8*lam01 - 12*lam20 +                                     &
         2*(8*lam11 - 9*lam30)*Sl**2 + (2*lam02 + 3*lam21)*Sl**4 +              &
         lam12*Sl**6)*Ss**2 + 2*Sqrt(2.Q+00)*Sl**2*                                  &
       (-48*lam01 + 72*lam20 + 72*(2*lam11 - 3*lam30)*Sl**2 +                   &
         2*(4*lam02 + 81*lam21 - 81*lam40)*Sl**4 +                              &
         3*(8*lam12 + 9*lam31)*Sl**6)*Ss**4 +                                   &
      9*(32*ck + 3*Sl**5*(40*lam21 + 66*lam40 - 9*lam31*Sl**2 +                 &
            6*lam50*Sl**2))*Ss**5 +                                             &
      4*Sqrt(2.Q+00)*(48*lam01 - 72*lam20 + 12*(-4*lam11 + 9*lam30)*Sl**2 +          &
         2*(4*lam02 + 27*lam21 - 81*lam40)*Sl**4 -                              &
         9*(4*lam12 - 6*lam31 + 9*lam50)*Sl**6)*Ss**6 +                         &
      54*(19*lam31 + 6*lam50)*Sl**5*Ss**7 -                                     &
      24*Sqrt(2.Q+00)*(4*lam11 + Sl**2*                                              &
          (4*lam02 - 6*lam21 + (-8*lam12 + 9*lam31)*Sl**2))*Ss**8 +             &
      32*Sqrt(2.Q+00)*(2*lam02 - 3*lam12*Sl**2)*Ss**10)/(486.Q+00*Sl**5*Ss**5)
  hd40s99=(2*Sqrt(2.Q+00)*Sl**12*(-5*lam02 + 12*lam12*Ss**2) -                       &
      2*Sqrt(2.Q+00)*Sl**10*(15*lam11 - 4*(7*lam02 + 9*lam21)*Ss**2 +                &
         18*(4*lam12 + 3*lam31)*Ss**4) +                                        &
      32*Sqrt(2.Q+00)*Sl**2*Ss**6*(24*lam01 - 36*lam20 +                             &
         9*(5*lam11 - 12*lam30)*Ss**2 + 4*(7*lam02 - 9*lam21)*Ss**4 +           &
         24*lam12*Ss**6) - 4*Sqrt(2.Q+00)*Sl**8*                                     &
       (30*lam01 - 45*lam20 + 9*(-11*lam11 + 12*lam30)*Ss**2 +                  &
         2*(11*lam02 + 90*lam21 - 81*lam40)*Ss**4 -                             &
         6*(8*lam12 + 63*lam31 - 54*lam50)*Ss**6) +                             &
      9*Sl**7*(5*Sqrt(2.Q+00)*ck + 6*(80*lam31 + 123*lam50)*Ss**7) -                 &
      320*Ss**7*(9*ck + Sqrt(2.Q+00)*Ss*                                             &
          (6*lam01 - 9*lam20 - 3*lam11*Ss**2 + 2*lam02*Ss**4)) -                &
      32*Sqrt(2.Q+00)*Sl**4*Ss**4*(-6*lam01 + 9*lam20 +                              &
         Ss**2*(33*lam11 - 54*lam30 +                                           &
            (11*lam02 + 18*lam21 - 81*lam40)*Ss**2 +                            &
            9*(4*lam12 - 3*lam31)*Ss**4)) +                                     &
      16*Sqrt(2.Q+00)*Sl**6*Ss**2*(12*lam01 - 18*lam20 +                             &
         Ss**2*(-39*lam11 + 54*lam30 +                                          &
            2*(2*lam02 + 54*lam21 - 81*lam40)*Ss**2 +                           &
            3*(8*lam12 + 9*lam31 - 54*lam50)*Ss**4)))/(1458.Q+00*Sl**7*Ss**7)
  hd50s99=(-5*(63*ck*(Sqrt(2.Q+00)*Sl**9 - 128*Ss**9) +                              &
        2*Sqrt(2.Q+00)*(Sl**2 - 2*Ss**2)**2*                                         &
         (-7*Sl**6*(12*lam01 - 18*lam20 + 3*lam11*Sl**2 + lam02*Sl**4) +        &
           Sl**4*(-216*lam01 + 324*lam20 +                                      &
              6*(28*lam11 - 45*lam30)*Sl**2 +                                   &
              5*(2*lam02 + 9*lam21)*Sl**4 + 15*lam12*Sl**6)*Ss**2 -             &
           2*Sl**2*(216*lam01 - 324*lam20 +                                     &
              108*(-2*lam11 + 3*lam30)*Sl**2 +                                  &
              (-4*lam02 + 99*lam21 - 162*lam40)*Sl**4 +                         &
              3*(4*lam12 + 9*lam31)*Sl**6)*Ss**4 -                              &
           4*(168*lam01 - 252*lam20 + 6*(-32*lam11 + 45*lam30)*Sl**2 +          &
              (-4*lam02 + 117*lam21 - 162*lam40)*Sl**4 +                        &
              3*(2*lam12 - 18*lam31 + 27*lam50)*Sl**6)*Ss**6 +                  &
           8*(42*lam11 + 5*(2*lam02 - 9*lam21)*Sl**2 +                          &
              3*(-4*lam12 + 9*lam31)*Sl**4)*Ss**8 +                             &
           16*(-14*lam02 + 15*lam12*Sl**2)*Ss**10)))/(4374.Q+00*Sl**9*Ss**9)
  hd01s99=(246*lam11 + 6*lam12*Sl**4 +                                          &
      (2*Sqrt(2.Q+00)*Sl**3*(-2*lam02 + 3*lam12*Ss**2))/Ss +                         &
      2*Ss**2*(95*lam02 + 18*lam21 + 24*lam12*Ss**2) +                          &
      (6*(-2*lam01*Ss*(2*Sqrt(2.Q+00)*Sl + 23*Ss) +                                  &
           3*ck*(2*Sl + Sqrt(2.Q+00)*Ss) +                                           &
           6*(Sqrt(2.Q+00)*Sl - 2*Ss)*Ss*(lam20 + 2*lam11*Ss**2)))/                  &
       (Ss**2*(Sl**2 - 2*Ss**2)) +                                              &
      Sl**2*(-5*lam02 + 18*(lam21 - 2*lam12*Ss**2)) -                           &
      (12*Sqrt(2.Q+00)*Sl*(lam11 + Ss**2*(-7*lam02 + 3*lam21 + lam12*Ss**2)))/       &
       Ss + (-36*ck + 4*Sqrt(2.Q+00)*Ss*                                             &
          (-6*lam01 + 9*lam20 + 3*lam11*Ss**2 - 2*lam02*Ss**4))/(Sl*Ss**2))/54.Q+00
  hd11s99=  (3*(-13*lam12 + 18*lam31)*Sl**2 +                                   &
      (2*Sqrt(2.Q+00)*Sl**3*(2*lam02 - 3*lam12*Ss**2))/Ss**3 +                       &
      6*(60*lam02 + 147*lam21 + (103*lam12 + 18*lam31)*Ss**2) +                 &
      (12*Sqrt(2.Q+00)*Sl*(lam11 +                                                   &
           Ss**2*(5*lam02 - 6*lam21 + 3*(7*lam12 - 3*lam31)*Ss**2)))/Ss**3      &
+ (72*ck + 8*Sqrt(2.Q+00)*Ss*(6*lam01 - 9*lam20 - 3*lam11*Ss**2 +                    &
            2*lam02*Ss**4))/(Sl**3*Ss**2) +                                     &
      (36*ck - 12*Sqrt(2.Q+00)*Ss*(-4*lam01 + 6*lam20 +                              &
            (13*lam11 - 9*lam30)*Ss**2 + 3*(-4*lam02 + lam21)*Ss**4 +           &
            4*lam12*Ss**6))/(Sl*Ss**4) +                                        &
      (18*(ck*(2*Sl + Sqrt(2.Q+00)*Ss) +                                             &
           2*Ss**3*(lam11*(-2*Sqrt(2.Q+00)*Sl + 31*Ss) -                             &
              3*(Sqrt(2.Q+00)*Sl - 2*Ss)*(lam30 + 2*lam21*Ss**2))))/                 &
       (-(Sl**2*Ss**4) + 2*Ss**6))/162.Q+00
  hd21s99=(9*(40*lam12 + 57*lam31) +                                            &
      (Sqrt(2.Q+00)*Sl**3*(-2*lam02 + 3*lam12*Ss**2))/Ss**5 -                        &
      (6*Sqrt(2.Q+00)*Sl*(lam11 + Ss**2*                                             &
            (lam02 - 3*lam21 + (-11*lam12 + 9*lam31)*Ss**2)))/Ss**5 +           &
      (-72*ck - 8*Sqrt(2.Q+00)*Ss*(6*lam01 - 9*lam20 - 3*lam11*Ss**2 +               &
            2*lam02*Ss**4))/(Sl**5*Ss**2) +                                     &
      (-36*ck + 4*Sqrt(2.Q+00)*Ss*(-6*lam01 + 9*lam20 +                              &
            18*(lam11 - lam30)*Ss**2 - 10*lam02*Ss**4 + 6*lam12*Ss**6))/        &
       (Sl**3*Ss**4) + (-18*ck +                                                &
         2*Sqrt(2.Q+00)*Ss*(-12*lam01 + 18*lam20 +                                   &
            3*(7*lam11 - 12*lam30)*Ss**2 +                                      &
            (22*lam02 - 72*lam21 + 27*lam40)*Ss**4 +                            &
            3*(22*lam12 - 9*lam31)*Ss**6))/(Sl*Ss**6) +                         &
      (9*(ck*(2*Sl + Sqrt(2.Q+00)*Ss) +                                              &
           6*Ss**5*(lam21*(2*Sqrt(2.Q+00)*Sl - 13*Ss) +                              &
              (Sqrt(2.Q+00)*Sl - 2*Ss)*(lam40 + 2*lam31*Ss**2))))/                   &
       (Ss**6*(Sl**2 - 2*Ss**2)))/81.Q+00
  hd31s99=(5*Sqrt(2.Q+00)*Sl**7*                                                     &
       (-9*ck + 24*lam01*Sl - 36*lam20*Sl + 6*lam11*Sl**3 +                     &
         2*lam02*Sl**5)*Ss - Sqrt(2.Q+00)*Sl**6*                                     &
       (120*lam01 - 180*lam20 + 54*(5*lam11 - 6*lam30)*Sl**2 +                  &
         2*(7*lam02 + 36*lam21)*Sl**4 + 15*lam12*Sl**6)*Ss**3 -                 &
      2*Sqrt(2.Q+00)*Sl**4*(24*lam01 - 36*lam20 +                                    &
         6*(-16*lam11 + 27*lam30)*Sl**2 +                                       &
         (34*lam02 - 171*lam21 + 162*lam40)*Sl**4 +                             &
         18*(lam12 - 3*lam31)*Sl**6)*Ss**5 +                                    &
      2*Sqrt(2.Q+00)*Sl**2*(48*lam01 - 72*lam20 - 48*lam11*Sl**2 +                   &
         Sl**4*(-8*lam02 + 54*(lam21 + 3*lam40) +                               &
            3*(88*lam12 - 63*lam31 + 54*lam50)*Sl**2))*Ss**7 -                  &
      18*(80*ck + 3*(47*lam31 + 6*lam50)*Sl**7)*Ss**8 -                         &
      4*Sqrt(2.Q+00)*(240*lam01 - 360*lam20 +                                        &
         108*(-2*lam11 + 3*lam30)*Sl**2 -                                       &
         2*(56*lam02 - 117*lam21 + 81*lam40)*Sl**4 +                            &
         3*(92*lam12 - 126*lam31 + 27*lam50)*Sl**6)*Ss**9 +                     &
      8*Sqrt(2.Q+00)*(60*lam11 - 2*(14*lam02 + 9*lam21)*Sl**2 +                      &
         27*(2*lam12 - lam31)*Sl**4)*Ss**11 +                                   &
      64*Sqrt(2.Q+00)*(-5*lam02 + 6*lam12*Sl**2)*Ss**13)/                            &
    (243.Q+00*Ss**8*(Sl**9 - 2*Sl**7*Ss**2))
  hd02s99=(219*lam12 - (4*Sqrt(2.Q+00)*Sl*(lam02 + 3*lam12*Ss**2))/Ss**3 +           &
      (18*(-2*Sqrt(2.Q+00)*lam01*Sl*Ss + 3*Sqrt(2.Q+00)*lam20*Sl*Ss -                     &
           2*lam11*(Sqrt(2.Q+00)*Sl - 8*Ss)*Ss**3 + ck*(5*Sl + Sqrt(2.Q+00)*Ss)))/        &
       (Ss**4*(Sl**2 - 2*Ss**2)**2) +                                           &
      (36*(-2*lam01*Ss*(2*Sqrt(2.Q+00)*Sl + 23*Ss) +                                 &
           3*ck*(2*Sl + Sqrt(2.Q+00)*Ss) +                                           &
           6*(Sqrt(2.Q+00)*Sl - 2*Ss)*Ss*(lam20 + 2*lam11*Ss**2)))/                  &
       (Ss**2*(-Sl**2 + 2*Ss**2)**3) +                                          &
      (-18*ck + 2*Sqrt(2.Q+00)*Ss*(-6*lam01 + 9*lam20 + 3*lam11*Ss**2 -              &
            2*lam02*Ss**4))/(Sl**3*Ss**4) +                                     &
      (3*(12*ck*Sl + 15*Sqrt(2.Q+00)*lam11*Sl*Ss**3 +                                &
           4*Ss**5*(lam02*(5*Sqrt(2.Q+00)*Sl + 44*Ss) -                              &
              3*(Sqrt(2.Q+00)*Sl - 2*Ss)*(lam21 + 2*lam12*Ss**2))))/                 &
       (-(Sl**2*Ss**6) + 2*Ss**8) +                                             &
      (36*ck + Sqrt(2.Q+00)*Ss**3*(33*lam11 +                                        &
            4*Ss**2*(-17*lam02 + 9*lam21 + 3*lam12*Ss**2)))/(Sl*Ss**6))/27.Q+00
  hd12s99=(2*(-2*Sqrt(2.Q+00)*lam02*Sl**12 +                                         &
        12*(26*lam12 + 3*lam31)*Sl**9*Ss**5 -                                   &
        4*Sqrt(2.Q+00)*Sl**4*Ss**4*                                                  &
         (18*lam01 - 27*lam20 + 3*(-19*lam11 + 9*lam30)*Ss**2 +                 &
           4*(5*lam02 - 9*lam21)*Ss**4 + 6*(13*lam12 - 6*lam31)*Ss**6) +        &
        2*Sqrt(2.Q+00)*Sl**8*(-12*lam01 + 18*lam20 +                                 &
           18*(lam11 - lam30)*Ss**2 + (37*lam02 - 18*lam21)*Ss**4 +             &
           9*(-17*lam12 + 10*lam31)*Ss**6) -                                    &
        12*Sl**5*Ss**2*(3*Sqrt(2.Q+00)*ck + 3*(31*lam11 + 6*lam30)*Ss**3 +           &
           12*lam21*Ss**5 - 4*(26*lam12 + 3*lam31)*Ss**7) +                     &
        16*Ss**7*(9*ck + Sqrt(2.Q+00)*Ss*                                            &
            (6*lam01 - 9*lam20 - 3*lam11*Ss**2 + 2*lam02*Ss**4)) -              &
        8*Sl**2*Ss**5*(27*ck +                                                  &
           Sqrt(2.Q+00)*Ss*(12*lam01 - 18*lam20 + 9*(lam11 - lam30)*Ss**2 +          &
              (-16*lam02 + 9*lam21)*Ss**4 + 6*lam12*Ss**6)) +                   &
        3*Sl**7*(3*Sqrt(2.Q+00)*ck -                                                 &
           16*(3*lam21*Ss**5 + (26*lam12 + 3*lam31)*Ss**7)) -                   &
        2*Sqrt(2.Q+00)*Sl**10*(3*lam11 +                                             &
           Ss**2*(4*lam02 - 9*(lam21 + (3*lam12 - 2*lam31)*Ss**2))) +           &
        6*Sqrt(2.Q+00)*Sl**6*Ss**2*(16*lam01 - 24*lam20 +                            &
           Ss**2*(-11*lam11 + 36*lam30 +                                        &
              Ss**2*(-16*lam02 + 27*lam21 + 94*lam12*Ss**2 -                    &
                 48*lam31*Ss**2)))))/(27.Q+00*Sl**5*Ss**5*(-Sl**2 + 2*Ss**2)**3)
  hd10s19=(3*Sqrt(2.Q+00)*lam12*Sl**7*Ss + 2*Sl**6*(lam02 + 3*lam12*Ss**2) -         &
      4*Sqrt(2.Q+00)*Sl**5*Ss*(lam02 - 9*lam21 + 6*lam12*Ss**2) +                    &
      2*Sqrt(2.Q+00)*Sl**3*Ss*(93*lam11 + 54*lam30 + 8*lam02*Ss**2 +                 &
         30*lam12*Ss**4) + 4*Sl**2*                                             &
       (6*lam01 - 9*lam20 + 9*(2*lam11 - 3*lam30)*Ss**2 -                       &
         (2*lam02 + 9*lam21)*Ss**4 + 6*lam12*Ss**6) +                           &
      2*Sl**4*(3*lam11 + Ss**2*(-2*lam02 + 9*lam21 - 12*lam12*Ss**2)) +         &
      4*Ss*(9*Sqrt(2.Q+00)*ck + 2*Ss*                                                &
          (6*lam01 - 9*lam20 - 3*lam11*Ss**2 + 2*lam02*Ss**4)) -                &
      4*Sl*(9*ck + Sqrt(2.Q+00)*Ss*(12*lam01 - 18*lam20 +                            &
            3*(37*lam11 + 9*lam30)*Ss**2 + 4*(lam02 + 9*lam21)*Ss**4 +          &
            12*lam12*Ss**6)))/(324.Q+00*Sl*Ss)
  hd20s19=(-(Sl**3*(-18*ck + 12*lam01*Sl - 18*lam20*Sl +                        &
           3*lam11*Sl**3 + lam02*Sl**5)) +                                      &
      2*Sl**2*(24*lam01 - 36*lam20 + 9*(5*lam11 - 6*lam30)*Sl**2 +              &
         (4*lam02 + 9*lam21)*Sl**4 + 3*lam12*Sl**6)*Ss**2 +                     &
      3*Sqrt(2.Q+00)*(-12*ck + Sl**3*                                                &
          (-48*lam11 + 72*lam30 + 3*(47*lam21 + 18*lam40)*Sl**2 +               &
            2*(-2*lam12 + 9*lam31)*Sl**4))*Ss**3 -                              &
      3*(16*lam01 - 24*lam20 - 36*(lam11 - 2*lam30)*Sl**2 +                     &
         2*(4*lam02 - 18*lam21 + 27*lam40)*Sl**4 +                              &
         (4*lam12 - 9*lam31)*Sl**6)*Ss**4 +                                     &
      6*Sqrt(2.Q+00)*Sl**3*(-159*lam21 - 27*lam40 + 8*lam12*Sl**2)*Ss**5 +           &
      2*(12*lam11 + 4*(4*lam02 - 9*lam21)*Sl**2 -                               &
         3*(4*lam12 + 9*lam31)*Sl**4)*Ss**6 -                                   &
      24*Sqrt(2.Q+00)*(2*lam12 + 9*lam31)*Sl**3*Ss**7 -                              &
      16*(lam02 - 3*lam12*Sl**2)*Ss**8)/(486.Q+00*Sl**3*Ss**3)
  hd30s19=(Sl**5*(-18*ck + 12*lam01*Sl - 18*lam20*Sl + 3*lam11*Sl**3 +          &
         lam02*Sl**5) - 3*Sl**4*                                                &
       (8*lam01 - 12*lam20 + 2*(8*lam11 - 9*lam30)*Sl**2 +                      &
         (2*lam02 + 3*lam21)*Sl**4 + lam12*Sl**6)*Ss**2 +                       &
      Sl**2*(-48*lam01 + 72*lam20 + 72*(2*lam11 - 3*lam30)*Sl**2 +              &
         2*(4*lam02 + 81*lam21 - 81*lam40)*Sl**4 +                              &
         3*(8*lam12 + 9*lam31)*Sl**6)*Ss**4 +                                   &
      9*Sqrt(2.Q+00)*(8*ck + 3*Sl**5*                                                &
          (-8*lam21 + 12*lam40 + 3*(7*lam31 + 2*lam50)*Sl**2))*Ss**5 +          &
      2*(48*lam01 - 72*lam20 + 12*(-4*lam11 + 9*lam30)*Sl**2 +                  &
         2*(4*lam02 + 27*lam21 - 81*lam40)*Sl**4 -                              &
         9*(4*lam12 - 6*lam31 + 9*lam50)*Sl**6)*Ss**6 -                         &
      54*Sqrt(2.Q+00)*(23*lam31 + 3*lam50)*Sl**5*Ss**7 -                             &
      12*(4*lam11 + Sl**2*(4*lam02 - 6*lam21 +                                  &
            (-8*lam12 + 9*lam31)*Sl**2))*Ss**8 +                                &
      16*(2*lam02 - 3*lam12*Sl**2)*Ss**10)/(486.Q+00*Sl**5*Ss**5)
  hd40s19=(-5*Sl**7*(-18*ck + 12*lam01*Sl - 18*lam20*Sl +                       &
         3*lam11*Sl**3 + lam02*Sl**5) +                                         &
      2*Sl**6*(48*lam01 - 72*lam20 + 9*(11*lam11 - 12*lam30)*Sl**2 +            &
         2*(7*lam02 + 9*lam21)*Sl**4 + 6*lam12*Sl**6)*Ss**2 -                   &
      2*Sl**4*(-48*lam01 + 72*lam20 + 12*(13*lam11 - 18*lam30)*Sl**2 +          &
         2*(11*lam02 + 90*lam21 - 81*lam40)*Sl**4 +                             &
         9*(4*lam12 + 3*lam31)*Sl**6)*Ss**4 +                                   &
      4*Sl**2*(96*lam01 - 144*lam20 + 12*(-11*lam11 + 18*lam30)*Sl**2 +         &
         4*(2*lam02 + 54*lam21 - 81*lam40)*Sl**4 +                              &
         3*(8*lam12 + 63*lam31 - 54*lam50)*Sl**6)*Ss**6 -                       &
      144*Sqrt(2.Q+00)*(5*ck + 3*(2*lam31 - 3*lam50)*Sl**7)*Ss**7 -             &
      8*(120*lam01 - 180*lam20 + 18*(-5*lam11 + 12*lam30)*Sl**2 +               &
         2*(11*lam02 + 18*lam21 - 81*lam40)*Sl**4 -                             &
         3*(8*lam12 + 9*lam31 - 54*lam50)*Sl**6)*Ss**8 +                        &
      16*(30*lam11 + 4*(7*lam02 - 9*lam21)*Sl**2 +                              &
         9*(-4*lam12 + 3*lam31)*Sl**4)*Ss**10 +                                 &
      64*(-5*lam02 + 6*lam12*Sl**2)*Ss**12)/(1458.Q+00*Sl**7*Ss**7)
  hd50s19=(5*(-126*ck*(Sl**9 - 16*Sqrt(2.Q+00)*Ss**9) +                              &
        (Sl**2 - 2*Ss**2)**2*(7*Sl**6*                                          &
            (12*lam01 - 18*lam20 + 3*lam11*Sl**2 + lam02*Sl**4) -               &
           Sl**4*(-216*lam01 + 324*lam20 +                                      &
              6*(28*lam11 - 45*lam30)*Sl**2 +                                   &
              5*(2*lam02 + 9*lam21)*Sl**4 + 15*lam12*Sl**6)*Ss**2 +             &
           2*Sl**2*(216*lam01 - 324*lam20 +                                     &
              108*(-2*lam11 + 3*lam30)*Sl**2 +                                  &
              (-4*lam02 + 99*lam21 - 162*lam40)*Sl**4 +                         &
              3*(4*lam12 + 9*lam31)*Sl**6)*Ss**4 +                              &
           4*(168*lam01 - 252*lam20 + 6*(-32*lam11 + 45*lam30)*Sl**2 +          &
              (-4*lam02 + 117*lam21 - 162*lam40)*Sl**4 +                        &
              3*(2*lam12 - 18*lam31 + 27*lam50)*Sl**6)*Ss**6 -                  &
           8*(42*lam11 + 5*(2*lam02 - 9*lam21)*Sl**2 +                          &
              3*(-4*lam12 + 9*lam31)*Sl**4)*Ss**8 +                             &
           16*(14*lam02 - 15*lam12*Sl**2)*Ss**10)))/(4374.Q+00*Sl**9*Ss**9)
  hd01s19=(Sqrt(2.Q+00)*(-96*lam11 +                                                 &
         Sl**2*(37*lam02 + 18*lam21 + 6*lam12*Sl**2)) +                         &
      (36*ck - 24*lam01*Sl + 36*lam20*Sl - 6*lam11*Sl**3 - 2*lam02*Sl**5)/      &
       (Sl**2*Ss) + (3*(2*lam11 +                                               &
           Sl**2*(14*lam02 - 6*lam21 + lam12*Sl**2))*Ss)/Sl -                   &
      2*Sqrt(2.Q+00)*(55*lam02 + 9*lam21)*Ss**2 -                                    &
      (2*(2*lam02 + 3*lam12*Sl**2)*Ss**3)/Sl - 24*Sqrt(2.Q+00)*lam12*Ss**4 +         &
      (18*ck*(Sqrt(2.Q+00)*Sl + 4*Ss) +                                              &
         12*Sl*(lam01*(19*Sqrt(2.Q+00)*Sl - 2*Ss) +                                  &
            3*(lam20 + lam11*Sl**2)*(4*Sqrt(2.Q+00)*Sl + Ss)))/                      &
       (Sl**4 - 2*Sl**2*Ss**2))/54.Q+00
  hd11s19=(3*Sqrt(2.Q+00)*(-24*lam02 - 84*lam21 +                                    &
         (53*lam12 + 18*lam31)*Sl**2) +                                         &
      (2*(-18*ck + 12*lam01*Sl - 18*lam20*Sl + 3*lam11*Sl**3 +                  &
           lam02*Sl**5))/(Sl**2*Ss**3) -                                        &
      (3*(24*ck - 8*lam01*Sl + 12*lam20*Sl +                                    &
           2*(7*lam11 - 18*lam30)*Sl**3 + 2*(-5*lam02 + 6*lam21)*Sl**5 +        &
           lam12*Sl**7))/(Sl**4*Ss) -                                           &
      (6*(2*lam11 + 3*Sl**2*(-4*lam02 + lam21 +                                 &
              (-7*lam12 + 3*lam31)*Sl**2))*Ss)/Sl**3 -                          &
      6*Sqrt(2.Q+00)*(71*lam12 + 9*lam31)*Ss**2 +                                    &
      (8*(lam02 - 3*lam12*Sl**2)*Ss**3)/Sl**3 +                                 &
      (-36*ck*(Sqrt(2.Q+00)*Sl + 4*Ss) +                                             &
         36*Sl**3*(3*(lam30 + lam21*Sl**2)*(4*Sqrt(2.Q+00)*Sl + Ss) +           &
            lam11*(35*Sqrt(2.Q+00)*Sl + 2*Ss)))/(Sl**6 - 2*Sl**4*Ss**2))/162.Q+00
  hd21s19=(-2*Sl**5*(-18*ck + 12*lam01*Sl - 18*lam20*Sl +                       &
         3*lam11*Sl**3 + lam02*Sl**5) +                                         &
      Sl**4*(24*lam01 - 36*lam20 + 18*(3*lam11 - 4*lam30)*Sl**2 -               &
         2*(lam02 - 9*lam21)*Sl**4 + 3*lam12*Sl**6)*Ss**2 +                     &
      2*Sl**4*(-6*lam11 + 36*lam30 +                                            &
         Sl**2*(28*lam02 - 36*lam21 + 54*lam40 +                                &
            3*(10*lam12 - 9*lam31)*Sl**2))*Ss**4 +                              &
      36*Sqrt(2.Q+00)*(2*ck + Sl**5*                                                 &
          (51*lam21 + 12*lam40 - 4*lam12*Sl**2 + 6*lam31*Sl**2))*Ss**5 +        &
      2*(48*lam01 - 72*lam20 + 12*(-5*lam11 + 6*lam30)*Sl**2 -                  &
         2*(32*lam02 - 72*lam21 + 27*lam40)*Sl**4 + 81*lam31*Sl**6)*Ss**6       &
+ 144*Sqrt(2.Q+00)*(2*lam12 + 3*lam31)*Sl**5*Ss**7 +                                 &
      4*(-12*lam11 + 16*lam02*Sl**2 + 3*(-20*lam12 + 9*lam31)*Sl**4)*           &
       Ss**8 + 16*(2*lam02 - 3*lam12*Sl**2)*Ss**10)/                            &
    (162.Q+00*Sl**5*Ss**5*(Sl**2 - 2*Ss**2))
  hd31s19=(10*Sl**7*(-18*ck + 12*lam01*Sl - 18*lam20*Sl +                       &
         3*lam11*Sl**3 + lam02*Sl**5) -                                         &
      Sl**6*(120*lam01 - 180*lam20 + 54*(5*lam11 - 6*lam30)*Sl**2 +             &
         2*(7*lam02 + 36*lam21)*Sl**4 + 15*lam12*Sl**6)*Ss**2 -                 &
      2*Sl**4*(24*lam01 - 36*lam20 + 6*(-16*lam11 + 27*lam30)*Sl**2 +           &
         (34*lam02 - 171*lam21 + 162*lam40)*Sl**4 +                             &
         18*(lam12 - 3*lam31)*Sl**6)*Ss**4 +                                    &
      2*Sl**2*(48*lam01 - 72*lam20 - 48*lam11*Sl**2 +                           &
         Sl**4*(-8*lam02 + 54*(lam21 + 3*lam40) +                               &
            3*(88*lam12 - 63*lam31 + 54*lam50)*Sl**2))*Ss**6 +                  &
      36*Sqrt(2.Q+00)*(-20*ck + 3*(67*lam31 + 12*lam50)*Sl**7)*Ss**7 -               &
      4*(240*lam01 - 360*lam20 + 108*(-2*lam11 + 3*lam30)*Sl**2 -               &
         2*(56*lam02 - 117*lam21 + 81*lam40)*Sl**4 +                            &
         3*(92*lam12 - 126*lam31 + 27*lam50)*Sl**6)*Ss**8 +                     &
      8*(60*lam11 - 2*(14*lam02 + 9*lam21)*Sl**2 +                              &
         27*(2*lam12 - lam31)*Sl**4)*Ss**10 +                                   &
      64*(-5*lam02 + 6*lam12*Sl**2)*Ss**12)/                                    &
    (486.Q+00*Sl**7*Ss**7*(Sl**2 - 2*Ss**2))
  hd02s19=(2*(-48*Sqrt(2.Q+00)*lam12 +                                               &
        (18*ck - 12*lam01*Sl + 18*lam20*Sl - 3*lam11*Sl**3 - lam02*Sl**5)/      &
         (Sl**4*Ss**3) - (36*ck + 30*lam11*Sl**3 +                              &
           Sl**5*(32*lam02 - 18*lam21 + 3*lam12*Sl**2))/(Sl**6*Ss) +            &
        ((-lam02 + 3*lam12*Sl**2)*Ss)/Sl**3 +                                   &
        (72*Sqrt(2.Q+00)*Sl**6*(lam21 + lam12*Sl**2) +                               &
           6*lam02*Sl**5*(34*Sqrt(2.Q+00)*Sl - 5*Ss) -                               &
           9*(8*ck + 7*lam11*Sl**3 - 2*Sl**5*(lam21 + lam12*Sl**2))*Ss)/        &
         (Sl**8 - 2*Sl**6*Ss**2) -                                              &
        (9*(ck*(Sqrt(2.Q+00)*Sl + 16*Ss) +                                           &
             2*Sl*(-2*lam01*Ss + 3*lam20*Ss +                                   &
                lam11*Sl**2*(-8*Sqrt(2.Q+00)*Sl + Ss))))/                            &
         (Sl**4*(Sl**2 - 2*Ss**2)**2) +                                         &
        (-54*ck*(Sqrt(2.Q+00)*Sl + 4*Ss) -                                           &
           36*Sl*(lam01*(19*Sqrt(2.Q+00)*Sl - 2*Ss) +                                &
              3*(lam20 + lam11*Sl**2)*(4*Sqrt(2.Q+00)*Sl + Ss)))/                    &
         (Sl**2*(Sl**2 - 2*Ss**2)**3)))/27.Q+00
  hd12s19=(2*((Sl**3*(-18*ck + 12*lam01*Sl - 18*lam20*Sl +                      &
             3*lam11*Sl**3 + lam02*Sl**5))/Ss**5 +                              &
        (Sl*(-36*ck + Sl*(24*lam01 - 36*lam20 + 18*lam30*Sl**2 +                &
                (10*lam02 - 9*lam21)*Sl**4)))/Ss**3 +                           &
        (36*lam01 - 54*lam20 - 3*lam11*Sl**2 +                                  &
           (11*lam02 - 36*lam21)*Sl**4 + 9*(-3*lam12 + 2*lam31)*Sl**6)/Ss +     &
        Sl**2*(2*lam02 - 3*lam12*Sl**2)*Ss +                                    &
        (12*Sqrt(2.Q+00)*(25*lam12 + 6*lam31)*Sl**7 +                                &
           36*(2*lam01 - 3*lam20)*Ss - 63*lam21*Sl**4*Ss -                      &
           6*(lam12 - 3*lam31)*Sl**6*Ss)/(Sl**2 - 2*Ss**2) +                    &
        (18*Sl*(ck*(Sqrt(2.Q+00)*Sl + 8*Ss) +                                        &
             Sl**3*(lam21*Sl**2*(8*Sqrt(2.Q+00)*Sl - Ss) -                           &
                (2*lam11 + 3*lam30)*Ss)))/(Sl**2 - 2*Ss**2)**2 -                &
        (36*Sl**3*(-(ck*(Sqrt(2.Q+00)*Sl + 4*Ss)) +                                  &
             Sl**3*(3*(lam30 + lam21*Sl**2)*(4*Sqrt(2.Q+00)*Sl + Ss) +               &
                lam11*(35*Sqrt(2.Q+00)*Sl + 2*Ss))))/(Sl**2 - 2*Ss**2)**3))/         &
    (27.Q+00*Sl**7)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  A00=Sqrt(hd00s11**2 + 4*hd00s19**2 - 2*hd00s11*hd00s99 + hd00s99**2)
  A10=(hd00s11*hd10s11 - hd00s99*hd10s11 + 4*hd00s19*hd10s19 -                  &
      hd00s11*hd10s99 + hd00s99*hd10s99)/A00
  A20=(-A10**2 + hd10s11**2 + 4*hd10s19**2 - 2*hd10s11*hd10s99 +                &
      hd10s99**2 + hd00s11*hd20s11 - hd00s99*hd20s11 + 4*hd00s19*hd20s19 -      &
      hd00s11*hd20s99 + hd00s99*hd20s99)/A00
  A30=(-3*A10*A20 + 3*hd10s11*hd20s11 - 3*hd10s99*hd20s11 +                     &
      12*hd10s19*hd20s19 - 3*hd10s11*hd20s99 + 3*hd10s99*hd20s99 +              &
      hd00s11*hd30s11 - hd00s99*hd30s11 + 4*hd00s19*hd30s19 -                   &
      hd00s11*hd30s99 + hd00s99*hd30s99)/A00
  A40=(-3*A20**2 - 4*A10*A30 + 3*hd20s11**2 + 12*hd20s19**2 -                   &
      6*hd20s11*hd20s99 + 3*hd20s99**2 + 4*hd10s11*hd30s11 -                    &
      4*hd10s99*hd30s11 + 16*hd10s19*hd30s19 - 4*hd10s11*hd30s99 +              &
      4*hd10s99*hd30s99 + hd00s11*hd40s11 - hd00s99*hd40s11 +                   &
      4*hd00s19*hd40s19 - hd00s11*hd40s99 + hd00s99*hd40s99)/A00
  A50=(-10*A20*A30 - 5*A10*A40 + 10*hd20s11*hd30s11 -                           &
      10*hd20s99*hd30s11 + 40*hd20s19*hd30s19 - 10*hd20s11*hd30s99 +            &
      10*hd20s99*hd30s99 + 5*hd10s11*hd40s11 - 5*hd10s99*hd40s11 +              &
      20*hd10s19*hd40s19 - 5*hd10s11*hd40s99 + 5*hd10s99*hd40s99 +              &
      hd00s11*hd50s11 - hd00s99*hd50s11 + 4*hd00s19*hd50s19 -                   &
      hd00s11*hd50s99 + hd00s99*hd50s99)/A00
  A01=(hd00s11*hd01s11 - hd00s99*hd01s11 + 4*hd00s19*hd01s19 -                  &
      hd00s11*hd01s99 + hd00s99*hd01s99)/A00
  A11=(-(A01*A10) + hd01s11*hd10s11 - hd01s99*hd10s11 +                         &
      4*hd01s19*hd10s19 - hd01s11*hd10s99 + hd01s99*hd10s99 +                   &
      hd00s11*hd11s11 - hd00s99*hd11s11 + 4*hd00s19*hd11s19 -                   &
      hd00s11*hd11s99 + hd00s99*hd11s99)/A00
  A21=(-2*A10*A11 - A01*A20 + 2*hd10s11*hd11s11 - 2*hd10s99*hd11s11 +           &
      8*hd10s19*hd11s19 - 2*hd10s11*hd11s99 + 2*hd10s99*hd11s99 +               &
      hd01s11*hd20s11 - hd01s99*hd20s11 + 4*hd01s19*hd20s19 -                   &
      hd01s11*hd20s99 + hd01s99*hd20s99 + hd00s11*hd21s11 -                     &
      hd00s99*hd21s11 + 4*hd00s19*hd21s19 - hd00s11*hd21s99 +                   &
      hd00s99*hd21s99)/A00
  A31= (-3*A11*A20 - 3*A10*A21 - A01*A30 + 3*hd11s11*hd20s11 -                  &
      3*hd11s99*hd20s11 + 12*hd11s19*hd20s19 - 3*hd11s11*hd20s99 +              &
      3*hd11s99*hd20s99 + 3*hd10s11*hd21s11 - 3*hd10s99*hd21s11 +               &
      12*hd10s19*hd21s19 - 3*hd10s11*hd21s99 + 3*hd10s99*hd21s99 +              &
      hd01s11*hd30s11 - hd01s99*hd30s11 + 4*hd01s19*hd30s19 -                   &
      hd01s11*hd30s99 + hd01s99*hd30s99 + hd00s11*hd31s11 -                     &
      hd00s99*hd31s11 + 4*hd00s19*hd31s19 - hd00s11*hd31s99 +                   &
      hd00s99*hd31s99)/A00
  A02=(-A01**2 + hd01s11**2 + 4*hd01s19**2 - 2*hd01s11*hd01s99 + hd01s99**2 +   &
      hd00s11*hd02s11 - hd00s99*hd02s11 + 4*hd00s19*hd02s19 -                   &
      hd00s11*hd02s99 + hd00s99*hd02s99)/A00
  A12=(-(A02*A10) - 2*A01*A11 + hd02s11*hd10s11 - hd02s99*hd10s11 +             &
      4*hd02s19*hd10s19 - hd02s11*hd10s99 + hd02s99*hd10s99 +                   &
      2*hd01s11*hd11s11 - 2*hd01s99*hd11s11 + 8*hd01s19*hd11s19 -               &
      2*hd01s11*hd11s99 + 2*hd01s99*hd11s99 + hd00s11*hd12s11 -                 &
      hd00s99*hd12s11 + 4*hd00s19*hd12s19 - hd00s11*hd12s99 +                   &
      hd00s99*hd12s99)/A00 
  B00=Sqrt(hd00p11**2 + 4*hd00p19**2 - 2*hd00p11*hd00p99 + hd00p99**2)
  B10=(hd00p11*hd10p11 - hd00p99*hd10p11 + 4*hd00p19*hd10p19 -                  &
      hd00p11*hd10p99 + hd00p99*hd10p99)/B00
  B20=(-B10**2 + hd10p11**2 + 4*hd10p19**2 - 2*hd10p11*hd10p99 +                &
      hd10p99**2 + hd00p11*hd20p11 - hd00p99*hd20p11 + 4*hd00p19*hd20p19 -      &
      hd00p11*hd20p99 + hd00p99*hd20p99)/B00
  B30=(-3*B10*B20 + 3*hd10p11*hd20p11 - 3*hd10p99*hd20p11 +                     &
      12*hd10p19*hd20p19 - 3*hd10p11*hd20p99 + 3*hd10p99*hd20p99 +              &
      hd00p11*hd30p11 - hd00p99*hd30p11 + 4*hd00p19*hd30p19 -                   &
      hd00p11*hd30p99 + hd00p99*hd30p99)/B00
  B40=(-3*B20**2 - 4*B10*B30 + 3*hd20p11**2 + 12*hd20p19**2 -                   &
      6*hd20p11*hd20p99 + 3*hd20p99**2 + 4*hd10p11*hd30p11 -                    &
      4*hd10p99*hd30p11 + 16*hd10p19*hd30p19 - 4*hd10p11*hd30p99 +              &
      4*hd10p99*hd30p99 + hd00p11*hd40p11 - hd00p99*hd40p11 +                   &
      4*hd00p19*hd40p19 - hd00p11*hd40p99 + hd00p99*hd40p99)/B00
  B50=(-10*B20*B30 - 5*B10*B40 + 10*hd20p11*hd30p11 -                           &
      10*hd20p99*hd30p11 + 40*hd20p19*hd30p19 - 10*hd20p11*hd30p99 +            &
      10*hd20p99*hd30p99 + 5*hd10p11*hd40p11 - 5*hd10p99*hd40p11 +              &
      20*hd10p19*hd40p19 - 5*hd10p11*hd40p99 + 5*hd10p99*hd40p99 +              &
      hd00p11*hd50p11 - hd00p99*hd50p11 + 4*hd00p19*hd50p19 -                   &
      hd00p11*hd50p99 + hd00p99*hd50p99)/B00
  B01=(hd00p11*hd01p11 - hd00p99*hd01p11 + 4*hd00p19*hd01p19 -                  &
      hd00p11*hd01p99 + hd00p99*hd01p99)/B00
  B11=(-(B01*B10) + hd01p11*hd10p11 - hd01p99*hd10p11 +                         &
      4*hd01p19*hd10p19 - hd01p11*hd10p99 + hd01p99*hd10p99 +                   &
      hd00p11*hd11p11 - hd00p99*hd11p11 + 4*hd00p19*hd11p19 -                   &
      hd00p11*hd11p99 + hd00p99*hd11p99)/B00
  B21=(-2*B10*B11 - B01*B20 + 2*hd10p11*hd11p11 - 2*hd10p99*hd11p11 +           &
      8*hd10p19*hd11p19 - 2*hd10p11*hd11p99 + 2*hd10p99*hd11p99 +               &
      hd01p11*hd20p11 - hd01p99*hd20p11 + 4*hd01p19*hd20p19 -                   &
      hd01p11*hd20p99 + hd01p99*hd20p99 + hd00p11*hd21p11 -                     &
      hd00p99*hd21p11 + 4*hd00p19*hd21p19 - hd00p11*hd21p99 +                   &
      hd00p99*hd21p99)/B00
  B31= (-3*B11*B20 - 3*B10*B21 - B01*B30 + 3*hd11p11*hd20p11 -                  &
      3*hd11p99*hd20p11 + 12*hd11p19*hd20p19 - 3*hd11p11*hd20p99 +              &
      3*hd11p99*hd20p99 + 3*hd10p11*hd21p11 - 3*hd10p99*hd21p11 +               &
      12*hd10p19*hd21p19 - 3*hd10p11*hd21p99 + 3*hd10p99*hd21p99 +              &
      hd01p11*hd30p11 - hd01p99*hd30p11 + 4*hd01p19*hd30p19 -                   &
      hd01p11*hd30p99 + hd01p99*hd30p99 + hd00p11*hd31p11 -                     &
      hd00p99*hd31p11 + 4*hd00p19*hd31p19 - hd00p11*hd31p99 +                   &
      hd00p99*hd31p99)/B00
  B02=(-B01**2 + hd01p11**2 + 4*hd01p19**2 - 2*hd01p11*hd01p99 + hd01p99**2 +   &
      hd00p11*hd02p11 - hd00p99*hd02p11 + 4*hd00p19*hd02p19 -                   &
      hd00p11*hd02p99 + hd00p99*hd02p99)/B00
  B12=(-(B02*B10) - 2*B01*B11 + hd02p11*hd10p11 - hd02p99*hd10p11 +             &
      4*hd02p19*hd10p19 - hd02p11*hd10p99 + hd02p99*hd10p99 +                   &
      2*hd01p11*hd11p11 - 2*hd01p99*hd11p11 + 8*hd01p19*hd11p19 -               &
      2*hd01p11*hd11p99 + 2*hd01p99*hd11p99 + hd00p11*hd12p11 -                 &
      hd00p99*hd12p11 + 4*hd00p19*hd12p19 - hd00p11*hd12p99 +                   &
      hd00p99*hd12p99)/B00 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  mboson2(1)=(a00 + hd00s11 + hd00s99)/2.Q+00
  mboson2(2:4)=hd00s22
  mboson2(5:8)=hd00s55
  mboson2(9)=(-a00 + hd00s11 + hd00s99)/2.Q+00
  mboson2(10)=(-b00 + hd00p11 + hd00p99)/2.Q+00
  mboson2(11:13)=hd00p22
  mboson2(14:17)=hd00p55
  mboson2(18)=(b00 + hd00p11 + hd00p99)/2.Q+00
  mbo2d10rho(1)=(a10 + hd10s11 + hd10s99)/2.Q+00
  mbo2d10rho(2:4)=hd10s22
  mbo2d10rho(5:8)=hd10s55
  mbo2d10rho(9)=(-a10 + hd10s11 + hd10s99)/2.Q+00
  mbo2d10rho(10)=(-b10 + hd10p11 + hd10p99)/2.Q+00
  mbo2d10rho(11:13)=hd10p22
  mbo2d10rho(14:17)=hd10p55
  mbo2d10rho(18)=(b10 + hd10p11 + hd10p99)/2.Q+00
  mbo2d20rho(1)=(a20 + hd20s11 + hd20s99)/2.Q+00
  mbo2d20rho(2:4)=hd20s22
  mbo2d20rho(5:8)=hd20s55
  mbo2d20rho(9)=(-a20 + hd20s11 + hd20s99)/2.Q+00
  mbo2d20rho(10)=(-b20 + hd20p11 + hd20p99)/2.Q+00
  mbo2d20rho(11:13)=hd20p22
  mbo2d20rho(14:17)=hd20p55
  mbo2d20rho(18)=(b20 + hd20p11 + hd20p99)/2.Q+00
  mbo2d30rho(1)=(a30 + hd30s11 + hd30s99)/2.Q+00
  mbo2d30rho(2:4)=hd30s22
  mbo2d30rho(5:8)=hd30s55
  mbo2d30rho(9)=(-a30 + hd30s11 + hd30s99)/2.Q+00
  mbo2d30rho(10)=(-b30 + hd30p11 + hd30p99)/2.Q+00
  mbo2d30rho(11:13)=hd30p22
  mbo2d30rho(14:17)=hd30p55
  mbo2d30rho(18)=(b30 + hd30p11 + hd30p99)/2.Q+00
  mbo2d40rho(1)=(a40 + hd40s11 + hd40s99)/2.Q+00
  mbo2d40rho(2:4)=hd40s22
  mbo2d40rho(5:8)=hd40s55
  mbo2d40rho(9)=(-a40 + hd40s11 + hd40s99)/2.Q+00
  mbo2d40rho(10)=(-b40 + hd40p11 + hd40p99)/2.Q+00
  mbo2d40rho(11:13)=hd40p22
  mbo2d40rho(14:17)=hd40p55
  mbo2d40rho(18)=(b40 + hd40p11 + hd40p99)/2.Q+00
  mbo2d50rho(1)=(a50 + hd50s11 + hd50s99)/2.Q+00
  mbo2d50rho(2:4)=hd50s22
  mbo2d50rho(5:8)=hd50s55
  mbo2d50rho(9)=(-a50 + hd50s11 + hd50s99)/2.Q+00
  mbo2d50rho(10)=(-b50 + hd50p11 + hd50p99)/2.Q+00
  mbo2d50rho(11:13)=hd50p22
  mbo2d50rho(14:17)=hd50p55
  mbo2d50rho(18)=(b50 + hd50p11 + hd50p99)/2.Q+00
  mbo2d01rho(1)=(a01 + hd01s11 + hd01s99)/2.Q+00
  mbo2d01rho(2:4)=hd01s22
  mbo2d01rho(5:8)=hd01s55
  mbo2d01rho(9)=(-a01 + hd01s11 + hd01s99)/2.Q+00
  mbo2d01rho(10)=(-b01 + hd01p11 + hd01p99)/2.Q+00
  mbo2d01rho(11:13)=hd01p22
  mbo2d01rho(14:17)=hd01p55
  mbo2d01rho(18)=(b01 + hd01p11 + hd01p99)/2.Q+00
  mbo2d11rho(1)=(a11 + hd11s11 + hd11s99)/2.Q+00
  mbo2d11rho(2:4)=hd11s22
  mbo2d11rho(5:8)=hd11s55
  mbo2d11rho(9)=(-a11 + hd11s11 + hd11s99)/2.Q+00
  mbo2d11rho(10)=(-b11 + hd11p11 + hd11p99)/2.Q+00
  mbo2d11rho(11:13)=hd11p22
  mbo2d11rho(14:17)=hd11p55
  mbo2d11rho(18)=(b11 + hd11p11 + hd11p99)/2.Q+00
  mbo2d21rho(1)=(a21 + hd21s11 + hd21s99)/2.Q+00
  mbo2d21rho(2:4)=hd21s22
  mbo2d21rho(5:8)=hd21s55
  mbo2d21rho(9)=(-a21 + hd21s11 + hd21s99)/2.Q+00
  mbo2d21rho(10)=(-b21 + hd21p11 + hd21p99)/2.Q+00
  mbo2d21rho(11:13)=hd21p22
  mbo2d21rho(14:17)=hd21p55
  mbo2d21rho(18)=(b21 + hd21p11 + hd21p99)/2.Q+00
  mbo2d31rho(1)=(a31 + hd31s11 + hd31s99)/2.Q+00
  mbo2d31rho(2:4)=hd31s22
  mbo2d31rho(5:8)=hd31s55
  mbo2d31rho(9)=(-a31 + hd31s11 + hd31s99)/2.Q+00
  mbo2d31rho(10)=(-b31 + hd31p11 + hd31p99)/2.Q+00
  mbo2d31rho(11:13)=hd31p22
  mbo2d31rho(14:17)=hd31p55
  mbo2d31rho(18)=(b31 + hd31p11 + hd31p99)/2.Q+00
  mbo2d02rho(1)=(a02 + hd02s11 + hd02s99)/2.Q+00
  mbo2d02rho(2:4)=hd02s22
  mbo2d02rho(5:8)=hd02s55
  mbo2d02rho(9)=(-a02 + hd02s11 + hd02s99)/2.Q+00
  mbo2d02rho(10)=(-b02 + hd02p11 + hd02p99)/2.Q+00
  mbo2d02rho(11:13)=hd02p22
  mbo2d02rho(14:17)=hd02p55
  mbo2d02rho(18)=(b02 + hd02p11 + hd02p99)/2.Q+00
  mbo2d12rho(1)=(a12 + hd12s11 + hd12s99)/2.Q+00
  mbo2d12rho(2:4)=hd12s22
  mbo2d12rho(5:8)=hd12s55
  mbo2d12rho(9)=(-a12 + hd12s11 + hd12s99)/2.Q+00
  mbo2d12rho(10)=(-b12 + hd12p11 + hd12p99)/2.Q+00
  mbo2d12rho(11:13)=hd12p22
  mbo2d12rho(14:17)=hd12p55
  mbo2d12rho(18)=(b12 + hd12p11 + hd12p99)/2.Q+00    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dmbo2drho(:,1)=mboson2
  dmbo2drho(:,2)=mbo2d10rho
  dmbo2drho(:,3)=mbo2d20rho
  dmbo2drho(:,4)=mbo2d30rho
  dmbo2drho(:,5)=mbo2d40rho
  dmbo2drho(:,6)=mbo2d50rho
  dmbo2drho(:,7)=mbo2d01rho
  dmbo2drho(:,8)=mbo2d11rho
  dmbo2drho(:,9)=mbo2d21rho
  dmbo2drho(:,10)=mbo2d31rho
  dmbo2drho(:,11)=mbo2d02rho
  dmbo2drho(:,12)=mbo2d12rho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  thetas=ATan((2*hd00s19)/(hd00s11 - hd00s99))/2.Q+00
  thetap=ATan((2*hd00p19)/(hd00p11 - hd00p99))/2.Q+00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !scalar,psendoscalar
  cos2phi(1)=0.5Q+00+((hd00s11-hd00s99)/Sqrt(4*hd00s19**2+(hd00s11-hd00s99)**2) &
   /2.Q+00+Sqrt(2.Q+00)*2*hd00s19/Sqrt(4*hd00s19**2+(hd00s11-hd00s99)**2))/3.Q+00
  cos2phi(2)=0.5Q+00+((hd00p99-hd00p11)/Sqrt(4*hd00p19**2+(hd00p11-hd00p99)**2) &
   /2.Q+00+Sqrt(2.Q+00)*2*hd00p19/Sqrt(4*hd00p19**2+(hd00p11-hd00p99)**2))/3.Q+00
  sin2phi(1)=0.5Q+00-((hd00s11-hd00s99)/Sqrt(4*hd00s19**2+(hd00s11-hd00s99)**2) &
   /2.Q+00+Sqrt(2.Q+00)*2*hd00s19/Sqrt(4*hd00s19**2+(hd00s11-hd00s99)**2))/3.Q+00
  sin2phi(2)=0.5Q+00-((hd00p99-hd00p11)/Sqrt(4*hd00p19**2+(hd00p11-hd00p99)**2) &
   /2.Q+00+Sqrt(2.Q+00)*2*hd00p19/Sqrt(4*hd00p19**2+(hd00p11-hd00p99)**2))/3.Q+00

end
