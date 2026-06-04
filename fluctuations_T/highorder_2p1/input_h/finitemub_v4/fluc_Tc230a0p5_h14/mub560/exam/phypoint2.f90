subroutine phypoint2(Nflow,yflow,rho0,mboson,mfermion,Vall,fpifk)
!calculate the physical point rho0,rho20
!wenrui wrote 2017.10.16
!wenrui changed 2017.11.10

  implicit none
  integer Nflow
  real(16) yflow(Nflow)
  integer N_str(4)
  !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck
  real(16) rho,rho2
  real(16) rho0(2)
  real(16) mboson2(8)
  real(16) mboson(8),mfermion(2)
  real(16) Vall 
  !total effective potential
  real(16) fpifk(2)
  !fpi fK
  real(16) lam00,lam10,lam20,lam30,lam40,lam50,lam60,lam70
  real(16) lam01,lam11,lam21,lam31,lam41,lam51
  real(16) lam02,lam12,lam22,lam32
  real(16) lam03,lam13
  real(16) hlk,hsk
  real(16) Z_pi,Z_K,Z_l,Z_s
  real(16) ck,jl,js
  real(16) kappa1,kappa2,Sl,Ss
  real(16) hd00s11,hd00s99,hd00s19,hd00s22,hd00s55
  real(16) hd00p11,hd00p99,hd00p19,hd00p22,hd00p55
  integer nn
  parameter(nn=2)
  real(16) x(nn)
  logical check
  external gapEq

  real*16 :: h_mod,hs_mod
  REAL*16 :: h_const,lamd10_const,lamd20_const,lamd01_const,Tcglue_const,hs_const
  REAL*16 :: SLSS_IR(2)

  common /strucFun/ N_str
  common /gapPara/ lam10,lam20,lam30,lam40,lam50,lam60,lam70,lam01,lam11,       &
                   lam21,lam31,lam41,lam51,lam02,lam12,lam22,lam32,lam03,       &
                   lam13,ck,kappa1,kappa2,jl,js

  COMMON /IR_RESULT/ SLSS_IR
  common /ini_const/ h_const,lamd10_const,lamd20_const,lamd01_const,Tcglue_const,hs_const

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)

  lam10=yflow(1)
  lam20=yflow(2)
  lam30=yflow(3)
  lam40=yflow(4)
  lam50=yflow(5)
  lam60=yflow(6)
  lam70=yflow(7)
  lam01=yflow(8)
  lam11=yflow(9)
  lam21=yflow(10)
  lam31=yflow(11)
  lam41=yflow(12)
  lam51=yflow(13)
  lam02=yflow(14)
  lam12=yflow(15)
  lam22=yflow(16)
  lam32=yflow(17)
  lam03=yflow(18)
  lam13=yflow(19)
  lam00=yflow(Nv+1)
  hlk=yflow((Nv+1)+1)
  hsk=yflow((Nv+1)+2)
  Z_pi=yflow((Nv+1)+(Nh+2)+1)
  Z_K=yflow((Nv+1)+(Nh+2)+2)
  Z_l=yflow((Nv+1)+(Nh+2)+3)
  Z_s=yflow((Nv+1)+(Nh+2)+4)
  ck=yflow((Nv+1)+(Nh+2)+Nz+1)
  Sl=yflow((Nv+1)+(Nh+2)+Nz+2)
  Ss=yflow((Nv+1)+(Nh+2)+Nz+3)
  jl=yflow((Nv+1)+(Nh+2)+Nz+4)
  js=yflow((Nv+1)+(Nh+2)+Nz+5)

  kappa1=(Sl**2 + Ss**2)/2.Q+00
  kappa2=(Sl**2 - 2*Ss**2)**2/24.Q+00

  rho0(1)=kappa1
  rho0(2)=kappa2

  SLSS_IR(1)=SL
  SLSS_IR(2)=SS

  fpifk(1)=Sl
  !fpi=sigmal
  fpifk(2)=(Sl+Sqrt(2.Q+00)*Ss)/2.Q+00
  !fK=(sigmal+Sqrt(2.Q+00)*sigmas)/2.

  rho=rho0(1)
  rho2=rho0(2)

  hd00s11=lam10 + (lam30*(kappa1 - rho)**2)/2.Q+00+ (lam50*(kappa1 - rho)**4)/24.Q+00+ &
    lam20*(-kappa1 + rho) + (lam40*(-kappa1 + rho)**3)/6.Q+00 -                     &
    (2*ck*(Sqrt(4*rho - 2*Sqrt(6.Q+00)*Sqrt(rho2)) +                                &
         Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2))))/(3.Q+00*Sqrt(3.Q+00)) +                      &
    lam21*(kappa1 - rho)*(kappa2 - rho2) +                                      &
    (lam12*(-kappa1 + rho)*(kappa2 - rho2)**2)/2.Q+00+                             &
    lam11*(-kappa2 + rho2) + (lam31*(kappa1 - rho)**2*(-kappa2 + rho2))/        &
     2.Q+00 + (2*(Sqrt(4*rho - 2*Sqrt(6.Q+00)*Sqrt(rho2)) -                             &
          2*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))**2*                                &
       (lam01 + (lam21*(kappa1 - rho)**2)/2.Q+00+ lam11*(-kappa1 + rho) +          &
         (lam31*(-kappa1 + rho)**3)/6.Q+00 +                                        &
         lam12*(kappa1 - rho)*(kappa2 - rho2) + lam02*(-kappa2 + rho2)))/       &
     27.Q+00 + ((2*Sqrt(2*rho - Sqrt(6.Q+00)*Sqrt(rho2)) +                              &
         Sqrt(2.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*                             &
       (((6*lam20 - (kappa1 - rho)*                                             &
               (6*lam30 + (kappa1 - rho)*                                       &
                  (-3*lam40 + kappa1*lam50 - lam50*rho)))*                      &
            (2*Sqrt(6*rho - 3*Sqrt(6.Q+00)*Sqrt(rho2)) +                            &
              Sqrt(6.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2))))/3.Q+00 +                   &
         ((Sqrt(2*rho - Sqrt(6.Q+00)*Sqrt(rho2)) +                                  &
              Sqrt(2.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*                        &
            (Sqrt(12*rho - 6*Sqrt(6.Q+00)*Sqrt(rho2)) -                             &
               2*Sqrt(3.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))**2*                  &
            (2*lam11 - 2*kappa2*lam12 - 2*kappa1*lam21 +                        &
              kappa1**2*lam31 + 2*lam21*rho - 2*kappa1*lam31*rho +              &
              lam31*rho**2 + 2*lam12*rho2))/(9.Q+00*Sqrt(3.Q+00))))/(18.Q+00*Sqrt(3.Q+00)) +    &
    ((Sqrt(4*rho - 2*Sqrt(6.Q+00)*Sqrt(rho2)) -                                     &
          2*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))**2*                                &
       (Sqrt(2*rho - Sqrt(6.Q+00)*Sqrt(rho2)) +                                     &
         Sqrt(2.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*                             &
       ((2*(lam02 + lam12*(-kappa1 + rho))*                                     &
            (Sqrt(2*rho - Sqrt(6.Q+00)*Sqrt(rho2)) +                                &
              Sqrt(2.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*                        &
            (Sqrt(12*rho - 6*Sqrt(6.Q+00)*Sqrt(rho2)) -                             &
               2*Sqrt(3.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))**2)/(9.Q+00*Sqrt(3.Q+00))    &
+ 2*Sqrt(3.Q+00)*(2*Sqrt(2*rho - Sqrt(6.Q+00)*Sqrt(rho2)) +                             &
            Sqrt(2.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*                          &
          (lam11 + (lam31*(kappa1 - rho)**2)/2.Q+00+ lam21*(-kappa1 + rho) +       &
            lam12*(-kappa2 + rho2))))/(162.Q+00*Sqrt(3.Q+00))

  hd00s22=lam10 + (lam30*(kappa1 - rho)**2)/2.Q+00+                                &
    (lam50*(kappa1 - rho)**4)/24.Q+00 + lam20*(-kappa1 + rho) +                     &
    (lam40*(-kappa1 + rho)**3)/6.Q+00 +                                             &
    (ck*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))/Sqrt(3.Q+00) +                             &
    lam21*(kappa1 - rho)*(kappa2 - rho2) +                                      &
    (lam12*(-kappa1 + rho)*(kappa2 - rho2)**2)/2.Q+00+                             &
    lam11*(-kappa2 + rho2) + (lam31*(kappa1 - rho)**2*(-kappa2 + rho2))/        &
     2.Q+00 + ((8*rho - 6*Sqrt(6.Q+00)*Sqrt(rho2))*                                     &
       (lam01 + (lam21*(kappa1 - rho)**2)/2.Q+00+ lam11*(-kappa1 + rho) +          &
         (lam31*(-kappa1 + rho)**3)/6.Q+00 +                                        &
         lam12*(kappa1 - rho)*(kappa2 - rho2) + lam02*(-kappa2 + rho2)))/6.Q+00

  hd00s55=(ck*Sqrt(12*rho - 6*Sqrt(6.Q+00)*Sqrt(rho2)) +                            &
      2*(2*rho + Sqrt(4*rho - 2*Sqrt(6.Q+00)*Sqrt(rho2))*                           &
          Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)) + Sqrt(6.Q+00)*Sqrt(rho2))*               &
       (lam01 + (lam21*(kappa1 - rho)**2)/2.Q+00+ lam11*(-kappa1 + rho) +          &
         (lam31*(-kappa1 + rho)**3)/6.Q+00 +                                        &
         lam12*(kappa1 - rho)*(kappa2 - rho2) + lam02*(-kappa2 + rho2)) +       &
      6*(lam10 + (lam30*(kappa1 - rho)**2)/2.Q+00+                                 &
         (lam50*(kappa1 - rho)**4)/24.Q+00 + lam20*(-kappa1 + rho) +                &
         (lam40*(-kappa1 + rho)**3)/6.Q+00 +                                        &
         lam21*(kappa1 - rho)*(kappa2 - rho2) +                                 &
         (lam12*(-kappa1 + rho)*(kappa2 - rho2)**2)/2.Q+00+                        &
         lam11*(-kappa2 + rho2) +                                               &
         (lam31*(kappa1 - rho)**2*(-kappa2 + rho2))/2.Q+00))/6.Q+00

  hd00s99=(12*ck*(2*Sqrt(12*rho - 6*Sqrt(6.Q+00)*Sqrt(rho2)) -                      &
         Sqrt(3.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2))) +                            &
      (2*(20*rho + 8*Sqrt(4*rho - 2*Sqrt(6.Q+00)*Sqrt(rho2))*                       &
            Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)) + 23*Sqrt(6.Q+00)*Sqrt(rho2))*          &
         (6*lam01 - 6*kappa1*lam11 + 3*kappa1**2*lam21 -                        &
           kappa1**3*lam31 + 6*lam11*rho - 6*kappa1*lam21*rho +                 &
           3*kappa1**2*lam31*rho + 3*lam21*rho**2 -                             &
           3*kappa1*lam31*rho**2 + lam31*rho**3 -                               &
           6*kappa2*(lam02 - kappa1*lam12 + lam12*rho) + 6*lam02*rho2 -         &
           6*kappa1*lam12*rho2 + 6*lam12*rho*rho2))/3.Q+00 +                        &
      108*(lam10 + (lam30*(kappa1 - rho)**2)/2.Q+00+                               &
         (lam50*(kappa1 - rho)**4)/24.Q+00 + lam20*(-kappa1 + rho) +                &
         (lam40*(-kappa1 + rho)**3)/6.Q+00 +                                        &
         lam21*(kappa1 - rho)*(kappa2 - rho2) +                                 &
         (lam12*(-kappa1 + rho)*(kappa2 - rho2)**2)/2.Q+00+                        &
         lam11*(-kappa2 + rho2) +                                               &
         (lam31*(kappa1 - rho)**2*(-kappa2 + rho2))/2.Q+00) +                       &
      2*Sqrt(3.Q+00)*(Sqrt(4*rho - 2*Sqrt(6.Q+00)*Sqrt(rho2)) -                         &
         2*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*                                    &
       (((6*lam20 - (kappa1 - rho)*                                             &
               (6*lam30 + (kappa1 - rho)*                                       &
                  (-3*lam40 + kappa1*lam50 - lam50*rho)))*                      &
            (Sqrt(12*rho - 6*Sqrt(6.Q+00)*Sqrt(rho2)) -                             &
              2*Sqrt(3.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2))))/3.Q+00 -                 &
         Sqrt(2/3.Q+00)*                                              &
          (Sqrt(12*rho - 6*Sqrt(6.Q+00)*Sqrt(rho2)) +                               &
            4*Sqrt(3.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*Sqrt(rho2)*             &
          (2*lam11 - 2*kappa2*lam12 - 2*kappa1*lam21 + kappa1**2*lam31 +        &
            2*lam21*rho - 2*kappa1*lam31*rho + lam31*rho**2 + 2*lam12*rho2      &
)) - 2*Sqrt(2.Q+00)*(Sqrt(4*rho - 2*Sqrt(6.Q+00)*Sqrt(rho2)) +                          &
         4*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*Sqrt(rho2)*                         &
       (-2*Sqrt(2/3.Q+00)*(lam02 + lam12*(-kappa1 + rho))*            &
          (Sqrt(12*rho - 6*Sqrt(6.Q+00)*Sqrt(rho2)) +                               &
            4*Sqrt(3.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*Sqrt(rho2) +            &
         2*Sqrt(3.Q+00)*(Sqrt(4*rho - 2*Sqrt(6.Q+00)*Sqrt(rho2)) -                      &
            2*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*                                 &
          (lam11 + (lam31*(kappa1 - rho)**2)/2.Q+00+ lam21*(-kappa1 + rho) +       &
            lam12*(-kappa2 + rho2))))/108.Q+00

  hd00s19=(12*Sqrt(3.Q+00)*ck*(Sqrt(2*rho - Sqrt(6.Q+00)*Sqrt(rho2)) -                  &
         Sqrt(2.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2))) -                            &
      (4*(-Sqrt(4*rho - 2*Sqrt(6.Q+00)*Sqrt(rho2)) +                                &
           2*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*                                  &
         (5*Sqrt(6*rho - 3*Sqrt(6.Q+00)*Sqrt(rho2)) +                               &
           7*Sqrt(6.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*                         &
         (lam01 + (lam21*(kappa1 - rho)**2)/2.Q+00+ lam11*(-kappa1 + rho) +        &
           (lam31*(-kappa1 + rho)**3)/6.Q+00 +                                      &
           lam12*(kappa1 - rho)*(kappa2 - rho2) + lam02*(-kappa2 + rho2))       &
)/Sqrt(3.Q+00) + 2*Sqrt(3.Q+00)*(2*Sqrt(2*rho - Sqrt(6.Q+00)*Sqrt(rho2)) +                  &
         Sqrt(2.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*                             &
       (((6*lam20 - (kappa1 - rho)*                                             &
               (6*lam30 + (kappa1 - rho)*                                       &
                  (-3*lam40 + kappa1*lam50 - lam50*rho)))*                      &
            (Sqrt(12*rho - 6*Sqrt(6.Q+00)*Sqrt(rho2)) -                             &
              2*Sqrt(3.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2))))/3.Q+00 -                 &
         Sqrt(2/3.Q+00)*                                              &
          (Sqrt(12*rho - 6*Sqrt(6.Q+00)*Sqrt(rho2)) +                               &
            4*Sqrt(3.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*Sqrt(rho2)*             &
          (2*lam11 - 2*kappa2*lam12 - 2*kappa1*lam21 + kappa1**2*lam31 +        &
            2*lam21*rho - 2*kappa1*lam31*rho + lam31*rho**2 + 2*lam12*rho2      &
)) + (2*(Sqrt(4*rho - 2*Sqrt(6.Q+00)*Sqrt(rho2)) -                                  &
            2*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))**2*                              &
         (Sqrt(2*rho - Sqrt(6.Q+00)*Sqrt(rho2)) +                                   &
           Sqrt(2.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*                           &
         (-2*Sqrt(2/3.Q+00)*(lam02 + lam12*(-kappa1 + rho))*          &
            (Sqrt(12*rho - 6*Sqrt(6.Q+00)*Sqrt(rho2)) +                             &
              4*Sqrt(3.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*Sqrt(rho2) +          &
           2*Sqrt(3.Q+00)*(Sqrt(4*rho - 2*Sqrt(6.Q+00)*Sqrt(rho2)) -                    &
              2*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))*                               &
            (lam11 + (lam31*(kappa1 - rho)**2)/2.Q+00+                             &
              lam21*(-kappa1 + rho) + lam12*(-kappa2 + rho2))))/                &
       (3.Q+00*Sqrt(3.Q+00)))/108.Q+00

  hd00p11=lam10 + (lam30*(kappa1 - rho)**2)/2.Q+00+ (lam50*(kappa1 - rho)**4)/24.Q+00 +&
    lam20*(-kappa1 + rho) + (lam40*(-kappa1 + rho)**3)/6.Q+00 +                     &
    (2*ck*(Sqrt(4*rho - 2*Sqrt(6.Q+00)*Sqrt(rho2)) +                                &
         Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2))))/(3.Q+00*Sqrt(3.Q+00)) +                      &
    lam21*(kappa1 - rho)*(kappa2 - rho2) +                                      &
    (lam12*(-kappa1 + rho)*(kappa2 - rho2)**2)/2.Q+00+                             &
    lam11*(-kappa2 + rho2) + (lam31*(kappa1 - rho)**2*(-kappa2 + rho2))/2.Q+00

  hd00p22=lam10 + (lam30*(kappa1 - rho)**2)/2.Q+00+                                &
    (lam50*(kappa1 - rho)**4)/24.Q+00 + lam20*(-kappa1 + rho) +                     &
    (lam40*(-kappa1 + rho)**3)/6.Q+00 -                                             &
    (ck*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))/Sqrt(3.Q+00) +                             &
    lam21*(kappa1 - rho)*(kappa2 - rho2) +                                      &
    (lam12*(-kappa1 + rho)*(kappa2 - rho2)**2)/2.Q+00+                             &
    lam11*(-kappa2 + rho2) + (lam31*(kappa1 - rho)**2*(-kappa2 + rho2))/        &
     2.Q+00 - Sqrt(2/3.Q+00)*Sqrt(rho2)*                                  &
     (lam01 + (lam21*(kappa1 - rho)**2)/2.Q+00+ lam11*(-kappa1 + rho) +            &
       (lam31*(-kappa1 + rho)**3)/6.Q+00 +                                          &
       lam12*(kappa1 - rho)*(kappa2 - rho2) + lam02*(-kappa2 + rho2))

  hd00p55=(-(ck*Sqrt(12*rho - 6*Sqrt(6.Q+00)*Sqrt(rho2))) +                         &
      (4*rho - 2*Sqrt(4*rho - 2*Sqrt(6.Q+00)*Sqrt(rho2))*                           &
          Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)) + 2*Sqrt(6.Q+00)*Sqrt(rho2))*             &
       (lam01 + (lam21*(kappa1 - rho)**2)/2.Q+00+ lam11*(-kappa1 + rho) +          &
         (lam31*(-kappa1 + rho)**3)/6.Q+00 +                                        &
         lam12*(kappa1 - rho)*(kappa2 - rho2) + lam02*(-kappa2 + rho2)) +       &
      6*(lam10 + (lam30*(kappa1 - rho)**2)/2.Q+00+                                 &
         (lam50*(kappa1 - rho)**4)/24.Q+00 + lam20*(-kappa1 + rho) +                &
         (lam40*(-kappa1 + rho)**3)/6.Q+00 +                                        &
         lam21*(kappa1 - rho)*(kappa2 - rho2) +                                 &
         (lam12*(-kappa1 + rho)*(kappa2 - rho2)**2)/2.Q+00+                        &
         lam11*(-kappa2 + rho2) +                                               &
        (lam31*(kappa1 - rho)**2*(-kappa2 + rho2))/2.Q+00))/6.Q+00

  hd00p99=(-4*Sqrt(2/3.Q+00)*ck*                                      &
       Sqrt(2*rho - Sqrt(6.Q+00)*Sqrt(rho2)) +                                      &
      (2*ck*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2)))/Sqrt(3.Q+00) +                         &
      2*Sqrt(6.Q+00)*Sqrt(rho2)*(lam01 + (lam21*(kappa1 - rho)**2)/2.Q+00+             &
         lam11*(-kappa1 + rho) + (lam31*(-kappa1 + rho)**3)/6.Q+00 +                &
         lam12*(kappa1 - rho)*(kappa2 - rho2) + lam02*(-kappa2 + rho2)) +       &
      6*(lam10 + (lam30*(kappa1 - rho)**2)/2.Q+00+                                 &
         (lam50*(kappa1 - rho)**4)/24.Q+00 + lam20*(-kappa1 + rho) +                &
         (lam40*(-kappa1 + rho)**3)/6.Q+00 +                                        &
         lam21*(kappa1 - rho)*(kappa2 - rho2) +                                 &
         (lam12*(-kappa1 + rho)*(kappa2 - rho2)**2)/2.Q+00+                        &
         lam11*(-kappa2 + rho2) +                                               &
         (lam31*(kappa1 - rho)**2*(-kappa2 + rho2))/2.Q+00))/6.Q+00

  hd00p19=-(ck*(Sqrt(2*rho - Sqrt(6.Q+00)*Sqrt(rho2)) -                             &
          Sqrt(2.Q+00)*Sqrt(rho + Sqrt(6.Q+00)*Sqrt(rho2))) +                           &
       Sqrt(rho2)*(6*lam01 - 6*kappa1*lam11 + 3*kappa1**2*lam21 -               &
          kappa1**3*lam31 + 6*lam11*rho - 6*kappa1*lam21*rho +                  &
          3*kappa1**2*lam31*rho + 3*lam21*rho**2 - 3*kappa1*lam31*rho**2 +      &
          lam31*rho**3 - 6*kappa2*(lam02 - kappa1*lam12 + lam12*rho) +          &
          6*lam02*rho2 - 6*kappa1*lam12*rho2 + 6*lam12*rho*rho2))/              &
    (3.Q+00*Sqrt(3.Q+00))

  mboson2(1)=hd00s99*Cos(ATan((2*hd00s19)/(hd00s11 - hd00s99))/2.Q+00)**2 -         &
    2*hd00s19*Cos(ATan((2*hd00s19)/(hd00s11 - hd00s99))/2.Q+00)*                    &
     Sin(ATan((2*hd00s19)/(hd00s11 - hd00s99))/2.Q+00) +                            &
    hd00s11*Sin(ATan((2*hd00s19)/(hd00s11 - hd00s99))/2.Q+00)**2
  mboson2(2)=hd00s22
  mboson2(3)=hd00s55
  mboson2(4)=hd00s11*Cos(ATan((2*hd00s19)/(hd00s11 - hd00s99))/2.Q+00)**2 +         &
    2*hd00s19*Cos(ATan((2*hd00s19)/(hd00s11 - hd00s99))/2.Q+00)*                    &
     Sin(ATan((2*hd00s19)/(hd00s11 - hd00s99))/2.Q+00) +                            &
    hd00s99*Sin(ATan((2*hd00s19)/(hd00s11 - hd00s99))/2.Q+00)**2
  mboson2(5)=hd00p99*Cos(ATan((2*hd00p19)/(hd00p11 - hd00p99))/2.Q+00)**2 -         &
    2*hd00p19*Cos(ATan((2*hd00p19)/(hd00p11 - hd00p99))/2.Q+00)*                    &
     Sin(ATan((2*hd00p19)/(hd00p11 - hd00p99))/2.Q+00) +                            &
    hd00p11*Sin(ATan((2*hd00p19)/(hd00p11 - hd00p99))/2.Q+00)**2
  mboson2(6)=hd00p22
  mboson2(7)=hd00p55
  mboson2(8)=hd00p11*Cos(ATan((2*hd00p19)/(hd00p11 - hd00p99))/2.Q+00)**2 +         &
    2*hd00p19*Cos(ATan((2*hd00p19)/(hd00p11 - hd00p99))/2.Q+00)*                    &
     Sin(ATan((2*hd00p19)/(hd00p11 - hd00p99))/2.Q+00) +                            &
    hd00p99*Sin(ATan((2*hd00p19)/(hd00p11 - hd00p99))/2.Q+00)**2

  mboson=sqrt(mboson2)

  call fityukawa(1.Q+0,h_mod,hs_mod)
  hlk=h_mod*h_const
!  hsk=h_mod*h_const
  hsk=hs_mod*hs_const

  mfermion(1)=hlk*Sqrt(2/3.Q+00)*Sqrt(2*rho - Sqrt(6.Q+00)*Sqrt(rho2))/2.Q+00
  mfermion(2)=hsk*Sqrt(2/3.Q+00)*Sqrt(rho+Sqrt(6.Q+00)*Sqrt(rho2))        &
              /Sqrt(2.Q+00)

  Vall= lam00+lam10*(rho-kappa1)+lam20/2.Q+00*(rho-kappa1)**2                       &
       +lam30/6.Q+00*(rho-kappa1)**3+lam40/24.Q+00*(rho-kappa1)**4                      &
       +lam50/120.Q+00*(rho-kappa1)**5+lam60/720.Q+00*(rho-kappa1)**6                   &
       +lam70/5040.Q+00*(rho-kappa1)**7+lam01*(rho2-kappa2)                         &
       +lam11*(rho-kappa1)*(rho2-kappa2)+lam21/2.Q+00*(rho-kappa1)**2*(rho2-kappa2) &
       +lam31/6.Q+00*(rho-kappa1)**3*(rho2-kappa2)                                  &
       +lam41/24.Q+00*(rho-kappa1)**4*(rho2-kappa2)                                 &
       +lam51/120.Q+00*(rho-kappa1)**5*(rho2-kappa2)                                &
       +lam02/2.Q+00*(rho2-kappa2)**2+lam12/2.Q+00*(rho-kappa1)*(rho2-kappa2)**2        &
       +lam22/4.Q+00*(rho-kappa1)**2*(rho2-kappa2)**2                               &
       +lam32/12.Q+00*(rho-kappa1)**3*(rho2-kappa2)**2                              &
       +lam03/6.Q+00*(rho2-kappa2)**3+lam13/6.Q+00*(rho-kappa1)*(rho2-kappa2)**3        &
       -ck*(Sl**2.Q+00*Ss)/2.Q+00/Sqrt(2.Q+00)-jl*Sl-js*Ss

end
