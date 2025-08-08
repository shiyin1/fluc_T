MODULE BF_POL_MOD
  CONTAINS
  SUBROUTINE BF_POL(MP2,MS2,MF2,K,L_CON,LB_CON,BF_POL_COM)
  USE TMU_MOD
  USE MU0P0_COM
  USE FNF_MOD
  IMPLICIT NONE
  
  REAL(16),INTENT(IN) :: MP2,MS2,MF2,K,L_CON,LB_CON
  COMPLEX(16) :: dB2F1AdlC,dB2F1AdlbC,dB1F1AdlC,dB1F1AdlbC,dB1F2AdlC,dB1F2AdlbC,dB2F1SdlC,dB2F1SdlbC,dB2F1PdlC,dB2F1PdlbC
  REAL(16) :: dB2F1Adl,dB2F1Adlb,dB1F1Adl,dB1F1Adlb,dB1F2Adl,dB1F2Adlb,dB2F1Sdl,dB2F1Sdlb,dB2F1Pdl,dB2F1Pdlb  
  REAL(16),INTENT(OUT),DIMENSION(10) :: BF_POL_COM
  REAL(16) :: nfd1lf,nfd1lbf,nfd1la,nfd1lba 
  REAL(16) :: nfd1xd1lf,nfd1xd1lbf,nfd1xd1la,nfd1xd1lba 
  REAL(16) :: nf0,nf1,nf2,L,LB,XFF,XFA

  L=L_CON
  LB=LB_CON
     XFF=K*SQRT(1.Q+00 + MF2) - MU
     CALL FNF012(XFF,T,L,LB,NF0,NF1,NF2)
     
     nfd1lf=-(nf2*(-1 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
     nfd1lbf=-0.5*(nf1*(-2 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
     nfd1xd1lf=-0.25*(nf2*(4 + 36*(2 + l*lb)*nf0**2 + 9*lb**2*nf1**2 - 36*l*nf2 + &
               72*l**2*nf2**2 + 6*nf0*(-10 - 6*l*lb + 18*lb*nf1 + 9*l*lb**2*nf1 + &
               24*l*nf2 + 18*l**2*lb*nf2 + 6*lb**2*nf2) +                         &
               6*lb*(-2*nf1 + 9*l*nf1*nf2 - 12*nf2**2)))/T
     nfd1xd1lbf=(6*lb*(1 - 3*nf0)*nf1**2 + 24*lb*(2 - 3*nf0)*nf0*nf2 +            &
                18*l**2*nf0*(-5 + 8*nf0 + 9*lb*nf1)*nf2 + 108*l**3*nf0*nf2**2 -   & 
                nf1*(4 - 21*nf0 + 18*nf0**2 + 72*lb**2*nf0*nf2) +                 &
                3*l*(12*nf0**3 + 6*nf0**2*(-3 + 5*lb*nf1) + 2*nf1*nf2 +           & 
                3*nf0*(2 - 7*lb*nf1 + 6*lb**2*nf1**2 - 2*nf1*nf2 - 8*lb*nf2**2)))/&
                (2.*T)
           
  L=LB_CON
  LB=L_CON
     XFA=K*SQRT(1.Q+00 + MF2) + MU
     CALL FNF012(XFA,T,L,LB,NF0,NF1,NF2)

     nfd1lba=-(nf2*(-1 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
     nfd1la=-0.5*(nf1*(-2 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
     nfd1xd1lba=-0.25*(nf2*(4 + 36*(2 + l*lb)*nf0**2 + 9*lb**2*nf1**2 - 36*l*nf2 + &
               72*l**2*nf2**2 + 6*nf0*(-10 - 6*l*lb + 18*lb*nf1 + 9*l*lb**2*nf1 + &
               24*l*nf2 + 18*l**2*lb*nf2 + 6*lb**2*nf2) +                         &
               6*lb*(-2*nf1 + 9*l*nf1*nf2 - 12*nf2**2)))/T
     nfd1xd1la=(6*lb*(1 - 3*nf0)*nf1**2 + 24*lb*(2 - 3*nf0)*nf0*nf2 +            &
               18*l**2*nf0*(-5 + 8*nf0 + 9*lb*nf1)*nf2 + 108*l**3*nf0*nf2**2 -   & 
               nf1*(4 - 21*nf0 + 18*nf0**2 + 72*lb**2*nf0*nf2) +                 &
               3*l*(12*nf0**3 + 6*nf0**2*(-3 + 5*lb*nf1) + 2*nf1*nf2 +           & 
               3*nf0*(2 - 7*lb*nf1 + 6*lb**2*nf1**2 - 2*nf1*nf2 - 8*lb*nf2**2)))/&
               (2.*T)
   
     dB1F1AdlC=(k**2*(nfd1la/(Sqrt(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) + (0,1)*p0)**2)) &
               + nfd1lf/(Sqrt(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) + (0,1)*p0)**2))))/2.                     
   
     dB1F1Adl=real(dB1F1AdlC)

     dB1F1AdlbC=(k**2*(nfd1lba/(Sqrt(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) + (0,1)*p0)**2)) &
                + nfd1lbf/(Sqrt(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) + (0,1)*p0)**2))))/2.

     dB1F1Adlb=real(dB1F1AdlbC)

     dB1F2AdlC=-0.5*(k**2*(-0.5*nfd1la/((1 + mf2)**1.5*(-k**2 + (-(k*Sqrt(1 + mf2)) + (0,1)*p0)**2)) &
               + (k*nfd1xd1la)/(2.*(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) + (0,1)*p0)**2))           &
               - nfd1lf/(2.*(1 + mf2)**1.5*(-k**2 + (k*Sqrt(1 + mf2) + (0,1)*p0)**2))                &
               + (k*nfd1xd1lf)/(2.*(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) + (0,1)*p0)**2))              &
               + (k*nfd1la*(-(k*Sqrt(1 + mf2)) + (0,1)*p0))/((1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2))  & 
               + (0,1)*p0)**2)**2) - (k*nfd1lf*(k*Sqrt(1 + mf2) + (0,1)*p0))/((1 + mf2)*(-k**2       &
               + (k*Sqrt(1 + mf2) + (0,1)*p0)**2)**2)))

     dB1F2Adl=real(dB1F2AdlC)

     dB1F2AdlbC=-0.5*(k**2*(-0.5*nfd1lba/((1 + mf2)**1.5*(-k**2 + (-(k*Sqrt(1 + mf2)) + (0,1)*p0)**2))& 
                + (k*nfd1xd1lba)/(2.*(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) + (0,1)*p0)**2))          &
                - nfd1lbf/(2.*(1 + mf2)**1.5*(-k**2 + (k*Sqrt(1 + mf2) + (0,1)*p0)**2))               &
                + (k*nfd1xd1lbf)/(2.*(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) + (0,1)*p0)**2))             & 
                + (k*nfd1lba*(-(k*Sqrt(1 + mf2)) + (0,1)*p0))/((1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) &
                + (0,1)*p0)**2)**2) - (k*nfd1lbf*(k*Sqrt(1 + mf2) + (0,1)*p0))/((1 + mf2)*(-k**2      &
                + (k*Sqrt(1 + mf2) + (0,1)*p0)**2)**2)))

     dB1F2Adlb=real(dB1F2AdlbC)

     dB2F1AdlC=-0.5*(k**2*((k**2*nfd1la)/(Sqrt(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) + (0,1)*p0)**2)**2)& 
               + (k**2*nfd1lf)/(Sqrt(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) + (0,1)*p0)**2)**2)))

     dB2F1Adl=real(dB2F1AdlC)

     dB2F1AdlbC=-0.5*(k**2*((k**2*nfd1lba)/(Sqrt(1 + mf2)*(-k**2 + (-(k*Sqrt(1 + mf2)) + (0,1)*p0)**2)**2)&
                + (k**2*nfd1lbf)/(Sqrt(1 + mf2)*(-k**2 + (k*Sqrt(1 + mf2) + (0,1)*p0)**2)**2)))

     dB2F1Adlb=real(dB2F1AdlbC)

     dB2F1PdlC=-0.5*(k**2*((k**2*nfd1la)/(Sqrt(1 + mf2)*(-(k**2*(1 + mp2)) + (-(k*Sqrt(1 + mf2)) + (0,1)*p0)**2)**2)& 
               + (k**2*nfd1lf)/(Sqrt(1 + mf2)*(-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) + (0,1)*p0)**2)**2)))

     dB2F1Pdl=real(dB2F1PdlC)

     dB2F1PdlbC=-0.5*(k**2*((k**2*nfd1lba)/(Sqrt(1 + mf2)*(-(k**2*(1 + mp2)) + (-(k*Sqrt(1 + mf2)) + (0,1)*p0)**2)**2)&
                + (k**2*nfd1lbf)/(Sqrt(1 + mf2)*(-(k**2*(1 + mp2)) + (k*Sqrt(1 + mf2) + (0,1)*p0)**2)**2)))

     dB2F1Pdlb=real(dB2F1PdlbC)

     dB2F1SdlC=-0.5*(k**2*((k**2*nfd1la)/(Sqrt(1 + mf2)*(-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) + (0,1)*p0)**2)**2) &
               + (k**2*nfd1lf)/(Sqrt(1 + mf2)*(-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) + (0,1)*p0)**2)**2)))

     dB2F1Sdl=real(dB2F1SdlC)

     dB2F1SdlbC=-0.5*(k**2*((k**2*nfd1lba)/(Sqrt(1 + mf2)*(-(k**2*(1 + ms2)) + (-(k*Sqrt(1 + mf2)) + (0,1)*p0)**2)**2)&
                + (k**2*nfd1lbf)/(Sqrt(1 + mf2)*(-(k**2*(1 + ms2)) + (k*Sqrt(1 + mf2) + (0,1)*p0)**2)**2)))

     dB2F1Sdlb=real(dB2F1SdlbC)
     
     BF_POL_COM(1)=dB1F1Adl
     BF_POL_COM(2)=dB1F1AdlB
     BF_POL_COM(3)=dB1F2Adl
     BF_POL_COM(4)=dB1F2AdlB
     BF_POL_COM(5)=dB2F1Adl
     BF_POL_COM(6)=dB2F1AdlB
     BF_POL_COM(7)=dB2F1Pdl
     BF_POL_COM(8)=dB2F1PdlB
     BF_POL_COM(9)=dB2F1Sdl
     BF_POL_COM(10)=dB2F1SdlB

  END SUBROUTINE BF_POL
  END MODULE BF_POL_MOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE DETAPSIDLLB_MOD
  CONTAINS
  SUBROUTINE DETAPSIDLLB(MP2,MS2,MF2,L,LB,K,H,G,ETAPHI,ETAPSI,ETAPSID1L,ETAPSID1LB)
  USE TMU_MOD
  USE CONST_SAVE , ONLY : NF,NC,V3,C2NC,PI
  USE BF_POL_MOD
  IMPLICIT NONE

  REAL(16),INTENT(OUT) :: ETAPSID1L,ETAPSID1LB
  REAL(16),INTENT(IN) :: MP2,MS2,MF2,K,L,LB,ETAPSI,ETAPHI,H,G
  REAL(16) :: ETAA
  REAL(16) :: dB2F1Adl,dB2F1Adlb,dB1F1Adl,dB1F1Adlb
  REAL(16) :: dB1F2Adl,dB1F2Adlb,dB2F1Sdl,dB2F1Sdlb,dB2F1Pdl,dB2F1Pdlb
  REAL(16),DIMENSION(10) :: BF_POL_COM

  ETAA=0.Q+00

  CALL BF_Pol(MP2,MS2,MF2,K,L,LB,BF_POL_COM)

  dB1F1Adl=BF_POL_COM(1)
  dB1F1AdlB=BF_POL_COM(2)
  dB1F2Adl=BF_POL_COM(3)
  dB1F2AdlB=BF_POL_COM(4)
  dB2F1Adl=BF_POL_COM(5)
  dB2F1AdlB=BF_POL_COM(6)
  dB2F1Pdl=BF_POL_COM(7)
  dB2F1PdlB=BF_POL_COM(8)
  dB2F1Sdl=BF_POL_COM(9)
  dB2F1SdlB=BF_POL_COM(10)
  
  ETAPSID1L=(C2NC*(2*dB2F1Adl*(4.Q+00 - ETAA) + 3*(dB1F1Adl - 2*dB1F2Adl)*(3.Q+00 - ETAPSI))*G**2*V3)/12.Q+00 &
            + ((4.Q+00 - ETAPSI)*H**2*(dB2F1Sdl + dB2F1Pdl*(-1.Q+00 + NF**2))*V3)/(12.Q+00*NF)
            
  ETAPSID1LB=(C2NC*(2*dB2F1Adlb*(4.Q+00 - ETAA) + 3*(dB1F1Adlb - 2*dB1F2Adlb)*(3.Q+00 - ETAPSI))*G**2*V3)/12.Q+00 &
             + ((4.Q+00 - ETAPSI)*H**2*(dB2F1Sdlb + dB2F1Pdlb*(-1.Q+00 + NF**2))*V3)/(12.Q+00*NF)

  END SUBROUTINE DETAPSIDLLB
  END MODULE DETAPSIDLLB_MOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE F_POL_MOD
    CONTAINS
    SUBROUTINE F_POL(MF2,K,L_CON,LB_CON,f2ad1l,f2ad1lb,f3ad1l,f3ad1lb)
    USE TMU_MOD
    USE FNF_MOD
    IMPLICIT NONE
  
    REAL(16),INTENT(IN) :: MF2,K,L_CON,LB_CON
    REAL(16),INTENT(OUT) :: f2ad1l,f2ad1lb,f3ad1l,f3ad1lb
    REAL(16) :: nfd1lf,nfd1lbf,nfd1la,nfd1lba 
    REAL(16) :: nfd1xd1lf,nfd1xd1lbf,nfd1xd1la,nfd1xd1lba,XFF,XFA
    REAL(16) :: nf0,nf1,nf2,L,LB,NFD2XD1LF,NFD2XD1LBF,NFD2XD1LA,NFD2XD1LBA
   
    L=L_CON
    LB=LB_CON
    XFF=K*SQRT(1.Q+00 + MF2) - MU
    CALL FNF012(XFF,T,L,LB,NF0,NF1,NF2)
   
    nfd1lf=-(nf2*(-1 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
    nfd1lbf=-0.5*(nf1*(-2 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
    nfd1xd1lf=-0.25*(nf2*(4 + 36*(2 + l*lb)*nf0**2 + 9*lb**2*nf1**2 - 36*l*nf2 + &
              72*l**2*nf2**2 + 6*nf0*(-10 - 6*l*lb + 18*lb*nf1 + 9*l*lb**2*nf1 + &
              24*l*nf2 + 18*l**2*lb*nf2 + 6*lb**2*nf2) +                         &
              6*lb*(-2*nf1 + 9*l*nf1*nf2 - 12*nf2**2)))/T
    nfd1xd1lbf=(6*lb*(1 - 3*nf0)*nf1**2 + 24*lb*(2 - 3*nf0)*nf0*nf2 +            &
               18*l**2*nf0*(-5 + 8*nf0 + 9*lb*nf1)*nf2 + 108*l**3*nf0*nf2**2 -   & 
               nf1*(4 - 21*nf0 + 18*nf0**2 + 72*lb**2*nf0*nf2) +                 &
               3*l*(12*nf0**3 + 6*nf0**2*(-3 + 5*lb*nf1) + 2*nf1*nf2 +           & 
               3*nf0*(2 - 7*lb*nf1 + 6*lb**2*nf1**2 - 2*nf1*nf2 - 8*lb*nf2**2)))/&
               (2.*T)
    nfd2xd1lf=-0.5*(nf2*(16*nf0**3 + 15*lb*nf0**2*nf1 - 24*(8*l - 3*lb**2)*nf0**2*nf2 +  &
              24*(6*l**2 - 20*lb - 9*l*lb**2)*nf0*nf2**2 +                               &
              3*(31 + 54*l*lb - 9*lb**3)*nf0**2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2) -      &
              3*(26*l + 9*l**2*lb - 69*lb**2)*nf0*nf2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2) +& 
              3*(3*l**2 - 23*lb)*nf2**2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2) +              &
              30*nf0*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2)**2 +                              &
              6*l*nf2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2)**2 +                             &
              (-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2)**3/4.))/T**2
    nfd2xd1lbf=-0.5*(nf1*(nf0**3 - 6*lb*nf0**2*nf1 - 3*(23*l - 3*lb**2)*nf0**2*nf2 +     &
               3*(69*l**2 - 26*lb - 9*l*lb**2)*nf0*nf2**2 +                              &
               30*nf0**2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2) -                             &
               6*(20*l + 9*l**2*lb - 6*lb**2)*nf0*nf2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2) +& 
               6*(3*l**2 - 8*lb)*nf2**2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2) -              &
               (3*(-31 + 9*l**3 - 54*l*lb)*nf0*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2)**2)/4. -& 
               (15*l*nf2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2)**2)/4. +                      &
               (-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2)**3))/T**2
  
    L=LB_CON
    LB=L_CON
    XFA=K*SQRT(1.Q+00 + MF2) + MU
    CALL FNF012(XFA,T,L,LB,NF0,NF1,NF2)
  
    nfd1lba=-(nf2*(-1 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
    nfd1la=-0.5*(nf1*(-2 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
    nfd1xd1lba=-0.25*(nf2*(4 + 36*(2 + l*lb)*nf0**2 + 9*lb**2*nf1**2 - 36*l*nf2 + &
              72*l**2*nf2**2 + 6*nf0*(-10 - 6*l*lb + 18*lb*nf1 + 9*l*lb**2*nf1 + &
              24*l*nf2 + 18*l**2*lb*nf2 + 6*lb**2*nf2) +                         &
              6*lb*(-2*nf1 + 9*l*nf1*nf2 - 12*nf2**2)))/T
    nfd1xd1la=(6*lb*(1 - 3*nf0)*nf1**2 + 24*lb*(2 - 3*nf0)*nf0*nf2 +           &
              18*l**2*nf0*(-5 + 8*nf0 + 9*lb*nf1)*nf2 + 108*l**3*nf0*nf2**2 -   & 
              nf1*(4 - 21*nf0 + 18*nf0**2 + 72*lb**2*nf0*nf2) +                 &
              3*l*(12*nf0**3 + 6*nf0**2*(-3 + 5*lb*nf1) + 2*nf1*nf2 +           & 
              3*nf0*(2 - 7*lb*nf1 + 6*lb**2*nf1**2 - 2*nf1*nf2 - 8*lb*nf2**2)))/&
              (2.*T)
    nfd2xd1lba=-0.5*(nf2*(16*nf0**3 + 15*lb*nf0**2*nf1 - 24*(8*l - 3*lb**2)*nf0**2*nf2 +  &
              24*(6*l**2 - 20*lb - 9*l*lb**2)*nf0*nf2**2 +                               &
              3*(31 + 54*l*lb - 9*lb**3)*nf0**2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2) -      &
              3*(26*l + 9*l**2*lb - 69*lb**2)*nf0*nf2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2) +& 
              3*(3*l**2 - 23*lb)*nf2**2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2) +              &
              30*nf0*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2)**2 +                              &
              6*l*nf2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2)**2 +                             &
              (-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2)**3/4.))/T**2
    nfd2xd1la=-0.5*(nf1*(nf0**3 - 6*lb*nf0**2*nf1 - 3*(23*l - 3*lb**2)*nf0**2*nf2 +     &
               3*(69*l**2 - 26*lb - 9*l*lb**2)*nf0*nf2**2 +                              &
               30*nf0**2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2) -                             &
               6*(20*l + 9*l**2*lb - 6*lb**2)*nf0*nf2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2) +& 
               6*(3*l**2 - 8*lb)*nf2**2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2) -              &
               (3*(-31 + 9*l**3 - 54*l*lb)*nf0*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2)**2)/4. -& 
               (15*l*nf2*(-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2)**2)/4. +                      &
               (-2 + 2*nf0 + 3*lb*nf1 + 6*l*nf2)**3))/T**2
  
    f2ad1l=-0.25*(-nfd1la - nfd1lf)/(1 + mf2)**1.5 + (-0.5*(k*nfd1xd1la)/Sqrt(1 + mf2)     &
           - (k*nfd1xd1lf)/(2.*Sqrt(1 + mf2)))/(2.*Sqrt(1 + mf2))
    f2ad1lb=-0.25*(-nfd1lba - nfd1lbf)/(1 + mf2)**1.5 + (-0.5*(k*nfd1xd1lba)/Sqrt(1 + mf2) &
            - (k*nfd1xd1lbf)/(2.*Sqrt(1 + mf2)))/(2.*Sqrt(1 + mf2))
    f3ad1l=((-3*(-nfd1la - nfd1lf))/(8.*(1 + mf2)**2.5) + (-0.5*(k*nfd1xd1la)/Sqrt(1 + mf2)            &
           - (k*nfd1xd1lf)/(2.*Sqrt(1 + mf2)))/(2.*(1 + mf2)**1.5) - ((k*nfd1xd1la)/(4.*(1 + mf2)**1.5)& 
           + (k*nfd1xd1lf)/(4.*(1 + mf2)**1.5) - (k**2*nfd2xd1la)/(4.*(1 + mf2))                       &
           - (k**2*nfd2xd1lf)/(4.*(1 + mf2)))/(2.*Sqrt(1 + mf2)))/2.
    f3ad1lb=((-3*(-nfd1lba - nfd1lbf))/(8.*(1 + mf2)**2.5) + (-0.5*(k*nfd1xd1lba)/Sqrt(1 + mf2)           &
            - (k*nfd1xd1lbf)/(2.*Sqrt(1 + mf2)))/(2.*(1 + mf2)**1.5) - ((k*nfd1xd1lba)/(4.*(1 + mf2)**1.5)& 
            + (k*nfd1xd1lbf)/(4.*(1 + mf2)**1.5) - (k**2*nfd2xd1lba)/(4.*(1 + mf2))                       &
            - (k**2*nfd2xd1lbf)/(4.*(1 + mf2)))/(2.*Sqrt(1 + mf2)))/2.
  
    END SUBROUTINE F_POL
    END MODULE F_POL_MOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE DETAPHIDLLB_MOD
  CONTAINS
  SUBROUTINE DETAPHIDLLB(MF2,L,LB,K,H,ETAPSI,ETAPHID1L,ETAPHID1LB)
  USE TMU_MOD
  USE CONST_SAVE , ONLY : NF,NC,V3,C2NC,PI
  USE F_POL_MOD
  IMPLICIT NONE
  
  REAL(16),INTENT(IN) :: MF2,K,L,LB,ETAPSI,H
  REAL(16),INTENT(OUT) :: ETAPHID1L,ETAPHID1LB
  REAL(16) f2ad1l,f2ad1lb,f3ad1l,f3ad1lb

  CALL F_POL(MF2,K,L,LB,f2ad1l,f2ad1lb,f3ad1l,f3ad1lb)
    
  ETAPHID1L=(-4*(-2.Q+00 + ETAPSI)*f3ad1l + (-3.Q+00 + 2*ETAPSI)*f2ad1l*H**2*NC)/(6.*PI**2)
              
  ETAPHID1LB=(-4*(-2.Q+00 + ETAPSI)*f3ad1lb + (-3.Q+00 + 2*ETAPSI)*f2ad1lb*H**2*NC)/(6.*PI**2)
  
  END SUBROUTINE DETAPHIDLLB
  END MODULE DETAPHIDLLB_MOD