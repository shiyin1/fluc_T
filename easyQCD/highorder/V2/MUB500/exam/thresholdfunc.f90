!Rui Wen 2021.02.01
MODULE F_THR_MOD
CONTAINS
SUBROUTINE F_THR(MF2,K,F2A,F3A,NFFDNX,NFADNX)
!CALCULATING THE RIGHT HAND SIDE OF DIFFERENTIAL EQUATIONS

  IMPLICIT NONE

  REAL(16),INTENT(IN) :: MF2,K
  REAL(16),INTENT(IN) :: NFFDNX(:),NFADNX(:)
  REAL(16) :: NFF,NFD1XF,NFD2XF,NFA,NFD1XA,NFD2XA
  REAL(16),INTENT(OUT) :: F2A,F3A

  NFF   =NFFDNX(1)
  NFD1XF=NFFDNX(2)
  NFD2XF=NFFDNX(3)

  NFA   =NFADNX(1)
  NFD1XA=NFADNX(2)
  NFD2XA=NFADNX(3)

  F2A=(1 - NFA + K*SQRT(1 + MF2)*NFD1XA + K*SQRT(1 + MF2)*NFD1XF - NFF)/        &
  (4.Q+0*(1 + MF2)**1.5)

  F3A=(-(K**2*(1 + MF2)*(NFD2XA + NFD2XF)) +                                    &
    3*(1 - NFA + K*SQRT(1 + MF2)*NFD1XA + K*SQRT(1 + MF2)*NFD1XF - NFF))/       &
  (16.Q+0*(1 + MF2)**2.5)

END SUBROUTINE F_THR

SUBROUTINE FT_THR(MF2,K,F2AT,F3AT,NFFDNX,NFADNX)

  IMPLICIT NONE

  REAL(16),INTENT(IN) :: MF2,K
  REAL(16),INTENT(IN) :: NFFDNX(6),NFADNX(6)
  REAL(16) :: NFF,NFD1XF,NFD2XF,NFA,NFD1XA,NFD2XA
  REAL(16),INTENT(OUT) :: F2AT,F3AT

  NFF   =NFFDNX(1)
  NFD1XF=NFFDNX(2)
  NFD2XF=NFFDNX(3)

  NFA   =NFADNX(1)
  NFD1XA=NFADNX(2)
  NFD2XA=NFADNX(3)

  F2AT=(-NFA + K*SQRT(1 + MF2)*(NFD1XA + NFD1XF) - NFF)/(4.Q+0*(1 + MF2)**1.5)

  F3AT=(-3*NFA - K*(-3*SQRT(1 + MF2)*NFD1XA - 3*SQRT(1 + MF2)*NFD1XF +           &
       K*(1 + MF2)*(NFD2XA + NFD2XF)) - 3*NFF)/(16.Q+0*(1 + MF2)**2.5)

END SUBROUTINE FT_THR
END MODULE F_THR_MOD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BB_THR_MOD
CONTAINS
SUBROUTINE BB_THR(MB2A,MB2B,K,FINVEBADN,FINVEBBDN,NBBADN,NBBBDN,B2B2)
!CALCULATING THE RIGHT HAND SIDE OF DIFFERENTIAL EQUATIONS

  IMPLICIT NONE

  REAL(16) MB2A,MB2B,K
  REAL(16) NBBA,NBBB,NBBAD1,NBBBD1
  REAL(16) FINVEBADN(2),FINVEBBDN(2),NBBADN(2),NBBBDN(2)
  REAL(16) FINVEBA,FINVEBB,FINVEBAD1,FINVEBBD1
  REAL(16) FABMBBABBD1AD0B,FABMBBABBD0AD1B,FABMBBABBD1AD1B
  REAL(16) B2B2

  NBBA  =NBBADN(1)
  NBBAD1=NBBADN(2)
  NBBB  =NBBBDN(1)
  NBBBD1=NBBBDN(2)
  
  FINVEBA  =FINVEBADN(1)
  FINVEBAD1=FINVEBADN(2)
  FINVEBB  =FINVEBBDN(1)
  FINVEBBD1=FINVEBBDN(2)

  FABMBBABBD1AD0B=-(1/(K**2*(MB2A - MB2B)**2))
  FABMBBABBD0AD1B=1/(K**2*(MB2A - MB2B)**2)
  FABMBBABBD1AD1B=-2/(K**2*(MB2A - MB2B)**3)

  B2B2=-(K**3*(FABMBBABBD0AD1B*(FINVEBAD1 + 2*FINVEBAD1*NBBA + 2*FINVEBA*NBBAD1) + &
       FABMBBABBD1AD1B*(FINVEBA + 2*FINVEBA*NBBA - FINVEBB*(1 + 2*NBBB)) -      &
       FABMBBABBD1AD0B*(FINVEBBD1 + 2*FINVEBBD1*NBBB + 2*FINVEBB*NBBBD1)))/2.Q+0

END SUBROUTINE BB_THR
END MODULE BB_THR_MOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BF_THR_MOD
CONTAINS
SUBROUTINE BF_MES_THR(MP2,MS2,MF2,T,K,NBPIDNX,NBSGDNX,NFFDNX,NFADNX,BF_MES_THR_COM)
!CALCULATING THE RIGHT HAND SIDE OF DIFFERENTIAL EQUATIONS

  USE MU0P0_COM
  IMPLICIT NONE

  REAL(16) MP2,MS2,MF2,T,K
  REAL(16),INTENT(OUT),DIMENSION(11) :: BF_MES_THR_COM
  REAL(16) NBPIDNX(6),NBSGDNX(6)
  REAL(16) NFFDNX(6),NFADNX(6)
  REAL(16) NFF,NFD1XF,NFD2XF,NFD3XF,NFD4XF,NFD5XF,NFA,NFD1XA,NFD2XA,NFD3XA,NFD4XA,NFD5XA
  REAL(16) NBPION,NBD1XPION,NBD2XPION,NBD3XPION,NBD4XPION,NBD5XPION
  REAL(16) NBSIGMA,NBD1XSIGMA,NBD2XSIGMA,NBD3XSIGMA,NBD4XSIGMA,NBD5XSIGMA

  REAL(16) B1F2P,B1F2S,B1F3P,B1F3S,B2F1AP,B2F1AS,B2F2AP,B2F2AS,B3F1AP,B2F3AP,B3F2AP
  COMPLEX(16) B1F2PC,B1F2SC,B1F3PC,B1F3SC,B2F1APC,B2F1ASC,B2F2APC,B2F2ASC,B3F1APC,B2F3APC,B3F2APC

  NBPION   =NBPIDNX(1)
  NBD1XPION=NBPIDNX(2)
  NBD2XPION=NBPIDNX(3)
  NBD3XPION=NBPIDNX(4)
  NBD4XPION=NBPIDNX(5)
  NBD5XPION=NBPIDNX(6)

  NBSIGMA   =NBSGDNX(1)
  NBD1XSIGMA=NBSGDNX(2)
  NBD2XSIGMA=NBSGDNX(3)
  NBD3XSIGMA=NBSGDNX(4)
  NBD4XSIGMA=NBSGDNX(5)
  NBD5XSIGMA=NBSGDNX(6)

  NFF   =NFFDNX(1)
  NFD1XF=NFFDNX(2)
  NFD2XF=NFFDNX(3)
  NFD3XF=NFFDNX(4)
  NFD4XF=NFFDNX(5)
  NFD5XF=NFFDNX(6)

  NFA   =NFADNX(1)
  NFD1XA=NFADNX(2)
  NFD2XA=NFADNX(3)
  NFD3XA=NFADNX(4)
  NFD4XA=NFADNX(5)
  NFD5XA=NFADNX(6)

  B1F2PC=(K**2*NFA)/(4.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                     &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -                      &
  (K**3*NFD1XA)/                                                                &
   (4.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -                      &
  (K**3*NFD1XF)/                                                                &
   (4.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
       (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +                         &
  (K**2*NFF)/(4.*(1 + MF2)**1.5*                                                &
     (-(K**2*(1 + MP2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +      &
  (K**4*NBPION)/                                                                &
   (2.*SQRT(1 + MP2)*(-(K**2*(1 + MF2)) +                                       &
        (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**2) +                  &
  (K**4*NBPION)/                                                                &
   (2.*SQRT(1 + MP2)*(-(K**2*(1 + MF2)) +                                       &
        (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**2) -                     &
  (K**3*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                      &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                  &
  (K**3*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                         &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) -                     &
  K**2/(4.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                  &
       (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)) +                        &
  K**4/(2.*SQRT(1 + MP2)*(-(K**2*(1 + MF2)) +                                   &
        (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)**2) -                 &
  (K**3*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                            &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2)                      

  B1F2P=REAL(B1F2PC)

  B1F2SC=(K**2*NFA)/(4.*(1 + MF2)**1.5*(-(K**2*(1 + MS2)) +                     &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -                      &
  (K**3*NFD1XA)/                                                                &
   (4.*(1 + MF2)*(-(K**2*(1 + MS2)) +                                           &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -                      &
  (K**3*NFD1XF)/                                                                &
   (4.*(1 + MF2)*(-(K**2*(1 + MS2)) +                                           &
       (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +                         &
  (K**2*NFF)/(4.*(1 + MF2)**1.5*                                                &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +      &
  (K**4*NBSIGMA)/                                                               &
   (2.*SQRT(1 + MS2)*(-(K**2*(1 + MF2)) +                                       &
        (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0)**2)**2) +                  &
  (K**4*NBSIGMA)/                                                               &
   (2.*SQRT(1 + MS2)*(-(K**2*(1 + MF2)) +                                       &
        (K*SQRT(1 + MS2) - MU0 + (0,1)*P0)**2)**2) -                     &
  (K**3*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                      &
   (2.*(1 + MF2)*(-(K**2*(1 + MS2)) +                                           &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                  &
  (K**3*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                         &
   (2.*(1 + MF2)*(-(K**2*(1 + MS2)) +                                           &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) -                     &
  K**2/(4.*(1 + MF2)**1.5*(-(K**2*(1 + MS2)) +                                  &
       (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)) +                        &
  K**4/(2.*SQRT(1 + MS2)*(-(K**2*(1 + MF2)) +                                   &
        (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0C)**2)**2) -                 &
  (K**3*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                            &
   (2.*(1 + MF2)*(-(K**2*(1 + MS2)) +                                           &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2)                      

  B1F2S=REAL(B1F2SC)

  B1F3PC=(-(K**4*NFA)/(4.*(1 + MF2)**1.5*                                       &
       (-(K**2*(1 + MP2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
           2)**2) + (3*K**2*NFA)/                                               &
     (8.*(1 + MF2)**2.5*(-(K**2*(1 + MP2)) +                                    &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -                    &
    (3*K**3*NFD1XA)/                                                            &
     (8.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) +                    &
    (K**4*NFD2XA)/                                                              &
     (8.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -                    &
    (K**4*NFF)/(4.*(1 + MF2)**1.5*                                              &
       (-(K**2*(1 + MP2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**     &
        2) - (3*K**3*NFD1XF)/                                                   &
     (8.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +                       &
    (K**4*NFD2XF)/                                                              &
     (8.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +                       &
    (3*K**2*NFF)/                                                               &
     (8.*(1 + MF2)**2.5*(-(K**2*(1 + MP2)) +                                    &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) -                       &
    (K**6*NBPION)/                                                              &
     (SQRT(1 + MP2)*(-(K**2*(1 + MF2)) +                                        &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (K**6*NBPION)/                                                              &
     (SQRT(1 + MP2)*(-(K**2*(1 + MF2)) +                                        &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**3) -                   &
    (3*K**3*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                  &
     (4.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (K**4*NFD1XA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                 &
     (2.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (K**4*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)/                 &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (K**4*NFD1XF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                    &
     (2.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) +                   &
    (3*K**3*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                     &
     (4.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) +                   &
    (K**4*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)/                    &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) +                   &
    K**4/(4.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) -                  &
    (3*K**2)/(8.*(1 + MF2)**2.5*                                                &
       (-(K**2*(1 + MP2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2))     &
- K**6/(SQRT(1 + MP2)*(-(K**2*(1 + MF2)) +                                      &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)**3) -               &
    (3*K**3*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                        &
     (4.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) -                  &
    (K**4*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)/                       &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3))/2.                

  B1F3P=REAL(B1F3PC)

  B1F3SC=(-(K**4*NFA)/(4.*(1 + MF2)**1.5*                                        &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
           2)**2) + (3*K**2*NFA)/                                               &
     (8.*(1 + MF2)**2.5*(-(K**2*(1 + MS2)) +                                    &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -                    &
    (3*K**3*NFD1XA)/                                                            &
     (8.*(1 + MF2)**2*(-(K**2*(1 + MS2)) +                                      &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) +                    &
    (K**4*NFD2XA)/                                                              &
     (8.*(1 + MF2)**1.5*(-(K**2*(1 + MS2)) +                                    &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -                    &
    (K**4*NFF)/(4.*(1 + MF2)**1.5*                                              &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**     &
        2) - (3*K**3*NFD1XF)/                                                   &
     (8.*(1 + MF2)**2*(-(K**2*(1 + MS2)) +                                      &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +                       &
    (K**4*NFD2XF)/                                                              &
     (8.*(1 + MF2)**1.5*(-(K**2*(1 + MS2)) +                                    &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +                       &
    (3*K**2*NFF)/                                                               &
     (8.*(1 + MF2)**2.5*(-(K**2*(1 + MS2)) +                                    &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) -                       &
    (K**6*NBSIGMA)/                                                             &
     (SQRT(1 + MS2)*(-(K**2*(1 + MF2)) +                                        &
          (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (K**6*NBSIGMA)/                                                             &
     (SQRT(1 + MS2)*(-(K**2*(1 + MF2)) +                                        &
          (K*SQRT(1 + MS2) - MU0 + (0,1)*P0)**2)**3) -                   &
    (3*K**3*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                  &
     (4.*(1 + MF2)**2*(-(K**2*(1 + MS2)) +                                      &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (K**4*NFD1XA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                 &
     (2.*(1 + MF2)**1.5*(-(K**2*(1 + MS2)) +                                    &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (K**4*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)/                 &
     ((1 + MF2)**1.5*(-(K**2*(1 + MS2)) +                                       &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (K**4*NFD1XF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                    &
     (2.*(1 + MF2)**1.5*(-(K**2*(1 + MS2)) +                                    &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) +                   &
    (3*K**3*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                     &
     (4.*(1 + MF2)**2*(-(K**2*(1 + MS2)) +                                      &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) +                   &
    (K**4*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)/                    &
     ((1 + MF2)**1.5*(-(K**2*(1 + MS2)) +                                       &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) +                   &
    K**4/(4.*(1 + MF2)**1.5*(-(K**2*(1 + MS2)) +                                &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) -                  &
    (3*K**2)/(8.*(1 + MF2)**2.5*                                                &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2))     &
- K**6/(SQRT(1 + MS2)*(-(K**2*(1 + MF2)) +                                      &
          (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0C)**2)**3) -               &
    (3*K**3*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                        &
     (4.*(1 + MF2)**2*(-(K**2*(1 + MS2)) +                                      &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) -                  &
    (K**4*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)/                       &
     ((1 + MF2)**1.5*(-(K**2*(1 + MS2)) +                                       &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3))/2.

  B1F3S=REAL(B1F3SC)


  B2F1APC=-(K**4*NFA)/(2.*SQRT(1 + MF2)*(-(K**2*(1 + MP2)) +                    &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) -                  &
  (K**4*NFF)/(2.*SQRT(1 + MF2)*                                                 &
     (-(K**2*(1 + MP2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2)     &
+ (K**3*NBD1XPION)/                                                             &
   (4.*(1 + MP2)*(-(K**2*(1 + MF2)) +                                           &
       (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)) -                      &
  (K**2*NBPION)/                                                                &
   (4.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                      &
       (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)) +                      &
  (K**3*NBD1XPION)/                                                             &
   (4.*(1 + MP2)*(-(K**2*(1 + MF2)) +                                           &
       (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)) -                         &
  (K**2*NBPION)/                                                                &
   (4.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                      &
       (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)) +                         &
  (K**3*NBPION*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0))/                   &
   (2.*(1 + MP2)*(-(K**2*(1 + MF2)) +                                           &
        (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**2) -                  &
  (K**3*NBPION*(K*SQRT(1 + MP2) - MU0 + (0,1)*P0))/                      &
   (2.*(1 + MP2)*(-(K**2*(1 + MF2)) +                                           &
        (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**2) +                     &
  K**4/(2.*SQRT(1 + MF2)*(-(K**2*(1 + MP2)) +                                   &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) -                    &
  K**2/(4.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                  &
       (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)) +                     &
  (K**3*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C))/                         &
   (2.*(1 + MP2)*(-(K**2*(1 + MF2)) +                                           &
        (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)**2)

  B2F1AP=REAL(B2F1APC)


  B2F1ASC=-(K**4*NFA)/(2.*SQRT(1 + MF2)*(-(K**2*(1 + MS2)) +                    &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) -                  &
  (K**4*NFF)/(2.*SQRT(1 + MF2)*                                                 &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2)     &
+ (K**3*NBD1XSIGMA)/                                                            &
   (4.*(1 + MS2)*(-(K**2*(1 + MF2)) +                                           &
       (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0)**2)) -                      &
  (K**2*NBSIGMA)/                                                               &
   (4.*(1 + MS2)**1.5*(-(K**2*(1 + MF2)) +                                      &
       (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0)**2)) +                      &
  (K**3*NBD1XSIGMA)/                                                            &
   (4.*(1 + MS2)*(-(K**2*(1 + MF2)) +                                           &
       (K*SQRT(1 + MS2) - MU0 + (0,1)*P0)**2)) -                         &
  (K**2*NBSIGMA)/                                                               &
   (4.*(1 + MS2)**1.5*(-(K**2*(1 + MF2)) +                                      &
       (K*SQRT(1 + MS2) - MU0 + (0,1)*P0)**2)) +                         &
  (K**3*NBSIGMA*(-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0))/                  &
   (2.*(1 + MS2)*(-(K**2*(1 + MF2)) +                                           &
        (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0)**2)**2) -                  &
  (K**3*NBSIGMA*(K*SQRT(1 + MS2) - MU0 + (0,1)*P0))/                     &
   (2.*(1 + MS2)*(-(K**2*(1 + MF2)) +                                           &
        (K*SQRT(1 + MS2) - MU0 + (0,1)*P0)**2)**2) +                     &
  K**4/(2.*SQRT(1 + MF2)*(-(K**2*(1 + MS2)) +                                   &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) -                    &
  K**2/(4.*(1 + MS2)**1.5*(-(K**2*(1 + MF2)) +                                  &
       (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0C)**2)) +                     &
  (K**3*(-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0C))/                         &
   (2.*(1 + MS2)*(-(K**2*(1 + MF2)) +                                           &
        (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0C)**2)**2)

  B2F1AS=REAL(B2F1ASC)


  B2F2APC=-(K**4*NFA)/(4.*(1 + MF2)**1.5*                                       &
     (-(K**2*(1 + MP2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**    &
      2) + (K**5*NFD1XA)/                                                       &
   (4.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                  &
  (K**5*NFD1XF)/                                                                &
   (4.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) -                     &
  (K**4*NFF)/(4.*(1 + MF2)**1.5*                                                &
     (-(K**2*(1 + MP2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2)     &
- (K**5*NBD1XPION)/                                                             &
   (4.*(1 + MP2)*(-(K**2*(1 + MF2)) +                                           &
        (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**2) +                  &
  (K**4*NBPION)/                                                                &
   (4.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                      &
        (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**2) -                  &
  (K**5*NBD1XPION)/                                                             &
   (4.*(1 + MP2)*(-(K**2*(1 + MF2)) +                                           &
        (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**2) +                     &
  (K**4*NBPION)/                                                                &
   (4.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                      &
        (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**2) +                     &
  (K**5*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                      &
   ((1 + MF2)*(-(K**2*(1 + MP2)) +                                              &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) -                  &
  (K**5*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                         &
   ((1 + MF2)*(-(K**2*(1 + MP2)) +                                              &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) -                     &
  (K**5*NBPION*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0))/                   &
   ((1 + MP2)*(-(K**2*(1 + MF2)) +                                              &
        (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**3) +                  &
  (K**5*NBPION*(K*SQRT(1 + MP2) - MU0 + (0,1)*P0))/                      &
   ((1 + MP2)*(-(K**2*(1 + MF2)) +                                              &
        (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**3) +                     &
  K**4/(4.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                  &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) +                    &
  K**4/(4.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                  &
        (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)**2) +                 &
  (K**5*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                            &
   ((1 + MF2)*(-(K**2*(1 + MP2)) +                                              &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3) -                    &
  (K**5*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C))/                         &
   ((1 + MP2)*(-(K**2*(1 + MF2)) +                                              &
        (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)**3)

  B2F2AP=REAL(B2F2APC)

  B2F2ASC=-(K**4*NFA)/(4.*(1 + MF2)**1.5*                                       &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**    &
      2) + (K**5*NFD1XA)/                                                       &
   (4.*(1 + MF2)*(-(K**2*(1 + MS2)) +                                           &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                  &
  (K**5*NFD1XF)/                                                                &
   (4.*(1 + MF2)*(-(K**2*(1 + MS2)) +                                           &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) -                     &
  (K**4*NFF)/(4.*(1 + MF2)**1.5*                                                &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2)     &
- (K**5*NBD1XSIGMA)/                                                            &
   (4.*(1 + MS2)*(-(K**2*(1 + MF2)) +                                           &
        (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0)**2)**2) +                  &
  (K**4*NBSIGMA)/                                                               &
   (4.*(1 + MS2)**1.5*(-(K**2*(1 + MF2)) +                                      &
        (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0)**2)**2) -                  &
  (K**5*NBD1XSIGMA)/                                                            &
   (4.*(1 + MS2)*(-(K**2*(1 + MF2)) +                                           &
        (K*SQRT(1 + MS2) - MU0 + (0,1)*P0)**2)**2) +                     &
  (K**4*NBSIGMA)/                                                               &
   (4.*(1 + MS2)**1.5*(-(K**2*(1 + MF2)) +                                      &
        (K*SQRT(1 + MS2) - MU0 + (0,1)*P0)**2)**2) +                     &
  (K**5*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                      &
   ((1 + MF2)*(-(K**2*(1 + MS2)) +                                              &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) -                  &
  (K**5*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                         &
   ((1 + MF2)*(-(K**2*(1 + MS2)) +                                              &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) -                     &
  (K**5*NBSIGMA*(-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0))/                  &
   ((1 + MS2)*(-(K**2*(1 + MF2)) +                                              &
        (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0)**2)**3) +                  &
  (K**5*NBSIGMA*(K*SQRT(1 + MS2) - MU0 + (0,1)*P0))/                     &
   ((1 + MS2)*(-(K**2*(1 + MF2)) +                                              &
        (K*SQRT(1 + MS2) - MU0 + (0,1)*P0)**2)**3) +                     &
  K**4/(4.*(1 + MF2)**1.5*(-(K**2*(1 + MS2)) +                                  &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) +                    &
  K**4/(4.*(1 + MS2)**1.5*(-(K**2*(1 + MF2)) +                                  &
        (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0C)**2)**2) +                 &
  (K**5*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                            &
   ((1 + MF2)*(-(K**2*(1 + MS2)) +                                              &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3) -                    &
  (K**5*(-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0C))/                         &
   ((1 + MS2)*(-(K**2*(1 + MF2)) +                                              &
        (-(K*SQRT(1 + MS2)) - MU0 + (0,1)*P0C)**2)**3)

  B2F2AS=REAL(B2F2ASC)


  B3F1APC=((K**6*NFA)/(SQRT(1 + MF2)*(-(K**2*(1 + MP2)) +                       &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) +                &
    (K**6*NFF)/(SQRT(1 + MF2)*(-(K**2*(1 + MP2)) +                              &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) +                   &
    (K**4*NBPION)/                                                              &
     (4.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                    &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (3*K**3*NBD1XPION)/                                                         &
     (8.*(1 + MP2)**2*(-(K**2*(1 + MF2)) +                                      &
         (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)) -                    &
    (K**4*NBD2XPION)/                                                           &
     (8.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                    &
         (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)) -                    &
    (3*K**2*NBPION)/                                                            &
     (8.*(1 + MP2)**2.5*(-(K**2*(1 + MF2)) +                                    &
         (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)) +                    &
    (K**4*NBPION)/                                                              &
     (4.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                    &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**2) +                   &
    (3*K**3*NBD1XPION)/                                                         &
     (8.*(1 + MP2)**2*(-(K**2*(1 + MF2)) +                                      &
         (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)) -                       &
    (K**4*NBD2XPION)/                                                           &
     (8.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                    &
         (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)) -                       &
    (3*K**2*NBPION)/                                                            &
     (8.*(1 + MP2)**2.5*(-(K**2*(1 + MF2)) +                                    &
         (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)) -                       &
    (K**4*NBD1XPION*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0))/              &
     (2.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                    &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (3*K**3*NBPION*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0))/               &
     (4.*(1 + MP2)**2*(-(K**2*(1 + MF2)) +                                      &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**2) -                &
    (K**4*NBPION*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)/              &
     ((1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                       &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**3) +                &
    (K**4*NBD1XPION*(K*SQRT(1 + MP2) - MU0 + (0,1)*P0))/                 &
     (2.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                    &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**2) -                   &
    (3*K**3*NBPION*(K*SQRT(1 + MP2) - MU0 + (0,1)*P0))/                  &
     (4.*(1 + MP2)**2*(-(K**2*(1 + MF2)) +                                      &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**2) -                   &
    (K**4*NBPION*(K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)/                 &
     ((1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                       &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**3) -                   &
    K**6/(SQRT(1 + MF2)*(-(K**2*(1 + MP2)) +                                    &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3) +                  &
    K**4/(4.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)**2) -               &
    (3*K**2)/(8.*(1 + MP2)**2.5*                                                &
       (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**     &
          2)) + (3*K**3*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C))/         &
     (4.*(1 + MP2)**2*(-(K**2*(1 + MF2)) +                                      &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)**2) -               &
    (K**4*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)/                    &
     ((1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                       &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)**3))/2.

  B3F1AP=REAL(B3F1APC)
  
  B2F3APC=((K**6*NFA)/(2.*(1 + MF2)**1.5*                                       &
       (-(K**2*(1 + MP2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
           2)**3) - (3*K**4*NFA)/                                               &
     (8.*(1 + MF2)**2.5*(-(K**2*(1 + MP2)) +                                    &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (3*K**5*NFD1XA)/                                                            &
     (8.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) -                &
    (K**6*NFD2XA)/                                                              &
     (8.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (K**6*NFF)/(2.*(1 + MF2)**1.5*                                              &
       (-(K**2*(1 + MP2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**     &
        3) + (3*K**5*NFD1XF)/                                                   &
     (8.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) -                   &
    (K**6*NFD2XF)/                                                              &
     (8.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) -                   &
    (3*K**4*NFF)/                                                               &
     (8.*(1 + MF2)**2.5*(-(K**2*(1 + MP2)) +                                    &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) +                   &
    (K**7*NBD1XPION)/                                                           &
     (2.*(1 + MP2)*(-(K**2*(1 + MF2)) +                                         &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (K**6*NBPION)/                                                              &
     (2.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                    &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**3) +                &
    (K**7*NBD1XPION)/                                                           &
     (2.*(1 + MP2)*(-(K**2*(1 + MF2)) +                                         &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**3) -                   &
    (K**6*NBPION)/                                                              &
     (2.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                    &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**3) +                   &
    (3*K**5*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                  &
     (2.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (K**6*NFD1XA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                 &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (3*K**6*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)/               &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**4) +                &
    (K**6*NFD1XF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                    &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) -                   &
    (3*K**5*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                     &
     (2.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) -                   &
    (3*K**6*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)/                  &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**4) +                   &
    (3*K**7*NBPION*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0))/               &
     ((1 + MP2)*(-(K**2*(1 + MF2)) +                                            &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**4) -                &
    (3*K**7*NBPION*(K*SQRT(1 + MP2) - MU0 + (0,1)*P0))/                  &
     ((1 + MP2)*(-(K**2*(1 + MF2)) +                                            &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**4) -                   &
    K**6/(2.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3) +                  &
    (3*K**4)/(8.*(1 + MF2)**2.5*                                                &
       (-(K**2*(1 + MP2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**        &
           2)**2) - K**6/                                                       &
     (2.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                    &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)**3) +               &
    (3*K**5*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                        &
     (2.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3) +                  &
    (3*K**6*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)/                     &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**4) +                  &
    (3*K**7*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C))/                     &
     ((1 + MP2)*(-(K**2*(1 + MF2)) +                                            &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)**4))/2.

  B2F3AP=REAL(B2F3APC)
  
  B3F2APC=((K**6*NFA)/(2.*(1 + MF2)**1.5*                                       &
       (-(K**2*(1 + MP2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
           2)**3) - (K**7*NFD1XA)/                                              &
     (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                         &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (K**7*NFD1XF)/                                                              &
     (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                         &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) +                   &
    (K**6*NFF)/(2.*(1 + MF2)**1.5*                                              &
       (-(K**2*(1 + MP2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**     &
        3) - (K**6*NBPION)/                                                     &
     (2.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                    &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (3*K**5*NBD1XPION)/                                                         &
     (8.*(1 + MP2)**2*(-(K**2*(1 + MF2)) +                                      &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (K**6*NBD2XPION)/                                                           &
     (8.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                    &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (3*K**4*NBPION)/                                                            &
     (8.*(1 + MP2)**2.5*(-(K**2*(1 + MF2)) +                                    &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**2) -                &
    (K**6*NBPION)/                                                              &
     (2.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                    &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**3) -                   &
    (3*K**5*NBD1XPION)/                                                         &
     (8.*(1 + MP2)**2*(-(K**2*(1 + MF2)) +                                      &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**2) +                   &
    (K**6*NBD2XPION)/                                                           &
     (8.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                    &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**2) +                   &
    (3*K**4*NBPION)/                                                            &
     (8.*(1 + MP2)**2.5*(-(K**2*(1 + MF2)) +                                    &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**2) -                   &
    (3*K**7*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                  &
     ((1 + MF2)*(-(K**2*(1 + MP2)) +                                            &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**4) +                &
    (3*K**7*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                     &
     ((1 + MF2)*(-(K**2*(1 + MP2)) +                                            &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**4) +                   &
    (K**6*NBD1XPION*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0))/              &
     ((1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                       &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (3*K**5*NBPION*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0))/               &
     (2.*(1 + MP2)**2*(-(K**2*(1 + MF2)) +                                      &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**3) +                &
    (3*K**6*NBPION*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)/            &
     ((1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                       &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0)**2)**4) -                &
    (K**6*NBD1XPION*(K*SQRT(1 + MP2) - MU0 + (0,1)*P0))/                 &
     ((1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                       &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**3) +                   &
    (3*K**5*NBPION*(K*SQRT(1 + MP2) - MU0 + (0,1)*P0))/                  &
     (2.*(1 + MP2)**2*(-(K**2*(1 + MF2)) +                                      &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**3) +                   &
    (3*K**6*NBPION*(K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)/               &
     ((1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                       &
          (K*SQRT(1 + MP2) - MU0 + (0,1)*P0)**2)**4) -                   &
    K**6/(2.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3) -                  &
    K**6/(2.*(1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)**3) +               &
    (3*K**4)/(8.*(1 + MP2)**2.5*                                                &
       (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) - MU0 +                         &
             (0,1)*P0C)**2)**2) -                                        &
    (3*K**7*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                        &
     ((1 + MF2)*(-(K**2*(1 + MP2)) +                                            &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**4) -                  &
    (3*K**5*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C))/                     &
     (2.*(1 + MP2)**2*(-(K**2*(1 + MF2)) +                                      &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)**3) +               &
    (3*K**6*(-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)/                  &
     ((1 + MP2)**1.5*(-(K**2*(1 + MF2)) +                                       &
          (-(K*SQRT(1 + MP2)) - MU0 + (0,1)*P0C)**2)**4))/2.Q+0

  B3F2AP=REAL(B3F2APC)
  BF_MES_THR_COM(1)=B1F2P
  BF_MES_THR_COM(2)=B1F2S
  BF_MES_THR_COM(3)=B1F3P
  BF_MES_THR_COM(4)=B1F3S
  BF_MES_THR_COM(5)=B2F1AP
  BF_MES_THR_COM(6)=B2F1AS
  BF_MES_THR_COM(7)=B2F2AP
  BF_MES_THR_COM(8)=B2F2AS
  BF_MES_THR_COM(9)=B3F1AP
  BF_MES_THR_COM(10)=B2F3AP
  BF_MES_THR_COM(11)=B3F2AP
END SUBROUTINE BF_MES_THR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BBF_THR(MP2,MS2,MF2,T,K,NBPIDNX,NBSGDNX,NFFDNX,NFADNX,BBF_THR_COM)
!CALCULATING THE RIGHT HAND SIDE OF DIFFERENTIAL EQUATIONS

  USE MU0P0_COM
  IMPLICIT NONE

  REAL(16) MP2,MS2,MF2,T,K
  REAL(16) NBPIDNX(6),NBSGDNX(6)
  REAL(16) NFFDNX(6),NFADNX(6)
  REAL(16) NBPION,NBD1XPION,NBD2XPION,NBD3XPION,NBD4XPION,NBD5XPION,&
          NBSIGMA,NBD1XSIGMA,NBD2XSIGMA,NBD3XSIGMA,NBD4XSIGMA,NBD5XSIGMA
  REAL(16) NFF,NFD1XF,NFD2XF,NFD3XF,NFD4XF,NFD5XF,NFA,NFD1XA,NFD2XA,NFD3XA,NFD4XA,NFD5XA
  
  COMPLEX(16) B2B1F1APSC,B2B1F1ASPC,B1B1F2APSC,B1B1F3APSC,B2B1F2APSC,B2B1F2ASPC,B2F3APC,B2F3ASC,B3F2APC,B3F2ASC
  REAL(16) B2B1F1APS,B2B1F1ASP,B1B1F2APS,B1B1F3APS,B2B1F2APS,B2B1F2ASP

  REAL(16),INTENT(OUT),DIMENSION(6) :: BBF_THR_COM

  NBPION   =NBPIDNX(1)
  NBD1XPION=NBPIDNX(2)
  NBD2XPION=NBPIDNX(3)
  NBD3XPION=NBPIDNX(4)
  NBD4XPION=NBPIDNX(5)
  NBD5XPION=NBPIDNX(6)

  NBSIGMA   =NBSGDNX(1)
  NBD1XSIGMA=NBSGDNX(2)
  NBD2XSIGMA=NBSGDNX(3)
  NBD3XSIGMA=NBSGDNX(4)
  NBD4XSIGMA=NBSGDNX(5)
  NBD5XSIGMA=NBSGDNX(6)

  NFF   =NFFDNX(1)
  NFD1XF=NFFDNX(2)
  NFD2XF=NFFDNX(3)
  NFD3XF=NFFDNX(4)
  NFD4XF=NFFDNX(5)
  NFD5XF=NFFDNX(6)

  NFA   =NFADNX(1)
  NFD1XA=NFADNX(2)
  NFD2XA=NFADNX(3)
  NFD3XA=NFADNX(4)
  NFD4XA=NFADNX(5)
  NFD5XA=NFADNX(6)


  B2B1F1APSC=-(K**3*NBD1XPION)/                                                 &
   (4.*(1 + MP2)*(MP2 - MS2)*(-(K**2*(1 + MF2)) +                               &
       (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0)**2)) +                      &
  (K**2*NBPION)/                                                                &
   (2.*SQRT(1 + MP2)*(MP2 - MS2)**2*                                            &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0)**2))     &
+ (K**2*NBPION)/                                                                &
   (4.*(1 + MP2)**1.5*(MP2 - MS2)*                                              &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0)**2))     &
- (K**3*NBD1XPION)/                                                             &
   (4.*(1 + MP2)*(MP2 - MS2)*(-(K**2*(1 + MF2)) +                               &
       (K*SQRT(1 + MP2) + MU0 - (0,1)*P0)**2)) +                         &
  (K**2*NBPION)/                                                                &
   (2.*SQRT(1 + MP2)*(MP2 - MS2)**2*                                            &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MP2) + MU0 - (0,1)*P0)**2)) +      &
  (K**2*NBPION)/                                                                &
   (4.*(1 + MP2)**1.5*(MP2 - MS2)*                                              &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MP2) + MU0 - (0,1)*P0)**2)) -      &
  (K**2*NBSIGMA)/                                                               &
   (2.*SQRT(1 + MS2)*(-MP2 + MS2)**2*                                           &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0)**2))     &
- (K**2*NBSIGMA)/                                                               &
   (2.*SQRT(1 + MS2)*(-MP2 + MS2)**2*                                           &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MS2) + MU0 - (0,1)*P0)**2)) +      &
  (K**6*NFA)/(2.*SQRT(1 + MF2)*                                                 &
     (-(K**2*(1 + MP2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**        &
         2)**2*(-(K**2*(1 + MS2)) +                                             &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) +                      &
  (K**6*NFF)/(2.*SQRT(1 + MF2)*                                                 &
     (-(K**2*(1 + MP2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2*     &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) -      &
  (K**3*NBPION*(-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0))/                   &
   (2.*(1 + MP2)*(MP2 - MS2)*(-(K**2*(1 + MF2)) +                               &
        (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0)**2)**2) +                  &
  (K**3*NBPION*(K*SQRT(1 + MP2) + MU0 - (0,1)*P0))/                      &
   (2.*(1 + MP2)*(MP2 - MS2)*(-(K**2*(1 + MF2)) +                               &
        (K*SQRT(1 + MP2) + MU0 - (0,1)*P0)**2)**2) +                     &
  K**2/(2.*SQRT(1 + MP2)*(MP2 - MS2)**2*                                        &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0C)**2))    &
+ K**2/(4.*(1 + MP2)**1.5*(MP2 - MS2)*                                          &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0C)**2))    &
- K**2/(2.*SQRT(1 + MS2)*(-MP2 + MS2)**2*                                       &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0C)**2))    &
- K**6/(2.*SQRT(1 + MF2)*(-(K**2*(1 + MP2)) +                                   &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)**2*                   &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2))    &
- (K**3*(-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0C))/                         &
   (2.*(1 + MP2)*(MP2 - MS2)*(-(K**2*(1 + MF2)) +                               &
        (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0C)**2)**2)                   

  B2B1F1APS=REAL(B2B1F1APSC)


  B2B1F1ASPC=-(K**2*NBPION)/(2.*SQRT(1 + MP2)*(MP2 - MS2)**2*                   &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0)**2))     &
- (K**2*NBPION)/                                                                &
   (2.*SQRT(1 + MP2)*(MP2 - MS2)**2*                                            &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MP2) + MU0 - (0,1)*P0)**2)) -      &
  (K**3*NBD1XSIGMA)/                                                            &
   (4.*(1 + MS2)*(-MP2 + MS2)*(-(K**2*(1 + MF2)) +                              &
       (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0)**2)) +                      &
  (K**2*NBSIGMA)/                                                               &
   (2.*SQRT(1 + MS2)*(-MP2 + MS2)**2*                                           &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0)**2))     &
+ (K**2*NBSIGMA)/                                                               &
   (4.*(1 + MS2)**1.5*(-MP2 + MS2)*                                             &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0)**2))     &
- (K**3*NBD1XSIGMA)/                                                            &
   (4.*(1 + MS2)*(-MP2 + MS2)*(-(K**2*(1 + MF2)) +                              &
       (K*SQRT(1 + MS2) + MU0 - (0,1)*P0)**2)) +                         &
  (K**2*NBSIGMA)/                                                               &
   (2.*SQRT(1 + MS2)*(-MP2 + MS2)**2*                                           &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MS2) + MU0 - (0,1)*P0)**2)) +      &
  (K**2*NBSIGMA)/                                                               &
   (4.*(1 + MS2)**1.5*(-MP2 + MS2)*                                             &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MS2) + MU0 - (0,1)*P0)**2)) +      &
  (K**6*NFA)/(2.*SQRT(1 + MF2)*                                                 &
     (-(K**2*(1 + MP2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*     &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**    &
      2) + (K**6*NFF)/                                                          &
   (2.*SQRT(1 + MF2)*(-(K**2*(1 + MP2)) +                                       &
       (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*                           &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2)     &
- (K**3*NBSIGMA*(-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0))/                  &
   (2.*(1 + MS2)*(-MP2 + MS2)*(-(K**2*(1 + MF2)) +                              &
        (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0)**2)**2) +                  &
  (K**3*NBSIGMA*(K*SQRT(1 + MS2) + MU0 - (0,1)*P0))/                     &
   (2.*(1 + MS2)*(-MP2 + MS2)*(-(K**2*(1 + MF2)) +                              &
        (K*SQRT(1 + MS2) + MU0 - (0,1)*P0)**2)**2) -                     &
  K**2/(2.*SQRT(1 + MP2)*(MP2 - MS2)**2*                                        &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0C)**2))    &
+ K**2/(2.*SQRT(1 + MS2)*(-MP2 + MS2)**2*                                       &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0C)**2))    &
+ K**2/(4.*(1 + MS2)**1.5*(-MP2 + MS2)*                                         &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0C)**2))    &
- K**6/(2.*SQRT(1 + MF2)*(-(K**2*(1 + MP2)) +                                   &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)*                       &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**       &
         2)**2) - (K**3*(-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0C))/         &
   (2.*(1 + MS2)*(-MP2 + MS2)*(-(K**2*(1 + MF2)) +                              &
        (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0C)**2)**2)                   

  B2B1F1ASP=REAL(B2B1F1ASPC)

  B1B1F2APSC=-(K**4*NBPION)/(2.*SQRT(1 + MP2)*(MP2 - MS2)*                      &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0)**2)**    &
      2) - (K**4*NBPION)/                                                       &
   (2.*SQRT(1 + MP2)*(MP2 - MS2)*                                               &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MP2) + MU0 - (0,1)*P0)**2)**2)     &
- (K**4*NBSIGMA)/                                                               &
   (2.*SQRT(1 + MS2)*(-MP2 + MS2)*                                              &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0)**2)**    &
      2) - (K**4*NBSIGMA)/                                                      &
   (2.*SQRT(1 + MS2)*(-MP2 + MS2)*                                              &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MS2) + MU0 - (0,1)*P0)**2)**2)     &
- (K**4*NFA)/(4.*(1 + MF2)**1.5*                                                &
     (-(K**2*(1 + MP2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*     &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2))     &
+ (K**5*NFD1XA)/                                                                &
   (4.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*                        &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2))     &
+ (K**5*NFD1XF)/                                                                &
   (4.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
       (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*                           &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) -      &
  (K**4*NFF)/(4.*(1 + MF2)**1.5*                                                &
     (-(K**2*(1 + MP2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*        &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +      &
  (K**5*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                      &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*                        &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**    &
      2) + (K**5*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/             &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2*                    &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2))     &
- (K**5*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                         &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
       (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*                           &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2)     &
- (K**5*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                         &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2*                       &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) -      &
  K**4/(2.*SQRT(1 + MP2)*(MP2 - MS2)*                                           &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0C)**       &
         2)**2) - K**4/                                                         &
   (2.*SQRT(1 + MS2)*(-MP2 + MS2)*                                              &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0C)**       &
         2)**2) + K**4/                                                         &
   (4.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                      &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)*                       &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2))    &
- (K**5*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C))/                         &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)*                       &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**       &
         2)**2) - (K**5*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C))/         &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)**2*                   &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2))    

  B1B1F2APS=REAL(B1B1F2APSC)
  

  B1B1F3APSC=((K**6*NBPION)/(SQRT(1 + MP2)*(MP2 - MS2)*                         &
       (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0)**      &
           2)**3) + (K**6*NBPION)/                                              &
     (SQRT(1 + MP2)*(MP2 - MS2)*                                                &
       (-(K**2*(1 + MF2)) + (K*SQRT(1 + MP2) + MU0 - (0,1)*P0)**2)**     &
        3) + (K**6*NBSIGMA)/                                                    &
     (SQRT(1 + MS2)*(-MP2 + MS2)*                                               &
       (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0)**      &
           2)**3) + (K**6*NBSIGMA)/                                             &
     (SQRT(1 + MS2)*(-MP2 + MS2)*                                               &
       (-(K**2*(1 + MF2)) + (K*SQRT(1 + MS2) + MU0 - (0,1)*P0)**2)**     &
        3) + (K**6*NFA)/                                                        &
     (4.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*                      &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
           2)**2) + (K**6*NFA)/                                                 &
     (4.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2*                  &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
          2)) - (3*K**4*NFA)/                                                   &
     (8.*(1 + MF2)**2.5*(-(K**2*(1 + MP2)) +                                    &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*                      &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
          2)) + (3*K**5*NFD1XA)/                                                &
     (8.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*                      &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
          2)) - (K**6*NFD2XA)/                                                  &
     (8.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*                      &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
          2)) + (K**6*NFF)/                                                     &
     (4.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*                         &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**     &
        2) + (K**6*NFF)/                                                        &
     (4.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2*                     &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2))      &
+ (3*K**5*NFD1XF)/                                                              &
     (8.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*                         &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2))      &
- (K**6*NFD2XF)/                                                                &
     (8.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*                         &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2))      &
- (3*K**4*NFF)/(8.*(1 + MF2)**2.5*                                              &
       (-(K**2*(1 + MP2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*      &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2))      &
+ (3*K**5*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                    &
     (4.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*                      &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
           2)**2) - (K**6*NFD1XA*                                               &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                           &
     (2.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*                      &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
           2)**2) + (3*K**5*NFA*                                                &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                           &
     (4.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2*                  &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
          2)) - (K**6*NFD1XA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/     &
     (2.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2*                  &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
          2)) - (K**6*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)/     &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*                      &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
           2)**3) - (K**6*NFA*                                                  &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)/                        &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2*                  &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
           2)**2) - (K**6*NFA*                                                  &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)/                        &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3*                  &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
          2)) + (K**6*NFD1XF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/        &
     (2.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*                         &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**     &
        2) - (3*K**5*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/            &
     (4.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*                         &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**     &
        2) + (K**6*NFD1XF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/           &
     (2.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                    &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2*                     &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2))      &
- (3*K**5*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                       &
     (4.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2*                     &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2))      &
- (K**6*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)/                      &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*                         &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**     &
        3) - (K**6*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)/           &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2*                     &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**     &
        2) - (K**6*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)/           &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3*                     &
       (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2))      &
+ K**6/(SQRT(1 + MP2)*(MP2 - MS2)*                                              &
       (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 -                         &
             (0,1)*P0C)**2)**3) +                                        &
    K**6/(SQRT(1 + MS2)*(-MP2 + MS2)*                                           &
       (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 -                         &
             (0,1)*P0C)**2)**3) -                                        &
    K**6/(4.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)*                     &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 +                         &
             (0,1)*P0C)**2)**2) -                                        &
    K**6/(4.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)**2*                 &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**     &
          2)) + (3*K**4)/                                                       &
     (8.*(1 + MF2)**2.5*(-(K**2*(1 + MP2)) +                                    &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)*                     &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**     &
          2)) - (3*K**5*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C))/         &
     (4.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)*                     &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 +                         &
             (0,1)*P0C)**2)**2) -                                        &
    (3*K**5*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C))/                     &
     (4.*(1 + MF2)**2*(-(K**2*(1 + MP2)) +                                      &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)**2*                 &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**     &
          2)) + (K**6*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)/        &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)*                     &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 +                         &
             (0,1)*P0C)**2)**3) +                                        &
    (K**6*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)/                    &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)**2*                 &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 +                         &
             (0,1)*P0C)**2)**2) +                                        &
    (K**6*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)/                    &
     ((1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                       &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)**3*                 &
       (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**     &
          2)))/2.Q+0

  B1B1F3APS=REAL(B1B1F3APSC)

  B2B1F2APSC=(K**5*NBD1XPION)/                                                  &
   (4.*(1 + MP2)*(MP2 - MS2)*(-(K**2*(1 + MF2)) +                               &
        (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0)**2)**2) -                  &
  (K**4*NBPION)/                                                                &
   (2.*SQRT(1 + MP2)*(MP2 - MS2)**2*                                            &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0)**2)**    &
      2) - (K**4*NBPION)/                                                       &
   (4.*(1 + MP2)**1.5*(MP2 - MS2)*                                              &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0)**2)**    &
      2) + (K**5*NBD1XPION)/                                                    &
   (4.*(1 + MP2)*(MP2 - MS2)*(-(K**2*(1 + MF2)) +                               &
        (K*SQRT(1 + MP2) + MU0 - (0,1)*P0)**2)**2) -                     &
  (K**4*NBPION)/                                                                &
   (2.*SQRT(1 + MP2)*(MP2 - MS2)**2*                                            &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MP2) + MU0 - (0,1)*P0)**2)**2)     &
- (K**4*NBPION)/                                                                &
   (4.*(1 + MP2)**1.5*(MP2 - MS2)*                                              &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MP2) + MU0 - (0,1)*P0)**2)**2)     &
+ (K**4*NBSIGMA)/                                                               &
   (2.*SQRT(1 + MS2)*(-MP2 + MS2)**2*                                           &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0)**2)**    &
      2) + (K**4*NBSIGMA)/                                                      &
   (2.*SQRT(1 + MS2)*(-MP2 + MS2)**2*                                           &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MS2) + MU0 - (0,1)*P0)**2)**2)     &
+ (K**6*NFA)/(4.*(1 + MF2)**1.5*                                                &
     (-(K**2*(1 + MP2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**        &
         2)**2*(-(K**2*(1 + MS2)) +                                             &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -                      &
  (K**7*NFD1XA)/                                                                &
   (4.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2*                    &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2))     &
- (K**7*NFD1XF)/                                                                &
   (4.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2*                       &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +      &
  (K**6*NFF)/(4.*(1 + MF2)**1.5*                                                &
     (-(K**2*(1 + MP2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2*     &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +      &
  (K**5*NBPION*(-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0))/                   &
   ((1 + MP2)*(MP2 - MS2)*(-(K**2*(1 + MF2)) +                                  &
        (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0)**2)**3) -                  &
  (K**5*NBPION*(K*SQRT(1 + MP2) + MU0 - (0,1)*P0))/                      &
   ((1 + MP2)*(MP2 - MS2)*(-(K**2*(1 + MF2)) +                                  &
        (K*SQRT(1 + MP2) + MU0 - (0,1)*P0)**2)**3) -                     &
  (K**7*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                      &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2*                    &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**    &
      2) - (K**7*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/             &
   ((1 + MF2)*(-(K**2*(1 + MP2)) +                                              &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3*                    &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2))     &
+ (K**7*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                         &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2*                       &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2)     &
+ (K**7*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                         &
   ((1 + MF2)*(-(K**2*(1 + MP2)) +                                              &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3*                       &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) -      &
  K**4/(2.*SQRT(1 + MP2)*(MP2 - MS2)**2*                                        &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0C)**       &
         2)**2) - K**4/                                                         &
   (4.*(1 + MP2)**1.5*(MP2 - MS2)*                                              &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0C)**       &
         2)**2) + K**4/                                                         &
   (2.*SQRT(1 + MS2)*(-MP2 + MS2)**2*                                           &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0C)**       &
         2)**2) - K**6/                                                         &
   (4.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                      &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)**2*                   &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2))    &
+ (K**5*(-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0C))/                         &
   ((1 + MP2)*(MP2 - MS2)*(-(K**2*(1 + MF2)) +                                  &
        (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0C)**2)**3) +                 &
  (K**7*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C))/                         &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)**2*                   &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**       &
         2)**2) + (K**7*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C))/         &
   ((1 + MF2)*(-(K**2*(1 + MP2)) +                                              &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)**3*                   &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2))

  B2B1F2APS=REAL(B2B1F2APSC)

  B2B1F2ASPC=(K**4*NBPION)/(2.*SQRT(1 + MP2)*(MP2 - MS2)**2*                    &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0)**2)**    &
      2) + (K**4*NBPION)/                                                       &
   (2.*SQRT(1 + MP2)*(MP2 - MS2)**2*                                            &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MP2) + MU0 - (0,1)*P0)**2)**2)     &
+ (K**5*NBD1XSIGMA)/                                                            &
   (4.*(1 + MS2)*(-MP2 + MS2)*(-(K**2*(1 + MF2)) +                              &
        (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0)**2)**2) -                  &
  (K**4*NBSIGMA)/                                                               &
   (2.*SQRT(1 + MS2)*(-MP2 + MS2)**2*                                           &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0)**2)**    &
      2) - (K**4*NBSIGMA)/                                                      &
   (4.*(1 + MS2)**1.5*(-MP2 + MS2)*                                             &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0)**2)**    &
      2) + (K**5*NBD1XSIGMA)/                                                   &
   (4.*(1 + MS2)*(-MP2 + MS2)*(-(K**2*(1 + MF2)) +                              &
        (K*SQRT(1 + MS2) + MU0 - (0,1)*P0)**2)**2) -                     &
  (K**4*NBSIGMA)/                                                               &
   (2.*SQRT(1 + MS2)*(-MP2 + MS2)**2*                                           &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MS2) + MU0 - (0,1)*P0)**2)**2)     &
- (K**4*NBSIGMA)/                                                               &
   (4.*(1 + MS2)**1.5*(-MP2 + MS2)*                                             &
     (-(K**2*(1 + MF2)) + (K*SQRT(1 + MS2) + MU0 - (0,1)*P0)**2)**2)     &
+ (K**6*NFA)/(4.*(1 + MF2)**1.5*                                                &
     (-(K**2*(1 + MP2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*     &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**    &
      2) - (K**7*NFD1XA)/                                                       &
   (4.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*                        &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**    &
      2) - (K**7*NFD1XF)/                                                       &
   (4.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
       (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*                           &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2)     &
+ (K**6*NFF)/(4.*(1 + MF2)**1.5*                                                &
     (-(K**2*(1 + MP2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*        &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2)     &
+ (K**5*NBSIGMA*(-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0))/                  &
   ((1 + MS2)*(-MP2 + MS2)*(-(K**2*(1 + MF2)) +                                 &
        (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0)**2)**3) -                  &
  (K**5*NBSIGMA*(K*SQRT(1 + MS2) + MU0 - (0,1)*P0))/                     &
   ((1 + MS2)*(-MP2 + MS2)*(-(K**2*(1 + MF2)) +                                 &
        (K*SQRT(1 + MS2) + MU0 - (0,1)*P0)**2)**3) -                     &
  (K**7*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                      &
   ((1 + MF2)*(-(K**2*(1 + MP2)) +                                              &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)*                        &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**    &
      3) - (K**7*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/             &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2*                    &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**    &
      2) + (K**7*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                &
   ((1 + MF2)*(-(K**2*(1 + MP2)) +                                              &
       (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)*                           &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3)     &
+ (K**7*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                         &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2*                       &
     (-(K**2*(1 + MS2)) + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2)     &
+ K**4/(2.*SQRT(1 + MP2)*(MP2 - MS2)**2*                                        &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MP2)) + MU0 - (0,1)*P0C)**       &
         2)**2) - K**4/                                                         &
   (2.*SQRT(1 + MS2)*(-MP2 + MS2)**2*                                           &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0C)**       &
         2)**2) - K**4/                                                         &
   (4.*(1 + MS2)**1.5*(-MP2 + MS2)*                                             &
     (-(K**2*(1 + MF2)) + (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0C)**       &
         2)**2) - K**6/                                                         &
   (4.*(1 + MF2)**1.5*(-(K**2*(1 + MP2)) +                                      &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)*                       &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**       &
         2)**2) + (K**5*(-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0C))/         &
   ((1 + MS2)*(-MP2 + MS2)*(-(K**2*(1 + MF2)) +                                 &
        (-(K*SQRT(1 + MS2)) + MU0 - (0,1)*P0C)**2)**3) +                 &
  (K**7*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C))/                         &
   ((1 + MF2)*(-(K**2*(1 + MP2)) +                                              &
       (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)*                       &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**       &
         2)**3) + (K**7*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C))/         &
   (2.*(1 + MF2)*(-(K**2*(1 + MP2)) +                                           &
        (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)**2*                   &
     (-(K**2*(1 + MS2)) + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0C)**2)**   &
      2)

  B2B1F2ASP=REAL(B2B1F2ASPC)

  BBF_THR_COM(1)=B2B1F1APS
  BBF_THR_COM(2)=B2B1F1ASP
  BBF_THR_COM(3)=B1B1F2APS
  BBF_THR_COM(4)=B1B1F3APS
  BBF_THR_COM(5)=B2B1F2APS
  BBF_THR_COM(6)=B2B1F2ASP
END SUBROUTINE BBF_THR

SUBROUTINE BF_A_THR(MF2,T,K,NBGLDNX,NFFDNX,NFADNX,BF_A_THR_COM)
!CALCULATING THE RIGHT HAND SIDE OF DIFFERENTIAL EQUATIONS

  USE MU0P0_COM
  IMPLICIT NONE

  REAL(16) MF2,T,K
  REAL(16) NBGLDNX(6)
  REAL(16) NBGLUON,NBD1XGLUON,NBD2XGLUON,NBD3XGLUON,NBD4XGLUON,NBD5XGLUON
  REAL(16) NFFDNX(6),NFADNX(6)
  REAL(16) NFF,NFD1XF,NFD2XF,NFD3XF,NFD4XF,NFD5XF,NFA,NFD1XA,NFD2XA,NFD3XA,NFD4XA,NFD5XA

  REAL(16) B1F1A,B2F1AA,B1F2A,B2F2AA,B1F3A,B3F1AA,B3F2AA,B2F3AA
  COMPLEX(16) B1F1AC,B2F1AAC,B1F2AC,B2F2AAC,B1F3AC,B3F1AAC,B3F2AAC,B2F3AAC

  REAL(16) :: BF_A_THR_COM(8)

  NBGLUON   =NBGLDNX(1)
  NBD1XGLUON=NBGLDNX(2)
  NBD2XGLUON=NBGLDNX(3)
  NBD3XGLUON=NBGLDNX(4)
  NBD4XGLUON=NBGLDNX(5)
  NBD5XGLUON=NBGLDNX(5)
  
  NFF   =NFFDNX(1)
  NFD1XF=NFFDNX(2)
  NFD2XF=NFFDNX(3)
  NFD3XF=NFFDNX(4)
  NFD4XF=NFFDNX(5)
  NFD5XF=NFFDNX(6)

  NFA   =NFADNX(1)
  NFD1XA=NFADNX(2)
  NFD2XA=NFADNX(3)
  NFD3XA=NFADNX(4)
  NFD4XA=NFADNX(5)
  NFD5XA=NFADNX(6)

  B1F1AC=-(K**2*NBGLUON)/(2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)) - &
  (K**2*NBGLUON)/(2.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)) +    &
  (K**2*NFA)/(2.*SQRT(1 + MF2)*                                                 &
     (-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) +               &
  (K**2*NFF)/(2.*SQRT(1 + MF2)*                                                 &
     (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) -                  &
  K**2/(2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)) -            &
  K**2/(2.*SQRT(1 + MF2)*(-K**2 +                                               &
       (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2))                          

  B1F1A=REAL(B1F1AC)


  B2F1AAC=(K**3*NBD1XGLUON)/                                                    &
   (4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)) -                 &
  (K**2*NBGLUON)/                                                               &
   (4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)) +                 &
  (K**3*NBD1XGLUON)/                                                            &
   (4.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)) -                  &
  (K**2*NBGLUON)/(4.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)) -    &
  (K**4*NFA)/(2.*SQRT(1 + MF2)*                                                 &
     (-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) -            &
  (K**4*NFF)/(2.*SQRT(1 + MF2)*                                                 &
     (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) +               &
  (K**3*NBGLUON*(-K - MU0 + (0,1)*P0))/                                  &
   (2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) -              &
  (K**3*NBGLUON*(K - MU0 + (0,1)*P0))/                                   &
   (2.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) -               &
  K**2/(4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)) +            &
  K**4/(2.*SQRT(1 + MF2)*(-K**2 +                                               &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) +                    &
  (K**3*(-K - MU0 + (0,1)*P0C))/                                         &
   (2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**2)               

  B2F1AA=REAL(B2F1AAC)
  

  B1F2AC=(K**4*NBGLUON)/(2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)** &
      2) + (K**4*NBGLUON)/                                                      &
   (2.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) +               &
  (K**2*NFA)/(4.*(1 + MF2)**1.5*                                                &
     (-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -               &
  (K**3*NFD1XA)/                                                                &
   (4.*(1 + MF2)*(-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2))     &
- (K**3*NFD1XF)/                                                                &
   (4.*(1 + MF2)*(-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +      &
  (K**2*NFF)/(4.*(1 + MF2)**1.5*                                                &
     (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) -                  &
  (K**3*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                      &
   (2.*(1 + MF2)*(-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**    &
      2) + (K**3*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                &
   (2.*(1 + MF2)*(-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2)     &
+ K**4/(2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**2) -         &
  K**2/(4.*(1 + MF2)**1.5*(-K**2 +                                              &
       (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)) -                        &
  (K**3*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                            &
   (2.*(1 + MF2)*(-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2)    

  B1F2A=REAL(B1F2AC)

  B2F2AAC=-(K**5*NBD1XGLUON)/                                                   &
   (4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) +              &
  (K**4*NBGLUON)/                                                               &
   (4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) -              &
  (K**5*NBD1XGLUON)/                                                            &
   (4.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) +               &
  (K**4*NBGLUON)/                                                               &
   (4.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) -               &
  (K**4*NFA)/(4.*(1 + MF2)**1.5*                                                &
     (-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +            &
  (K**5*NFD1XA)/                                                                &
   (4.*(1 + MF2)*(-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**    &
      2) + (K**5*NFD1XF)/                                                       &
   (4.*(1 + MF2)*(-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2)     &
- (K**4*NFF)/(4.*(1 + MF2)**1.5*                                                &
     (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) -               &
  (K**5*NBGLUON*(-K - MU0 + (0,1)*P0))/                                  &
   (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**3 +                   &
  (K**5*NBGLUON*(K - MU0 + (0,1)*P0))/                                   &
   (-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**3 +                    &
  (K**5*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                      &
   ((1 + MF2)*(-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3)     &
- (K**5*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                         &
   ((1 + MF2)*(-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) +      &
  K**4/(4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**2) +         &
  K**4/(4.*(1 + MF2)**1.5*(-K**2 +                                              &
        (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) -                    &
  (K**5*(-K - MU0 + (0,1)*P0C))/                                         &
   (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**3 +                  &
  (K**5*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                            &
   ((1 + MF2)*(-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3)       

  B2F2AA=REAL(B2F2AAC)


  B1F3AC=(-((K**6*NBGLUON)/                                                     &
       (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**3) -              &
    (K**6*NBGLUON)/                                                             &
     (-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**3 -                  &
    (K**4*NFA)/(4.*(1 + MF2)**1.5*                                              &
       (-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +          &
    (3*K**2*NFA)/                                                               &
     (8.*(1 + MF2)**2.5*(-K**2 +                                                &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -                    &
    (3*K**3*NFD1XA)/                                                            &
     (8.*(1 + MF2)**2*(-K**2 +                                                  &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) +                    &
    (K**4*NFD2XA)/                                                              &
     (8.*(1 + MF2)**1.5*(-K**2 +                                                &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -                    &
    (K**4*NFF)/(4.*(1 + MF2)**1.5*                                              &
       (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) -             &
    (3*K**3*NFD1XF)/                                                            &
     (8.*(1 + MF2)**2*(-K**2 +                                                  &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +                       &
    (K**4*NFD2XF)/                                                              &
     (8.*(1 + MF2)**1.5*(-K**2 +                                                &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +                       &
    (3*K**2*NFF)/                                                               &
     (8.*(1 + MF2)**2.5*(-K**2 +                                                &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) -                       &
    (3*K**3*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                  &
     (4.*(1 + MF2)**2*(-K**2 +                                                  &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (K**4*NFD1XA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                 &
     (2.*(1 + MF2)**1.5*(-K**2 +                                                &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (K**4*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)/                 &
     ((1 + MF2)**1.5*(-K**2 +                                                   &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (K**4*NFD1XF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                    &
     (2.*(1 + MF2)**1.5*(-K**2 +                                                &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) +                   &
    (3*K**3*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                     &
     (4.*(1 + MF2)**2*(-K**2 +                                                  &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) +                   &
    (K**4*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)/                    &
     ((1 + MF2)**1.5*(-K**2 +                                                   &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) -                   &
    K**6/(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**3 +            &
    K**4/(4.*(1 + MF2)**1.5*(-K**2 +                                            &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) -                  &
    (3*K**2)/(8.*(1 + MF2)**2.5*                                                &
       (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)) -               &
    (3*K**3*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                        &
     (4.*(1 + MF2)**2*(-K**2 +                                                  &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) -                  &
    (K**4*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)/                       &
     ((1 + MF2)**1.5*(-K**2 +                                                   &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3))/2.                

  B1F3A=REAL(B1F3AC)

  B3F1AAC=((K**4*NBGLUON)/                                                      &
     (4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) +            &
    (3*K**3*NBD1XGLUON)/                                                        &
     (8.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)) -               &
    (K**4*NBD2XGLUON)/                                                          &
     (8.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)) -               &
    (3*K**2*NBGLUON)/                                                           &
     (8.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)) +               &
    (K**4*NBGLUON)/                                                             &
     (4.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) +             &
    (3*K**3*NBD1XGLUON)/                                                        &
     (8.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)) -                &
    (K**4*NBD2XGLUON)/                                                          &
     (8.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)) -                &
    (3*K**2*NBGLUON)/                                                           &
     (8.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)) +                &
    (K**6*NFA)/(SQRT(1 + MF2)*(-K**2 +                                          &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) +                &
    (K**6*NFF)/(SQRT(1 + MF2)*(-K**2 +                                          &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) -                   &
    (K**4*NBD1XGLUON*(-K - MU0 + (0,1)*P0))/                             &
     (2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) +            &
    (3*K**3*NBGLUON*(-K - MU0 + (0,1)*P0))/                              &
     (4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) -            &
    (K**4*NBGLUON*(-K - MU0 + (0,1)*P0)**2)/                             &
     (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**3 +                 &
    (K**4*NBD1XGLUON*(K - MU0 + (0,1)*P0))/                              &
     (2.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) -             &
    (3*K**3*NBGLUON*(K - MU0 + (0,1)*P0))/                               &
     (4.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) -             &
    (K**4*NBGLUON*(K - MU0 + (0,1)*P0)**2)/                              &
     (-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**3 +                  &
    K**4/(4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**2) -       &
    (3*K**2)/(8.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)) -      &
    K**6/(SQRT(1 + MF2)*(-K**2 +                                                &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3) +                  &
    (3*K**3*(-K - MU0 + (0,1)*P0C))/                                     &
     (4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**2) -           &
    (K**4*(-K - MU0 + (0,1)*P0C)**2)/                                    &
     (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**3)/2.              

  B3F1AA=REAL(B3F1AAC)

  B3F2AAC=(-(K**6*NBGLUON)/                                                     &
     (2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**3) -            &
    (3*K**5*NBD1XGLUON)/                                                        &
     (8.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) +            &
    (K**6*NBD2XGLUON)/                                                          &
     (8.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) +            &
    (3*K**4*NBGLUON)/                                                           &
     (8.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) -            &
    (K**6*NBGLUON)/                                                             &
     (2.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**3) -             &
    (3*K**5*NBD1XGLUON)/                                                        &
     (8.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) +             &
    (K**6*NBD2XGLUON)/                                                          &
     (8.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) +             &
    (3*K**4*NBGLUON)/                                                           &
     (8.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) +             &
    (K**6*NFA)/(2.*(1 + MF2)**1.5*                                              &
       (-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) -          &
    (K**7*NFD1XA)/                                                              &
     (2.*(1 + MF2)*(-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**      &
           2)**3) - (K**7*NFD1XF)/                                              &
     (2.*(1 + MF2)*(-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**     &
        3) + (K**6*NFF)/                                                        &
     (2.*(1 + MF2)**1.5*(-K**2 +                                                &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) +                   &
    (K**6*NBD1XGLUON*(-K - MU0 + (0,1)*P0))/                             &
     (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**3 -                 &
    (3*K**5*NBGLUON*(-K - MU0 + (0,1)*P0))/                              &
     (2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**3) +            &
    (3*K**6*NBGLUON*(-K - MU0 + (0,1)*P0)**2)/                           &
     (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**4 -                 &
    (K**6*NBD1XGLUON*(K - MU0 + (0,1)*P0))/                              &
     (-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**3 +                  &
    (3*K**5*NBGLUON*(K - MU0 + (0,1)*P0))/                               &
     (2.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**3) +             &
    (3*K**6*NBGLUON*(K - MU0 + (0,1)*P0)**2)/                            &
     (-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**4 -                  &
    (3*K**7*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                  &
     ((1 + MF2)*(-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**     &
        4) + (3*K**7*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/            &
     ((1 + MF2)*(-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**4)      &
- K**6/(2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**3) +         &
    (3*K**4)/(8.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**       &
        2) - K**6/                                                              &
     (2.*(1 + MF2)**1.5*(-K**2 +                                                &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3) -                  &
    (3*K**5*(-K - MU0 + (0,1)*P0C))/                                     &
     (2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**3) +           &
    (3*K**6*(-K - MU0 + (0,1)*P0C)**2)/                                  &
     (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**4 -                &
    (3*K**7*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                        &
     ((1 + MF2)*(-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**4))/   &
  2.                                                                            

  B3F2AA=REAL(B3F2AAC)

  B2F3AAC=((K**7*NBD1XGLUON)/                                                   &
     (2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**3) -            &
    (K**6*NBGLUON)/                                                             &
     (2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**3) +            &
    (K**7*NBD1XGLUON)/                                                          &
     (2.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**3) -             &
    (K**6*NBGLUON)/                                                             &
     (2.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**3) +             &
    (K**6*NFA)/(2.*(1 + MF2)**1.5*                                              &
       (-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) -          &
    (3*K**4*NFA)/                                                               &
     (8.*(1 + MF2)**2.5*(-K**2 +                                                &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (3*K**5*NFD1XA)/                                                            &
     (8.*(1 + MF2)**2*(-K**2 +                                                  &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) -                &
    (K**6*NFD2XA)/                                                              &
     (8.*(1 + MF2)**1.5*(-K**2 +                                                &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (K**6*NFF)/(2.*(1 + MF2)**1.5*                                              &
       (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) +             &
    (3*K**5*NFD1XF)/                                                            &
     (8.*(1 + MF2)**2*(-K**2 +                                                  &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) -                   &
    (K**6*NFD2XF)/                                                              &
     (8.*(1 + MF2)**1.5*(-K**2 +                                                &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) -                   &
    (3*K**4*NFF)/                                                               &
     (8.*(1 + MF2)**2.5*(-K**2 +                                                &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) +                   &
    (3*K**7*NBGLUON*(-K - MU0 + (0,1)*P0))/                              &
     (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**4 -                 &
    (3*K**7*NBGLUON*(K - MU0 + (0,1)*P0))/                               &
     (-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**4 +                  &
    (3*K**5*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                  &
     (2.*(1 + MF2)**2*(-K**2 +                                                  &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (K**6*NFD1XA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                 &
     ((1 + MF2)**1.5*(-K**2 +                                                   &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (3*K**6*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)/               &
     ((1 + MF2)**1.5*(-K**2 +                                                   &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**4) +                &
    (K**6*NFD1XF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                    &
     ((1 + MF2)**1.5*(-K**2 +                                                   &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) -                   &
    (3*K**5*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                     &
     (2.*(1 + MF2)**2*(-K**2 +                                                  &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) -                   &
    (3*K**6*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)/                  &
     ((1 + MF2)**1.5*(-K**2 +                                                   &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**4) -                   &
    K**6/(2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**3) -       &
    K**6/(2.*(1 + MF2)**1.5*(-K**2 +                                            &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3) +                  &
    (3*K**4)/(8.*(1 + MF2)**2.5*                                                &
       (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) +            &
    (3*K**7*(-K - MU0 + (0,1)*P0C))/                                     &
     (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**4 +                &
    (3*K**5*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                        &
     (2.*(1 + MF2)**2*(-K**2 +                                                  &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3) +                  &
    (3*K**6*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)/                     &
     ((1 + MF2)**1.5*(-K**2 +                                                   &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**4))/2.                

  B2F3AA=REAL(B2F3AAC)

  BF_A_THR_COM(1)= B1F1A
  BF_A_THR_COM(2)= B2F1AA
  BF_A_THR_COM(3)= B1F2A
  BF_A_THR_COM(4)= B2F2AA
  BF_A_THR_COM(5)= B1F3A
  BF_A_THR_COM(6)= B3F1AA
  BF_A_THR_COM(7)= B3F2AA
  BF_A_THR_COM(8)= B2F3AA
END SUBROUTINE BF_A_THR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BF_AT_THR(MF2,T,K,NBGLDNX,NFFDNX,NFADNX,BF_AT_THR_COM)
!CALCULATING THE RIGHT HAND SIDE OF DIFFERENTIAL EQUATIONS

  USE MU0P0_COM
  IMPLICIT NONE

  REAL(16) MF2,T,K
  REAL(16) NBGLDNX(6)
  REAL(16) NBGLUON,NBD1XGLUON,NBD2XGLUON,NBD3XGLUON,NBD4XGLUON,NBD5XGLUON
  REAL(16) NFFDNX(6),NFADNX(6)
  REAL(16) NFF,NFD1XF,NFD2XF,NFD3XF,NFD4XF,NFD5XF,NFA,NFD1XA,NFD2XA,NFD3XA,NFD4XA,NFD5XA

  REAL(16) B1F1AT,B2F1AAT,B1F2AT,B2F2AAT,B1F3AT,B3F1AAT
  COMPLEX(16) B1F1AC,B2F1AAC,B1F2AC,B2F2AAC,B1F3AC,B3F1AAC

  REAL(16) :: BF_AT_THR_COM(6)

  NBGLUON   =NBGLDNX(1)
  NBD1XGLUON=NBGLDNX(2)
  NBD2XGLUON=NBGLDNX(3)
  NBD3XGLUON=NBGLDNX(4)
  NBD4XGLUON=NBGLDNX(5)
  NBD5XGLUON=NBGLDNX(5)

  NFF   =NFFDNX(1)
  NFD1XF=NFFDNX(2)
  NFD2XF=NFFDNX(3)
  NFD3XF=NFFDNX(4)
  NFD4XF=NFFDNX(5)
  NFD5XF=NFFDNX(6)

  NFA   =NFADNX(1)
  NFD1XA=NFADNX(2)
  NFD2XA=NFADNX(3)
  NFD3XA=NFADNX(4)
  NFD4XA=NFADNX(5)
  NFD5XA=NFADNX(6)

  B1F1AC=-(K**2*NBGLUON)/(2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)) - &
  (K**2*NBGLUON)/(2.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)) +    &
  (K**2*NFA)/(2.*SQRT(1 + MF2)*                                                 &
     (-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) +               &
  (K**2*NFF)/(2.*SQRT(1 + MF2)*                                                 &
     (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) 

  B1F1AT=REAL(B1F1AC)

  B2F1AAC=(K**3*NBD1XGLUON)/                                                    &
   (4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)) -                 &
  (K**2*NBGLUON)/                                                               &
   (4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)) +                 &
  (K**3*NBD1XGLUON)/                                                            &
   (4.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)) -                  &
  (K**2*NBGLUON)/(4.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)) -    &
  (K**4*NFA)/(2.*SQRT(1 + MF2)*                                                 &
     (-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) -            &
  (K**4*NFF)/(2.*SQRT(1 + MF2)*                                                 &
     (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) +               &
  (K**3*NBGLUON*(-K - MU0 + (0,1)*P0))/                                  &
   (2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) -              &
  (K**3*NBGLUON*(K - MU0 + (0,1)*P0))/                                   &
   (2.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) 

  B2F1AAT=REAL(B2F1AAC)

  B1F2AC=(K**4*NBGLUON)/(2.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)** &
      2) + (K**4*NBGLUON)/                                                      &
   (2.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) +               &
  (K**2*NFA)/(4.*(1 + MF2)**1.5*                                                &
     (-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -               &
  (K**3*NFD1XA)/                                                                &
   (4.*(1 + MF2)*(-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2))     &
- (K**3*NFD1XF)/                                                                &
   (4.*(1 + MF2)*(-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +      &
  (K**2*NFF)/(4.*(1 + MF2)**1.5*                                                &
     (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) -                  &
  (K**3*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                      &
   (2.*(1 + MF2)*(-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**    &
      2) + (K**3*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                &
   (2.*(1 + MF2)*(-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) 

  B1F2AT=REAL(B1F2AC)

  B2F2AAC=-(K**5*NBD1XGLUON)/                                                   &
   (4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) +              &
  (K**4*NBGLUON)/                                                               &
   (4.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) -              &
  (K**5*NBD1XGLUON)/                                                            &
   (4.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) +               &
  (K**4*NBGLUON)/                                                               &
   (4.*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) -               &
  (K**4*NFA)/(4.*(1 + MF2)**1.5*                                                &
     (-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +            &
  (K**5*NFD1XA)/                                                                &
   (4.*(1 + MF2)*(-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**    &
      2) + (K**5*NFD1XF)/                                                       &
   (4.*(1 + MF2)*(-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2)     &
- (K**4*NFF)/(4.*(1 + MF2)**1.5*                                                &
     (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) -               &
  (K**5*NBGLUON*(-K - MU0 + (0,1)*P0))/                                  &
   (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**3 +                   &
  (K**5*NBGLUON*(K - MU0 + (0,1)*P0))/                                   &
   (-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**3 +                    &
  (K**5*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                      &
   ((1 + MF2)*(-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3)     &
- (K**5*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                         &
   ((1 + MF2)*(-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) 

  B2F2AAT=REAL(B2F2AAC)

  B1F3AC=(-((K**6*NBGLUON)/                                                     &
       (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**3) -              &
    (K**6*NBGLUON)/                                                             &
     (-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**3 -                  &
    (K**4*NFA)/(4.*(1 + MF2)**1.5*                                              &
       (-K**2 + (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +          &
    (3*K**2*NFA)/                                                               &
     (8.*(1 + MF2)**2.5*(-K**2 +                                                &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -                    &
    (3*K**3*NFD1XA)/                                                            &
     (8.*(1 + MF2)**2*(-K**2 +                                                  &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) +                    &
    (K**4*NFD2XA)/                                                              &
     (8.*(1 + MF2)**1.5*(-K**2 +                                                &
         (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)) -                    &
    (K**4*NFF)/(4.*(1 + MF2)**1.5*                                              &
       (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) -             &
    (3*K**3*NFD1XF)/                                                            &
     (8.*(1 + MF2)**2*(-K**2 +                                                  &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +                       &
    (K**4*NFD2XF)/                                                              &
     (8.*(1 + MF2)**1.5*(-K**2 +                                                &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) +                       &
    (3*K**2*NFF)/                                                               &
     (8.*(1 + MF2)**2.5*(-K**2 +                                                &
         (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)) -                       &
    (3*K**3*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                  &
     (4.*(1 + MF2)**2*(-K**2 +                                                  &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (K**4*NFD1XA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0))/                 &
     (2.*(1 + MF2)**1.5*(-K**2 +                                                &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**2) +                &
    (K**4*NFA*(-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)/                 &
     ((1 + MF2)**1.5*(-K**2 +                                                   &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) -                &
    (K**4*NFD1XF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                    &
     (2.*(1 + MF2)**1.5*(-K**2 +                                                &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) +                   &
    (3*K**3*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0))/                     &
     (4.*(1 + MF2)**2*(-K**2 +                                                  &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**2) +                   &
    (K**4*NFF*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)/                    &
     ((1 + MF2)**1.5*(-K**2 +                                                   &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) -                   &
    (0.)*K**6/(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**3 +            &
    (0.)*K**4/(4.*(1 + MF2)**1.5*(-K**2 +                                            &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) -                  &
    (0.)*(3*K**2)/(8.*(1 + MF2)**2.5*                                                &
       (-K**2 + (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)) -               &
    (0.)*(3*K**3*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C))/                        &
     (4.*(1 + MF2)**2*(-K**2 +                                                  &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**2) -                  &
    (0.)*(K**4*(K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)/                       &
     ((1 + MF2)**1.5*(-K**2 +                                                   &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3))/2.Q+00                

  B1F3AT=REAL(B1F3AC)

  B3F1AAC=((K**4*NBGLUON)/                                                      &
     (4.Q+0*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) +            &
    (3*K**3*NBD1XGLUON)/                                                        &
     (8.Q+0*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)) -               &
    (K**4*NBD2XGLUON)/                                                          &
     (8.*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)) -               &
    (3*K**2*NBGLUON)/                                                           &
     (8.Q+0*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)) +               &
    (K**4*NBGLUON)/                                                             &
     (4.Q+0*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) +             &
    (3*K**3*NBD1XGLUON)/                                                        &
     (8.Q+0*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)) -                &
    (K**4*NBD2XGLUON)/                                                          &
     (8.Q+0*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)) -                &
    (3*K**2*NBGLUON)/                                                           &
     (8.Q+0*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)) +                &
    (K**6*NFA)/(SQRT(1 + MF2)*(-K**2 +                                          &
          (-(K*SQRT(1 + MF2)) - MU0 + (0,1)*P0)**2)**3) +                &
    (K**6*NFF)/(SQRT(1 + MF2)*(-K**2 +                                          &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0)**2)**3) -                   &
    (K**4*NBD1XGLUON*(-K - MU0 + (0,1)*P0))/                             &
     (2.Q+0*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) +            &
    (3*K**3*NBGLUON*(-K - MU0 + (0,1)*P0))/                              &
     (4.Q+0*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**2) -            &
    (K**4*NBGLUON*(-K - MU0 + (0,1)*P0)**2)/                             &
     (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0)**2)**3 +                 &
    (K**4*NBD1XGLUON*(K - MU0 + (0,1)*P0))/                              &
     (2.Q+0*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) -             &
    (3*K**3*NBGLUON*(K - MU0 + (0,1)*P0))/                               &
     (4.Q+0*(-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**2) -             &
    (K**4*NBGLUON*(K - MU0 + (0,1)*P0)**2)/                              &
     (-(K**2*(1 + MF2)) + (K - MU0 + (0,1)*P0)**2)**3 +                  &
    (0.)*K**4/(4.Q+0*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**2) -       &
    (0.)*(3*K**2)/(8.Q+0*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)) -      &
    (0.)*K**6/(SQRT(1 + MF2)*(-K**2 +                                                &
          (K*SQRT(1 + MF2) - MU0 + (0,1)*P0C)**2)**3) +                  &
    (0.)*(3*K**3*(-K - MU0 + (0,1)*P0C))/                                     &
     (4.Q+0*(-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**2) -           &
    (0.Q+0)*(K**4*(-K - MU0 + (0,1)*P0C)**2)/                                    &
     (-(K**2*(1 + MF2)) + (-K - MU0 + (0,1)*P0C)**2)**3)/2.Q+0              

  B3F1AAT=REAL(B3F1AAC)

  BF_AT_THR_COM(1)=B1F1AT
  BF_AT_THR_COM(2)=B2F1AAT
  BF_AT_THR_COM(3)=B1F2AT
  BF_AT_THR_COM(4)=B2F2AAT
  BF_AT_THR_COM(5)=B1F3AT
  BF_AT_THR_COM(6)=B3F1AAT

END SUBROUTINE BF_AT_THR
END MODULE BF_THR_MOD
