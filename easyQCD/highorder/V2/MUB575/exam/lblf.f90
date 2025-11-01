MODULE LTM_MOD
CONTAINS 
SUBROUTINE LBTEM(K,DMBO2DRHO,NBDNX,ETAPHI,LBT)
!INPUT DIMENSIONLESS MASSES


  IMPLICIT NONE
  INTEGER(2),PARAMETER :: N=2
  !THE NUMBER OF MESONS
  REAL(16),INTENT(IN) :: K
  REAL(16) DMBO2DRHO(N,6),NBDNX(N,6)
  REAL(16) LBT(N,6)
  REAL(16) MB2D0RHO(N),MB2D1RHO(N),MB2D2RHO(N),MB2D3RHO(N),MB2D4RHO(N),MB2D5RHO(N)
  REAL(16) NBD0X(N),NBD1X(N),NBD2X(N),NBD3X(N),NBD4X(N),NBD5X(N)
  REAL(16) ETAPHI
  REAL(16) LBT0(N),LBT1(N),LBT2(N),LBT3(N),LBT4(N),LBT5(N)

  MB2D0RHO=DMBO2DRHO(:,1)
  MB2D1RHO=DMBO2DRHO(:,2)
  MB2D2RHO=DMBO2DRHO(:,3)
  MB2D3RHO=DMBO2DRHO(:,4)
  MB2D4RHO=DMBO2DRHO(:,5)
  MB2D5RHO=DMBO2DRHO(:,6)

  NBD0X=NBDNX(:,1)
  NBD1X=NBDNX(:,2)
  NBD2X=NBDNX(:,3)
  NBD3X=NBDNX(:,4)
  NBD4X=NBDNX(:,5)
  NBD5X=NBDNX(:,6)

  LBT0=(2*(1 - ETAPHI/5.Q+0)*(0.5Q+0 + NBD0X))/(3.Q+0*SQRT(1 + MB2D0RHO))
  LBT1=-((1 - ETAPHI/5.Q+0)*MB2D1RHO*(0.5Q+0 + NBD0X))/                              &
     (3.Q+0*(1 + MB2D0RHO)**1.5) +                                                      &
    ((1 - ETAPHI/5.Q+0)*K*MB2D1RHO*NBD1X)/(3.Q+0*(1 + MB2D0RHO))
  LBT2=((1 - ETAPHI/5.Q+0)*MB2D1RHO**2*(0.5Q+0 + NBD0X))/                            &
     (2.Q+0*(1 + MB2D0RHO)**2.5) -                                                      &
    ((1 - ETAPHI/5.Q+0)*MB2D2RHO*(0.5Q+0 + NBD0X))/(3.Q+0*(1 + MB2D0RHO)**1.5) -             &
    ((1 - ETAPHI/5.Q+0)*K*MB2D1RHO**2*NBD1X)/(2.Q+0*(1 + MB2D0RHO)**2) +                   &
    ((1 - ETAPHI/5.Q+0)*K*MB2D2RHO*NBD1X)/(3.Q+0*(1 + MB2D0RHO)) +                         &
    ((1 - ETAPHI/5.Q+0)*K**2*MB2D1RHO**2*NBD2X)/(6.Q+0*(1 + MB2D0RHO)**1.5)
  LBT3=(-5*(1 - ETAPHI/5.Q+0)*MB2D1RHO**3*(0.5Q+0 + NBD0X))/                         &
     (4.Q+0*(1 + MB2D0RHO)**3.5) +                                                      &
    (3*(1 - ETAPHI/5.Q+0)*MB2D1RHO*MB2D2RHO*(0.5Q+0 + NBD0X))/                       &
     (2.Q+0*(1 + MB2D0RHO)**2.5) -                                                      &
    ((1 - ETAPHI/5.Q+0)*MB2D3RHO*(0.5Q+0 + NBD0X))/(3.Q+0*(1 + MB2D0RHO)**1.5) +             &
    (5*(1 - ETAPHI/5.Q+0)*K*MB2D1RHO**3*NBD1X)/(4.Q+0*(1 + MB2D0RHO)**3) -                 &
    (3*(1 - ETAPHI/5.Q+0)*K*MB2D1RHO*MB2D2RHO*NBD1X)/(2.Q+0*(1 + MB2D0RHO)**2) +           &
    ((1 - ETAPHI/5.Q+0)*K*MB2D3RHO*NBD1X)/(3.Q+0*(1 + MB2D0RHO)) -                         &
    ((1 - ETAPHI/5.Q+0)*K**2*MB2D1RHO**3*NBD2X)/(2.Q+0*(1 + MB2D0RHO)**2.5) +              &
    ((1 - ETAPHI/5.Q+0)*K**2*MB2D1RHO*MB2D2RHO*NBD2X)/                             &
     (2.Q+0*(1 + MB2D0RHO)**1.5) +                                                      &
    ((1 - ETAPHI/5.Q+0)*K**3*MB2D1RHO**3*NBD3X)/(12.Q+0*(1 + MB2D0RHO)**2)
  LBT4=(35*(1 - ETAPHI/5.Q+0)*MB2D1RHO**4*(0.5Q+0 + NBD0X))/                         &
     (8.Q+0*(1 + MB2D0RHO)**4.5) -                                                      &
    (15*(1 - ETAPHI/5.Q+0)*MB2D1RHO**2*MB2D2RHO*(0.5Q+0 + NBD0X))/                   &
     (2.Q+0*(1 + MB2D0RHO)**3.5) +                                                      &
    (3*(1 - ETAPHI/5.Q+0)*MB2D2RHO**2*(0.5Q+0 + NBD0X))/                             &
     (2.Q+0*(1 + MB2D0RHO)**2.5) +                                                      &
    (2*(1 - ETAPHI/5.Q+0)*MB2D1RHO*MB2D3RHO*(0.5Q+0 + NBD0X))/                       &
     (1 + MB2D0RHO)**2.5 - ((1 - ETAPHI/5.Q+0)*MB2D4RHO*(0.5Q+0 + NBD0X))/                &
     (3.Q+0*(1 + MB2D0RHO)**1.5) -                                                      &
    (35*(1 - ETAPHI/5.Q+0)*K*MB2D1RHO**4*NBD1X)/(8.Q+0*(1 + MB2D0RHO)**4) +                &
    (15*(1 - ETAPHI/5.Q+0)*K*MB2D1RHO**2*MB2D2RHO*NBD1X)/                          &
     (2.Q+0*(1 + MB2D0RHO)**3) - (3*(1 - ETAPHI/5.Q+0)*K*MB2D2RHO**2*NBD1X)/               &
     (2.Q+0*(1 + MB2D0RHO)**2) - (2*(1 - ETAPHI/5.Q+0)*K*MB2D1RHO*MB2D3RHO*                &
       NBD1X)/(1 + MB2D0RHO)**2 +                                                    &
    ((1 - ETAPHI/5.Q+0)*K*MB2D4RHO*NBD1X)/(3.Q+0*(1 + MB2D0RHO)) +                         &
    (15*(1 - ETAPHI/5.Q+0)*K**2*MB2D1RHO**4*NBD2X)/(8.Q+0*(1 + MB2D0RHO)**3.5) -           &
    (3*(1 - ETAPHI/5.Q+0)*K**2*MB2D1RHO**2*MB2D2RHO*NBD2X)/                        &
     (1 + MB2D0RHO)**2.5 + ((1 - ETAPHI/5.Q+0)*K**2*MB2D2RHO**2*NBD2X)/                 &
     (2.Q+0*(1 + MB2D0RHO)**1.5) +                                                      &
    (2*(1 - ETAPHI/5.Q+0)*K**2*MB2D1RHO*MB2D3RHO*NBD2X)/                           &
     (3.Q+0*(1 + MB2D0RHO)**1.5) -                                                      &
    (5*(1 - ETAPHI/5.Q+0)*K**3*MB2D1RHO**4*NBD3X)/(12.Q+0*(1 + MB2D0RHO)**3) +             &
    ((1 - ETAPHI/5.Q+0)*K**3*MB2D1RHO**2*MB2D2RHO*NBD3X)/                          &
     (2.Q+0*(1 + MB2D0RHO)**2) + ((1 - ETAPHI/5.Q+0)*K**4*MB2D1RHO**4*NBD4X)/              &
     (24.Q+0*(1 + MB2D0RHO)**2.5)
  LBT5=(-315*(1 - ETAPHI/5.Q+0)*MB2D1RHO**5*(0.5Q+0 + NBD0X))/                       &
     (16.Q+0*(1 + MB2D0RHO)**5.5) +                                                     &
    (175*(1 - ETAPHI/5.Q+0)*MB2D1RHO**3*MB2D2RHO*(0.5Q+0 + NBD0X))/                  &
     (4.Q+0*(1 + MB2D0RHO)**4.5) -                                                      &
    (75*(1 - ETAPHI/5.Q+0)*MB2D1RHO*MB2D2RHO**2*(0.5Q+0 + NBD0X))/                   &
     (4.Q+0*(1 + MB2D0RHO)**3.5) -                                                      &
    (25*(1 - ETAPHI/5.Q+0)*MB2D1RHO**2*MB2D3RHO*(0.5Q+0 + NBD0X))/                   &
     (2.Q+0*(1 + MB2D0RHO)**3.5) +                                                      &
    (5*(1 - ETAPHI/5.Q+0)*MB2D2RHO*MB2D3RHO*(0.5Q+0 + NBD0X))/                       &
     (1 + MB2D0RHO)**2.5 + (5*(1 - ETAPHI/5.Q+0)*MB2D1RHO*MB2D4RHO*                     &
       (0.5Q+0 + NBD0X))/(2.Q+0*(1 + MB2D0RHO)**2.5) -                                    &
    ((1 - ETAPHI/5.Q+0)*MB2D5RHO*(0.5Q+0 + NBD0X))/(3.Q+0*(1 + MB2D0RHO)**1.5) +             &
    (315*(1 - ETAPHI/5.Q+0)*K*MB2D1RHO**5*NBD1X)/(16.Q+0*(1 + MB2D0RHO)**5) -              &
    (175*(1 - ETAPHI/5.Q+0)*K*MB2D1RHO**3*MB2D2RHO*NBD1X)/                         &
     (4.Q+0*(1 + MB2D0RHO)**4) + (75*(1 - ETAPHI/5.Q+0)*K*MB2D1RHO*MB2D2RHO**2*            &
       NBD1X)/(4.Q+0*(1 + MB2D0RHO)**3) +                                               &
    (25*(1 - ETAPHI/5.Q+0)*K*MB2D1RHO**2*MB2D3RHO*NBD1X)/                          &
     (2.Q+0*(1 + MB2D0RHO)**3) - (5*(1 - ETAPHI/5.Q+0)*K*MB2D2RHO*MB2D3RHO*NBD1X)/         &
     (1 + MB2D0RHO)**2 - (5*(1 - ETAPHI/5.Q+0)*K*MB2D1RHO*MB2D4RHO*NBD1X)/              &
     (2.Q+0*(1 + MB2D0RHO)**2) + ((1 - ETAPHI/5.Q+0)*K*MB2D5RHO*NBD1X)/                    &
     (3.Q+0*(1 + MB2D0RHO)) - (35*(1 - ETAPHI/5.Q+0)*K**2*MB2D1RHO**5*NBD2X)/              &
     (4.Q+0*(1 + MB2D0RHO)**4.5) +                                                      &
    (75*(1 - ETAPHI/5.Q+0)*K**2*MB2D1RHO**3*MB2D2RHO*NBD2X)/                       &
     (4.Q+0*(1 + MB2D0RHO)**3.5) -                                                      &
    (15*(1 - ETAPHI/5.Q+0)*K**2*MB2D1RHO*MB2D2RHO**2*NBD2X)/                       &
     (2.Q+0*(1 + MB2D0RHO)**2.5) -                                                      &
    (5*(1 - ETAPHI/5.Q+0)*K**2*MB2D1RHO**2*MB2D3RHO*NBD2X)/                        &
     (1 + MB2D0RHO)**2.5 + (5*(1 - ETAPHI/5.Q+0)*K**2*MB2D2RHO*MB2D3RHO*NBD2X)/         &
     (3.Q+0*(1 + MB2D0RHO)**1.5) +                                                      &
    (5*(1 - ETAPHI/5.Q+0)*K**2*MB2D1RHO*MB2D4RHO*NBD2X)/                           &
     (6.Q+0*(1 + MB2D0RHO)**1.5) +                                                      &
    (35*(1 - ETAPHI/5.Q+0)*K**3*MB2D1RHO**5*NBD3X)/(16.Q+0*(1 + MB2D0RHO)**4) -            &
    (25*(1 - ETAPHI/5.Q+0)*K**3*MB2D1RHO**3*MB2D2RHO*NBD3X)/                       &
     (6.Q+0*(1 + MB2D0RHO)**3) + (5*(1 - ETAPHI/5.Q+0)*K**3*MB2D1RHO*MB2D2RHO**2*          &
       NBD3X)/(4.Q+0*(1 + MB2D0RHO)**2) +                                               &
    (5*(1 - ETAPHI/5.Q+0)*K**3*MB2D1RHO**2*MB2D3RHO*NBD3X)/                        &
     (6.Q+0*(1 + MB2D0RHO)**2) - (5*(1 - ETAPHI/5.Q+0)*K**4*MB2D1RHO**5*NBD4X)/            &
     (16.Q+0*(1 + MB2D0RHO)**3.5) +                                                     &
    (5*(1 - ETAPHI/5.Q+0)*K**4*MB2D1RHO**3*MB2D2RHO*NBD4X)/                        &
     (12.Q+0*(1 + MB2D0RHO)**2.5) +                                                     &
    ((1 - ETAPHI/5.Q+0)*K**5*MB2D1RHO**5*NBD5X)/(48.Q+0*(1 + MB2D0RHO)**3)

  LBT(:,1)=LBT0
  LBT(:,2)=LBT1
  LBT(:,3)=LBT2
  LBT(:,4)=LBT3
  LBT(:,5)=LBT4
  LBT(:,6)=LBT5
END SUBROUTINE LBTEM

SUBROUTINE LFTEM_STRANGE(K,MF2,NFF,NFA,ETAPSI,LFS_0)
  IMPLICIT NONE
  REAL(16),INTENT(IN) :: K
  REAL(16),INTENT(IN) :: MF2,NFF,NFA
  REAL(16),INTENT(IN) :: ETAPSI
  REAL(16),INTENT(OUT) :: LFS_0
  LFS_0=((1 - ETAPSI/4.Q+0)*(1 - NFA - NFF))/(3.Q+0*SQRT(1 + MF2))
END SUBROUTINE LFTEM_STRANGE

SUBROUTINE LFTEM(K,DMFE2DRHO,NFFDNX,NFADNX,ETAPSI,LFT)

  IMPLICIT NONE
  REAL(16),INTENT(IN) :: K
  REAL(16),INTENT(IN) :: DMFE2DRHO(:),NFFDNX(:),NFADNX(:)
  REAL(16),INTENT(IN) :: ETAPSI
  REAL(16),INTENT(OUT) :: LFT(:)
  REAL(16) MF2,MF2D1RHO
  REAL(16) NFF,NFD1XF,NFD2XF,NFD3XF,NFD4XF,NFD5XF
  REAL(16) NFA,NFD1XA,NFD2XA,NFD3XA,NFD4XA,NFD5XA

  MF2=DMFE2DRHO(1)
  MF2D1RHO=DMFE2DRHO(2)

  NFF=NFFDNX(1)
  NFD1XF=NFFDNX(2)
  NFD2XF=NFFDNX(3)
  NFD3XF=NFFDNX(4)
  NFD4XF=NFFDNX(5)
  NFD5XF=NFFDNX(6)

  NFA=NFADNX(1)
  NFD1XA=NFADNX(2)
  NFD2XA=NFADNX(3)
  NFD3XA=NFADNX(4)
  NFD4XA=NFADNX(5)
  NFD5XA=NFADNX(6)

  LFT(1)=((1 - ETAPSI/4.Q+0)*(1 - NFA - NFF))/(3.Q+0*SQRT(1 + MF2))
  LFT(2)=((1 - ETAPSI/4.Q+0)*(-(K*MF2D1RHO*NFD1XA)/                                   &
          (2.Q+0*SQRT(1 + MF2)) -                                                  &
         (K*MF2D1RHO*NFD1XF)/(2.Q+0*SQRT(1 + MF2))))/                              &
     (3.Q+0*SQRT(1 + MF2)) -                                                       &
    ((1 - ETAPSI/4.Q+0)*MF2D1RHO*(1 - NFA - NFF))/(6.Q+0*(1 + MF2)**1.5)
  LFT(3)=-((1 - ETAPSI/4.Q+0)*MF2D1RHO*                                               &
        (-(K*MF2D1RHO*NFD1XA)/(2.Q+0*SQRT(1 + MF2)) -                              &
          (K*MF2D1RHO*NFD1XF)/(2.Q+0*SQRT(1 + MF2))))/                             &
     (3.Q+0*(1 + MF2)**1.5) +                                                      &
    ((1 - ETAPSI/4.Q+0)*((K*MF2D1RHO**2*NFD1XA)/(4.Q+0*(1 + MF2)**1.5) +              &
         (K*MF2D1RHO**2*NFD1XF)/(4.Q+0*(1 + MF2)**1.5) -                           &
         (K**2*MF2D1RHO**2*NFD2XA)/(4.Q+0*(1 + MF2)) -                             &
         (K**2*MF2D1RHO**2*NFD2XF)/(4.Q+0*(1 + MF2))))/                            &
     (3.Q+0*SQRT(1 + MF2)) +                                                       &
    ((1 - ETAPSI/4.Q+0)*MF2D1RHO**2*(1 - NFA - NFF))/(4.Q+0*(1 + MF2)**2.5)
  LFT(4)=(3*(1 - ETAPSI/4.Q+0)*MF2D1RHO**2*                                           &
       (-(K*MF2D1RHO*NFD1XA)/(2.Q+0*SQRT(1 + MF2)) -                               &
         (K*MF2D1RHO*NFD1XF)/(2.Q+0*SQRT(1 + MF2))))/                              &
     (4.Q+0*(1 + MF2)**2.5) -                                                      &
    ((1 - ETAPSI/4.Q+0)*MF2D1RHO*                                                  &
       ((K*MF2D1RHO**2*NFD1XA)/(4.Q+0*(1 + MF2)**1.5) +                            &
         (K*MF2D1RHO**2*NFD1XF)/(4.Q+0*(1 + MF2)**1.5) -                           &
         (K**2*MF2D1RHO**2*NFD2XA)/(4.Q+0*(1 + MF2)) -                             &
         (K**2*MF2D1RHO**2*NFD2XF)/(4.Q+0*(1 + MF2))))/                            &
     (2.Q+0*(1 + MF2)**1.5) +                                                      &
    ((1 - ETAPSI/4.Q+0)*((-3*K*MF2D1RHO**3*NFD1XA)/                                &
          (8.Q+0*(1 + MF2)**2.5) -                                                 &
         (3*K*MF2D1RHO**3*NFD1XF)/(8.Q+0*(1 + MF2)**2.5) +                         &
         (3*K**2*MF2D1RHO**3*NFD2XA)/(8.Q+0*(1 + MF2)**2) +                        &
         (3*K**2*MF2D1RHO**3*NFD2XF)/(8.Q+0*(1 + MF2)**2) -                        &
         (K**3*MF2D1RHO**3*NFD3XA)/(8.Q+0*(1 + MF2)**1.5) -                        &
         (K**3*MF2D1RHO**3*NFD3XF)/(8.Q+0*(1 + MF2)**1.5)))/                       &
     (3.Q+0*SQRT(1 + MF2)) -                                                       &
    (5*(1 - ETAPSI/4.Q+0)*MF2D1RHO**3*(1 - NFA - NFF))/(8.Q+0*(1 + MF2)**3.5)
  LFT(5)=(-5*(1 - ETAPSI/4.Q+0)*MF2D1RHO**3*                                          &
       (-(K*MF2D1RHO*NFD1XA)/(2.Q+0*SQRT(1 + MF2)) -                               &
         (K*MF2D1RHO*NFD1XF)/(2.Q+0*SQRT(1 + MF2))))/                              &
     (2.Q+0*(1 + MF2)**3.5) +                                                      &
    (3*(1 - ETAPSI/4.Q+0)*MF2D1RHO**2*                                             &
       ((K*MF2D1RHO**2*NFD1XA)/(4.Q+0*(1 + MF2)**1.5) +                            &
         (K*MF2D1RHO**2*NFD1XF)/(4.Q+0*(1 + MF2)**1.5) -                           &
         (K**2*MF2D1RHO**2*NFD2XA)/(4.Q+0*(1 + MF2)) -                             &
         (K**2*MF2D1RHO**2*NFD2XF)/(4.Q+0*(1 + MF2))))/                            &
     (2.Q+0*(1 + MF2)**2.5) -                                                      &
    (2*(1 - ETAPSI/4.Q+0)*MF2D1RHO*                                                &
       ((-3*K*MF2D1RHO**3*NFD1XA)/(8.Q+0*(1 + MF2)**2.5) -                         &
         (3*K*MF2D1RHO**3*NFD1XF)/(8.Q+0*(1 + MF2)**2.5) +                         &
         (3*K**2*MF2D1RHO**3*NFD2XA)/(8.Q+0*(1 + MF2)**2) +                        &
         (3*K**2*MF2D1RHO**3*NFD2XF)/(8.Q+0*(1 + MF2)**2) -                        &
         (K**3*MF2D1RHO**3*NFD3XA)/(8.Q+0*(1 + MF2)**1.5) -                        &
         (K**3*MF2D1RHO**3*NFD3XF)/(8.Q+0*(1 + MF2)**1.5)))/                       &
     (3.Q+0*(1 + MF2)**1.5) +                                                      &
    ((1 - ETAPSI/4.Q+0)*((15*K*MF2D1RHO**4*NFD1XA)/                                &
          (16.Q+0*(1 + MF2)**3.5) +                                                &
         (15*K*MF2D1RHO**4*NFD1XF)/(16.Q+0*(1 + MF2)**3.5) -                       &
         (15*K**2*MF2D1RHO**4*NFD2XA)/(16.Q+0*(1 + MF2)**3) -                      &
         (15*K**2*MF2D1RHO**4*NFD2XF)/(16.Q+0*(1 + MF2)**3) +                      &
         (3*K**3*MF2D1RHO**4*NFD3XA)/(8.Q+0*(1 + MF2)**2.5) +                      &
         (3*K**3*MF2D1RHO**4*NFD3XF)/(8.Q+0*(1 + MF2)**2.5) -                      &
         (K**4*MF2D1RHO**4*NFD4XA)/(16.Q+0*(1 + MF2)**2) -                         &
         (K**4*MF2D1RHO**4*NFD4XF)/(16.Q+0*(1 + MF2)**2)))/                        &
     (3.Q+0*SQRT(1 + MF2)) +                                                       &
    (35*(1 - ETAPSI/4.Q+0)*MF2D1RHO**4*(1 - NFA - NFF))/                           &
     (16.Q+0*(1 + MF2)**4.5)
  LFT(6)=(175*(1 - ETAPSI/4.Q+0)*MF2D1RHO**4*                                  &
       (-(K*MF2D1RHO*NFD1XA)/(2.Q+0*SQRT(1 + MF2)) -                               &
         (K*MF2D1RHO*NFD1XF)/(2.Q+0*SQRT(1 + MF2))))/                              &
     (16.Q+0*(1 + MF2)**4.5) -                                                     &
    (25*(1 - ETAPSI/4.Q+0)*MF2D1RHO**3*                                            &
       ((K*MF2D1RHO**2*NFD1XA)/(4.Q+0*(1 + MF2)**1.5) +                            &
         (K*MF2D1RHO**2*NFD1XF)/(4.Q+0*(1 + MF2)**1.5) -                           &
         (K**2*MF2D1RHO**2*NFD2XA)/(4.Q+0*(1 + MF2)) -                             &
         (K**2*MF2D1RHO**2*NFD2XF)/(4.Q+0*(1 + MF2))))/                            &
     (4.Q+0*(1 + MF2)**3.5) +                                                      &
    (5*(1 - ETAPSI/4.Q+0)*MF2D1RHO**2*                                             &
       ((-3*K*MF2D1RHO**3*NFD1XA)/(8.Q+0*(1 + MF2)**2.5) -                         &
         (3*K*MF2D1RHO**3*NFD1XF)/(8.Q+0*(1 + MF2)**2.5) +                         &
         (3*K**2*MF2D1RHO**3*NFD2XA)/(8.Q+0*(1 + MF2)**2) +                        &
         (3*K**2*MF2D1RHO**3*NFD2XF)/(8.Q+0*(1 + MF2)**2) -                        &
         (K**3*MF2D1RHO**3*NFD3XA)/(8.Q+0*(1 + MF2)**1.5) -                        &
         (K**3*MF2D1RHO**3*NFD3XF)/(8.Q+0*(1 + MF2)**1.5)))/                       &
     (2.Q+0*(1 + MF2)**2.5) -                                                      &
    (5*(1 - ETAPSI/4.Q+0)*MF2D1RHO*                                                &
       ((15*K*MF2D1RHO**4*NFD1XA)/(16.Q+0*(1 + MF2)**3.5) +                        &
         (15*K*MF2D1RHO**4*NFD1XF)/(16.Q+0*(1 + MF2)**3.5) -                       &
         (15*K**2*MF2D1RHO**4*NFD2XA)/(16.Q+0*(1 + MF2)**3) -                      &
         (15*K**2*MF2D1RHO**4*NFD2XF)/(16.Q+0*(1 + MF2)**3) +                      &
         (3*K**3*MF2D1RHO**4*NFD3XA)/(8.Q+0*(1 + MF2)**2.5) +                      &
         (3*K**3*MF2D1RHO**4*NFD3XF)/(8.Q+0*(1 + MF2)**2.5) -                      &
         (K**4*MF2D1RHO**4*NFD4XA)/(16.Q+0*(1 + MF2)**2) -                         &
         (K**4*MF2D1RHO**4*NFD4XF)/(16.Q+0*(1 + MF2)**2)))/                        &
     (6.Q+0*(1 + MF2)**1.5) +                                                      &
    ((1 - ETAPSI/4.Q+0)*((-105*K*MF2D1RHO**5*NFD1XA)/                              &
          (32.Q+0*(1 + MF2)**4.5) -                                                &
         (105*K*MF2D1RHO**5*NFD1XF)/(32.Q+0*(1 + MF2)**4.5) +                      &
         (105*K**2*MF2D1RHO**5*NFD2XA)/(32.Q+0*(1 + MF2)**4) +                     &
         (105*K**2*MF2D1RHO**5*NFD2XF)/(32.Q+0*(1 + MF2)**4) -                     &
         (45*K**3*MF2D1RHO**5*NFD3XA)/(32.Q+0*(1 + MF2)**3.5) -                    &
         (45*K**3*MF2D1RHO**5*NFD3XF)/(32.Q+0*(1 + MF2)**3.5) +                    &
         (5*K**4*MF2D1RHO**5*NFD4XA)/(16.Q+0*(1 + MF2)**3) +                       &
         (5*K**4*MF2D1RHO**5*NFD4XF)/(16.Q+0*(1 + MF2)**3) -                       &
         (K**5*MF2D1RHO**5*NFD5XA)/(32.Q+0*(1 + MF2)**2.5) -                       &
         (K**5*MF2D1RHO**5*NFD5XF)/(32.Q+0*(1 + MF2)**2.5)))/                      &
     (3.Q+0*SQRT(1 + MF2)) -                                                       &
    (315*(1 - ETAPSI/4.Q+0)*MF2D1RHO**5*(1 - NFA - NFF))/                          &
     (32.Q+0*(1 + MF2)**5.5)

END SUBROUTINE LFTEM
END MODULE LTM_MOD
