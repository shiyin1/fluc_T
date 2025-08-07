MODULE FNF_MOD
CONTAINS
SUBROUTINE FNF012(X,T,L,LB,NF0,NF1,NF2)
  !FERMION DISTRIBUTION FUNCTION

  IMPLICIT NONE

  REAL(16),INTENT(IN) :: X,T,L,LB
  REAL(16),INTENT(OUT) :: NF0,NF1,NF2
  REAL(16) :: OVER

  OVER=X/T
  IF(OVER > 60.Q+0)THEN
    NF0=0.Q+0
    NF1=0.Q+0
    NF2=0.Q+0
  ELSE
    NF0=1.Q+0/(1.Q+0+3.Q+0*LB*EXP(X/T)+3.Q+0*L*EXP(2.Q+0*X/T)+EXP(3.Q+0*X/T))
    NF1=(2*EXP(X/T))/(1.Q+0+3.Q+0*LB*EXP(X/T)+3.Q+0*L*EXP(2.Q+0*X/T)+EXP(3.Q+0*X/T))
    NF2=(EXP(2.Q+0*X/T))/(1.Q+0+3.Q+0*LB*EXP(X/T)+3.Q+0*L*EXP(2.Q+0*X/T)+EXP(3.Q+0*X/T))
  END IF
END SUBROUTINE FNF012
END MODULE FNF_MOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE NFDX_MOD
CONTAINS
SUBROUTINE NFDX(T,MU,L,LB,NF0,NF1,NF2,NFDNX)
!CALCULATE THE DERIVATIVES OF FERMI DISTRIBUTION FUNCTIONS

  IMPLICIT NONE

  REAL(16),INTENT(IN) :: T,MU
  REAL(16),INTENT(IN) :: L,LB,NF0,NF1,NF2
  REAL(16),INTENT(OUT) :: NFDNX(6)

  NFDNX(1)=NF0 + LB*NF1 + L*NF2
  NFDNX(2)=(-2*LB*NF1 - L*NF2 + 3*(NF0**2 + (LB*NF1 + L*NF2)**2 +               &
       NF0*(-1 + 2*LB*NF1 + 2*L*NF2)))/T
  NFDNX(3)=(18*NF0**3 + 27*NF0**2*(-1 + 2*LB*NF1 + 2*L*NF2) +                      &
    (-1 + 3*LB*NF1 + 3*L*NF2)*(6*LB**2*NF1**2 + 4*LB*NF1*(-1 + 3*L*NF2) +       &
       L*NF2*(-1 + 6*L*NF2)) +                                                  &
    9*NF0*(1 + 6*LB**2*NF1**2 + 2*L*NF2*(-2 + 3*L*NF2) +                        &
       LB*NF1*(-5 + 12*L*NF2)))/T**2
  NFDNX(4)=(162*NF0**4 + 162*LB**4*NF1**4 + 324*NF0**3*(-1 + 2*LB*NF1 + 2*L*NF2) +   &
    216*LB**3*NF1**3*(-1 + 3*L*NF2) +                                           &
    8*LB*NF1*(-1 + 3*L*NF2)*(1 + 9*L*NF2*(-1 + 3*L*NF2)) +                      &
    L*NF2*(-1 + 3*L*NF2)*(1 + 18*L*NF2*(-1 + 3*L*NF2)) +                        &
    12*LB**2*NF1**2*(7 + 9*L*NF2*(-5 + 9*L*NF2)) +                              &
    27*NF0**2*(7 + 36*LB**2*NF1**2 + 4*L*NF2*(-7 + 9*L*NF2) +                   &
       8*LB*NF1*(-4 + 9*L*NF2)) +                                               &
    3*NF0*(-9 + 216*LB**3*NF1**3 + 36*LB**2*NF1**2*(-7 + 18*L*NF2) +            &
       8*LB*NF1*(11 + 27*L*NF2*(-2 + 3*L*NF2)) +                                &
       2*L*NF2*(29 + 18*L*NF2*(-5 + 6*L*NF2))))/T**3
  NFDNX(5)=(1944*NF0**5 + 4860*NF0**4*(-1 + 2*LB*NF1 + 2*L*NF2) +                  &
    810*NF0**3*(5 + 24*LB**2*NF1**2 + 4*L*NF2*(-5 + 6*L*NF2) +                  &
       2*LB*NF1*(-11 + 24*L*NF2)) +                                             &
    45*NF0**2*(-27 + 432*LB**3*NF1**3 + 108*LB**2*NF1**2*(-5 + 12*L*NF2) +      &
       2*L*NF2*(83 + 216*L*NF2*(-1 + L*NF2)) +                                  &
       2*LB*NF1*(107 + 162*L*NF2*(-3 + 4*L*NF2))) +                             &
    (-1 + 3*LB*NF1 + 3*L*NF2)*(648*LB**4*NF1**4 +                               &
       864*LB**3*NF1**3*(-1 + 3*L*NF2) +                                        &
       L*NF2*(-1 + 6*L*NF2)*(1 + 36*L*NF2*(-1 + 3*L*NF2)) +                     &
       8*LB*NF1*(-1 + 3*L*NF2)*(2 + 27*L*NF2*(-1 + 4*L*NF2)) +                  &
       12*LB**2*NF1**2*(26 + 9*L*NF2*(-19 + 36*L*NF2))) +                       &
    3*NF0*(27 + 3240*LB**4*NF1**4 + 1620*LB**3*NF1**3*(-3 + 8*L*NF2) +          &
       60*LB**2*NF1**2*(41 + 108*L*NF2*(-2 + 3*L*NF2)) +                        &
       10*L*NF2*(-26 + 3*L*NF2*(43 + 108*L*NF2*(-1 + L*NF2))) +                 &
       5*LB*NF1*(-95 + 12*L*NF2*(61 + 27*L*NF2*(-7 + 8*L*NF2)))))/T**4
  NFDNX(6)=(29160*NF0**6 + 29160*LB**6*NF1**6 +                                    &
    87480*NF0**5*(-1 + 2*LB*NF1 + 2*L*NF2) +                                    &
    58320*LB**5*NF1**5*(-1 + 3*L*NF2) +                                         &
    9*NF0*(-27 + 2*LB*NF1*(412 +                                                &
          45*LB*NF1*(-77 + 12*LB*NF1*(22 + 3*LB*NF1*(-11 + 6*LB*NF1)))) +       &
       374*L*NF2 + 900*L*LB*NF1*(-2 + 3*LB*NF1)*                                &
        (5 + 12*LB*NF1*(-2 + 3*LB*NF1))*NF2 +                                   &
       540*L**2*(-5 + 9*LB*NF1*(9 + 4*LB*NF1*(-9 + 10*LB*NF1)))*NF2**2 +        &
       2160*L**3*(5 + 6*LB*NF1*(-8 + 15*LB*NF1))*NF2**3 +                       &
       3240*L**4*(-7 + 30*LB*NF1)*NF2**4 + 19440*L**5*NF2**5) +                 &
    L*NF2*(-1 + 3*L*NF2)*(1 + 90*L*NF2*(1 - 6*L*NF2)**2*(-1 + 3*L*NF2)) +       &
    3240*LB**4*NF1**4*(13 + 27*L*NF2*(-3 + 5*L*NF2)) +                          &
    7290*NF0**4*(13 + 60*LB**2*NF1**2 + 4*L*NF2*(-13 + 15*L*NF2) +              &
       8*LB*NF1*(-7 + 15*L*NF2)) +                                              &
    12960*LB**3*NF1**3*(-1 + 3*L*NF2)*(1 + L*NF2*(-7 + 15*L*NF2)) +             &
    4*LB*NF1*(-1 + 3*L*NF2)*(8 +                                                &
       45*L*NF2*(-1 + 3*L*NF2)*(5 + 36*L*NF2*(-1 + 3*L*NF2))) +                 &
    1620*NF0**3*(-27 + 360*LB**3*NF1**3 +                                       &
       36*LB**2*NF1**2*(-13 + 30*L*NF2) +                                       &
       LB*NF1*(197 + 216*L*NF2*(-4 + 5*L*NF2)) +                                &
       4*L*NF2*(41 + 9*L*NF2*(-11 + 10*L*NF2))) +                               &
    27*NF0**2*(279 + 16200*LB**4*NF1**4 +                                       &
       12960*LB**3*NF1**3*(-2 + 5*L*NF2) +                                      &
       90*LB**2*NF1**2*(163 + 72*L*NF2*(-11 + 15*L*NF2)) +                      &
       20*LB*NF1*(-172 + 9*L*NF2*(133 + 360*L*NF2*(-1 + L*NF2))) +              &
       20*L*NF2*(-119 + 9*L*NF2*(53 + 18*L*NF2*(-6 + 5*L*NF2)))) +              &
    6*LB**2*NF1**2*(248 + 45*L*NF2*                                             &
        (-85 + 9*L*NF2*(59 + 12*L*NF2*(-14 + 15*L*NF2)))))/T**5

END SUBROUTINE NFDX

SUBROUTINE NFDX_COM(MF2,T,MU,L,LB,K,NFFFDNX,NFAFDNX)
!CALCULATING THE RIGHT HAND SIDE OF DIFFERENTIAL EQUATIONS

  USE FNF_MOD
  IMPLICIT NONE

  REAL(16) MF2,T,MU,L,LB,K
  REAL(16) X,NF0,NF1,NF2
  REAL(16) NFFFDNX(6),NFAFDNX(6)

  X=K*SQRT(1.Q+0 + MF2) - MU
  CALL FNF012(X,T,L,LB,NF0,NF1,NF2)
  CALL NFDX(T,MU,L,LB,NF0,NF1,NF2,NFFFDNX)

  X=K*SQRT(1.Q+0 + MF2) + MU
  CALL FNF012(X,T,LB,L,NF0,NF1,NF2)
  CALL NFDX(T,MU,LB,L,NF0,NF1,NF2,NFAFDNX)

END SUBROUTINE NFDX_COM
END MODULE NFDX_MOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE FNB_MOD
CONTAINS
SUBROUTINE FNB(X,T,NB)
  !BOSON DISTRIBUTION FUNCTION
  IMPLICIT NONE

  REAL(16),INTENT(IN) :: X, T
  REAL(16),INTENT(OUT) :: NB
  REAL(16) :: OVER

  OVER=X/T
  IF(OVER>60.Q+0)THEN
    NB=0.Q+0
  ELSEIF(OVER<1.Q-9)THEN
    NB=1.Q+0/OVER
  ELSE
    NB=1.Q+0/(EXP(OVER)-1.Q+0)
  END IF
END SUBROUTINE FNB
END MODULE FNB_MOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE NBDX_MOD
CONTAINS
SUBROUTINE NBDX_COM(MB2,T,K,NBDX,NBBDN,FINVEBDN)
!CALCULATING THE RIGHT HAND SIDE OF DIFFERENTIAL EQUATIONS

  USE FNB_MOD
  IMPLICIT NONE

  REAL(16) MB2,T,K
  REAL(16) NB
  REAL(16) NBDX(6)
  REAL(16) NBBDN(2)
  REAL(16) FINVEBDN(2)

  CALL FNB(K*SQRT(1.Q+0 + MB2),T,NB)

  NBDX(1)=NB
  NBDX(2)=-((NB*(1 + NB))/T)
  NBDX(3)=(NB*(1 + NB)*(1 + 2*NB))/T**2
  NBDX(4)=-((NB*(1 + NB)*(1 + 6*NB*(1 + NB)))/T**3)
  NBDX(5)=(NB*(1 + NB)*(1 + 2*NB)*(1 + 12*NB*(1 + NB)))/T**4
  NBDX(6)=-((NB*(1 + NB)*(1 + 30*NB*(1 + NB)*(1 + 2*NB)**2))/T**5)

  NBBDN(1)=NBDX(1)
  NBBDN(2)=(K*NBDX(2))/(2.Q+0*SQRT(1.Q+0 + MB2))

  FINVEBDN(1)=1/(K*SQRT(1.Q+0 + MB2))
  FINVEBDN(2)=-1/(2.Q+0*K*(1.Q+0 + MB2)**1.5)

END SUBROUTINE NBDX_COM
END MODULE NBDX_MOD
