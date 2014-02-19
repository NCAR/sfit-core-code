!-----------------------------------------------------------------------------
!    Copyright (c) 2013-2014 NDACC/IRWG
!    This file is part of sfit.
!
!    sfit is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.
!
!    sfit is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with sfit.  If not, see <http://www.gnu.org/licenses/>
!-----------------------------------------------------------------------------

      MODULE SYNSPEC

      USE PARAMS
      USE TRANSMIS
      USE SOLAR
      USE BANDPARAM

      IMPLICIT NONE

      COMPLEX(DBLE_COMPLEX), DIMENSION(:), ALLOCATABLE :: IMGG
      INTEGER, DIMENSION(MAXBND)      :: MPT
      INTEGER, DIMENSION(MAXBND)      :: MFFT
      INTEGER, DIMENSION(MAXBND)      :: LOWFIL
      INTEGER, DIMENSION(MAXBND)      :: HIFILL
      INTEGER, DIMENSION(MAXBND)      :: NSTZ1
      INTEGER, DIMENSION(MAXBND)      :: NSTZ2
      INTEGER, DIMENSION(MAXBND)      :: NSTART

      REAL(DOUBLE), DIMENSION(500)    :: APD
      REAL(DOUBLE), DIMENSION(MAXEAP) :: EAPF
      REAL(DOUBLE), DIMENSION(MAXEAP) :: EAPF0
      REAL(DOUBLE), DIMENSION(MAXEAP) :: EAPX
      REAL(DOUBLE), DIMENSION(MAXEAP) :: EPHSF
      REAL(DOUBLE), DIMENSION(MAXEAP) :: EPHSF0
      REAL(DOUBLE), DIMENSION(MAXEAP) :: EPHSX

      LOGICAL      :: F_EPHASE=.FALSE., F_EAPOD=.FALSE.
      INTEGER      :: IEPHS, JEPHS, NEPHS, IEAP, JEAP, NEAP
      INTEGER      :: NEAPRT = 0, NEPHSRT = 0
      REAL(DOUBLE) :: EPHSPAR, SEPHSPAR, EAPPAR, SEAPPAR, AMP

      INTEGER, DIMENSION(MAXBND)      :: IENORM, IAP
      REAL(DOUBLE), DIMENSION(MMAX)   :: TOBS, TOBS_ORIG
      REAL(DOUBLE) :: RMS

      REAL(DOUBLE) :: C0, C1, C2, C4

      CONTAINS

!----------------------------------------------------------------------
      SUBROUTINE FSPEC2(IBAND, MONONE, PHI)


! --- APPLY APODIZATION AND PHASE ERROR TO INTERFEROGRAM
! --- COMPUTE  COMPLEX SPECTRUM BY TAKING INVERSE FFT
! --- REVISED  SEPT 24, 1996  CPR

      INTEGER, INTENT(IN)         :: IBAND, MONONE
      REAL(DOUBLE), INTENT(INOUT) :: PHI

      INTEGER, DIMENSION(:),ALLOCATABLE :: IW
      INTEGER      :: NZ1, NZ2, NMON, LF, I, MINUS, NAERR
      REAL(DOUBLE) :: WMEAN, DX, PHI0, X, AP, FACTOR, T2, T4, T6, T8
      !REAL(DOUBLE) , EXTERNAL :: EPHS, EAPDZ, APDZ

!      REAL*4 TOBS

!  --- RETRIEVE NORTON AND BEER APODIZING COEFFICIENTS
      IF (IAP(IBAND)>=0 .AND. IAP(IBAND)<4) THEN
         CALL BOBAPD(IBAND)
      ENDIF

      WMEAN = (WSTOP(IBAND)-WSTART(IBAND))/2.D0 + WSTART(IBAND)

! --- APODIZE AND FIELD OF VIEW CORRECTION

      NZ1 = NSTZ1(IBAND)
      NZ2 = NSTZ2(IBAND)
      NMON = NM(IBAND)
      LF = LOWFIL(IBAND)

! --- PHASE ERROR
! --- ADD PHASE FUNCTION (SAME FOR ALL BANDPASSES AND SPECTRA) TO PHASE
! --- CONSTANTS RETRIEVED FOR EACH BANDPASS AND SPECTRUM

      DX = 1.D0/(MPT(IBAND)*DN(IBAND))
      !print*, iband, dx, mpt(iband), dn(iband)
      PHI0 = PHI
      DO I = 1, NZ1
         !X = I/(MPT(IBAND)*DN(IBAND))
         X = REAL(I,KIND=8)*DX
         !print*, i, nz1, x
         IF (IEPHS > 0) PHI = PHI0 + EPHS(X,DX,PMAX(IBAND))
         IMGG(I+1) = IMGG(I+1)*EXP(CMPLX(0.D0,(-PHI),KIND = 8))
         IMGG(MPT(IBAND)+1-I) = CONJG(IMGG(I+1))
      END DO

! --- CENTER FRINGE
      X = 0.D0
      IF (IEPHS > 0) PHI = PHI0 + EPHS(X,DX,PMAX(IBAND))
      IMGG(1) = IMGG(1)*EXP(CMPLX(0.D0,(-PHI),KIND = 8))

! --- APODIZE AND FIELD OF VIEW CORRECTION

      DO I = 1, NZ1
         X = I/(MPT(IBAND)*DN(IBAND))
! --- MODEL EMPIRICAL APODIZATION
         AP = 1.D0
         IF (IEAP > 0) AP = EAPDZ(X,DX,PMAX(IBAND))
!  --- CONSTRAIN TO BE POSITIVE
         AP = DMAX1(0.D0,AP)
         FACTOR = WMEAN*X*OMEGA(IBAND)/2.D0
         IF (FACTOR <= 0.7D0) THEN
            T2 = FACTOR*FACTOR
            T4 = T2*T2
            T6 = T4*T2
            T8 = T4*T4
            FACTOR = 1.D0 - T2/6.D0 + T4/120.D0 - T6/5040.D0 + T8/362880.D0
         ELSE
            FACTOR = SIN(FACTOR)/FACTOR
         ENDIF
         IMGG(I+1)            = IMGG(I+1)*FACTOR*AP*APDZ(IBAND,X)
         IMGG(MPT(IBAND)+1-I) = IMGG(MPT(IBAND)+1-I)*FACTOR*AP*APDZ(IBAND,X)
      END DO

      IMGG(NZ1:NZ2) = 0.D0

      ALLOCATE (IW(MFFT(IBAND)+1), STAT=NAERR)

      IF (NAERR /= 0) THEN
          WRITE(16, *) 'FSPECT2: COULD NOT ALLOCATE IW ARRAY'
          WRITE(16, *) 'ERROR NUMBER = ', NAERR
          WRITE(00, *) 'FSPECT2: COULD NOT ALLOCATE IW ARRAY'
          WRITE(00, *) 'ERROR NUMBER = ', NAERR
          CALL SHUTDOWN
          STOP 4
      ENDIF

! ---COMPLEX INVERSE TRANSFORM
      MINUS = -MFFT(IBAND)

      CALL FFT (IMGG, MINUS, IW)
!
! --- ADD OFFSET AND STORE IN TCONV ARRAY
      TCONV(MONONE:NMON-1+MONONE) = IMGG(1+LF:NMON+LF) + AMP
!print *, 'iband', iband
!print *, TCONV(MONONE:NMON-1+MONONE)

      DEALLOCATE(IW)

      RETURN

      END SUBROUTINE FSPEC2


!----------------------------------------------------------------------
      SUBROUTINE FSPEC1(IBAND, MONONE, MXONE)

! --- TRUNCATE FFT OF SPECTRUM TO SIMULATE FOURIER SPECTROMETER
! --- UNAPODIZED OUTPUT

      INTEGER, INTENT(IN)               :: IBAND, MONONE, MXONE
      INTEGER                           :: HP, LF, NMON, MPTBND, MFFTBD, NAERR
      INTEGER, DIMENSION(:),ALLOCATABLE :: IW

! --- SAVE AMPLITUDE-- ASSUME SAMPLING INTO WINGS
! --- FILL COMPLEX ARRAY IMGG

      ALLOCATE (IW(MFFT(IBAND)+1), STAT=NAERR)
      IF (NAERR /= 0) THEN
          WRITE(16, *) 'FSPECT1: COULD NOT ALLOCATE IW ARRAY'
          WRITE(16, *) 'ERROR NUMBER = ', NAERR
          WRITE(00, *) 'FSPECT1: COULD NOT ALLOCATE IW ARRAY'
          WRITE(00, *) 'ERROR NUMBER = ', NAERR
          CALL SHUTDOWN
          STOP 4
      ENDIF

      LF        = LOWFIL(IBAND)
      NMON      = NM(IBAND)
      MPTBND    = MPT(IBAND)
      MFFTBD    = MFFT(IBAND)
      IMGG(:LF) = 0.D0
      HP        = HIFILL(IBAND) + 1
      IMGG(HP:MPTBND) = 0.D0

      IF( .NOT. IFCO )THEN
         AMP = TCALC(2,MONONE)
         IMGG(1+LF:NMON+LF) = TCALC(2,MONONE:NMON-1+MONONE) - AMP
      ELSE
         AMP = TCALC(2,MONONE)*TCO(MXONE)
         IMGG(1+LF:NMON+LF) = TCALC(2,MONONE:NMON-1+MONONE)*TCO(MXONE:NMON-1+MXONE) - AMP
      ENDIF

!  --- COMPUTE FORWARD TRANSFORM
      CALL FFT (IMGG, MFFTBD, IW)

      DEALLOCATE(IW)

      RETURN

      END SUBROUTINE FSPEC1


!----------------------------------------------------------------------
      REAL(KIND(0.0D0)) FUNCTION EPHS (X, DX, PMA)

      REAL(DOUBLE), INTENT(IN) :: X, DX, PMA
      INTEGER                  :: I, IMAX
      REAL(DOUBLE)             :: XP

      XP = X/PMA
      EPHS = 0.D0

      IF( F_EPHASE )THEN
         SELECT CASE (IEPHS)

         CASE (1, 4)
   ! --- INTERPOLATE THE PHASE FUNCTION TO PATH DIFFERENCE X
            IF( X .LT. (EPHSX(1) - DX) .OR. X .GT. (EPHSX(JEPHS) + DX) )THEN
               !PRINT *, X, EPHSX(1), DX, EPHSX(JEPHS)
               WRITE (16, *) 'PATH DIFFERENCE', X, ' IS OUT OF RANGE OF EPHSX'
               WRITE (00, *) 'PATH DIFFERENCE', X, ' IS OUT OF RANGE OF EPHSX'
               CALL SHUTDOWN
               STOP 4
            ENDIF
            IMAX = JEPHS
            DO I = 2, JEPHS
               IF (EPHSX(I) <= X) CYCLE
               IMAX = I
               EXIT
            END DO
            !print*, imax, x, jephs, dx, EPHSX(IMAX), EPHSX(jephs)
            EPHS = EPHSF(IMAX-1) + (X - EPHSX(IMAX-1))*(EPHSF(IMAX)-EPHSF(IMAX-1))&
               /(EPHSX(IMAX)-EPHSX(IMAX-1))

         CASE (2)
   ! --- EPHS IS A POLYNOMIAL WITH NEPHS TERMS
            EPHS = 0.D0
            DO I = 1, NEPHS
               EPHS = EPHS + (EPHSF(I)-1.D0)*XP**I
            END DO

         CASE DEFAULT
            EPHS = 0.D0

         END SELECT

      ENDIF

      RETURN

      END FUNCTION EPHS


!----------------------------------------------------------------------
      REAL(KIND(0.0D0)) FUNCTION EAPDZ (X, DX, PMA)

      REAL(DOUBLE), INTENT(IN) :: X, DX, PMA
      INTEGER                  :: I, IMAX=0
      REAL(DOUBLE)             :: XP

      EAPDZ = 0.0D0
      XP = X/PMA


      IF( F_EAPOD )THEN

         SELECT CASE (IEAP)
         CASE (1, 4)
   ! --- INTERPOLATE THE APODIZATION FUNCTION TO PATH DIFFERENCE X
            IF (X<EAPX(1) - DX .OR. X>EAPX(JEAP)+DX) THEN
               WRITE(16, *) 'PATH DIFFERENCE', X, ' IS OUT OF RANGE OF EAPX'
               WRITE(00, *) 'PATH DIFFERENCE', X, ' IS OUT OF RANGE OF EAPX'
               CALL SHUTDOWN
               STOP 4
            ENDIF
            DO I = 2, JEAP
               IF (EAPX(I) <= X) CYCLE
               IMAX = I
               EXIT
            END DO
            EAPDZ = EAPF(IMAX-1) + (X - EAPX(IMAX-1))*(EAPF(IMAX)-EAPF(IMAX-1))/(&
               EAPX(IMAX)-EAPX(IMAX-1))

         CASE (2)
   ! --- EAPDZ IS A POLYNOMIAL WITH NEAP TERMS
            EAPDZ = 1.D0
            DO I = 1, NEAP
               EAPDZ = EAPDZ + (EAPF(I)-1.D0)*XP**I
            END DO

         CASE (3)
   ! --- EAPDZ IS A FOURIER SERIES WITH NEAP FREQUENCIES
            EAPDZ = 1.D0
            DO I = 1, NEAP
               EAPDZ = EAPDZ + (1.D0 - EAPF(I+1))*SIN(2.D0*PI*I*EAPF(1)*XP) + (&
                  1.D0 - EAPF(I+2))*COS(2.D0*PI*I*EAPF(1)*XP)
            END DO

         CASE DEFAULT
            WRITE(16, *) ' EAPDZ.F : ERROR IEAP OUT OF RANGE (1-4) : ', IEAP
            WRITE(00, *) ' EAPDZ.F : ERROR IEAP OUT OF RANGE (1-4) : ', IEAP
            CALL SHUTDOWN
            STOP 4
         END SELECT

      ENDIF

      RETURN
      END FUNCTION EAPDZ


!----------------------------------------------------------------------
      SUBROUTINE BOBAPD(IBAND)

      INTEGER :: IAP1, IBAND

!  --- RETRIEVE NORTON AND BEER APODIZING FUNCTIONS - REVISED FEB 22, 1990
      IAP1 = IAP(IBAND) + 1

      SELECT CASE (IAP1)
!  IAP=0 (BOX CAR APODIZATION)
      CASE DEFAULT
         C0 = 1.D0
         C1 = 0.D0
         C2 = 0.D0
         C4 = 0.D0
         RETURN
!  IAP=1  (WEAK APODIZATION)
      CASE (2)
         C0 = 0.384093D0
         C1 = -0.087577D0
         C2 = 0.703484D0
         C4 = 0.D0
         RETURN
!  IAP=2 (MODERATE APODIZATION)
      CASE (3)
         C0 = 0.152442D0
         C1 = -0.136176D0
         C2 = 0.983734D0
         C4 = 0.D0
         RETURN
!  IAP=3 (STRONG APODIZATION)
      CASE (4)
         C0 = 0.045335D0
         C1 = 0.D0
         C2 = 0.554883D0
         C4 = 0.399782D0
         RETURN
!  IAP=4 (DENVER DATA)
      CASE (5)
         C0 = 0.D0
         C1 = 0.D0
         C2 = 1.D0
         C4 = 0.D0
         RETURN
      END SELECT

      END SUBROUTINE BOBAPD


!----------------------------------------------------------------------
      REAL(KIND(0.0D0)) FUNCTION APDZ (IBAND,X)

!      CALCULATE APPLIED APODIZING FUNCTION FOR SFIT2- REVISED:  JULY 3, 1997
!      APODIZING FUNCTIONS 0-3 FROM
!        REF: NORTON AND BEER, J. OPT. SOC. AMER. 66, 259 (1976).
!           REVISED COEFFICIENTS FROM J. OPT. SOC. AMER. 67, 419 (1977)

      INTEGER, INTENT(IN) :: IBAND
      REAL(DOUBLE)        :: X, W, PART

      APDZ = 0.0D0
      IF (IAP(IBAND) >= 0 .AND. IAP(IBAND) <= 9) THEN
         IF (IAP(IBAND) <= 4) THEN
!  --- NORTON AND BEER APODZING FUNCTIONS
            IF (X >= PMAX(IBAND)) THEN
               APDZ = 0.0D0
               RETURN
            ENDIF
            W = 1.0D0 - (X/PMAX(IBAND))**2
            APDZ = C0 + C1*W + C2*W*W + C4*W*W*W*W
            RETURN
         ENDIF
         IF (IAP(IBAND) <= 5) THEN
!  --- IAP=5 (TRIANGLE APODIZING FUNCTION)
            APDZ = 1.0D0 - X/PMAX(IBAND)
            RETURN
         ENDIF
         IF (IAP(IBAND) <= 6) THEN
!  --- IAP=6 (HAPP-GENZEL)
            APDZ = 0.54D0 + 0.46D0*COS(3.14159D0*X/PMAX(IBAND))
            RETURN
         ENDIF
!  --- IAP=7 (KPNO ATMOSPHERIC SPECTRA)
         IF (IAP(IBAND) > 7) GO TO 104
         APDZ = 1.D0
         PART = X/PMAX(IBAND)
         IF (PART <= 0.9D0) RETURN
         IF (PART <= 1.0D0) THEN
!  --- APPLY COSIGN SQUARE TO TAIL OF INTERFEROGRAM
            PART = 10.D0*(PART - 0.9D0)
            PART = 1.570795D0*PART
            APDZ = COS(PART)
            APDZ = APDZ*APDZ
            RETURN
         ENDIF
   15    CONTINUE
         APDZ = 0.D0
         RETURN
!  --- IAP=8 (KPNO ATMOSPHERIC SPECTRA)
  104    CONTINUE
         IF (IAP(IBAND) <= 8) THEN
            APDZ = 1.D0
            PART = X/PMAX(IBAND)
            IF (PART <= 0.95D0) RETURN
            IF (PART > 1.0D0) GO TO 15
!  --- APPLY COSIGN SQUARE TO TAIL OF INTERFEROGRAM
            PART = 10.D0*(PART - 0.95D0)
            PART = 1.570795D0*PART
            APDZ = COS(PART)
            APDZ = APDZ*APDZ
            RETURN
         ENDIF
!  IAP=9  HAMMING FUNCTION
         PART = X/PMAX(IBAND)
         APDZ = 0.53856D0 + 0.46144D0*COS(3.141592654D0*PART)
         IF (PART > 1.D0) APDZ = 0.D0
         RETURN
      ENDIF

!      YOU HAVE ASKED FOR NON-EXISTANT APODIZATION
      WRITE(16, 11) IAP
      WRITE(00, 11) IAP
      CALL SHUTDOWN
      STOP 4

   11 FORMAT(' NO SUCH APODIZING FUNCTION - IAP =',I4)

      END FUNCTION APDZ


!----------------------------------------------------------------------
      SUBROUTINE FFT(Z, MVAL, IWK)

      INTEGER  :: MVAL
      INTEGER  :: IWK(*)
      COMPLEX(DBLE_COMPLEX)  :: Z(*)

      INTEGER :: K01, J2, M, MP, N, I, MM, KN, MK, KB, K0, K2, JJ, K, J, ISP, &
         JSP, K1, K3=0
      REAL(DOUBLE) :: A0, A1, A3, B1, B2, B3, CK, SK, SQ, ONE, ZERO, TWOPI, &
         SYGN, RAD, C1, C2=0.0, C3=0.0, S1, S2=0.0, S3=0.0, TEMP, XN, A2, B0
      REAL(DOUBLE), DIMENSION(2) :: Z0, Z1, Z2, Z3
      COMPLEX(DBLE_COMPLEX) :: ZA0, ZA1, ZA2, ZA3, AK2
!-----------------------------------------------
!     E2.4
!***********************************************************************FFT    4
! VERSION TO RUN ON IBM PC-AT AND CLONES
! USES DOUBLE PREISION  REV DATE MAY 9,1987
!
!   FUNCTION            - COMPUTE THE FAST FOURIER TRANSFORM, GIVEN A
!                           COMPLEX VECTOR OF LENGTH EQUAL TO A POWER
!                           OF TWO
!   USAGE               - CALL FFT (Z,M,IWK)
!   PARAMETERS  Z       - COMPLEX VECTOR OF LENGTH N=2**M
!                           WHICH CONTAINS ON INPUT THE
!                           DATA TO BE TRANSFORMED. ON
!                           OUTPUT,A CONTAINS THE FOURIER
!                           COEFFICIENTS.
!                M      - N = 2**M IS THE NUMBER OF DATA POINTS.
!                         M= +N FFT WILL PERFORM FOURIER
!                             TRANSFORM.
!                         M= -N FFT WILL PERFORM INVERSE
!                             TRANSFORM.
!                IWK    - WORK AREA VECTOR OF LENGTH M+1.
!   PRECISION           - SINGLE
!   LANGUAGE            - FORTRAN
!   LATEST REVISION     - APRIL 16, 1980
!-----------------------------------------------------------------------FFT   25
!
!     IMPLICIT REAL ( KIND=8 )  (A-H,O-Z)
      EQUIVALENCE (ZA0, Z0(1)), (ZA1, Z1(1)), (ZA2, Z2(1)), (ZA3, Z3(1)), (A0, &
!      EQUIVALENCE (ZA0, Z0), (ZA1, Z1), (ZA2, Z2), (ZA3, Z3), (A0, &
         Z0(1)), (B0, Z0(2)), (A1, Z1(1)), (B1, Z1(2)), (A2, Z2(1)), (B2, Z2(2)&
         ), (A3, Z3(1)), (B3, Z3(2))
      DATA SQ, SK, CK/ .70710678118655D0, .38268343236509D0, .92387953251129D0&
         /
      DATA TWOPI/ 6.2831853071796D0/
      DATA ZERO/ 0.0D0/
      DATA ONE/ 1.0D0/
!                                  SQ=SQRT2/2,SK=SIN(PI/8),CK=COS(PI/8)
!                                  TWOPI=2*PI
      SYGN = 1.0D0
      IF (MVAL < 0) SYGN = -1.0D0
      M = IABS(MVAL)
      MP = M + 1
      N = 2**M
!     CMPEPSILON used for comparison --  IF (SYGN .NE. 1.0D0) THEN
      IF (ABS(SYGN - 1.0) > ABS(SYGN + 1.0)*CMPEPSILON) THEN
         DO I = 1, N
!        Z(I) = DCONJG(Z(I))
            Z(I) = CONJG(Z(I))
         END DO
      ENDIF
      IWK(1) = 1
      MM = (M/2)*2
      KN = N + 1
!                                  INITIALIZE WORK VECTOR
      DO I = 2, MP
         IWK(I) = IWK(I-1) + IWK(I-1)
      END DO
      RAD = TWOPI/N
      MK = M - 4
      KB = 1
      IF (MM /= M) THEN
         K2 = KN
         K0 = IWK(MM+1) + KB
         K01 = K0
         J2 = MIN0(KB + 1,K01)
         DO K0 = K01, J2, -1
            K2 = K2 - 1
            AK2 = Z(K2)
            Z(K2) = Z(K0-1) - AK2
            Z(K0-1) = Z(K0-1) + AK2
         END DO
      ENDIF
      C1 = ONE
      S1 = ZERO
      JJ = 0
      K = MM - 1
      J = 4
      IF (K >= 1) GO TO 30
      GO TO 9005
   20 CONTINUE
      IF (IWK(J) > JJ) GO TO 25
      JJ = JJ - IWK(J)
      J = J - 1
      IF (IWK(J) > JJ) GO TO 25
      JJ = JJ - IWK(J)
      J = J - 1
      K = K + 2
      GO TO 20
   25 CONTINUE
      JJ = IWK(J) + JJ
      J = 4
   30 CONTINUE
      ISP = IWK(K)
      IF (JJ == 0) GO TO 40
!                                  RESET TRIGONOMETRIC PARAMETERS
      C2 = JJ*ISP*RAD
      C1 = COS(C2)
      S1 = SIN(C2)
   35 CONTINUE
      C2 = C1*C1 - S1*S1
      S2 = C1*(S1 + S1)
      C3 = C2*C1 - S2*S1
      S3 = C2*S1 + S2*C1
   40 CONTINUE
      JSP = ISP + KB
!                                  DETERMINE FOURIER COEFFICIENTS
!                                    IN GROUPS OF 4
      IF (S1 /= ZERO) THEN
         DO I = 1, ISP
            K0 = JSP - I
            K1 = K0 + ISP
            K2 = K1 + ISP
            K3 = K2 + ISP
            ZA0 = Z(K0)
            ZA1 = Z(K1)
            ZA2 = Z(K2)
            ZA3 = Z(K3)
            TEMP = A1
            A1 = A1*C1 - B1*S1
            B1 = TEMP*S1 + B1*C1
            TEMP = A2
            A2 = A2*C2 - B2*S2
            B2 = TEMP*S2 + B2*C2
            TEMP = A3
            A3 = A3*C3 - B3*S3
            B3 = TEMP*S3 + B3*C3
            TEMP = A0 + A2
            A2 = A0 - A2
            A0 = TEMP
            TEMP = A1 + A3
            A3 = A1 - A3
            A1 = TEMP
            TEMP = B0 + B2
            B2 = B0 - B2
            B0 = TEMP
            TEMP = B1 + B3
            B3 = B1 - B3
            B1 = TEMP
!        Z(K0) = DCMPLX(A0+A1,B0+B1)
!        Z(K1) = DCMPLX(A0-A1,B0-B1)
!        Z(K2) = DCMPLX(A2-B3,B2+A3)
!        Z(K3) = DCMPLX(A2+B3,B2-A3)
            Z(K0) = CMPLX(A0 + A1,B0 + B1,KIND = 8)
            Z(K1) = CMPLX(A0 - A1,B0 - B1,KIND = 8)
            Z(K2) = CMPLX(A2 - B3,B2 + A3,KIND = 8)
            Z(K3) = CMPLX(A2 + B3,B2 - A3,KIND = 8)
         END DO
      ELSE
         DO I = 1, ISP
            K0 = JSP - I
            K1 = K0 + ISP
            K2 = K1 + ISP
            K3 = K2 + ISP
            ZA0 = Z(K0)
            ZA1 = Z(K1)
            ZA2 = Z(K2)
            ZA3 = Z(K3)
            TEMP = A0 + A2
            A2 = A0 - A2
            A0 = TEMP
            TEMP = A1 + A3
            A3 = A1 - A3
            A1 = TEMP
            TEMP = B0 + B2
            B2 = B0 - B2
            B0 = TEMP
            TEMP = B1 + B3
            B3 = B1 - B3
            B1 = TEMP
!        Z(K0) = DCMPLX(A0+A1,B0+B1)
!        Z(K1) = DCMPLX(A0-A1,B0-B1)
!        Z(K2) = DCMPLX(A2-B3,B2+A3)
!        Z(K3) = DCMPLX(A2+B3,B2-A3)
            Z(K0) = CMPLX(A0 + A1,B0 + B1,KIND = 8)
            Z(K1) = CMPLX(A0 - A1,B0 - B1,KIND = 8)
            Z(K2) = CMPLX(A2 - B3,B2 + A3,KIND = 8)
            Z(K3) = CMPLX(A2 + B3,B2 - A3,KIND = 8)
         END DO
      ENDIF
      IF (K <= 1) GO TO 55
      K = K - 2
      GO TO 30
   55 CONTINUE
      KB = K3 + ISP
!                                  CHECK FOR COMPLETION OF FINAL
!                                    ITERATION
      IF (KN <= KB) GO TO 9005
      IF (J /= 1) GO TO 60
      K = 3
      J = MK
      GO TO 20
   60 CONTINUE
      J = J - 1
      C2 = C1
      IF (J /= 2) GO TO 65
      C1 = C1*CK + S1*SK
      S1 = S1*CK - C2*SK
      GO TO 35
   65 CONTINUE
      C1 = (C1 - S1)*SQ
      S1 = (C2 + S1)*SQ
      GO TO 35
 9005 CONTINUE
!     CMPEPSILON used for comparison -- IF (SYGN .NE. 1.D0)  THEN
      IF (ABS(SYGN - 1.D0) > ABS(SYGN + 1.D0)*CMPEPSILON) THEN
         XN = N
         DO I = 1, N
!     Z(I)=DCONJG(Z(I))/XN
            Z(I) = CONJG(Z(I))/XN
         END DO
      ENDIF

      CALL QXZ136 (Z, M, IWK)

      RETURN

      END SUBROUTINE FFT


!----------------------------------------------------------------------
      SUBROUTINE QXZ136(Z, M, IWK)

!-FFRDR2--------S-------LIBRARY 3---------------------------------------
!
!***********************************************************************
!
!   FUNCTION            - THIS SUBROUTINE PERMUTES A COMPLEX DATA VECTOR
!                           IN REVERSE BINARY ORDER TO NORMAL ORDER. THE
!                           ROUTINE CAN ALSO BE USED TO PERMUTE A COM-
!                           PLEX DATA VECTOR IN NORMAL ORDER TO REVERSE
!                           BINARY ORDER SINCE THE PERMUTATION IS SYM-
!                           METRIC.
!   USAGE               - CALL QXZ136(Z,M,IWK)
!   PARAMETERS  Z       - COMPLEX VECTOR OF LENGTH N=2**M WHICH
!                           CONTAINS ON INPUT THE DATA TO BE
!                           PERMUTED. ON OUTPUT, Z CONTAINS THE
!                           PERMUTED DATA VECTOR.
!                M      - N=2**M IS THE NUMBER OF DATA POINTS.
!                IWK    - WORK AREA VECTOR OF LENGTH M+1
!   PRECISION           - SINGLE
!   LANGUAGE            - FORTRAN
!   LATEST REVISION     - MARCH 16, 1973
!-----------------------------------------------------------------

      INTEGER :: M
      INTEGER :: IWK(*)
      COMPLEX(DBLE_COMPLEX) :: Z(*)

      INTEGER :: MP, JJ, I, N2, N4, N8=0, LM, NN, J, JK, K
      COMPLEX(DBLE_COMPLEX) :: TEMP

      IF (M > 1) THEN
         MP = M + 1
         JJ = 1
!                                  INITIALIZE WORK VECTOR
         IWK(1) = 1
         DO I = 2, MP
            IWK(I) = IWK(I-1)*2
         END DO
         N4 = IWK(MP-2)
         IF (M > 2) N8 = IWK(MP-3)
         N2 = IWK(MP-1)
         LM = N2
         NN = IWK(MP) + 1
         MP = MP - 4
!                                  DETERMINE INDICES AND SWITCH A*S
         J = 2
   10    CONTINUE
         JK = JJ + N2
         TEMP = Z(J)
         Z(J) = Z(JK)
         Z(JK) = TEMP
         J = J + 1
         IF (JJ <= N4) THEN
            JJ = JJ + N4
         ELSE
            JJ = JJ - N4
            IF (JJ <= N8) THEN
               JJ = JJ + N8
            ELSE
               JJ = JJ - N8
               K = MP
               DO WHILE(IWK(K) < JJ)
                  JJ = JJ - IWK(K)
                  K = K - 1
               END DO
               JJ = IWK(K) + JJ
            ENDIF
         ENDIF
         IF (JJ > J) THEN
            K = NN - J
            JK = NN - JJ
            TEMP = Z(J)
            Z(J) = Z(JJ)
            Z(JJ) = TEMP
            TEMP = Z(K)
            Z(K) = Z(JK)
            Z(JK) = TEMP
         ENDIF
         J = J + 1
!                                  CYCLE REPEATED UNTIL LIMITING NUMBER
!                                    OF CHANGES IS ACHIEVED
         IF (J <= LM) GO TO 10
      ENDIF

      RETURN

      END SUBROUTINE QXZ136


      END MODULE SYNSPEC
