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

      MODULE XSECTIONS

      USE PARAMS
      USE RETVPARAM
      USE VIBFCN
      USE TIPS
      USE ISOTOPE
      USE LINEPARAM
      USE BANDPARAM
      USE MOLCPARAM
      USE VOIGT_SDV_LM
      USE LINESHAPE_PCQSDHC
      IMPLICIT NONE

      INTEGER :: NMONSM, NCROSS
      INTEGER, PARAMETER :: MOLMAXP1 = MOLMAX + 1
      REAL(DOUBLE), DIMENSION(:,:,:), ALLOCATABLE :: CROSS

      CONTAINS

!----------------------------------------------------------------------
      SUBROUTINE KROSSR( NR_LEVEL, ICOUNT )

! --- COMPUTE OPTICAL DEPTHS FOR LAYER NUM_LAYER
! --- COMPUTE INITIAL OPTICAL DEPTHS FOR ALL LAYERS IF NUM_LAYERS = 1
! --- REVISION DATE:  JULY 16, 1996
! --- DECLARE VARIABLES INTEGER(KIND=4) TO AVOID OVERFLOW
! --- IN PRESSURE-BROADENED LOWER LAYERS

      LOGICAL      :: PRTDEBUG = .FALSE.
      INTEGER      :: NR_LEVEL, K_START, K_END, ICOUNT !, NRESET=0
      INTEGER      :: JSTART, JSTOP, N1, K, I, J, INDXX, IBAND, LMIN, LMAX
      INTEGER      :: N, IMOL, NPOINT, MO, ISO
      !INTEGER      :: II, IJ
      REAL(DOUBLE) :: TXE, VIBFAC, STIMFC, SSL, ACOFB, SCOFB, ALOR, ADOP, &
                      AKZERO, YDUM, OPTMAX, XDUM, AKV, OPTCEN, DELLOR, WLIN, START, &
                      SSTOP, ANUZ, QT, QTSTDTEMP, GI, SSLOLD, BETAP, GZ, LMTVAL
      REAL(DOUBLE) :: AKV_R, AKV_I, G2, LM, S0=0.0D0, S2=0.0D0
      REAL(DOUBLE) :: ANUVC = 0.0d0, ETA0=0.0D0

      REAL (DOUBLE), DIMENSION(4) :: SDVLM_PARAM ! PARAMETERS FOR SDV AND/OR LINEMIXING
                                                 ! CALCULATION ACCORDING TO BOONE

      !REAL(DOUBLE) , EXTERNAL :: VOIGT
      !REAL(DOUBLE) , EXTERNAL :: GALATRY
      !REAL(DOUBLE) , EXTERNAL :: BETAT

      INTERFACE
         REAL(DOUBLE) FUNCTION VOIGT (X, Y)
            USE PARAMS
            REAL(DOUBLE) , INTENT(IN) :: X
            REAL(DOUBLE) , INTENT(IN) :: Y
         END FUNCTION VOIGT
         REAL(DOUBLE) FUNCTION GALATRY (X, Y, Z)
            USE PARAMS, ONLY:  DOUBLE, DBLE_COMPLEX
            REAL(DOUBLE)  :: X
            REAL(DOUBLE)  :: Y
            REAL(DOUBLE) , INTENT(IN) :: Z
         END FUNCTION GALATRY
         REAL(DOUBLE) FUNCTION BETAT (MOL, T)
            USE PARAMS
            INTEGER , INTENT(IN) :: MOL
            REAL(DOUBLE) , INTENT(IN) :: T
            END FUNCTION BETAT
      END INTERFACE

      IF (ICOUNT.EQ.1) PRINT *, ' COMPUTING CROSS-SECTIONS...'

      GI = 0
      QT = 0.0D0
      QTSTDTEMP = 0.0D0
!print*,'kro'
!  --- INITIALIZE X-SECTION ARRAYS
      N1 = NRET + 1

      IF (NR_LEVEL.LT.0) THEN
         CROSS(:N1,:KMAX,:NCROSS) = 0.D0
         K_START = 1
         K_END = KMAX + NCELL !NPATH
      ELSE
         !PRINT*, 'LAYER : ', NR_LEVEL
         ! RUN ON TWO LEVELS - THE CURRENT AND THE PREVIOUS TO UN-PRETURB IT
         ! *** ASSUME PERTURBATION IN CONSECUTIVE LEVEL ORDER ****
         K_START = NR_LEVEL
         IF (NR_LEVEL .GT. 1 ) K_START = K_START - 1
         K_END = NR_LEVEL
         IF (K_END .GT. (KMAX + NCELL)) K_END=KMAX + NCELL ! final iteration in FM has NR_LEVEL=KMAX+1, and should undo the lowest level perturbation
         CROSS(:N1,K_START:K_END,:NCROSS) = 0.D0
      END IF
      INDXX = 0
!                    ------------ LOOP OVER BANDPASSES
      DO IBAND = 1, NBAND
!                    ------------ LOOP OVER LAYERS

         !print*, 'kro 1 ', k_start, k_end
         DO K = K_START, K_END
!            print *, 'icount', icount, 'level', K
!                    ------------ LOOP OVER SPECTRAL LINES

!print*, 'kro 2 ', k, T(K), P(K)

            LMIN = LINE1(IBAND)
            LMAX = LINE2(IBAND)
            DO N = LMIN, LMAX
!print*, 'kro 3 ',iband, k, n, lmin, lmax

!  --- SELECT DISTANCE FROM LINE CENTER FOR CALCULATIONS
               IF( IFMIX(LGAS(N)) == 0 ) CYCLE
               IMOL = LGAS(N)
               TXE  = RCONST2*ETWO(N)*(1.D0/T(K)-1.D0/STDTEMP)
!              ... DEFAULT MOL AND ISO ID
               MO  = ICODE(IMOL)
               ISO = ISCODE(IMOL)
!print*, 'kro 4 ', mo
!  --- CHECK FOR ISOTOPE SEPARATION
!  ---  eg TIPS IS NOT AWARE OF ON-THE-FLY ISOTOPE SEPARATION
               DO I=1, NISOSEP
                   IF((MO .EQ. NEWID(I)) .AND. (ISO .EQ. NEWISO(I)))THEN
                       MO  = OLDID(I)
                       ISO = OLDISO(I)
                   ENDIF
               ENDDO

               IF( USE_TIPS )THEN
                  CALL BD_TIPS_2017(MO, STDTEMP, ISO, GI, QTSTDTEMP)
                  CALL BD_TIPS_2017(MO, T(K), ISO, GI, QT)
               ELSE
                  QT = -1.0
               ENDIF

!              --- STIMULATED EMISSION CORRECTION TO LINE INTENSITY
               STIMFC = (1.D0 - EXP((-RCONST2*AZERO(N)/T(K))))/(1.D0 - EXP((-RCONST2*AZERO(N)/STDTEMP)))

               IF (PRTDEBUG) THEN
                  PRINT *
                  PRINT *,"IBAND,K,N,ICODE,ISCODE,MOL,ISO =", IBAND,K,N, ICODE(IMOL),ISCODE(IMOL), MO,ISO
                  PRINT *,"               T(K), QTSTD, QT = ", T(K), QTSTDTEMP, QT
                  VIBFAC = QV(IMOL,KMAX+1)/QV(IMOL,K)
                  SSLOLD = ST296(N)*(STDTEMP/T(K))*(STDTEMP/T(K))**TDEP(ICODE(IMOL))*VIBFAC*STIMFC*EXP((-TXE))
               ENDIF

               IF (QTSTDTEMP <= 0.0 .OR. QT <= 0.0) THEN
!              --- USE NON-TIPS METHOD VIBRATIONAL PARTITION FUNCTION IF SPECIES IS NOT INCLUDED IN TIPS
                  VIBFAC = QV(IMOL,KMAX+1)/QV(IMOL,K)
                  SSL = ST296(N)*(STDTEMP/T(K))*(STDTEMP/T(K))**TDEP(ICODE(IMOL))*VIBFAC*STIMFC*EXP((-TXE))
                  IF (PRTDEBUG) THEN
                     PRINT *,"         T RATIO: TIPS, VIBFCN = ", QTSTDTEMP/QT, (STDTEMP/T(K))**TDEP(ICODE(IMOL))*VIBFAC
                     PRINT *,"              USING VIBFCN SSL = ", SSL
                  ENDIF
               ELSE
!               --- USE TIPS
                  SSL = ST296(N)*(STDTEMP/T(K))*(QTSTDTEMP/QT)*STIMFC*EXP((-TXE))
                  IF (PRTDEBUG) THEN
                     PRINT *,"         T RATIO: TIPS, VIBFCN = ", QTSTDTEMP/QT, (STDTEMP/T(K))**TDEP(ICODE(IMOL))*VIBFAC
                     PRINT *,"                USING TIPS SSL = ", SSL
                  ENDIF
               ENDIF

               IF (PRTDEBUG) THEN
                  PRINT *, "                      TIPS SSL = ", SSL
                  PRINT *, "                    VIBFCN SSL = ", SSLOLD
               ENDIF

! --- PUT MO AND ISO BACK TO NON-ISOTOPE SEPARATION AND TO 49 FROM 7 FOR CIA
!               MO  = ICODE(IMOL)
!               ISO = ISCODE(IMOL)

               ! RESET SPEED DEPENDANCY PARAMETER
               G2 = 0.0D0
               ! --- ACCOUNT FOR O2
               IF( HFLAG(N,FCIA_FLAG) .OR. HFLAG(N,SCIA_FLAG) .OR. MO.EQ.7 )THEN
                  ACOFB = AAA(N) + (XGAS(IMOL,K)-0.21D0)*(SSS(N)-AAA(N))/0.79D0
                  ACOFB = ACOFB * P(K)
                  SCOFB = 0.0D0
               ELSEIF  ( LSM_SDV.and.HFLAG(N,SDV_FLAG) ) THEN
                  ACOFB = GAMMA0(N)*P(K)*(1.0D0 - XGAS(IMOL,K))
                  SCOFB = SSS(N)*P(K)*XGAS(IMOL,K)
                  G2 = GAMMA2(N)*P(k) * (1.0D0 - XGAS(IMOL,K))! not yet implemented: + SELF_GAMMA2(N)*P(K)*XGAS(IMOL,K)
               ELSE
                  ACOFB = AAA(N)*P(K)*(1.0D0 - XGAS(IMOL,K))
                  SCOFB = SSS(N)*P(K)*XGAS(IMOL,K)
               ENDIF

               ! PRESSURE SHIFT PARAMETER FOR pCqSD
               IF ( FPS.AND.(LSHAPEMODEL.EQ.4) ) THEN
                  IF ( .false..and.HFLAG(N,SDV_FLAG) ) THEN
                     S0 = SHIFT0(N)*P(K) * (1.0D0 - XGAS(IMOL,K))! + SELF_SHIFT0(N)*P(K) * XGAS(IMOL,K)          ! PRESSURE SHIFT IF SDV IS USED
                  ELSE
                     S0 = PSLIN(N)*P(K) * (1.0D0 - XGAS(IMOL,K))! + SELF_SHIFT0(N)*P(K) * XGAS(IMOL,K)        ! PRESSURE SHIFT IF SDV IS NOT USED
                  END IF
                  S2 = SHIFT2(N)*P(k)* (1.0D0 - XGAS(IMOL,K))    ! PRESSURE SHIFT OF GAMMA2
               ELSE
                  S0 = 0.0D0
                  S2 = 0.0D0
               END IF


               ALOR = (ACOFB + SCOFB)*(STDTEMP/T(K))**TDLIN(N)
               ADOP = RFACTOR*SQRT(T(K)/GMASS(N))*AZERO(N)

               IF( HFLAG(N,FCIA_FLAG) .OR. HFLAG(N,SCIA_FLAG) )ADOP=1.0D0

               AKZERO = SSL/ADOP*ALOGSQ/PISQ

               IF( HFLAG(N,FCIA_FLAG) )AKZERO=AKZERO*(1.0D0 - XGAS(IMOL,K))*PMASMX(K)
               IF( HFLAG(N,SCIA_FLAG) )AKZERO=AKZERO*XGAS(IMOL,K)*PMASMX(K)

               YDUM = ALOGSQ*ALOR/ADOP

               ! SPEED DEPENDENT VOIGT
               SDVLM_PARAM(1:4) = 0.0D0
               IF (HFLAG(N,SDV_FLAG)) THEN
                  SDVLM_PARAM(1) = GAMMA2(N)*P(K) ! ASYMMETRY FOR SDV (MIXING COEFFICIENT)
                  SDVLM_PARAM(2) = GAMMA0(N)      ! PRESSURE BROADENING FOR SDV
                  SDVLM_PARAM(3) = 0.0!ETA2(N)    ! PRESSURE SHIFT FOR SPEED DEPENDENT VOIGT
               END IF
               ! LINE MIXING PARAMETRIZATION ACCORDING TO HASE.  FOR CO2 ONLY?
               IF (HFLAG(N,LM_FLAG)) THEN
                  SDVLM_PARAM(2) = (ACOFB + SCOFB)   ! PRESSURE BROADENING FOR VOIGT
                  LMTVAL = (T(K) - 260.0D0) / 60.0D0 ! LM-REF TEMPERATURES: 200/260/320 K
                  SDVLM_PARAM(4) = YLM(N) * (1.0D0 + LMTVAL *(LMTK1(N) + LMTVAL * LMTK2(N)))
                  LM = SDVLM_PARAM(4)*P(K)
               END IF

! --- CHECK FOR ZERO PARAMETERS IN SVD_PARAM ESPECIALLY 3 (ETA2) MUST BE NON ZERO IF THIS IS NOT THE CASE,
! --- USE VOIGT LINESHAPE -- NOTE LINEMIXING PROMPTS USE OF SDVMIX OR VOIGTMIX !
               ! IF(( ABS(SDVLM_PARAM(3)) .LT. TINY(0.0D0) ) .and. &
               !      ( HFLAG(N,SDV_FLAG)) .and. &
               !      ( K .EQ. K_START ))THEN
               !    NRESET = NRESET + 1
               !    WRITE(0,100) AZERO(N), HFLAG(N,1:8), SDVLM_PARAM(1:4), NRESET
               !    HFLAG(N,SDV_FLAG) = .FALSE.
               ! ENDIF
               IF ((LSHAPEMODEL.eq.3).and.(HFLAG(N,SDV_FLAG) .OR. HFLAG(N,LM_FLAG)) )THEN
                  CALL SDV_MISC(ALOGSQ/ADOP, STDTEMP/T(K), P(K), TDLIN(N), SDVLM_PARAM)
               ENDIF
               DO I = 1, NRET + NCELL
                  IF (LGAS(N) /= IRET(I)) CYCLE
                  NPOINT = I
                  OPTMAX = PMASMX(K)*XORG(I,K)
                  !print*, 'kro xorg(ik) ',i,k, XORG(I,K), PMASMX(K), optmax
                  GO TO 349
               END DO
               NPOINT = NRET + 1
               OPTMAX = PMASMX(K)*XGAS(IMOL,K)
!  --- CALCULATE LINE CENTER OPTICAL DEPTH
  349          CONTINUE
               XDUM = 0.D0
               IF((LSHAPEMODEL.EQ.2).and.HFLAG(N,GALATRY_FLAG)) THEN
                  BETAP = BETA(N)*P(K)
                  BETAP = BETAP * BETAT(ICODE(IMOL),T(K))
                  GZ = ALOGSQ*BETAP/ADOP
                  AKV = AKZERO*GALATRY(XDUM,YDUM,GZ)
               ELSEIF((LSHAPEMODEL.EQ.3).and.HFLAG(N,SDV_FLAG).and.HFLAG(N,LM_FLAG)) THEN
                  CALL SDVMIX(XDUM*ADOP/ALOGSQ,AKV)
                  AKV = AKZERO * AKV
               ELSEIF((LSHAPEMODEL.EQ.3).AND.(HFLAG(N,LM_FLAG))) THEN
                  CALL VOIGTMIX(XDUM*ADOP/ALOGSQ,AKV)
                  AKV = AKZERO * AKV
               ELSEIF(LSHAPEMODEL.EQ.4) THEN
                  ! pCqSDHC MODEL (Tran)
                  call pCqSDHC(azero(N),ADOP,ALOR,G2,S0, S2, ANUVC, ETA0, azero(N),AKV_R,AKV_I)
                  AKV = SSL * (AKV_R + LM * AKV_I)! ALL OTHER PARTS OF AKZERO ARE ALREADY PART OF
                  ! AKV_R AND AKV_I
               ELSE
                  AKV = AKZERO * VOIGT(XDUM,YDUM)
               ENDIF
               OPTCEN = AKV*OPTMAX


!print*, 'kro 7 ',dist, DELLOR, SSL, ALOR, OPTMAX, TAUMIN
!IF (DELLOR > DELNU) print*, DELLOR

               !  --- EXTEND CALCULATIONS FURTHER INTO THE WINGS IF NECESSARY
               IF (ICOUNT.EQ.1) THEN
                  !  --- SKIP OVER LINE IF OPTICAL DEPTH AT LINE CENTER IS LESS
                  !  --- THAN TAUMIN
                  DIST(N,K) = -1.0D0
                  IF (OPTCEN.GE.TAUMIN) THEN
                     !print*, 'kro 6 here'
                     DIST(N,K) = DELNU
                     !  --- CALCULATE DISTANCE FROM LINE CENTER CORRESPONDING TO OPTICAL
                     !  --- DEPTH OF TAUMIN FOR A LORENTZ LINE
                     DELLOR = SQRT(SSL*ALOR*OPTMAX/(TAUMIN*PI))
                     IF (DELLOR > DELNU) DIST(N,K) = DELLOR
                  ELSE ! THIS NEEDS TO BE HERE, THE TEST 4 LINES DOWN DOES NOT SUFFICE, DONT KNOW WHY. -- MP
                     CYCLE
                  END IF
               ENDIF
               IF (DIST(N,K).LE.0.0D0) CYCLE
!                  print *, icount, k, N
!               end IF
               !  --- CORRECT POSITION FOR PRESSURE SHIFT
               WLIN = AZERO(N) + P(K)*PSLIN(N)
!print*, 'kro 8 azero ', AZERO(N), P(K), PSLIN(N)
!  --- IF NO PRESSURE SHIFT, pCqSDHC calculates it own pressure shift
               IF( (.NOT. FPS).OR.(LSHAPEMODEL.EQ.4) ) WLIN = AZERO(N)
               START = WLIN - DIST(N,K)
               SSTOP = WLIN + DIST(N,K)
               JSTART = FLOOR((START - WMON(IBAND))/DN(IBAND) + 1.00000001D0)
               JSTOP = FLOOR((SSTOP - WMON(IBAND))/DN(IBAND) + 1.00000001D0)
!print*, 'kro 9 ',WLIN, dist, START, sstop, jstart, jstop, WMON(IBAND)
               IF (JSTOP < 1) CYCLE
               IF (JSTART > NM(IBAND)) CYCLE
               JSTART = MAX0(1,JSTART)
               JSTOP = MIN0(NM(IBAND),JSTOP)
               DO J = JSTART, JSTOP
                  ANUZ = WMON(IBAND) + (J - 1)*DN(IBAND)
                  XDUM = ALOGSQ*(ANUZ - WLIN)/ADOP
                  IF ((LSHAPEMODEL.EQ.2).AND.HFLAG(N,GALATRY_FLAG)) THEN
                     XDUM = ABS(XDUM)
                     AKV  = AKZERO*GALATRY(XDUM,YDUM,GZ)
                  ELSEIF((LSHAPEMODEL.EQ.3).and.HFLAG(N,SDV_FLAG).and.HFLAG(N,LM_FLAG)) THEN
                     CALL SDVMIX(XDUM*ADOP/ALOGSQ,AKV)
                     AKV = AKZERO * AKV
                  ELSEIF((LSHAPEMODEL.EQ.3).AND.HFLAG(N,LM_FLAG)) THEN
                     CALL VOIGTMIX(XDUM*ADOP/ALOGSQ,AKV)
                     AKV = AKZERO * AKV
                  ELSEIF ((LSHAPEMODEL.EQ.4).OR.HFLAG(N,SDV_FLAG).OR.HFLAG(N,LM_FLAG)) THEN
                     ! pCqSDHC MODEL (Tran)
                     call pCqSDHC(WLIN,ADOP,ALOR,G2,S0, S2, ANUVC,ETA0,&
                          ANUZ,AKV_R,AKV_I)
                     AKV = SSL * (AKV_R + LM * AKV_I)! ALL OTHER PARTS OF AKZERO ARE ALREADY PART OF
                     ! AKV_R AND AKV_I
                  ELSE
                     XDUM = ABS(XDUM)
                     AKV  = AKZERO*VOIGT(XDUM,YDUM)
                  ENDIF
                  CROSS(NPOINT,K,J+INDXX) = CROSS(NPOINT,K,J+INDXX) + AKV*OPTMAX
                  !print*, npoint, indxx, j, k, CROSS(NPOINT,K,J+INDXX), AKV, OPTMAX, xdum, ydum
               ENDDO ! J
            ENDDO ! LINES
            !print*, 'kro 10 ',j-1, NPOINT, K, J, INDXX !, CROSS(NPOINT,K,J-1+INDXX)
         ENDDO ! LAYERS
         INDXX = INDXX + NM(IBAND)
      ENDDO !BANDS

      !WRITE (*, 1001)
! 1000 FORMAT('  Beginning cross section calculations...')
! 1001 FORMAT('  Cross section calculations completed.')

      RETURN

 !100  FORMAT(F12.5, 1X, 8L2, 4D15.7, " CLEAR SVD & LM FLAG : ", I5 )

      END SUBROUTINE KROSSR

      END MODULE XSECTIONS
