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
      USE GEN_LINESHAPE

      IMPLICIT NONE

      INTEGER :: NMONSM, NCROSS
      INTEGER, PARAMETER :: MOLMAXP1 = MOLMAX + 1
      REAL(DOUBLE), DIMENSION(:,:,:), ALLOCATABLE :: CROSS

      CONTAINS

!----------------------------------------------------------------------
      SUBROUTINE KROSSR( NR_LEVEL )

! --- COMPUTE OPTICAL DEPTHS FOR LAYER NUM_LAYER
! --- COMPUTE INITIAL OPTICAL DEPTHS FOR ALL LAYERS IF NUM_LAYERS = 1
! --- REVISION DATE:  JULY 16, 1996
! --- DECLARE VARIABLES INTEGER(KIND=4) TO AVOID OVERFLOW
! --- IN PRESSURE-BROADENED LOWER LAYERS

      LOGICAL      :: PRTDEBUG = .FALSE.
      INTEGER      :: NR_LEVEL, K_START, K_END, NRESET=0
      INTEGER      :: JSTART, JSTOP, N1, K, I, J, INDXX, IBAND, LMIN, LMAX
      INTEGER      :: N, IMOL, NPOINT, MO, ISO
      REAL(DOUBLE) :: DIST, TXE, VIBFAC, STIMFC, SSL, ACOFB, SCOFB, ALOR, ADOP, &
                      AKZERO, YDUM, OPTMAX, XDUM, AKV, OPTCEN, DELLOR, WLIN, START, &
                      SSTOP, ANUZ, QT, QTSTDTEMP, GI, SSLOLD, BETAP, GZ, LMTVAL
      REAL(DOUBLE) :: AKV_R, AKV_I, G2, LM, S0, S2, ANUVC, ETA

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

      GI = 0
      QT = 0.0D0
      QTSTDTEMP = 0.0D0

!  --- INITIALIZE X-SECTION ARRAYS
      N1 = NRET + 1

      IF (NR_LEVEL.LT.0) THEN
         CROSS(:N1,:KMAX,:NCROSS) = 0.D0
         K_START = 1
         K_END = KMAX
      ELSE
         !PRINT*, 'LAYER : ', NR_LEVEL
         K_START = NR_LEVEL
         K_END = NR_LEVEL
         CROSS(:N1,K_START:K_END,:NCROSS) = 0.D0
      END IF
      INDXX = 0
!                    ------------ LOOP OVER BANDPASSES
      DO IBAND = 1, NBAND
!                    ------------ LOOP OVER LAYERS
         DO K = K_START, K_END
!                    ------------ LOOP OVER SPECTRAL LINES
            LMIN = LINE1(IBAND)
            LMAX = LINE2(IBAND)
            DO N = LMIN, LMAX


!  --- SELECT DISTANCE FROM LINE CENTER FOR CALCULATIONS
               DIST = DELNU
               IF( IFMIX(LGAS(N)) == 0 ) CYCLE
               IMOL = LGAS(N)
               TXE  = RCONST2*ETWO(N)*(1.D0/T(K)-1.D0/STDTEMP)
!              ... DEFAULT MOL AND ISO ID
               MO  = ICODE(IMOL)
               ISO = ISCODE(IMOL)

!  --- CHECK FOR ISOTOPE SEPARATION
               DO I=1, NISOSEP
                   IF((MO .EQ. NEWID(I)) .AND. (ISO .EQ. NEWISO(I)))THEN
                       MO  = OLDID(I)
                       ISO = OLDISO(I)
                   ENDIF
               ENDDO

! --- CASE FOR CIA O2 (49/1 FCIA, 49/2 SCIA) DUMMY BACK TO O2 FOR PART. FCN.
!               IF( HFLAG(N,FCIA_FLAG) .OR. HFLAG(N,SCIA_FLAG) )THEN
!                  MO  = 7
!                  ISO = 1
!               ENDIF

               CALL BD_TIPS_2008_SFIT2(MO, STDTEMP, ISO, GI, QTSTDTEMP)
               CALL BD_TIPS_2008_SFIT2(MO, T(K), ISO, GI, QT)
               STIMFC = (1.D0 - EXP((-RCONST2*AZERO(N)/T(K))))/(1.D0 - EXP((-RCONST2*AZERO(N)/STDTEMP)))
!               VIBFAC = QV(IMOL,KMAX+1)/QV(IMOL,K)
!               SSLOLD = ST296(N)*(STDTEMP/T(K))*(STDTEMP/T(K))**TDEP(ICODE(IMOL))*VIBFAC*STIMFC*EXP((-TXE))

               IF (PRTDEBUG) THEN
                  PRINT *
                  PRINT *,"IBAND,K,N,ID1,ID2,QT,QTSTD,T(K)=", &
                           IBAND,K,N,ICODE(IMOL),ISCODE(IMOL),QT,QTSTDTEMP,T(K)
                  PRINT *,"   VIB PART RESULT =", &
                       (STDTEMP/T(K))*(STDTEMP/T(K))**TDEP(ICODE(IMOL))*VIBFAC
               ENDIF

               IF (QTSTDTEMP <= 0.0 .OR. QT <= 0.0) THEN
!                 ... USE NON-TIPS METHOD VIBRATIONAL PARTITION FUNCTION
                  VIBFAC = QV(IMOL,KMAX+1)/QV(IMOL,K)
!               --- STIMULATED EMISSION CORRECTION TO LINE INTENSITY
                  SSL = ST296(N)*(STDTEMP/T(K))*(STDTEMP/T(K))**TDEP(ICODE(IMOL))*VIBFAC*STIMFC*EXP((-TXE))
                  !SSL = SSLOLD
                  IF (PRTDEBUG) THEN
                     PRINT *,"   VIB PARTITION SSL =",SSL
                  ENDIF
               ELSE
!                 ... USE TIPS
                  SSL = ST296(N)*(STDTEMP/T(K))*(QTSTDTEMP/QT)*STIMFC*EXP((-TXE))
                  IF (PRTDEBUG) THEN
                     PRINT *,"   QTSTDTEMP/QT=",QTSTDTEMP/QT
                     PRINT *,"   TIPS SSL, VIB PARTITION SSL =",SSL,SSLOLD
                  ENDIF
               ENDIF

! --- PUT MO AND ISO BACK TO NON-ISOTOPE SEPARATION AND TO 49 FROM 7 FOR CIA
!               MO  = ICODE(IMOL)
!               ISO = ISCODE(IMOL)

! --- ACCOUNT FOR O2
               SDVLM_PARAM(1:4) = 0.0D0 ! STRUCTURE FOR INITIALISATION OF BOONE'S  
                                        ! SDV AND LM ROUTINES
               G2 = 0.0D0
               ! IF PCQSDHC LINEMODEL IS USED AND NO SPEED DEPENDENCY, S0 CONTAINS THE FREQUENCY SHIFT 
               IF ( FPS ) S0 = PSLIN(N)*P(K)
               IF( HFLAG(N,FCIA_FLAG) .OR. HFLAG(N,SCIA_FLAG) .OR. MO.EQ.7 )THEN
                  ACOFB = AAA(N) + (XGAS(IMOL,K)-0.21D0)*(SSS(N)-AAA(N))/0.79D0
                  ACOFB = ACOFB * P(K)
                  SCOFB = 0.0D0
               ELSEIF  ( LSM_SDV.and.LSHAPEMODEL.NE.3.and.HFLAG(N,SDV_FLAG) ) THEN
                  ACOFB = GAMMA0(N)*P(K)*(1.0D0 - XGAS(IMOL,K))
                  SCOFB = SSS(N)*P(K)*XGAS(IMOL,K)
                  G2 = GAMMA0(N)*GAMMA2(N)*P(k) ! Note HITRAN gives data by Devi, the function 
                                                ! implemented here (Tran) defines this parameter differently
                  IF ( FPS ) THEN
                     S0 = SHIFT0(N)*P(K)
                     S2 = SHIFT2(N)*P(k)*T(K)
                  ELSE
                     S0 = 0.0D0
                     S2 = 0.0D0
                  END IF
               ELSEIF ( LSHAPEMODEL.EQ.3.AND.HFLAG(N,SDV_FLAG) ) THEN
                  ! Boone's implementation of SDV
                  SDVLM_PARAM(1) = GAMMA2(N)*P(k) ! PRESSURE NARROWING FOR SDV
                  SDVLM_PARAM(2) = ACOFB + SCOFB  ! PRESSURE BROADENING FOR SDV
                  SDVLM_PARAM(3) = 0.0!ETA2(N)    ! PRESSURE SHIFT IN SDV
               ELSE
                  ACOFB = AAA(N)*P(K)*(1.0D0 - XGAS(IMOL,K))
                  SCOFB = SSS(N)*P(K)*XGAS(IMOL,K)
               END IF

               ANUVC = 0.0D0
               IF ( LSM_DICKE.AND.LSHAPEMODEL.NE.2.AND.HFLAG(N,GALATRY_FLAG) ) THEN 
                  ! Tran's implemeatation of Line shape
                  BETAP = BETA(N)*P(K)
                  ! RELATION GIVEN BY NGO ET.AL. 'AN ISOLATED LINE-SHAPE MODEL ...', JQRST, 2013
                  ANUVC = BETAP !3.0D0/4.0D0*BETAP
               END IF

               eta = 0.0d0
               if(LSM_CORR.and.LSHAPEMODEL.EQ.0) then
                  ! not yet implementated because no measured line parameters exist.
                  eta = 0.0d0
               end if


               ALOR = (ACOFB + SCOFB)*(STDTEMP/T(K))**TDLIN(N)
               ADOP = RFACTOR*SQRT(T(K)/GMASS(N))*AZERO(N)

               IF( HFLAG(N,FCIA_FLAG) .OR. HFLAG(N,SCIA_FLAG) )ADOP=1.0D0

               AKZERO = SSL/ADOP*ALOGSQ/PISQ

               LM = 0.0D0
               ! LINE MIXING PARAMETRIZATION ACCORDING TO FRANK HASE.
               IF (LSHAPEMODEL.NE.3.AND.HFLAG(N,LM_FLAG) ) THEN 
                     LMTVAL = (T(K) - 260.0D0) / 60.0D0 ! LM-REF TEMPERATURES: 200/260/320 K
                     LM = YLM(N) * (1.0D0 + LMTVAL *(LMTK1(N) + LMTVAL * LMTK2(N)))
                     LM = LM*P(K)
               ELSEIF (HFLAG(N,LM_FLAG)) THEN
                  SDVLM_PARAM(2) = (ACOFB + SCOFB)   ! PRESSURE BROADENING FOR BOONE'S VOIGT
                  LMTVAL = (T(K) - 260.0D0) / 60.0D0 ! LM-REF TEMPERATURES: 200/260/320 K
                  SDVLM_PARAM(4) = YLM(N) * (1.0D0 + LMTVAL *(LMTK1(N) + LMTVAL * LMTK2(N)))
                  !                     print *, AZERO(N), LMTK1(N), LMTK2(N), YLM(N), SDVLM_PARAM(4)*P(K)
               END IF



               IF( HFLAG(N,FCIA_FLAG) )AKZERO=AKZERO*(1.0D0 - XGAS(IMOL,K))*PMASMX(K)
               IF( HFLAG(N,SCIA_FLAG) )AKZERO=AKZERO*XGAS(IMOL,K)*PMASMX(K)

               YDUM = ALOGSQ*ALOR/ADOP

                  
! --- CHECK FOR ZERO PARAMETERS IN SVD_PARAM ESPECIALLY 1 (GAMMA2) MUST BE NON ZERO IF THIS IS NOT THE CASE,
! --- USE VOIGT LINESHAPE. THE MODEL 4 IS NOT AFFECTED BY THIS
               IF (LSHAPEMODEL.EQ.3) THEN
                  IF(( ABS(SDVLM_PARAM(1)) .LT. 1.0d-10) .and. &
                       ( HFLAG(N,SDV_FLAG)) .and. &
                       ( K .EQ. K_START ))THEN
                     NRESET = NRESET + 1
                     WRITE(0,100) AZERO(N), HFLAG(N,1:8), SDVLM_PARAM(1:4), NRESET
                     HFLAG(N,SDV_FLAG) = .FALSE.
                  ENDIF
                  IF ( HFLAG(N,SDV_FLAG) .OR. HFLAG(N,LM_FLAG) )THEN
                     CALL SDV_MISC(ALOGSQ/ADOP, STDTEMP/T(K), P(K), TDLIN(N), SDVLM_PARAM)
                  ENDIF
               END if
               DO I = 1, NRET
                  IF (LGAS(N) /= IRET(I)) CYCLE
                  NPOINT = I
                  OPTMAX = PMASMX(K)*XORG(I,K)
                  GO TO 349
               END DO
               NPOINT = NRET + 1
               OPTMAX = PMASMX(K)*XGAS(IMOL,K)
!  --- CALCULATE LINE CENTER OPTICAL DEPTH
  349          CONTINUE
               XDUM = 0.D0
               IF(LSHAPEMODEL.EQ.2.and.HFLAG(N,GALATRY_FLAG)) THEN
!               IF(HFLAG(N,GALATRY_FLAG)) THEN
                  BETAP = BETA(N)*P(K)
                  BETAP = BETAP * BETAT(ICODE(IMOL),T(K))
                  GZ = ALOGSQ*BETAP/ADOP
                  AKV = AKZERO*GALATRY(XDUM,YDUM,GZ)
               ELSEIF(LSHAPEMODEL.EQ.3.AND.HFLAG(N,SDV_FLAG)) THEN
                  ! SDV WITH OR WITHOUT LINE MIXING (Boone) 
                  CALL SDVMIX(XDUM*ADOP/ALOGSQ,AKV)
                  AKV = AKZERO * AKV
               ELSEIF(LSHAPEMODEL.EQ.3.AND.HFLAG(N,LM_FLAG)) THEN
                  ! VOIGT WITH LINE MIXING (Boone) 
                  CALL VOIGTMIX(XDUM*ADOP/ALOGSQ,AKV_R, AKV_I)
                  AKV = AKZERO * (AKV_R +SDVLM_PARAM(4)*AKV_I) 
               ELSEIF (LSHAPEMODEL.EQ.1) THEN
                  ! VOIGT WITHOUT LINEMIXING
                  AKV = AKZERO * VOIGT(XDUM,YDUM)
               ELSEIF (LSHAPEMODEL.EQ.0) THEN
                  ! pCqSDHC MODEL (Tran)
                  call pCqSDHC(azero(N),ADOP,ALOR,G2,S0, S2, ANUVC,ETA,&
                       azero(N),AKV_R,AKV_I)
                  AKV = SSL * (AKV_R + LM * AKV_I)! ALL OTHER PARTS OF AKZERO ARE ALREADY PART OF 
                                                  ! AKV_R AND AKV_I
               ELSE
                  write(*,*) 'Unknow line shape model: ', LSHAPEMODEL
                  write(16,*) 'Unknow line shape model: ', LSHAPEMODEL
               END IF
               OPTCEN = AKV*OPTMAX

!  --- SKIP OVER LINE IF OPTICAL DEPTH AT LINE CENTER IS LESS
!  --- THAN TAUMIN
               IF (OPTCEN < TAUMIN) CYCLE
!  --- CALCULATE DISTANCE FROM LINE CENTER CORRESPONDING TO OPTICAL
!  --- DEPTH OF TAUMIN FOR A LORENTZ LINE
               DELLOR = SQRT(SSL*ALOR*OPTMAX/(TAUMIN*PI))
!  --- EXTEND CALCULATIONS FURTHER INTO THE WINGS IF NECESSARY
               IF (DELLOR > DELNU) DIST = DELLOR
!  --- CORRECT POSITION FOR PRESSURE SHIFT 
               WLIN = AZERO(N) + P(K)*PSLIN(N)
!!  --- IF NO PRESSURE SHIFT
               IF( .NOT. FPS ) WLIN = AZERO(N)
!!  --- If using the pCqSDHC line shape, the pressure shifted is already included
               IF( LSHAPEMODEL.EQ.0 ) WLIN = AZERO(N)
               START = WLIN - DIST
               SSTOP = WLIN + DIST
               JSTART = FLOOR((START - WMON(IBAND))/DN(IBAND) + 1.00000001D0)
               JSTOP = FLOOR((SSTOP - WMON(IBAND))/DN(IBAND) + 1.00000001D0)
               IF (JSTOP < 1) CYCLE
               IF (JSTART > NM(IBAND)) CYCLE
               JSTART = MAX0(1,JSTART)
               JSTOP = MIN0(NM(IBAND),JSTOP)
               DO J = JSTART, JSTOP
                  ANUZ = WMON(IBAND) + (J - 1)*DN(IBAND)
                  XDUM = ALOGSQ*(ANUZ - WLIN)/ADOP
                  IF (LSHAPEMODEL.EQ.2.AND.HFLAG(N,GALATRY_FLAG)) THEN
!                  IF (HFLAG(N,GALATRY_FLAG)) THEN
                     XDUM = ABS(XDUM)
                     AKV  = AKZERO*GALATRY(XDUM,YDUM,GZ)
                  ELSEIF(LSHAPEMODEL.EQ.3.AND.HFLAG(N,SDV_FLAG)) THEN
                     ! SDV WITH OR WITHOUT LINE MIXING 
                     CALL SDVMIX(XDUM*ADOP/ALOGSQ,AKV)
                     AKV = AKZERO * AKV
                  ELSEIF(LSHAPEMODEL.EQ.3.AND.HFLAG(N,LM_FLAG)) THEN
                     ! VOIGT WITH LINE MIXING 
                     CALL VOIGTMIX(XDUM*ADOP/ALOGSQ,AKV_R, AKV_I)
                     AKV = AKZERO * (AKV_R +SDVLM_PARAM(4)*P(K)*AKV_I) 
                  ELSEif(LSHAPEMODEL.EQ.1) THEN
                     ! VOIGT WITHOUT LINEMIXING
                     XDUM = ABS(XDUM)
                     AKV  = AKZERO*VOIGT(XDUM,YDUM)
                  ELSEIF (LSHAPEMODEL.EQ.0) THEN
                     call pCqSDHC(WLIN,ADOP,ALOR,G2,S0,S2,ANUVC,ETA,&
                          ANUZ,AKV_R,AKV_I)
                     AKV = SSL*(AKV_R + LM*AKV_I)! ALL OTHER PARTS OF AKZERO ARE ALREADY PART OF 
                     ! AKV_R AND AKV_I
                  END IF
                  CROSS(NPOINT,K,J+INDXX) = CROSS(NPOINT,K,J+INDXX) + AKV*OPTMAX
                  !                  print*, j, k, CROSS(NPOINT,K,J+INDXX)
               ENDDO ! J
            ENDDO ! LINES
            !print*, j-1, CROSS(NPOINT,K,J-1+INDXX)
         ENDDO ! LAYERS
         INDXX = INDXX + NM(IBAND)
      ENDDO !BANDS

      !WRITE (*, 1001)
! 1000 FORMAT('  Beginning cross section calculations...')
! 1001 FORMAT('  Cross section calculations completed.')

      RETURN

 100  FORMAT(F12.5, 1X, 8L2, 4D15.7, " CLEAR SVD & LM FLAG : ", I5 )

      END SUBROUTINE KROSSR

      END MODULE XSECTIONS
