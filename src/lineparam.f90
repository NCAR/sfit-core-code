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

      MODULE LINEPARAM

! --- HITRAN LINE PARAMETERS
! --- LNMAX : MAXIMUM NUMBER OF LINES FROM HITRAN LIST
! --- ASCII IMPLEMENTATION OF HBIN NOT COMPLETE AFTER GALATRY

      USE PARAMS
      USE BANDPARAM
      USE RETVPARAM
      USE DATAFILES
      USE MOLCPARAM
      USE ISOTOPE
      USE RAYTRACE

      IMPLICIT NONE

      LOGICAL,      DIMENSION(:,:), ALLOCATABLE :: HFLAG   ! FLAGS FOR: GALATRY, FCIA, SCIA, LMIX, SVD, 6, 7, 8
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: ST296   ! INTENSITY AT STANDARD TEMPERATURE
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: AAA     ! AIR BROADENING HALFWIDTH
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: SSS     ! SELF BROADENING HALFWIDTH
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: AZERO   ! WAVENUMBER
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: ETWO    ! LOWER STATE ENERGY
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: GMASS   ! MOLECULAR MASS
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: PSLIN   ! PRESSURE SHIFT
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: TDLIN   ! COEFFICIENT TEMP. DEP. AIR H-W
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: BETA    ! BETA0 FOR GALATRY ROUTINE
                                                           ! ANUVC FOR PCQSDHV MODEL
      INTEGER,      DIMENSION(:),   ALLOCATABLE :: LGAS    ! LINE INDEX NUMBER
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: GAMMA0  ! GAMMA0 FOR SDV (BOONE, PCQSDHC) ROUTINE
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: GAMMA2  ! GAMMA2 FOR SDV (BOONE, PCQSDHC) ROUTINE
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: SELF_GAMMA2  ! GAMMA2 (SELF) FOR SDV (BOONE, PCQSDHC) ROUTINE
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: SHIFT0  ! PRESSURE SHIFT FOR LSHAPE (TRAN) ROUTINE
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: SELF_SHIFT0  ! PRESSURE SHIFT SELF FOR LSHAPE (TRAN) ROUTINE
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: SHIFT2  ! P-SHIFT OF GAMMA 2 FOR LSHAPE (TRAN) ROUTINE
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: ETA     ! CORRELATION FOR PCQSDHC MODEL
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: LMTK1, LMTK2, YLM ! LINE MIXING PARAMETERS


      INTEGER :: NR_SDVGAS, NR_LMGAS

      REAL(DOUBLE) :: AUXLMTK1, AUXLMTK2, AUXYLM ! LINE MIXING PARAMETERS FOR PARAMETRISATION
                                                 ! OF Y GIVEN BY HASE

      INTEGER, DIMENSION(MAXBND)     :: LINE1   ! BEGINNING LINE FOR BAND I
      INTEGER, DIMENSION(MAXBND)     :: LINE2   ! ENDING LINE FOR BAND I

      INTEGER                        :: UFLAG, UGAS
      REAL(DOUBLE)                   :: UVAL
      LOGICAL :: USE_LM = .FALSE.               ! SWITCH ON LINEMIXING FOR LINES WITH LINEMIXING
                                                ! PARAMETERS AND EITHER VOIGT OR SDV LINESHAPE
      INTEGER :: LSHAPEMODEL                    ! USER CHOICE OF LINE SHAPE MODEL
                                                ! 0 = CHOOSE LINEPARAM DEPENDING ON EXISTANCE OF PARAMETERS
                                                ! 1 = FORCE VOIGT FOR ALL LINES
                                                ! 2 = USE GALATRY FOR LINES WITH PARAMETERS, VOIGT ELSE
                                                ! 3 = USE Boone SDV FOR LINES WITH PARAMETERS

      ! if lshapemodel = 0, the following switches may be used to switch certain features
      ! by default all switched on for default line shape
      logical :: lsm_sdv = .false.              ! if TRUE, speed dependent Voigt is used
      logical :: lsm_dicke = .false.            ! if TRUE, Dicke narrowing is calculated
      logical :: lsm_corr = .false.            ! if TRUE, crosscorrelation between SDV and Dicke 
                                                ! narrowing is calculated

      INTEGER, PARAMETER   :: GALATRY_FLAG=1,FCIA_FLAG=2,SCIA_FLAG=3,SDV_FLAG=4,&
           CORR_FLAG=5,LM_1ST_FLAG=6, LM_FULL_FLAG=7

      TYPE, PUBLIC :: HITRANDATA
         INTEGER  :: MO              ! MOL ID
         INTEGER  :: IS              ! ISOTOPE ID #
         REAL(8)  :: NU              ! WAVENUMBER
         REAL(8)  :: SL              ! INTENSITY [CM-1/(MOLEC/CM-2)]
         REAL(4)  :: EA              ! EINSTEIN A COEFF
         REAL(4)  :: AH              ! AIR BROADENED HALFWIDTH [CM-1/ATM]
         REAL(4)  :: SH              ! SELF BROADENED HALFWIDTH [CM-1/ATM]
         REAL(8)  :: EL              ! LOWER STATE ENERGY [CM-1]
         REAL(4)  :: TX              ! TEMPERATURE EXPONENT
         REAL(4)  :: PS              ! PRESSURE SHIFT [CM-1]
         CHARACTER (LEN=60) :: QA    ! QUANTA DATA
         CHARACTER (LEN=18) :: ER    ! ERROR AND REF FLAGS
         CHARACTER (LEN=1)  :: LM    ! LINE MIXING FLAGS
         REAL(4) :: UW               ! UPPER STAT WT
         REAL(4) :: LW               ! LOWER STAT WT
         REAL(4) :: BT               ! GALATRY BETA0
         REAL(4) :: GAMMA0           ! GAMMA 0 FOREIGN
         REAL(4) :: GAMMA2           ! GAMMA 2 NARROWING PARAMETER FOREIGN
         REAL(4) :: SELF_GAMMA2      ! GAMMA 2 NARROWING PARAMETER SELF
         REAL(4) :: SHIFT0           ! PRESSURE SHIFT FOREIGN FOR  GEN LINESHAPE
         REAL(4) :: SELF_SHIFT0      ! PRESSURE SHIFT SELF FOR GEN LINESHAPE
         REAL(4) :: SHIFT2           ! PRESSURE SHIFT FOR GAMMA 2 IN GEN LINESHAPE
         REAL(4) :: ETA              ! ETA 
         REAL(4) :: LMTK1            ! LMTK1 for Line Mixing (1st order expansion, Rosenkranz model)
         REAL(4) :: LMTK2            ! LMTK2 for Line Mixing (1st order expansion, Rosenkranz model)
         REAL(4) :: YLM              ! YLM for Line Mixing (1st order expansion, Rosenkranz model)
         LOGICAL :: FLAG(8)          ! GAL, FCIA, SCIA, SVD, CORR, LM_1ST, LM_FULL, 8
      END TYPE HITRANDATA

      CONTAINS

!----------------------------------------------------------------------
      SUBROUTINE OPTLIN

      INTEGER      :: LINE, IBAND, ISO, MO, L, K, I, NUMLIN, TNBAND
      INTEGER      :: NLINES, NLINES_GALATRY, NLINES_1ST_LM, NLINES_FULL_LM, NLINES_SDV, NLINES_FCIA, NLINES_SCIA, NLINES_ETA
      REAL(8)      :: WAVLO, WAVHI, PS, TW, ELOWER, SW, AW, SLINE, WAVNUM=0.0, DENLIN, B0
      REAL(8)      :: G0, G2, SG2, S0, SS0, S2, ETA0, L1, L2, L3, TW5, TW6
      CHARACTER (LEN=200)        :: CHLINE
      INTEGER, DIMENSION(MAXGAS) :: NMOLINE
      LOGICAL                    :: HBIN = .TRUE., HF(8)

      TYPE (HITRANDATA)  :: HLP

!  --- EXTRACT LINES FROM TAPE14 FOR SPECTRAL BANDPASSES
!  --- LINES ARE ASSUMED TO BE IN ASCENDING WAVENUMBER ORDER

      IFMIX(:MAXGAS)=0
      NGAS = 0
      NLINES = 0
      NLINES_1ST_LM = 0
      NLINES_FULL_LM = 0
      NLINES_SDV = 0
      NLINES_ETA = 0
      NLINES_FCIA = 0
      NLINES_SCIA = 0
      NLINES_GALATRY = 0

      WRITE (*, 120) TRIM(TFILE(14))
      WRITE (16, 121) TRIM(TFILE(14))

      write(16,*) 'LINESHAPEMODEL = ', LSHAPEMODEL
      write(16,*) 'DICKE NARROWING ', LSM_DICKE
      write(16,*) 'SPEED DEPENDANCY OF PRESSURE BROADENING ', LSM_SDV
      write(16,*) 'CORRELATION OF PRESSURE AND TEMPERATURE BROADENING ', LSM_CORR
      ! Check for consistency in line parameters and line features
      IF ((LSHAPEMODEL == 1).AND.LSM_DICKE) THEN
         write (*,*) 'Line shape model Voigtfunction does not support Dicke narrowing'
         write (16,*) 'Line shape model Voigtfunction does not support Dicke narrowing'
         call shutdown
         stop 2
      end IF
      IF (((LSHAPEMODEL == 1).OR.(LSHAPEMODEL == 2)).AND.LSM_SDV) THEN
         write (*,*) 'Line shape model Voigtfunction does not support speed dependancy'
         write (16,*) 'Line shape model Voigtfunction does not support speed dependancy'
         call shutdown
         stop 2
      end IF
      IF (((LSHAPEMODEL == 1).OR.(LSHAPEMODEL == 2)).AND.LSM_CORR) THEN
         write (*,*) 'Line shape model does not support coorelation between pressure and temperature broadening'
         write (16,*) 'Line shape model does not support coorelation between pressure and temperature broadening'
         call shutdown
         stop 2
      end IF
      IF (((LSHAPEMODEL == 1).OR.(LSHAPEMODEL == 2)).AND.USE_LM) THEN
         write (*,*) 'Line shape model does not support line mixing'
         write (16,*) 'Line shape model does not support line mixing'
         call shutdown
         stop 2
      end IF
!
!  --- OPEN HITRANLINE DATA FILE
      IF( HBIN )THEN
         OPEN( UNIT=14, FILE=TRIM(TFILE(14)), STATUS='OLD', ERR=669, FORM='UNFORMATTED' )
         READ(14) TNBAND
         IF( TNBAND .NE. NBAND )THEN
            WRITE(16,150) TNBAND, NBAND, '  BAND NUMBER IN HITRAN LISTING AND RETRIEVAL MISMATCH'
            WRITE( 6,150) TNBAND, NBAND, '  BAND NUMBER IN HITRAN LISTING AND RETRIEVAL MISMATCH'
            !STOP ': BAND NUMBER IN HITRAN LISTING AND RETRIEVAL MISMATCH'
         ENDIF
         DO I=1, TNBAND
            READ(14) K, TW5, TW6
            IF(( (TW5-WAVE5(K)) .GT. 0.0) .OR. ((TW6-WAVE6(K)) .LT. 0.0 )) THEN
               !print*, TW5, WAVE5(K), TW6, WAVE6(K)
               WRITE(16,151) K, TW5-WAVE5(K), TW6-WAVE6(K),'  HITRAN LIST WINDOW AND RETRIEVAL WINDOW MISMATCH'
               WRITE( 6,151) K, TW5-WAVE5(K), TW6-WAVE6(K),'  HITRAN LIST WINDOW AND RETRIEVAL WINDOW MISMATCH'
               !STOP ': LINEPARAM HITRAN LIST & RETRIEVAL WINDOW MISMATCH'
            ENDIF
         ENDDO
      ELSE
         OPEN( UNIT=14, FILE=TRIM(TFILE(14)), STATUS='OLD', ERR=669 )
      ENDIF

      ALLOCATE( ST296(LNMAX), AAA(LNMAX), SSS(LNMAX), AZERO(LNMAX), ETWO(LNMAX), &
           GMASS(LNMAX), PSLIN(LNMAX), TDLIN(LNMAX), BETA(LNMAX), LGAS(LNMAX), &
           GAMMA0(LNMAX), GAMMA2(LNMAX), SELF_GAMMA2(LNMAX), SHIFT0(LNMAX), SELF_SHIFT0(LNMAX), SHIFT2(LNMAX), ETA(LNMAX), &
           LMTK1(LNMAX), LMTK2(LNMAX), YLM(LNMAX), HFLAG(LNMAX,8) )

      LINE = 0
 10   CONTINUE

!  --- READ HITRAN FORMAT LINE LIST (MODIFIED MOLECULE NUMBERS)
      IF( HBIN )THEN
         !READ(14, END=230) MO, ISO, WAVNUM, SLINE, EA, AH, SH, ELOWER, TX, P4, QA, ER, LM, UW, LW, BT
         READ(14, END=230) HLP
         MO     = HLP%MO
         ISO    = HLP%IS
         WAVNUM = HLP%NU
         SLINE  = HLP%SL
         ELOWER = HLP%EL
         AW     = REAL( HLP%AH, 8 )
         SW     = REAL( HLP%SH, 8 )
         TW     = REAL( HLP%TX, 8 )
         PS     = REAL( HLP%PS, 8 )
         B0     = REAL( HLP%BT, 8 )
         G0     = REAL( HLP%GAMMA0, 8 )
         G2     = REAL( HLP%GAMMA2, 8 )
         SG2     = REAL( HLP%SELF_GAMMA2, 8 )
         S0     = REAL( HLP%SHIFT0, 8 )
         SS0     = REAL( HLP%SELF_SHIFT0, 8 )
         S2     = REAL( HLP%SHIFT2, 8 )
         ETA0     = REAL( HLP%ETA, 8 )
         L1     = REAL( HLP%LMTK1, 8 )
         L2     = REAL( HLP%LMTK2, 8 )
         L3     = REAL( HLP%YLM, 8 )
         HF     = HLP%FLAG
         !WRITE(12,*), MO, ISO, WAVNUM, SLINE, AW, SW, ELOWER, TW, PS, B0
      ELSE
         READ (14, '(A)', END=230 ) CHLINE
         READ (CHLINE, 107, END=230) MO, ISO, WAVNUM, SLINE, AW, SW, ELOWER, TW, PS, B0
         !write(13,*), MO, ISO, WAVNUM, SLINE, AW, SW, ELOWER, TW, PS, B0
      ENDIF
      !print *, MO, ISO, WAVNUM, SLINE, EA, AH, SH, ELOWER, TX, P4, QA, ER, LM, UW, LW, BT

      LINE = LINE +1
      IF( LINE .EQ. 1 )WAVLO = WAVNUM

!  --- CHECK FOR MOL ID OUT OF BOUNDS
      IF (MO > MOLTOTAL) GO TO 665

!  --- IDENTIFY SPECIES AND STORE AFGL CODE IN ARRAY ICODE
      IF (NGAS .LE. 0) THEN
         NGAS          = 1
         ICODE(NGAS)   = MO        ! STORE HITRAN ID
         ISCODE(NGAS)  = ISO       ! STORE HITRAN ID2
         NMOLINE(NGAS) = 1         ! STORE # OF LINES FOR THIS MOLECULE
         GO TO 2268
      ENDIF

!  --- CHECK IF WE ALREADY HAVE A LINE OF THIS GAS (AND ITS PROFILE)
      DO L = 1, NGAS
         IF (ICODE(L) == MO) THEN
            NMOLINE(L) = NMOLINE(L) + 1
            GOTO 2269
         ENDIF
      END DO

!  --- NEW GAS - ADD TO LIST SORTED BY GAS
      NGAS          = NGAS + 1         ! TOTAL # OF GASES
      IF( NGAS > MAXGAS )GOTO 670
      ICODE(NGAS)   = MO               ! HITRAN ID OF THIS LINE
      ISCODE(NGAS)  = ISO              ! HITRAN ISOTOPE ID OF THIS GAS
      NMOLINE(NGAS) = 1                ! # OF LINES FOR THIS GAS

2268  CONTINUE

!  --- LOOP THROUGH MIXING RATIO LIST FOR PROFILE FOR THIS LINE
      DO I=1, NMOL
!print*, 'linep ',TRIM(HMOLS(I)) , '  ',  TRIM(NAME(ICODE(NGAS)))
         IF( TRIM(HMOLS(I)) .EQ. TRIM(NAME(ICODE(NGAS))))THEN
            XGAS(NGAS,:NPATH) = FXGAS(icode(ngas),:NPATH)
            IFMIX(NGAS) = 1
            IF( SUM(FXGAS(icode(ngas),:NPATH)) .LE. 0.0D0 ) IFMIX(NGAS) = 0
!print*, 'linep ',I, TRIM(NAME(ICODE(NGAS))), ngas, npath, ICODE(NGAS), ifmix(ngas), FXGAS(icode(ngas),:NPATH), XGAS(NGAS,:NPATH)
!  --- ADD LINE TO LIST
            NLINES = NLINES + 1
            LGAS(NLINES) = NGAS
            GO TO 2270
          ENDIF
      ENDDO

!  --- MIXING RATIO PROFILE NOT FOUND FOR NEW GAS--SET VMRS AND IFMIX TO ZERO
!  --- DON'T ADD LINE TO LINELIST
      WRITE (16, 88) NAME(ICODE(NGAS))
      XGAS(NGAS,:KMAX) = 0.D0
      IFMIX(NGAS) = 0
      GOTO 10

!  --- "OLD" GAS - TEST IF VMRS ARE ZERO - IF SO, DON'T ADD LINE
2269  CONTINUE
      IF (IFMIX(L) == 0) GOTO 10

!  --- NONZERO VMRS FOR GAS--ADD LINE TO LIST
      NLINES = NLINES + 1
      LGAS(NLINES) = L

2270  CONTINUE
      IF( NLINES .GT. LNMAX )GOTO 7
      AZERO(NLINES) = WAVNUM
      IF( UFLAG .EQ. 1 ) THEN
         IF( MO .EQ. UGAS ) THEN
            PRINT *, ' CHANGING LINE STRENGTH : ', NAME(UGAS), '  ',UGAS, '  ', UVAL, ' !!!!!!!!!!!!!!!!!!!!!!!!'
            SLINE = SLINE + UVAL*SLINE
         ENDIF
      ENDIF
      IF( UFLAG .EQ. 2 ) THEN
         IF( MO .EQ. UGAS ) THEN
            PRINT *, ' CHANGING Air-broadened .5-width : ', NAME(UGAS), '  ',UGAS, '  ', UVAL, ' !!!!!!!!!!!!!!!!!!!!!!!!'
            AW = AW + UVAL*AW
         ENDIF
      ENDIF
      ST296(NLINES)  = SLINE
      AAA(NLINES)    = AW
      SSS(NLINES)    = SW
      ETWO(NLINES)   = ELOWER
      TDLIN(NLINES)  = TW
      PSLIN(NLINES)  = PS
      BETA(NLINES)   = B0
      GAMMA0(NLINES) = G0
      GAMMA2(NLINES) = G2
      SELF_GAMMA2(NLINES) = SG2

      ! Default pressure shift is PSLIN
      if (S0.lt.tiny(S0)) then
         SHIFT0(NLINES)   = PS
      else
         SHIFT0(NLINES)   = S0
      end if
      SELF_SHIFT0(NLINES)   = SS0
      SHIFT2(NLINES)   = S2
      ETA(NLINES)   = ETA0
      LMTK1(NLINES)  = L1
      LMTK2(NLINES)  = L2
      YLM(NLINES)    = L3

      HFLAG(NLINES,1:8) = .FALSE.
      IF( LSM_DICKE.AND.HF(GALATRY_FLAG) .AND. ((LSHAPEMODEL == 0).OR.(LSHAPEMODEL==2))) THEN
         HFLAG(NLINES,GALATRY_FLAG) = .TRUE.
         NLINES_GALATRY = NLINES_GALATRY + 1
         !print *, HFLAG(NLINES,GALATRY_FLAG), nlines
      END IF
      IF( HF(FCIA_FLAG)) THEN
         HFLAG(NLINES,FCIA_FLAG) = .TRUE.
         NLINES_FCIA = NLINES_FCIA + 1
      END IF
      IF( HF(SCIA_FLAG)) THEN
         HFLAG(NLINES,SCIA_FLAG) = .TRUE.
         NLINES_SCIA = NLINES_SCIA + 1
      END IF
      IF( LSM_SDV.and.HF(SDV_FLAG).AND.(LSHAPEMODEL == 0)) THEN
         HFLAG(NLINES,SDV_FLAG) = .TRUE.
         NLINES_SDV = NLINES_SDV + 1
      END IF
      IF( LSM_CORR.and.HF(CORR_FLAG).AND.(LSHAPEMODEL == 0)) THEN
         HFLAG(NLINES,CORR_FLAG) = .TRUE.
         NLINES_ETA = NLINES_ETA + 1
      END IF
      IF( HF(LM_1ST_FLAG).AND.USE_LM.AND.((LSHAPEMODEL == 0))) THEN
         HFLAG(NLINES,LM_1ST_FLAG) = .TRUE.
         NLINES_1ST_LM = NLINES_1ST_LM + 1
      END IF
      IF( HF(LM_FULL_FLAG).AND.USE_LM.AND.((LSHAPEMODEL == 0))) THEN
         HFLAG(NLINES,LM_FULL_FLAG) = .TRUE.
         NLINES_FULL_LM = NLINES_FULL_LM + 1
      END IF

!  --- CHECK TO BE SURE ISOTOPIC MOLECULAR WEIGHT IS DEFINED
      IF( ISO .LT. 1 .OR. ISO .GT. NHIISO(ICODE(LGAS(NLINES))) )GOTO 668
      GMASS(NLINES) = XMASS(ISO,MO)

!  --- CONVERSION OF LINE STRENGTH TO CM-2ATM-1 UNITS
      ST296(NLINES) = ST296(NLINES)*SCHMIT*ZEROC/STDTEMP

! --- FOR O2 CONTINUUM (49/1 FCIA, 49/2 SCIA)
      IF( HFLAG(NLINES,FCIA_FLAG) .OR. HFLAG(NLINES,SCIA_FLAG) )THEN
!         AAA(NLINES) = 0.1D0
!         SSS(NLINES) = 0.1D0
          if (AW .le. 0.0) AAA(NLINES) = 0.1D0

          if (SW .le. 0.0) SSS(NLINES) = AAA(NLINES)
         ST296(NLINES) = ST296(NLINES)*SCHMIT*ZEROC/STDTEMP
      ENDIF

      GOTO 10

230   CONTINUE

      CLOSE(14)

      WAVHI = WAVNUM
      NUMLIN = NLINES
      WRITE (16, 886) NUMLIN
!  --- TEST TO BE SURE LINE LIST ISN'T EMPTY
      IF (NUMLIN <= 0) GO TO 667

      DO I = 1, NGAS
         IF (IFMIX(I) == 0) THEN
            WRITE (16, 700) NAME(ICODE(I))
         ELSE
            WRITE (16, 710) NAME(ICODE(I))
         ENDIF
      END DO

      WRITE (16, 891) NGAS
!  --- OUTPUT NUMBER OF GASES USED
      WRITE (16, 892) SUM( IFMIX )

!  --- PRINT ONLY THOSE PROFILES USED TO DETAIL
      DO I = 1, NGAS
         IF (IFMIX(I) == 0) CYCLE
         WRITE (16, 882) NAME(ICODE(I)), ICODE(I), NHIISO(ICODE(I)), NMOLINE(I)
         WRITE (16,883) ( XMASS( L, ICODE(I) ), L=1, NHIISO(ICODE(I)) )
         WRITE (16, 884) (XGAS(I,K),K=1,KMAX)
      END DO

      WRITE (16, 300) WAVLO, WAVHI

      DENLIN = NLINES/(WAVHI - WAVLO)
      WRITE (16, 887) DENLIN
      IF (DENLIN>500 .AND. DENLIN<1000.) THEN
         TAUMIN = 1.E-07
         WRITE (16, 888) TAUMIN
      ELSE IF (DENLIN >= 1000) THEN
         TAUMIN = 1.E-08
         WRITE (16, 889) TAUMIN
      ELSE
         WRITE (16, 890) TAUMIN
      ENDIF

!  --- DETERMINE RANGE OF LINES TO CONSIDER FOR EACH BANDPASS
      WRITE (16, 301)
      DO IBAND = 1, NBAND
         IF (NSCAN(IBAND) == 0) CYCLE
         NLINES = 0
         DO I = 1, NUMLIN
            IF (AZERO(I) < WAVE5(IBAND)) CYCLE
            IF (AZERO(I) >= WAVE6(IBAND)) CYCLE
            NLINES = NLINES + 1
            IF (NLINES == 1) LINE1(IBAND) = I
            LINE2(IBAND) = I
         END DO
         I = IBAND
         IF (NLINES == 0) GO TO 667
         WRITE (16, 302) IBAND, WAVE5(IBAND), LINE1(IBAND), WAVE6(IBAND), LINE2(IBAND)
      END DO

! ---REUSE HF -- GALATRY_FLAG=1,FCIA_FLAG=2,SCIA_FLAG=3,SDV_FLAG=4,CORR_FLAG=5,LM_1ST_FLAG=6,LM_FULL_FLAG=7
      HF(1:8) = .FALSE.
      DO I=1, NUMLIN
         IF( HFLAG(I,GALATRY_FLAG) ) HF(GALATRY_FLAG) = .TRUE.
         IF( HFLAG(I,FCIA_FLAG) )    HF(FCIA_FLAG) = .TRUE.
         IF( HFLAG(I,SCIA_FLAG) )    HF(SCIA_FLAG) = .TRUE.
         IF( HFLAG(I,SDV_FLAG) )     HF(SDV_FLAG) = .TRUE.
         IF( HFLAG(I,CORR_FLAG) )     HF(CORR_FLAG) = .TRUE.
         IF( HFLAG(I,LM_1ST_FLAG) )      HF(LM_1ST_FLAG) = .TRUE.
         IF( HFLAG(I,LM_FULL_FLAG) )      HF(LM_FULL_FLAG) = .TRUE.
      ENDDO

      !FORALL( I=1:NLINES, HFLAG(I,GALATRY_FLAG) ) HF(GALATRY_FLAG) = .TRUE.
      !FORALL( I=1:NLINES, HFLAG(I,FCIA_FLAG) ) HF(FCIA_FLAG) = .TRUE.
      !FORALL( I=1:NLINES, HFLAG(I,SCIA_FLAG) ) HF(SCIA_FLAG) = .TRUE.
      !FORALL( I=1:NLINES, HFLAG(I,SDV_FLAG) ) HF(SDV_FLAG) = .TRUE.
      !FORALL( I=1:NLINES, HFLAG(I,LM_FLAG) ) HF(LM_FLAG) = .TRUE.


      WRITE(*,131) HF(FCIA_FLAG) ,NLINES_FCIA
      WRITE(*,132) HF(SCIA_FLAG) ,NLINES_SCIA
      WRITE(*,130) HF(GALATRY_FLAG) ,NLINES_GALATRY
      WRITE(*,133) HF(SDV_FLAG) ,NLINES_SDV
      WRITE(*,134) HF(CORR_FLAG) ,NLINES_ETA
      WRITE(*,135) HF(LM_1ST_FLAG) ,NLINES_1ST_LM
      WRITE(*,136) HF(LM_FULL_FLAG) ,NLINES_FULL_LM

      WRITE(16,141) HF(FCIA_FLAG) ,NLINES_FCIA
      WRITE(16,142) HF(SCIA_FLAG) ,NLINES_SCIA
      WRITE(16,140) HF(GALATRY_FLAG) ,NLINES_GALATRY
      WRITE(16,143) HF(SDV_FLAG) ,NLINES_SDV
      WRITE(16,144) HF(CORR_FLAG) ,NLINES_ETA
      WRITE(16,145) HF(LM_1ST_FLAG) ,NLINES_1ST_LM
      WRITE(16,146) HF(LM_FULL_FLAG) ,NLINES_FULL_LM

      RETURN

    7 CONTINUE
      WRITE (16, 112) LNMAX
      WRITE ( 0, 112) LNMAX
      CALL SHUTDOWN
      STOP '2'
  665 CONTINUE
      WRITE (16, 110) MO
      WRITE ( 0, 110) MO
      CALL SHUTDOWN
      STOP '2'
  667 CONTINUE
      WRITE (16, 101) I
      WRITE ( 0, 101) I
      CALL SHUTDOWN
      STOP '2'
  668 CONTINUE
      WRITE (16, 102) LINE, ICODE(LGAS(NLINES)), ISO, NHIISO(ICODE(LGAS(NLINES)))
      WRITE (16, 103) AZERO(LINE)
      WRITE ( 0, 102) LINE, ICODE(LGAS(NLINES)), ISO, NHIISO(ICODE(LGAS(NLINES)))
      WRITE ( 0, 103) AZERO(LINE)
      CALL SHUTDOWN
      STOP 2
  669 CONTINUE
      WRITE (16, 113) TRIM(TFILE(14))
      WRITE ( 0, 113) TRIM(TFILE(14))
      CALL SHUTDOWN
      STOP '2'
  670 CONTINUE
      WRITE (16, 114) MAXGAS
      WRITE ( 0, 114) MAXGAS
      CALL SHUTDOWN
      STOP '2'


   88 FORMAT(A7,': OPTLIN: MIXING RATIO NOT FOUND TAPE12, SET TO ZERO FOR ALL LAYERS')
  101 FORMAT(' ABORT *** NO LINES FOR BANDPASS #',I3)
  102 FORMAT(/,' ABORT *** ISOTOPE CODE OUT OF RANGE',/,&
         ' ERROR MESSAGE FROM SUBROUTINE OPTLIN',/,' LINE # =',I5,&
         ' MOLECULE # =',I3,' ISOTOPE # =',I2,/,' MAXIMUM ISOTOPE # =',I2)
  103 FORMAT(/,' LINE POSITION =',F10.4)
  107 FORMAT(I2,I1,F12.6,1P,E10.3,10X,0P,F5.4,F5.4,F10.4,F4.2,F8.6,30X,30X,6X,12X,1X,14X,f10.5)
  110 FORMAT(' ABORT !!! OPTLIN - MOLECULE ID OUT OF BOUNDS: ',I4)
  112 FORMAT(/,'LINEPARAM: ABORT !!! NUMBER OF LINES EXCEEDS ',I8)
  113 FORMAT(/,' INPUT FILE OPEN ERROR UNIT 14-LINE FILE = "',A,'"')
  114 FORMAT(/,' ABORT !!! NUMBER OF MOLECULES EXCEEDS',I5)

  120 FORMAT( '   HITRAN FILE : ', A )
  121 FORMAT(/' HITRAN LINELIST FILE : ', A )
  130 FORMAT( '     GALATRY FLAG & LINES WITH GALATRY PARAMETERS FOUND       :',L3, I7)
  131 FORMAT( '     FCIA FLAG & FCIA LINES FOUND                             :',L3, I7)
  132 FORMAT( '     SCIA FLAG & SCIA LINES FOUND                             :',L3, I7)
  133 FORMAT( '     SDV FLAG & LINES WITH SDV PARAMETERS FOUND               :',L3, I7)
  134 FORMAT( '     CORR FLAG & LINES WITH CORR PARAMETERS FOUND             :',L3, I7)
  135 FORMAT( '     LMIX 1ST FLAG & LINES WITH LMIX 1ST PARAMETERS FOUND     :',L3, I7)
  136 FORMAT( '     LMIX FULL FLAG & LINES WITH LMIX FULL PARAMETERS FOUND   :',L3, I7)
  140 FORMAT(/' GALATRY FLAG & LINES WITH GALATRY PARAMETERS FOUND       :',L3, I7)
  141 FORMAT( ' FCIA FLAG & FCIA LINES FOUND                             :',L3, I7)
  142 FORMAT( ' SCIA FLAG & SCIA LINES FOUND                             :',L3, I7)
  143 FORMAT( ' SDV FLAG & LINES WITH SDV PARAMETERS FOUND               :',L3, I7)
  144 FORMAT( ' CORR FLAG & LINES WITH CORR PARAMETERS FOUND             :',L3, I7)
  145 FORMAT( ' LMIX 1ST FLAG & LINES WITH LMIX 1ST PARAMETERS FOUND     :',L3, I7)
  146 FORMAT( ' FULL LMIX FLAG & LINES WITH FULL LMIX PARAMETERS FOUND   :',L3, I7)
  150 FORMAT( 2I6, A)
  151 FORMAT( I6, 2F14.6, A)

  300 FORMAT(/,' FILE LIMITS FOR LINES=',2F12.6)
  301 FORMAT(/,' WAVENUMBER LIMITS AND INDICES FOR INCLUDED LINES',/,&
         ' BANDPASS   LOWER NU LIMIT  INDEX1     UPPER NU LIMIT  INDEX2')
  302 FORMAT(I6,F18.5,I9,F17.5,I11)
  700 FORMAT('           ZERO GAS : ',A,' NO EFFECT IN CALCULATION')
  710 FORMAT('          USING GAS : ',A)
  882 FORMAT(/,' MIXING RATIO PROFILE FOR: ',A7,' ID: ',I2,' NISO:',I4,' LINES:', I5)
  883 FORMAT( 10(F6.1,2X) )
  884 FORMAT(6(1P,E12.4))
  886 FORMAT(/,' TOTAL NUMBER OF LINES IN ANALYSIS LIST       : ',I6)
  887 FORMAT(/,' LINE DENSITY (LINES/WAVENUMBER)              : ',F7.1)
  888 FORMAT(' GT.500 LINES-WAVENUMBER-TAUMIN DECREASED TO  : ',1P,D9.2)
  889 FORMAT(' GT.1000 LINES-WAVENUMBER-TAUMIN DECREASED TO : ',1P,D9.2)
  890 FORMAT(' TAUMIN DEFINED AS                            : ',1P,D9.2)
  891 FORMAT(/,' NUMBER OF GASES FOUND IN MICROWINDOWS        : ',I2)
  892 FORMAT(' NUMBER OF GASES INCLUDED IN CALCULATION      : ',I2)

      RETURN

      END SUBROUTINE OPTLIN

      SUBROUTINE RELEASE_MEM_LP

      IF( ALLOCATED( ST296 ))THEN
         DEALLOCATE( ST296, AAA, SSS, AZERO, ETWO, GMASS, PSLIN, TDLIN, BETA, LGAS, &
              GAMMA0, GAMMA2, SELF_GAMMA2, SHIFT0, SHIFT2, LMTK1, LMTK2, YLM, HFLAG)

      ENDIF

      END SUBROUTINE RELEASE_MEM_LP

      END MODULE LINEPARAM
