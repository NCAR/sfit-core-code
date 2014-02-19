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

      MODULE INITIALIZE

      USE PARAMS
      USE RETVPARAM
      USE TRANSMIS
      USE MOLCPARAM
      USE XSECTIONS
      USE DATAFILES
      USE SYNSPEC
      USE LINEPARAM
      USE BANDPARAM
      USE RAYTRACE
      USE SOLAR
      USE WRITEOUT
      USE CHANNEL


      IMPLICIT NONE

      LOGICAL :: RETFLG

      INTEGER :: NSTNR = 0
      INTEGER :: ISOFLAG, IPFLAG
      INTEGER :: ISMIX, NVAR = 0, NFIT = 0

      REAL(DOUBLE) :: SPHS, PHS, SWSHFT, WSHFT, TOL, &
                      SBCKOFF, BCKOFF, SBCKCRV, BCKCRV, SBCKSL, BCKSL

      REAL(DOUBLE), DIMENSION(MMAX)    :: WWV
      REAL(DOUBLE), DIMENSION(MAXSNR)  :: WWV0, WWV1, GSTNR

      CHARACTER (LEN=14), DIMENSION(NMAX) :: PNAME
      CHARACTER (LEN=14), DIMENSION(NMAX) :: ORIG_PNAME
      REAL(DOUBLE), DIMENSION(NMAX)       :: PARM  = 0.0D0
      REAL(DOUBLE), DIMENSION(NMAX)       :: SPARM = 0.0D0
      CHARACTER (LEN=14), DIMENSION(5)    :: CPNAM
      REAL(DOUBLE), DIMENSION (5)         :: SCOPAR, COPAR

      CHARACTER(LEN=1024)                  :: S_KB_PROFILE_GASES
      CHARACTER(LEN=14), DIMENSION(MOLMAX) :: S_KB_PRF_GAS

      LOGICAL :: F_KB           = .FALSE.
      LOGICAL :: F_KB_PROFILE   = .FALSE.
      LOGICAL :: F_KB_TEMP      = .FALSE.
      LOGICAL :: F_KB_SLOPE     = .FALSE.
      LOGICAL :: F_KB_CURVATURE = .FALSE.
      LOGICAL :: F_KB_SOLSHFT   = .FALSE.
      LOGICAL :: F_KB_SOLSTRNTH = .FALSE.
      LOGICAL :: F_KB_PHASE     = .FALSE.
      LOGICAL :: F_KB_IFDIFF    = .FALSE.
      LOGICAL :: F_KB_WSHIFT    = .FALSE.
      LOGICAL :: F_KB_EAP       = .FALSE.
      LOGICAL :: F_KB_EPHS      = .FALSE.
      LOGICAL :: F_KB_ZSHIFT    = .FALSE.
      LOGICAL :: F_KB_LINE      = .FALSE.
      INTEGER :: I_KB_LINE_TYPE = 0
      LOGICAL :: F_KB_SZA       = .FALSE.
      LOGICAL :: F_KB_FOV       = .FALSE.
      LOGICAL :: F_KB_OPD       = .FALSE.

      CHARACTER(LEN=1024)                  :: S_KB_LINE_GASES
      CHARACTER(LEN=14), DIMENSION(MOLMAX) :: S_KB_LINE_GAS

      LOGICAL, DIMENSION(NMAX) :: IS_IN_KB = .true.


      CONTAINS

!----------------------------------------------------------------------
      SUBROUTINE SETUP1

!  --- SETUP CALCULATIONS FOR THE VARIOUS BANDPASSES

      INTEGER        :: IBAND, NEXTRA, NAERR, MAXMPT
      REAL(DOUBLE)   :: DWAVE, WAVE1, WAVE2, AN


!  --- READ DATA FOR EACH BANDPASS.  STORE IN ARRAY TOBS
!  --- ALSO TABULATE LIST OF SCANS BEING FIT IN EACH REGION
!        NSCAN(IBAND)    = NUMBER OF SCANS TO BE FIT FOR BANDPASS IBAND
!        ISCAN(IBAND,I)) = SPECTRUM INDEX FOR THE ITH SCAN FITTED

      CALL GETSPEC()

      WRITE(16,105) NATMOS

!  --- COMPUTE INTERVAL FOR MONOCHROMATIC CALCULATIONS FOR EACH BANDPASS
!         ADD ABOUT 10 TIMES THE RESOLUTION TO BOTH SIDES OF THE INTEGRATION
!         INTERVAL TO ALLOW FOR WAVELENGTH SHIFTS
      MAXMPT = 0
      WRITE (16, 200)
      DO IBAND = 1, NBAND

         IF (NSCAN(IBAND) == 0) CYCLE

         DWAVE  = 10.D0/PMAX(IBAND)
         NEXTRA = NINT( DWAVE/DN(IBAND))
         WAVE1  = WSTART(IBAND) - NEXTRA*DN(IBAND)
         WAVE2  = WSTOP(IBAND)  + NEXTRA*DN(IBAND)
         NSTART(IBAND) = NEXTRA + 1
         WMON(IBAND)   = WAVE1
         NM(IBAND)     = FLOOR((WAVE2 - WAVE1)/DN(IBAND) + 1.000000001D0)
!  --- INTERVAL FOR INPUT OF LINE DATA
         WAVE5(IBAND) = WAVE1 - DLINES
         WAVE6(IBAND) = WAVE2 + DLINES

!  --- ESTIMATE 2**M FOR FFT
         AN = NM(IBAND)
         MFFT(IBAND) = FLOOR(LOG(AN)/LOG(2.0D0)) + 1
         MPT(IBAND)  = 2**MFFT(IBAND)
         IF (MPT(IBAND) > MAXMPT) MAXMPT = MPT(IBAND)
         LOWFIL(IBAND) = (MPT(IBAND)-NM(IBAND))/2
         HIFILL(IBAND) = LOWFIL(IBAND) + NM(IBAND)
         NSTZ1(IBAND)  = FLOOR(DN(IBAND)*PMAX(IBAND)*MPT(IBAND)) + 1
         NSTZ2(IBAND)  = MPT(IBAND) - NSTZ1(IBAND)

         WRITE (16, 201) WSTART(IBAND), WSTOP(IBAND), NM(IBAND), MPT(IBAND), NSTZ1(IBAND), NSTZ2(IBAND)
      END DO

!  --- COUNT TOTAL NUMBER OF MONOCHROMATIC POINTS (NMONSM)
!  --- AND THE TOTAL NUMBER OF CROSS SECTION POINTS (NCROSS)
!  --- TO BE CALCULATED AND CHECK FOR ARRAY OVERFLOWS
      NMONSM = 0
      NCROSS = 0
      NCROSS = SUM(NM(:NBAND))
      NMONSM = DOT_PRODUCT(NM(:NBAND),NSCAN(:NBAND))

      ALLOCATE (CROSS(NRET+1,KMAX,NCROSS), STAT=NAERR)
      IF (NAERR /= 0) THEN
         WRITE (16, *) 'INITIALIZE: COULD NOT ALLOCATE CROSS ARRAY ERROR NUMBER = ', NAERR
         WRITE ( 0, *) 'INITIALIZE: COULD NOT ALLOCATE CROSS ARRAY ERROR NUMBER = ', NAERR
         CALL SHUTDOWN
         STOP 2
      ENDIF
      ALLOCATE (CROSS_FACMAS(NRET+1,KMAX,NMONSM), STAT=NAERR)
      IF (NAERR /= 0) THEN
         WRITE (16, *) 'INITIALIZE: COULD NOT ALLOCATE CROSS_FACMAS ARRAY ERROR NUMBER = ', NAERR
         WRITE ( 0, *) 'INITIALIZE: COULD NOT ALLOCATE CROSS_FACMAS ARRAY ERROR NUMBER = ', NAERR
         CALL SHUTDOWN
         STOP 2
      ENDIF
      ALLOCATE (TCO(NCROSS), STAT=NAERR)
      IF (NAERR /= 0) THEN
         WRITE (16, *) 'INITIALIZE: COULD NOT ALLOCATE TCO ARRAY ERROR NUMBER = ', NAERR
         WRITE ( 0, *) 'INITIALIZE: COULD NOT ALLOCATE TCO ARRAY ERROR NUMBER = ', NAERR
         CALL SHUTDOWN
         STOP 2
      ENDIF
      ALLOCATE (TCONV(NMONSM),STAT=NAERR)
      IF (NAERR /= 0) THEN
         WRITE (16, *) 'INITIALIZE: COULD NOT ALLOCATE TCONV ARRAY ERROR NUMBER = ', NAERR
         WRITE ( 0, *) 'INITIALIZE: COULD NOT ALLOCATE TCONV ARRAY ERROR NUMBER = ', NAERR
         CALL SHUTDOWN
         STOP 2
      ENDIF
      ALLOCATE (TCALC(2,NMONSM), STAT=NAERR)
      IF (NAERR /= 0) THEN
         WRITE (16, *) 'INITIALIZE: COULD NOT ALLOCATE TCALC ARRAY ERROR NUMBER = ', NAERR
         WRITE ( 0, *) 'INITIALIZE: COULD NOT ALLOCATE TCALC ARRAY ERROR NUMBER = ', NAERR
         CALL SHUTDOWN
         STOP 2
      ENDIF
      ALLOCATE (TCALC_I(2,NMONSM), STAT=NAERR)
      IF (NAERR /= 0) THEN
         WRITE (16, *) 'INITIALIZE: COULD NOT ALLOCATE TCALC_I ARRAY ERROR NUMBER = ', NAERR
         WRITE ( 0, *) 'INITIALIZE: COULD NOT ALLOCATE TCALC_I ARRAY ERROR NUMBER = ', NAERR
         CALL SHUTDOWN
         STOP 2
      ENDIF
      ALLOCATE (IMGG(MAXMPT), STAT=NAERR)
      IF (NAERR /= 0) THEN
         WRITE (16, *) 'INITIALIZE: COULD NOT ALLOCATE IMGG ARRAY ERROR NUMBER = ', NAERR
         WRITE ( 0, *) 'INITIALIZE: COULD NOT ALLOCATE IMGG ARRAY ERROR NUMBER = ', NAERR
         CALL SHUTDOWN
         STOP 2
      ENDIF
      ALLOCATE (TCALC_E(2,NMONSM, KMAX+1), STAT=NAERR)
      IF (NAERR /= 0) THEN
         WRITE (16, *) 'INITIALIZE: COULD NOT ALLOCATE TCALC_E ARRAY ERROR NUMBER = ', NAERR
         WRITE ( 0, *) 'INITIALIZE: COULD NOT ALLOCATE TCALC_E ARRAY ERROR NUMBER = ', NAERR
         CALL SHUTDOWN
         STOP 2
      ENDIF
      ALLOCATE (TCALC_S(2,NMONSM, KMAX), STAT=NAERR)
      IF (NAERR /= 0) THEN
         WRITE (16, *) 'INITIALIZE: COULD NOT ALLOCATE TCALC_S ARRAY ERROR NUMBER = ', NAERR
         WRITE ( 0, *) 'INITIALIZE: COULD NOT ALLOCATE TCALC_S ARRAY ERROR NUMBER = ', NAERR
         CALL SHUTDOWN
         STOP 2
      ENDIF

      WRITE (16, 212) NMONSM, NCROSS

      RETURN

  105 FORMAT(/,' TOTAL NUMBER OF SPECTRAL DATA POINTS TO FIT =',I6)
  200 FORMAT(/,'    WSTART     WSTOP    NMON    MPT    NSTZ1    NSTZ2')
  201 FORMAT(2F10.4,I8,I7,2I9)
  212 FORMAT(/,' NMONSM   =',I8,/,&
               ' NCROSS   =',I8)

      END SUBROUTINE SETUP1


! ------------------------------------------------------------------------------

      SUBROUTINE SETUP2

      CHARACTER (LEN=64)  :: ILSHEAD
      INTEGER             :: I, J, K
      REAL(DOUBLE)        :: DUMMY

! ---- STORE CO LINE LIST
      IF( IFCO )THEN
         WRITE (*, *) ' READING SOLAR LINE LIST FILE...'
         CALL SOLARFH( 0 )

! ---- IF NO SOLAR CO LINES FOUND, RESET IFCO=0
         IF (NCOLNS == 0) IFCO = .FALSE.

! --- DEFAULT VALUES TO ZERO IN CASE PRINT BELOW
         EAPX(:MAXEAP) = 0.0
         EAPF(:MAXEAP) = 0.0
         EPHSX(:MAXEAP) = 0.0
         EPHSF(:MAXEAP) = 0.0
      ENDIF

! --- IF USING EMPIRICAL APODIZATION, READ PARAMETERS
      IF( F_EAPOD )THEN
         WRITE (*, *) ' READING EMPIRICAL MODULATION PARAMETER FILE...'
         CALL FILEOPEN( 23, 3 )
         IF (IEAP == 4) THEN
            ! READS LINEFIT FORMAT ASSUMES 20 VALUES)
            JEAP = 20
            READ (23, '(A64)') ILSHEAD
            ! DON'T READ EPHS HERE SINCE THERE COULD BE TWO DIFFERENT FILES
            READ (23, *, ERR=120) (EAPX(I),EAPF(I),DUMMY,I=1,JEAP)
         ELSE
            ! ORIGINAL FILE FORMAT
            READ (23, *) JEAP
            READ (23, *) (EAPF(I),I=1,JEAP)
            WRITE (16, '(/A)') 'EMPIRICAL MODULATION FUNCTION COEFFICIENTS'
            WRITE (16, *) (EAPF(I),I=1,JEAP)
            IF (IEAP == 1) THEN
               READ (23, *) (EAPX(I),I=1,JEAP)
               WRITE (16, *) (EAPX(I),I=1,JEAP)
            ENDIF
         ENDIF
         CALL FILECLOSE(23, 2)
      ENDIF

! --- IF USING EMPIRICAL PHASE FUNCTION, READ PARAMETERS
      IF( F_EPHASE )THEN
         WRITE (*, *) ' READING EMPIRICAL PHASE FUNCTION FILE...'
         CALL FILEOPEN(24, 3)
         IF (IEPHS == 4) THEN
            ! READS LINEFIT FORMAT ASSUMES 20 VALUES)
            JEPHS = 20
            READ (24, '(A64)') ILSHEAD
            ! DON'T READ EAPF HERE SINCE THERE COULD BE TWO DIFFERENT FILES
            READ (24, *, ERR=125) (EPHSX(I),DUMMY,EPHSF(I),I=1,JEPHS)
         ELSE
            ! ORIGINAL FILE FORMAT
            READ (24, *) JEPHS
            READ (24, *) (EPHSF(I),I=1,JEPHS)
            WRITE (16, '(/A)') 'EMPIRICAL PHASE FUNCTION COEFFICIENTS'
            WRITE (16, *) (EPHSF(I),I=1,JEPHS)
            IF (IEPHS == 1) THEN
               READ (24, *) (EPHSX(I),I=1,JEPHS)
               WRITE (16, *) (EPHSX(I),I=1,JEPHS)
            ENDIF
         ENDIF
         CALL FILECLOSE(24,2)
      ENDIF
      IF (IEAP==4 .OR. IEPHS==4) THEN
         ! IF USING BOTH IEAP=4 OR IEPHS=4 PRINT AS FOUND IN ONE FILE
         !  WHICH ASSUMES EAPX IS THE SAME AS EPHSX
         WRITE (16, 202) JEAP
         WRITE (16, 203) (EPHSX(I),EAPF(I),EPHSF(I),I=1,JEAP)
      ENDIF

!  --- INPUT ATMOSPHERIC LINE DATA FROM TAPE14
      WRITE (*, *) ' READING ATMOSPHERIC LINE LIST FILE...'
      CALL OPTLIN

!  --- PRINT OUT T-DEPENDENCE OF HALFWIDTHS
!      WRITE(16,3661)
!      WRITE(16,3662) (NAME(ICODE(I)),THALF(ICODE(I)),I=1,NGAS)

! --- FOR EACH GAS, IDENTIFY IF IT IS A RETRIEVAL GAS.   IF IT IS,
!     STORE THE STARTING PROFILE IN ARRAY XORG
      DO J = 1, NRET
         DO I = 1, NGAS
            IF (IGAS(J) == ICODE(I)) GO TO 232
         END DO
         WRITE (16, 3663) J
         WRITE (00, 3663) J
         CALL SHUTDOWN
         STOP 2
  232    CONTINUE
         IRET(J) = I

! --- CHECK FOR A 0.0 IN AN INITIAL PROFILE TO BE RETRIEVED
         DO K = 1, KMAX
            IF (XGAS(I,K) <= 0.0D0) GO TO 107
            X(J,K)    = XGAS(I,K)
            XORG(J,K) = XGAS(I,K)
         END DO
      END DO

      IF( NRET .EQ. 0 )THEN
         X(1,:KMAX) = 1.0D0
         XORG(1,:KMAX) = 1.0D0
      ENDIF

      PLANCK_C1 = 2.0D0 * C_PLANCK * V_LIGHT ** 2.0D0 * 100.0D0 ** 4.0D0
      PLANCK_C2 = 100.0D0 * C_PLANCK * V_LIGHT / C_BOLTZ

      RETURN

  107 CONTINUE
      WRITE (16, 668) NAME(ICODE(I))
      WRITE (00, 668) NAME(ICODE(I))
      CALL SHUTDOWN
      STOP 2

  120 CONTINUE
      WRITE (16, 130) JEAP, TFILE(23)
      WRITE (16, 130) JEAP, TFILE(23)
      CALL SHUTDOWN
      STOP 2

  125 CONTINUE
      WRITE (16, 135) JEPHS, TFILE(24)
      WRITE (16, 135) JEPHS, TFILE(24)
      CALL SHUTDOWN
      STOP 2

  130 FORMAT(/,' ABORT -SETUP2- ERROR READING EAP FILE. ',I5,' VALUE REQUIRED',/&
         ,' FILENAME: "',A,'"')
  135 FORMAT(/,' ABORT -SETUP2- ERROR READING EPHS FILE. ',I5,' VALUE REQUIRED',&
         /,' FILENAME: "',A,'"')

  202 FORMAT(/,' TABULAR FORM OF ILS PARAMETERS, ASSUMING N= ',I3,/,&
         '  OPD    MODULATION       PHASE')
  203 FORMAT(F7.3,2ES12.4)
  668 FORMAT(/,' ABORT -SETUP2- ZERO VMR VALUE FOUND IN RETRIEVAL GAS PROFILE: ',A7)
! 3661 FORMAT(/,/,' TEMPERATURE DEPENDENCE OF HALFWIDTHS',/,' GAS        TDEP'/)
! 3662 FORMAT(1X,A7,F6.2)
 3663 FORMAT(/,' ABORT -SETUP2- NO LINES OR PROFILE FOR RETRIEVAL GAS #',I3)

      RETURN

      END SUBROUTINE SETUP2

! ------------------------------------------------------------------------------

      SUBROUTINE SETUP3( XSC_DETAIL, NR_LEVEL )

      LOGICAL, INTENT(IN) :: XSC_DETAIL
      INTEGER, INTENT(IN) :: NR_LEVEL

! --- IF NR_LEVEL = -1 CALCULATE CROSSSECTIONS FOR ALL ALTITUDE LEVELS,
!      ELSE ONLY FOR THE LEVEL: NR_LEVEL

! --- CALCULATE VIBRATIONAL PARTITION FUNCTION FOR ALL GASES AT EACH LAYER AND AT 296K.
      IF( NR_LEVEL .EQ. -1 )WRITE (*, *) ' CALCULATING PARTITION FUNCTIONS...'
      CALL QVIB( XSC_DETAIL )

! --- COMPUTE CROSS SECTIONS FOR RETRIEVAL AND BACKGROUND GASES
      IF( NR_LEVEL .EQ. -1 )WRITE (*, *) ' CALCULATING CROSS SECTIONS...'
      CALL KROSSR( NR_LEVEL )

      RETURN

      END SUBROUTINE SETUP3

!-------------------------------------------------------------------------------

      SUBROUTINE GETSPEC( )

      INTEGER            :: IBAND, JSCAN, YYYY, MO, DD, HH, MI
      INTEGER            :: MAXPT, NREF, NPFILE, I, J, NPTSB, ISPECKODE, ISZA
      REAL(DOUBLE)       :: R4AMP, WLIM1, WLIM2, WHI, WLOW, SMM, WAVE, TAVE, SPACE
      REAL(DOUBLE)       :: SZA1, ROE1, LAT1, LON1, SECS, BSNR
      CHARACTER (LEN=80) :: TITLE


      DATA MAXPT / MMAX /

! --- SUBROUTINE TO READ ASCII ATMOSPHERIC SPECTRAL DATA
! --- FIRST_CALL IS BY BAND - ONE SPACING FOR ALL SPECTRA IN A BAND

!  --- OPEN ASCII SPECTRAL DATA
      WRITE (*, *) ' READING ASCII SPECTRA FILE: ', TFILE(15)(1:LEN_TRIM(TFILE(15)))
      CALL FILEOPEN( 15, 3 )

      WRITE (6, *) 'NFIT  BAND  SCAN/BAND  SCAN_ID  SCAN_CODE    SPACING                   RANGE         SNR'

      NFITS  = 0
      NATMOS = 0
      JSCAN  = 0
      NSPEC  = 0
      NSCAN(:MAXSPE)         = 0
      ISCAN(:MAXBND,:MAXSPE) = 0
      ISPEC(:MAXSPE)         = 0

! --- LOOP OVER BANDS AND SAVE EACH FOUND SPECTRUM
! --- BANDS ARE DEFINED IN SFIT4 INPUT FILE
! --- ALL SPECTRA FOR A BAND MUST BE IN ORDER
! --- POINT SPACING FOR THE FIRST SPECTRA BLOCK IN A BAND DEFINES THE SPACING FOR THAT BAND
! --- HERE SZA1 IS ASTRONOMICAL SZA --- RAYTRACE HAS NOT BEEN RUN

      L3: DO IBAND = 1, NBAND

   19    CONTINUE
         READ(15, *, END=21) SZA1, ROE1, LAT1, LON1, BSNR
         READ(15, *, END=21) YYYY, MO, DD, HH, MI, SECS
         READ(15, 888) TITLE

! CHECK THAT ALL NUMBERS ARE FINITE
         IF (ISNAN(SZA1).OR.ISNAN(ROE1).OR.ISNAN(LAT1).OR.ISNAN(LON1).OR.ISNAN(BSNR).OR.&
              ISNAN(SECS)) THEN
            WRITE(16,*) "NAN DETECTED IN TAPE 15 (SPECTRUM)"
            WRITE(0,*) "NAN DETECTED IN TAPE 15 (SPECTRUM)"
            CALL SHUTDOWN
            STOP 2
         END IF
         GO TO 22

   21    CONTINUE
         REWIND(15)
         CYCLE L3

   22    CONTINUE
         READ (15, *) WLOW, WHI, SPACE, NPFILE
         IF (ISNAN(WLOW).OR.ISNAN(WHI).OR.ISNAN(SPACE))THEN
            WRITE(16,*) "NAN DETECTED IN TAPE 15 (SPECTRUM)"
            WRITE(0,*) "NAN DETECTED IN TAPE 15 (SPECTRUM)"
            CALL SHUTDOWN
            STOP 2
         END IF

         WLIM1 = WAVE3(IBAND)
         WLIM2 = WAVE4(IBAND)

! -- IF NOT THIS BAND THEN DUMMY READ BLOCK AND GET NEXT
         IF( WLIM1>WHI .OR. WLIM2<WLOW )THEN
            DO I = 1, NPFILE
               READ (15, *) R4AMP
            ENDDO
            GOTO 19
         ENDIF

         ISPECKODE = NINT(1000.*SZA1)
! --- SAVE INDEX NUMBER OF SCAN (SZA) JSCAN
         L4: DO J=1, MAXSPE
            IF( ISPEC(J) .EQ. 0 .OR. ISPEC(J) .EQ. ISPECKODE )THEN
               ISPEC(J) = ISPECKODE
               ISZA = J
               EXIT L4
            ENDIF
         ENDDO L4
         IF( ISPEC(ISZA) .NE. ISPECKODE )GOTO 50

! --- RUNNING NUMBER OF SPECTRA (SZA'S)
         IF( ISZA .GT. NSPEC )NSPEC = ISZA

! --- STORE PARAMETERS FOR RAYTRACING FOR THIS SPECTRUM
! --- REDUNDANT FOR SAY SECOND BAND FROM SAME SPECTRUM
         ASTANG(ISZA)    = SZA1
         REARTH(ISZA)    = ROE1
         REFLAT(ISZA)    = LAT1
         XVB   (ISZA)    = (WLIM1 + WLIM2) / 2.

! --- FIRST INDEX THIS SPECTRUM
         NREF                  = NATMOS + 1

! --- NUMBER OF SCANS PER BAND
         NSCAN(IBAND)          = NSCAN(IBAND) + 1
         JSCAN                 = NSCAN(IBAND)

! --- FIRST INDEX OF SPECTRA IN TOBS ARRAY
         ISCNDX(1,IBAND,JSCAN) = NATMOS+ 1

! --- SCAN INDEX FOR THIS BAND / SPECTRUM
         ISCAN(IBAND,JSCAN)    = ISZA

! --- SNR BY SCAN
         SCNSNR(IBAND,JSCAN)   = BSNR

! --- RUNNING NUMBER OF FITS
         NFITS                 = NFITS + 1

! --- LOCATION OF EACH SCAN / BAND
         LOCS(IBAND,JSCAN)%YYYY = YYYY
         LOCS(IBAND,JSCAN)%MO   = MO
         LOCS(IBAND,JSCAN)%DD   = DD
         LOCS(IBAND,JSCAN)%HH   = HH
         LOCS(IBAND,JSCAN)%MI   = MI
         LOCS(IBAND,JSCAN)%SECS = SECS
         LOCS(IBAND,JSCAN)%ELON = LON1
         LOCS(IBAND,JSCAN)%NLAT = LAT1

! --- CHECK THAT SPACING FOR THIS BAND IS THE SAME FOR ALL SPECTRA
         IF( JSCAN .EQ. 1 )THEN
            SPAC(IBAND) = SPACE
         ELSE
            IF( ABS(SPAC(IBAND)-SPACE) > SPACE/100000. )THEN
               WRITE(16,106)
               WRITE(16,*) "POINT SPACING MUST BE THE SAME FOR ALL SPECTRA IN BAND"
               WRITE(00,106)
               WRITE(00,*) "POINT SPACING MUST BE THE SAME FOR ALL SPECTRA IN BAND"
               CALL SHUTDOWN
               STOP 2
            ENDIF
         ENDIF

         STITLE(NFITS) = TITLE

         WRITE(6,108) NFITS, IBAND, NSCAN(IBAND), ISZA, ISPEC(ISZA), SPAC(IBAND), WAVE3(IBAND), &
               WAVE4(IBAND), SCNSNR(IBAND,JSCAN)

         !IF( IBAND .EQ. NBAND )WRITE (31, 10) TITLE
         WRITE (6, '(4X,A76)') TITLE

         NPTSB = 0
         SMM   = 0.D0
         L5: DO I = 1, NPFILE
            READ (15, *, END=20) R4AMP
            IF (ISNAN(R4AMP))THEN
               WRITE(16,*) "NAN DETECTED IN TAPE 15 (SPECTRUM)"
               WRITE(0,*) "NAN DETECTED IN TAPE 15 (SPECTRUM)"
               CALL SHUTDOWN
               STOP 2
            END IF
            WAVE         = WLOW + REAL((I - 1),8)*SPAC(IBAND)
            IF (WAVE<WLIM1 .OR. WAVE>WLIM2) CYCLE L5
            NPTSB        = NPTSB + 1
            NATMOS       = NATMOS + 1
            IF (NATMOS > MAXPT) GO TO 40
            WWV(NATMOS)  = WAVE
            TOBS(NATMOS) = R4AMP
            SMM          = SMM + TOBS(NATMOS)
            IF (NPTSB == 1) WSTART(IBAND) = WAVE
            WSTOP(IBAND) = WAVE
         ENDDO L5

   20    CONTINUE

         IF( NPTSB .EQ. 0 )GOTO 60
         ISCNDX(2,IBAND,JSCAN) = NATMOS

!  --- ADJUST POINT SPACING FOR MONOCHROMATIC CALCULATIONS
!         CHOOSE SPACING SO THAT SPECTRAL DATA POINT SPACING
!         IS A MULTIPLE OF THE MONOCHROMATIC POINT SPACING
         IF( NSCAN(IBAND) .EQ. 1 )THEN
            WRITE(16,109) IBAND
            WRITE(16,110) WAVFAC(IBAND), PMAX(IBAND), FOVDIA(IBAND), DN(IBAND), IAP(IBAND), &
                  (SCNSNR(IBAND,JSCAN),JSCAN=1,NSCAN(IBAND))
            NSPAC(IBAND) = FLOOR(SPAC(IBAND)/DN(IBAND) + 1.0000001D0)
            DN(IBAND)    = SPAC(IBAND)/NSPAC(IBAND)
            NPRIM(IBAND) = NPTSB
            WRITE(16,566) DN(IBAND), NSPAC(IBAND)
         ENDIF

         WRITE(16,*) ""
         WRITE(16,10) STITLE(NFITS)
         WRITE(16,12) ISPEC(ISZA)
         WRITE(16,11) WLOW, WHI, SPAC(IBAND), NPFILE

         WRITE(16,102) WSTART(IBAND), WSTOP(IBAND), NPTSB, NATMOS

! ---  NORMALIZE AMPLITUDES TO AVERAGE VALUE IF ABSORPTION MEASUREMENTS
         IF( IEMISSION .EQ. 0 .OR. IENORM(IBAND) .eq. 1) THEN
            TAVE              = SMM/REAL(NPTSB,8)
            TOBS(NREF:NATMOS) = TOBS(NREF:NATMOS)/TAVE
            !print *, 'tave ', tave, NREF, NATMOS
            !print *, maxval(TOBS(NREF:NATMOS))
            !TOBS(NREF:NATMOS) = TOBS(NREF:NATMOS)/ maxval(TOBS(NREF:NATMOS))
         END IF

! --- IF WE GET HERE WE NEED TO READ ANOTHER BLOCK IN T15ASC
         GOTO 19

      ENDDO L3

      WRITE(16,111) NSPEC, (ISPEC(I),I=1,NSPEC)

! --- ONLY USE ONE VALUE OF LATIDUTE FOR ALL RAYTRACES
      REF_LAT = REFLAT(1)

      CLOSE(15)
      RETURN

! --- BUFFER TRUNCATION
   40 CONTINUE
      WRITE (16, 101) MAXPT
      WRITE (16, 102) WSTART(IBAND), WSTOP(IBAND), NPTSB, NATMOS
      WRITE (00, 101) MAXPT
      WRITE (00, 102) WSTART(IBAND), WSTOP(IBAND), NPTSB, NATMOS
      CALL SHUTDOWN
      STOP 2

! --- TOO MANY SPECTRA
   50 CONTINUE
      WRITE(16,105) MAXSPE
      WRITE(00,105) MAXSPE
      CALL SHUTDOWN
      STOP 2

! --- NO POINTS FOUND
   60 CONTINUE
      WRITE (16, 103)
      WRITE (00, 103)
      CALL SHUTDOWN
      STOP 2

   10 FORMAT(1X,A80)
   11 FORMAT(  ' FIRST POINT (CM-1)                   : ',F12.4, /, &
               ' LAST POINT (CM-1)                    : ',F12.4, /, &
               ' POINT SPACING (CM-1)                 : ',F12.8, /, &
               ' NUMBER OF POINTS                     : ',I12)
   12 FORMAT(  ' SPECTRUM CODE                        : ',I12)
  102 FORMAT(  ' WSTART                               : ',F12.4, /, &
               ' WSTOP                                : ',F12.4, /, &
               ' NPRIME                               : ',I12, /, &
               ' NATMOS                               : ',I12)
  110 FORMAT(  ' WAVFAC                               : ',F12.7,/,&
               ' MAX OPD [CM-1]                       : ',F12.2,/,&
               ' FOV [MR]                             : ',F12.5,/,&
               ' REQUESTED POINT SPACING              : ',F12.7,/,&
               ' APODIZATION CODE                     : ',I12,/,  &
               ' SNR                                  : ',40F12.7)
  109 FORMAT(/,' BAND                                 : ',I12)
  111 FORMAT(/,' NUMBER OF UNIQUE SPECTRA             : ',I12, /,' SZA CODES:',/, 40I10)
  566 FORMAT(  ' MONOCHROMATIC SPACE                  : ',F12.8, /, &
               ' NSPAC                                : ',I12)
  108 FORMAT(I5,I6,I11,I9,I11,1X,F10.7, 2X, F10.3,' -',F10.3,F12.5)

  101 FORMAT(' GETSPEC: ABORT-SPECTRAL DATA ARRAY SIZE LIMIT EXCEEDED-MAX =',I5)
  103 FORMAT(/,/,5X,'GETSPEC: ABORT...NO POINTS'/,/)
  105 FORMAT(" GETSPEC: ABORT... ATTEMPT TO READ TOO MANY SPECTRA. MAX = ",I5)
  106 FORMAT(" GETSPEC: ABORT... POINT SPACING MUST BE THE SAME FOR ALL SPECTRA")

  888 FORMAT(A80)

      RETURN

      END SUBROUTINE GETSPEC

!-------------------------------------------------------------------------------

      SUBROUTINE FILSE( SED, NFIT )

! SFIT 4
! NO DEFAULT - INPUT ONE SNR FOR EACH BAND
! DON'T LIKE IT  - SHOULD BE IN T15 ASCII FILE AS A SNR PER BAND/SPECTRUM
! IS NOW MAR 2013 vP1.7

! --- FILLS SED VECTOR (VARIANCES OF MEASUREMENTS)
! --- USES DEFAULT S/N UNLESS A DIFFERENT VALUE IS SPECIFIED
! --- FOR A GIVEN WAVENUMBER INTERVAL
! --- MAXSNR=MAXIMUM NUMBER OF ALTERNATE S/N VALUES
! --- MMAX=MAXIMUM NUMBER OF SPECTRAL DATA POINTS

      INTEGER, INTENT(IN)       :: NFIT
      !REAL(DOUBLE), INTENT(IN)  :: SNR
      REAL(DOUBLE), INTENT(OUT) :: SED(NFIT)

      INTEGER :: I, K, IBAND, JSCAN, IW
      REAL(DOUBLE), ALLOCATABLE, DIMENSION(:) :: STNR

      ALLOCATE(STNR(NFIT))

      IF( NSTNR .NE. 0 )THEN
         WRITE(16,11) NSTNR
         DO I = 1, NSTNR
            WRITE(16,12) WWV0(I), WWV1(I), GSTNR(I)
         END DO
      ENDIF

      DO IW = 1, NFIT
         DO IBAND = 1, NBAND
            IF ((WWV(IW) .LT. WAVE3(IBAND)) .OR. (WWV(IW) .GT. WAVE4(IBAND))) CYCLE
            DO JSCAN=1, NSCAN(IBAND)
               IF(( IW.GE.ISCNDX(1,IBAND,JSCAN)) .AND. (IW.LE.ISCNDX(2,IBAND,JSCAN)) )THEN
                  IF ((IEMISSION .eq. 0) .OR. (IENORM(IBAND) .eq. 1)) THEN
                     STNR(IW) = 1.0D0 / SCNSNR(IBAND,JSCAN)
                  ELSE
                     STNR(IW) = SCNSNR(IBAND,JSCAN)
                  ENDIF
                  DO K = 1, NSTNR
                     IF ((WWV(IW) .LT. WWV0(K)) .OR. (WWV(IW) .GT. WWV1(K))) CYCLE
                     IF ((IEMISSION .EQ. 0) .OR.( IENORM(IBAND) .eq. 1)) THEN
                        STNR(IW) = 1.0D0 / GSTNR(K)
                     ELSE
                        STNR(IW) = GSTNR(IBAND)
                     ENDIF
                  ENDDO ! NSTNR
               ENDIF ! ISCNDX
            ENDDO ! NSCAN
         ENDDO ! NBAND
      ENDDO ! NFIT
      !DO I=1, NFIT
      !PRINT*, I, WWV(I), 1.0/STNR(I)
      !ENDDO

! --- MEAN SNR
      SNR = SUM( 1.0D0/STNR(:NFIT) )/REAL(NFIT,8)
      SED = STNR(:NFIT)*STNR(:NFIT)

!      IF (F_WRTSEFILE) THEN
!         CALL FILEOPEN( 67, 2 )
!         WRITE( 67, *) NFIT
!         WRITE( 67, *) (SED(I),I=1,g)
!         CALL FILECLOSE( 67, 1 )
!      ENDIF

      DEALLOCATE( STNR )

      RETURN

 11   FORMAT(/,"SNR DE-WEIGHTING INTERVALS : ", I5,/"     LOW WAVEN#     HIGH WAVE#        NEW SNR")
 12   FORMAT(2F15.5,F15.3)

    END SUBROUTINE FILSE


!-------------------------------------------------------------------------------

    SUBROUTINE INIT_PARM()

! --- INITIALIZE PARM VECTOR - FITTED PARAMETERS FOR FORWARD MODEL

      IMPLICIT NONE

      INTEGER :: I, KK, N

      NVAR = 0

      !  --- BACKGROUND FITTING
      IF( F_BACKG )THEN
         IF (NBACK == 2) THEN
            IF (NFITS > 0) THEN
               !  --- BACKGROUND SLOPE - NBACK=2
               do i = 1,nfits
                  WRITE(PNAME(NVAR+I), '(a10,i1)') 'BckGrdSlp_', i
               end do
               PARM(NVAR+1:NFITS) = BCKSL
               SPARM(NVAR+1:NFITS) = SBCKSL
               !  --- BACKGROUND CURVATURE - NBACK=3
               NVAR = NFITS
            ENDIF
         ELSE
            IF (NBACK == 3) THEN
               IF (NFITS > 0) THEN
                  !  --- BACKGROUND SLOPE - NBACK=2
                  do i = 1,nfits
                     write(PNAME(I*2-1+NVAR), '(a10,i1)') 'BckGrdSlp_', i
                  end do
                  PARM(NVAR+1:NFITS*2-1+NVAR:2) = BCKSL
                  SPARM(NVAR+1:NFITS*2-1+NVAR:2) = SBCKSL
                  !  --- BACKGROUND CURVATURE - NBACK=3
                  do i = 1,nfits
                     write(PNAME(I*2+NVAR), '(a10,i1)') 'BckGrdCur_', i
                  end do
                  PARM(NVAR+2:NFITS*2+NVAR:2) = BCKCRV
                  SPARM(NVAR+2:NFITS*2+NVAR:2) = SBCKCRV
                  NVAR = NFITS*2 + NVAR
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      NBKFIT = (NBACK - 1)*NFITS

      ! ---  WAVENUMBER SHIFT
      IF (F_WSHIFT) THEN

         !  --- SINGLE WAVELENGTH SHIFT FOR ALL BANDPASSES
         IF (ISPARM <= 1) THEN
            NVAR = NVAR + 1
            PNAME(NVAR) = 'SWNumShft'
            PARM(NVAR) = WSHFT
            SPARM(NVAR) = SWSHFT
            NSHIFT = 1
         ELSE
            !  --- INDEPENDENT WAVELENGTH SHIFT FOR EACH BANDPASS (ISPARM=2)
            !  --- OR INDEPENDENT WAVELENGTH SHIFT FOR EACH FITTING (ISPARM=3)
            IF (ISPARM == 2) N = NBAND
            IF (ISPARM == 3) N = NFITS
            IF (N > 0) THEN
               do i = 1,N
                  WRITE(PNAME(NVAR+I), '(a10,i1)') 'IWNumShft_', i
               end do
               PARM(NVAR+1:N+NVAR) = WSHFT
               SPARM(NVAR+1:N+NVAR) = SWSHFT
               NVAR = N + NVAR
            ENDIF
            NSHIFT = N
         ENDIF
      ELSE
         ISPARM = 0
      ENDIF

      !  --- FIT ZERO LEVELS IF IZERO=1
      !  --- TOTAL NUMBER OF ZERO LEVEL FITS=NZERO
      NZERO = 0
      DO I = 1, NBAND
         IF (F_ZSHIFT(I)) THEN
            IF (IZERO(I) .NE. 1 ) CYCLE
            N = NSCAN(I)
            IF (N > 0) THEN
               do kk = 1, n
                  write(PNAME(kk+NVAR),'(a8,i1)') 'ZeroLev_',i
               end do
               PARM(NVAR+1:N+NVAR) = ZSHIFT(I,1)
               SPARM(NVAR+1:N+NVAR) = SZERO(I)
               NVAR = N + NVAR
               NZERO = N + NZERO
            ENDIF
         ELSE
            IZERO(I) = 0
         END IF
      ENDDO

      !  --- SOLAR LINES INCLUSION
      NSOLAR = 0
      IF( IFCO )THEN
         DO I = 1, 5
            IF( .NOT. F_RTSOL(I) )CYCLE
            NVAR = NVAR + 1
            NSOLAR = NSOLAR + 1
            IF (NSOLAR .EQ. 1) NSOLAR1 = NVAR
            PNAME(NVAR) = CPNAM(I)
            PARM(NVAR)  = CIPARM(I)
            SPARM(NVAR) = SCPARM(I)
         END DO
      ENDIF

      !  --- EMPIRICAL APODIZATION
      IF( F_RTAPOD )THEN
         NEAPRT = NEAP
         IF (NEAPRT > 0) THEN
            EAPF0(:NEAPRT) = EAPF(:NEAPRT)
            PNAME(NVAR+1:NEAPRT+NVAR) = 'EmpApdFcn'
            PARM(NVAR+1:NEAPRT+NVAR) = EAPPAR
            SPARM(NVAR+1:NEAPRT+NVAR) = SEAPPAR
            NVAR = NEAPRT + NVAR
         ENDIF
      ENDIF

      !  --- EMPIRICAL PHASE FUNCTION
      IF( F_RTPHASE )THEN
         NEPHSRT = NEPHS
         IF (NEPHSRT > 0) THEN
            EPHSF0(:NEPHSRT) = EPHSF(:NEPHSRT)
            PNAME(NVAR+1:NEPHSRT+NVAR) = 'EmpPhsFnc'
            PARM(NVAR+1:NEPHSRT+NVAR) = EPHSPAR
            SPARM(NVAR+1:NEPHSRT+NVAR) = SEPHSPAR
            NVAR = NEPHSRT + NVAR
         ENDIF
      ENDIF

      !  --- DIFFERENTIAL WAVENUMBER SHIFT FOR RETRIEVAL GASES
      NDIFF = 0
      IF (IFDIFF .AND. NRET.GT.1) THEN
         do kk = 2, nret
            PNAME(kk+NVAR-1) = 'DWNumShft_'//trim(GAS(kk))
         end do
         PARM(NVAR+1:NRET-1+NVAR)  = WSHFT
         SPARM(NVAR+1:NRET-1+NVAR) = SWSHFT
         NDIFF = NRET - 1
         NVAR  = NRET - 1 + NVAR
      ENDIF

      !  --- SKIP IF NOT IFPHASE; FIT PHASE ERRORS IF IFPHASE=T
      !  --- TOTAL NUMBER OF PHASE ERROR FITS=NPHASE
      NPHASE = 0
      IF( IFPHASE )THEN
         DO I = 1, NBAND
            N = NSCAN(I)
            IF (N > 0) THEN
               do kk = 1, n
                  write(PNAME(kk+NVAR),'(a8,i1)') 'SPhsErr_',i
               end do
               PARM(NVAR+1:N+NVAR) = PHS
               SPARM(NVAR+1:N+NVAR) = SPHS
               NVAR = N + NVAR
               NPHASE = N + NPHASE
            ENDIF
         END DO
      ENDIF

      ! --- LINE PARAMETER RETRIEVAL
      ! --- LINE INTENSITY
      if (IFLINE /= 0) THEN
         IF ( NILINE /= 0 ) then
            PNAME(NVAR+1:NVAR+NILINE) = 'LineInt'
            PARM(NVAR+1:NVAR+NILINE)  = 0.0D0
            SPARM(NVAR+1:NVAR+NILINE) = 1.0D0
            NVAR = NVAR + NILINE
         end IF
         IF ( NPLINE /= 0 ) then
            PNAME(NVAR+1:NVAR+NILINE) = 'LinePAir'
            PARM(NVAR+1:NVAR+NPLINE)  = 0.0D0
            SPARM(NVAR+1:NVAR+NPLINE) = 1.0D0
            NVAR = NVAR + NPLINE
         end IF
         IF ( NTLINE /= 0 ) then
            PNAME(NVAR+1:NVAR+NTLINE) = 'LineTAir'
            PARM(NVAR+1:NVAR+NTLINE)  = 0.0D0
            SPARM(NVAR+1:NVAR+NTLINE) = 1.0D0
            NVAR = NVAR + NTLINE
         end IF
      end if

      if (ifsza /= 0) then
         PNAME(NVAR+1:NVAR+NSPEC) = 'SZA'
         PARM(NVAR+1:NVAR+NSPEC)  = 0.0D0
         SPARM(NVAR+1:NVAR+NSPEC) = 1.0D0
         NVAR = NVAR + NSPEC
      end IF

      do i = 1,nband
         if (iffov /= 0) then
            write(PNAME(NVAR+1:NVAR+2), '(a4,i1)'), 'FOV_', i
            PARM(NVAR+1:NVAR+2)  = 0.0D0
            SPARM(NVAR+1:NVAR+2) = 1.0D0
            NVAR = NVAR + 1
         end if
         if (ifopd /= 0) then
            write(PNAME(NVAR+1:NVAR+2), '(a4,i1)'), 'OPD_', i
            PARM(NVAR+1:NVAR+2)  = 0.0D0
            SPARM(NVAR+1:NVAR+2) = 1.0D0
            NVAR = NVAR + 1
         end if
      end do

      !  ---  RETRIEVAL GAS MIXING RATIOS
      !  ---  MIXING RATIO AT ISMIX +1
      ISMIX = NVAR
      DO KK = 1, NRET
         NGIDX(KK,1,0) = NVAR + 1
         IF( IFPRF(KK) )THEN
            !  --- RETRIEVING VERTICAL PROFILE
            N = NLAYERS
            PNAME(NVAR+1:N+NVAR) = NAME(IGAS(KK)) ! L->H ALT
            IF (ILOGRETRIEVAL(KK)/=0) THEN !MP
               ! WE WANT A LOGARITHMIC PROFILE
               PARM(NVAR+1:N+NVAR) = LOG(XORG(KK,:N))
            ELSE
               PARM(NVAR+1:N+NVAR) = COLSF(KK) !1.D0
            ENDIF
            SPARM(NVAR+1:N+NVAR) = SIG(:N,KK)
            NVAR = NVAR + N
         ELSE
            ! --- SCALING VERTICAL DISTRIBUTION
            NVAR = NVAR + 1
            PNAME(NVAR) = NAME(IGAS(KK))
            PARM(NVAR)  = COLSF(KK)
            SPARM(NVAR) = SCOLSF(KK)
         ENDIF
         NGIDX(KK,2,0) = NVAR
      END DO


      ! --- TEMPERATURE RETRIEVAL
      IF( IFTEMP )THEN
         NTEMP = NLAYERS
         PNAME(NVAR+1:NVAR+NTEMP) = 'TEMPERAT'
         PARM(NVAR+1:NVAR+NTEMP)  = 1.D0
         SPARM(NVAR+1:NVAR+NTEMP) = TSIGMA(1:NTEMP)
         NTEMP1 = NVAR + 1
         NVAR = NVAR + NTEMP
      ENDIF


      ! --- INSERT  CHANNEL PARAMETERS INTO STATE VECTOR PARM()
      CALL INSERT_CHANNEL_PARMS (NVAR, PARM, PNAME, SPARM)

      ! --- TEST FOR OVERFLOWS
      IF (NVAR > NMAX) GO TO 556
      NFIT = NATMOS
      IF (NFIT > MMAX) GO TO 557

      ! CHECK FOR THE FINITE VALUES OF SPARM IF THE SETUP IS FOR RETRIEVAL
      IF( RETFLG .AND. .NOT.(ALL(SPARM(:NVAR).GT.TINY(SPARM(1)))) )THEN
         WRITE(*,*) 'FOUND SIGMA VALUE EQUAL TO ZERO'
         WRITE(16,*) 'FOUND SIGMA VALUE EQUAL TO ZERO'
         DO I=1, NVAR
            WRITE(16,*) I, PNAME(I), PARM(I), SPARM(I)
            WRITE(00,*) I, PNAME(I), PARM(I), SPARM(I)
         ENDDO
         CALL SHUTDOWN
         STOP 2
      END IF

      RETURN

556   CONTINUE
      WRITE (16, 223)
      WRITE (16, *) 'NVAR =', NVAR, ' TO BIG, NMAX = ', NMAX
      WRITE (00, 223)
      WRITE (00, *) 'NVAR =', NVAR, ' TO BIG, NMAX = ', NMAX
      CALL SHUTDOWN
      STOP 2

557   CONTINUE
      WRITE (16, 224)
      WRITE (16, *) 'NFIT=', NFIT, 'TO BIG, MAX FIT PARMS =', MMAX
      WRITE (00, 224)
      WRITE (00, *) 'NFIT=', NFIT, 'TO BIG, MAX FIT PARMS =', MMAX
      CALL SHUTDOWN
      STOP 2

223   FORMAT(/,' ABORT--TOO MANY VARIABLES')
224   FORMAT(/,' ABORT--TOO MANY SPECTRAL DATA POINTS')


      END SUBROUTINE INIT_PARM


      SUBROUTINE FILSA( SA )

! --- CREATE SA FROM SPARM VECTOR
! --- FILL OFF AXIS VALUES FOR RETRIEVAL GASES AND TEMPERATURE

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(INOUT)   :: SA(:,:)
      LOGICAL                       :: FILOPEN = .FALSE.
      INTEGER                       :: I, J, KK, N, JROW, JCOL, INDXX
      REAL(DOUBLE)                  :: TSAHWD, DELZ = 0.0D0, RHO  = 0.0D0

!  --- FILL DIAGONAL ELEMENTS OF SA
      DO I = 1, NVAR
         SA(I,I) = SPARM(I)*SPARM(I)
      ENDDO

!  --- OFF DIAGONAL ELEMENTS OF SA MATRIX (A PRIORI COMPONENTS)
      INDXX = ISMIX
      DO KK = 1, NRET
         !print *, 'fill off diag ', kk, IFPRF(KK)
         N = 1
         IF( IFPRF(KK) ) THEN
            N = NLAYERS
            SELECT CASE ( IFOFF(KK) )
            CASE ( 1:3 )
!  --- FILL OFF DIAGONAL ELEMENTS OF SA
            DO I = 1, NLAYERS
              DO J = 1, NLAYERS
                IF (I == J) CYCLE
                JROW = I + INDXX
                JCOL = J + INDXX
                IF( ZBAR(I) < ZGMIN(KK) ) CYCLE
                IF( ZBAR(J) < ZGMIN(KK) ) CYCLE
                IF( ZBAR(I) > ZGMAX(KK) ) CYCLE
                IF( ZBAR(J) > ZGMAX(KK) ) CYCLE

                DELZ = ZBAR(I) - ZBAR(J)

                SELECT CASE ( IFOFF(KK) )
                CASE (1)       !gaussian
                  RHO = (ALOGSQ*DELZ/ZWID(KK))**2
                  RHO = MIN( RHO, 90.0D0 )
                  SA(JROW,JCOL) = SPARM(JROW)*SPARM(JCOL)*EXP(-RHO)
                CASE (2)       !exponential
                  RHO = ABS(ALOGSQ*DELZ/ZWID(KK))
                  RHO = MIN( RHO, 90.0D0 ) !664
                  SA(JROW,JCOL) = SPARM(JROW)*SPARM(JCOL)*EXP(-RHO)
                CASE (3)
                  WRITE(16,*) "IFOFF=3 NOT SUPPORTED"
                  WRITE(00,*) "IFOFF=3 NOT SUPPORTED"
                  CALL SHUTDOWN
                  STOP 2
                END SELECT

              END DO
            END DO
! --- READ IN FULL COVARIANCE FROM FILE
            CASE ( 4 )
               INQUIRE( UNIT=62, OPENED=FILOPEN )
               IF ( .NOT. FILOPEN )CALL FILEOPEN( 62, 3 )
               DO I = 1, NLAYERS
                  READ( 62,* ) (SA( I+INDXX, J+INDXX), J=1, N)
               END DO
            CASE ( 0 )
               PRINT *, ' PROFILE RETRIEVAL GAS: ', NAME(IGAS(KK)), ' NO OFF DIAGONAL VALUES SET.'
            END SELECT
         ENDIF
         INDXX = INDXX + N
      END DO

      CALL FILECLOSE( 62, 2 )


!  --- FILL OFF DIAGONAL ELEMENTS OF SA of T
      IF( IFTEMP )THEN
         INDXX = NTEMP1
            DO I = 1, NLAYERS
              DO J = 1, NLAYERS
                IF (I == J) CYCLE
                JROW = I + INDXX
                JCOL = J + INDXX
                DELZ = ZBAR(I) - ZBAR(J)
                TSAHWD = 20.0D0
                SELECT CASE ( 1 )
                CASE (1)       !gaussian
                  RHO = (ALOGSQ*DELZ/TSAHWD)**2
                  RHO = MIN( RHO, 90.0D0 )
                  SA(JROW,JCOL) = SPARM(JROW)*SPARM(JCOL)*EXP(-RHO)
                CASE (2)       !exponential
                  RHO = ABS(ALOGSQ*DELZ/TSAHWD)
                  RHO = MIN( RHO, 90.0D0 ) !664
                  SA(JROW,JCOL) = SPARM(JROW)*SPARM(JCOL)*EXP(-RHO)
                CASE (3)
                  WRITE(16,*) "IFOFF=3 NOT SUPPORTED"
                  WRITE(00,*) "IFOFF=3 NOT SUPPORTED"
                  CALL SHUTDOWN
                  STOP 2
                END SELECT
              END DO
            END DO
       ENDIF

!  --- WRITE OUT FULL SA MATRIX
      IF (F_WRTSA) THEN
         CALL FILEOPEN( 63, 1 )
         WRITE(63,*) TRIM(TAG), ' FULL INITIAL STATE VECTOR COVARIANCE N X N MATRIX'
         WRITE(63,*) NVAR, NVAR
         WRITE(63,260) ADJUSTR(PNAME(:NVAR))
         DO I=1,NVAR
            WRITE(63,261) (SA(I,J),J=1,NVAR)
         END DO
         CALL FILECLOSE( 63, 1 )
      ENDIF

      RETURN

  260 FORMAT( 2000( 12X, A14 ))
  261 FORMAT( 2000ES26.18 )


      END SUBROUTINE FILSA




!-------------------------------------------------------------------------------
      SUBROUTINE RELEASE_MEM_INT

      IF( ALLOCATED( CROSS )       )DEALLOCATE (CROSS)
      IF( ALLOCATED( CROSS_FACMAS ))DEALLOCATE (CROSS_FACMAS)
      IF( ALLOCATED( TCO )         )DEALLOCATE (TCO)
      IF( ALLOCATED( TCONV )       )DEALLOCATE (TCONV)
      IF( ALLOCATED( TCALC )       )DEALLOCATE (TCALC)
      IF( ALLOCATED( TCALC_I )     )DEALLOCATE (TCALC_I)
      IF( ALLOCATED( TCALC_E )     )DEALLOCATE (TCALC_E)
      IF( ALLOCATED( TCALC_S )     )DEALLOCATE (TCALC_S)
      IF( ALLOCATED( IMGG )        )DEALLOCATE (IMGG)

      RETURN

      END SUBROUTINE RELEASE_MEM_INT

      END MODULE INITIALIZE
