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

      MODULE READIN

      USE PARAMS
      USE RETVPARAM
      USE TRANSMIS
      USE MOLCPARAM
      USE XSECTIONS
      USE DATAFILES
      USE SYNSPEC
      USE LINEPARAM
      USE SOLAR
      USE BANDPARAM
      USE INITIALIZE
      USE OPT
      USE CHANNEL

      IMPLICIT NONE

      CONTAINS


!----------------------------------------------------------------------
      SUBROUTINE READCK1(NLEV, NEGFLAG)

      USE ISOTOPE

      IMPLICIT NONE

      INTEGER, INTENT(OUT)  :: NLEV
      INTEGER, INTENT(OUT)  :: NEGFLAG

      INTEGER   :: I, NRMAX, NPGAS, J, N

      NRMAX = MOLMAX
      NPGAS = 0

! -- CHECK THAT EVERY GAS IS ONLY ONCE IN THE RETRIEVAL LIST

      DO I=1, NRET
         DO J=1, NRET
            if ((I.ne.J).and.(GAS(I).EQ.GAS(J))) THEN
               WRITE (*, *) "GAS ", TRIM(GAS(I)), " DEFINED TWICE IN GAS...LIST"
               WRITE (16, *) "GAS ", TRIM(GAS(I)), " DEFINED TWICE IN GAS...LIST"
               exit
            end if
         end DO
      end DO
! --- DOUBLE CHECK CHECK THAT PROFILE REIEVALS ARE AHEAD OF COLUMNS IN LIST
      I=0
      DO J=1, NRET
        IF( IFPRF(J) ) I=I+1
      ENDDO
      DO J=1, I
! --- SHOULD NOT GET HERE DUE TO CHECKS IN BINPUT_PARSE...
         IF( .NOT. IFPRF(J) )THEN
            WRITE(16,*) 'PUT COLUMN RETREAVAL GAS AFTER LAST PROFILE GAS'
            WRITE(00,*) 'PUT COLUMN RETREAVAL GAS AFTER LAST PROFILE GAS'
            CALL SHUTDOWN
            STOP '2'
         ENDIF
      ENDDO

      WRITE (16, 402) NLAYERS
      WRITE (16, 400) NRMAX, NRET
      WRITE (16, 410) USEISO

      IF (NLAYERS .NE. NLEV) THEN
         WRITE (16, *) "NUMBER OF LAYERS FROM INPUT  ",NLAYERS," FOR GAS ",GAS(J)
         WRITE (16, *) "DOES NOT MATCH LAYERS FROM STATION.LAYERS FILE (USED IN RAYTRACING) ",NLEV
         WRITE (00, *) "NUMBER OF LAYERS FROM INPUT  ",NLAYERS," FOR GAS ",GAS(J)
         WRITE (00, *) "DOES NOT MATCH LAYERS FROM STATION.LAYERS FILE (USED IN RAYTRACING) ",NLEV
         CALL SHUTDOWN
         STOP '2'
      ENDIF

! --- SEE IF WE NEED TO SEPARATE OUT ISOTOPES
      IF ( USEISO ) CALL RDISOFILE( 16 )

      IF( NRET .LE. NRMAX .AND. NRET .GE. 1 )THEN
         DO J = 1, NRET
            WRITE (16, 600) J, GAS(J)
            !print *,J, GAS(J)
            DO I = 1, MOLTOTAL
!               write(*,*) gas(j), name(i)
               IF (GAS(J) == NAME(I)) GO TO 176
            END DO
            WRITE (16, 610) GAS(J)
            WRITE (00, 610) GAS(J)
            CALL SHUTDOWN
            STOP '2'

  176       CONTINUE
            IGAS(J) = I
            WRITE (16, 601) IFPRF(J)
            IF( .NOT. IFPRF(J) )THEN
! --- FOR COLUMN RETRIEVAL THE LOG FUNCTION IS SHUT OFF AUTOMATICALLY
!               ILOGRETRIEVAL(J) = 0
               WRITE (16, 401) COLSF(J), SCOLSF(J)
               !write(*,'(a, 2i5,a10,2f12.4)') ' readck1 ',nret, j, gas(j), COLSF(J), SCOLSF(J)
               IF( ABS(COLSF(J)) .LT. TINY(0.0D0) )THEN
                  WRITE (16, *) "APRIORI VMR SCALE FOR COLUMN GAS: ",GAS(J), " IS NOT SET IN SFIT4.CTL FILE"
                  WRITE ( 0, *) "APRIORI VMR SCALE FOR COLUMN GAS: ",GAS(J), " IS NOT SET IN SFIT4.CTL FILE"
                  CALL SHUTDOWN
                  STOP '2'
               ENDIF
               IF( ABS(SCOLSF(J)) .LT. TINY(0.0D0) )THEN
                  WRITE (16, *) "COLUMN SIGMA FOR GAS: ",GAS(J), " IS NOT SET IN SFIT4.CTL FILE"
                  WRITE ( 0, *) "COLUMN SIGMA FOR GAS: ",GAS(J), " IS NOT SET IN SFIT4.CTL FILE"
                  CALL SHUTDOWN
                  STOP '2'
               ENDIF

               CYCLE
            ENDIF

! --- THIS WILL BE A PROFILE RETRIEVAL
            NPGAS = NPGAS + 1
            IF (NPGAS > MAXPRF) GO TO 301

! --- DIAGONAL ELEMENTS OF SA OR FILENAME OF FULL COVARIANCE
            IF( .NOT. CORRELATE(J) )IFOFF(J) = 0
            IF( ABS(COLSF(J)) .LT. TINY(0.0D0) )THEN
               WRITE (16, *) "APRIORI VMR SCALE PROFILE FOR  ",GAS(J), " IS NOT SET IN SFIT4.CTL FILE"
               WRITE ( 0, *) "APRIORI VMR SCALE PROFILE FOR  ",GAS(J), " IS NOT SET IN SFIT4.CTL FILE"
               CALL SHUTDOWN
               STOP '2'
            ENDIF
            WRITE (16, 611) COLSF(J)
            SELECT CASE ( IFOFF(J) )
            CASE (0)    ! NO INTERLAYER CORRELATION
               WRITE (16, 613)
               WRITE (16, 612) (SIG(I,J),I=1,NLAYERS)
               WRITE (16, 403)
            CASE (1)    ! GAUSSIAN ILC (ORIGINAL)
               WRITE (16, 613)
               WRITE (16, 612) (SIG(I,J),I=1,NLAYERS)
               WRITE (16, 614) ZWID(J), ZGMIN(J), ZGMAX(J)
            CASE (2)    ! EXPONENTIAL ILC
               WRITE (16, 613)
               WRITE (16, 612) (SIG(I,J),I=1,NLAYERS)
               WRITE (16, 615) ZWID(J), ZGMIN(J), ZGMAX(J)
            !CASE (3)    ! NOT USED

            CASE (4)    ! READ IN FILE AS FULL SA ( SA.INPUT )
               WRITE (16, 617 ) N, N, TRIM( TFILE(62) )
               WRITE (16, 612) (SIG(I,J),I=1,NLAYERS)
               !SIG( 1:N, J ) = 0.0D0
            CASE (5)    ! READ IN FILE AS FULL SA INVERSE (SA.INPUT )
               WRITE (16, 618 ) N, N, TRIM( TFILE(62) )
               WRITE (16, 612) (SIG(I,J),I=1,NLAYERS)
               !SIG( 1:N, J ) = 0.0D0
            CASE DEFAULT
               WRITE(16,*) ' READCK1: FLAG IFOFF MUST BE ONE OF 0, 1, 2, 4, 5'
               WRITE(00,*) ' READCK1: FLAG IFOFF MUST BE ONE OF 0, 1, 2, 4, 5'
               CALL SHUTDOWN
               STOP '2'
            END SELECT

            WRITE(16,619) ILOGRETRIEVAL(J)
         END DO

         WRITE (16, 620) DELNU

         WRITE(16,*)''
         WRITE(16,*) ' LINE SHAPE MODEL:'
         SELECT CASE ( LSHAPEMODEL )  ! USER CHOICE OF LINE SHAPE MODEL
            CASE (0)
               WRITE (16,*) '  0 = CHOOSE MODEL DEPENDING ON EXISTANCE OF PARAMETERS'
            CASE (1)
               WRITE (16,*) '  1 = FORCE VOIGT FOR ALL LINES'
            CASE (2)
               WRITE (16,*) '  2 = USE GALATRY FOR LINES WITH PARAMETERS, VOIGT ELSE'
            CASE (3)
               WRITE (16,*) '  3 = USE SDV & LINE MIXING FOR LINES WITH PARAMETERS'
            CASE DEFAULT
               WRITE(16,*)' LINE SHAPE MODEL FLAG OUT OF RANGE MIST BE 0, 1, 2, 3)'
               WRITE(00,*)' LINE SHAPE MODEL FLAG OUT OF RANGE MIST BE 0, 1, 2, 3)'
               CALL SHUTDOWN
               STOP '2'
         END SELECT

         NEGFLAG = -1
         IF( ITRMAX .LT. 0 ) THEN
            ITRMAX = IABS(ITRMAX)
            NEGFLAG = 0
         ENDIF
         WRITE (16, 650) ITRMAX

      ELSE IF( NRET .EQ. 0 )THEN
         WRITE(16,630)
      ENDIF

      RETURN

      IF( CONVERGENCE .LT. 0.0 )THEN
         WRITE(16,651)
         WRITE(00,651)
         CALL SHUTDOWN
         STOP '2'
      ENDIF

!      WRITE (16, 605) NRMAX
!      CLOSE(16)
!      STOP

  301 CONTINUE
      WRITE (16, 606) NPGAS, MAXPRF
      WRITE (00, 606) NPGAS, MAXPRF
      CALL SHUTDOWN
      STOP '2'


  402 FORMAT(/,' NUMBER OF RETRIEVAL LAYERS IN INPUT VARIANCE VECTORS : ', I5 )
  400 FORMAT(  ' NUMBER OF RETRIEVAL GASES (MAX=',I2,')                   : ', I5 )
  410 FORMAT(  ' ISOTOPE SEPARATION FLAG                              : ', L5 )
  401 FORMAT(  ' COLUMN RETRIEVAL SCALE AND VARIANCE        : ',2F10.5)
  403 FORMAT(/,' OFF DIAGNOAL COEFFICIENTS SET TO ZERO')

  600 FORMAT(/,' RETRIEVAL GAS #      ',I2, '                    : ', A7)
  601 FORMAT(  ' PROFILE RETRIEVAL CODE                     : ',L5 )

!  605 FORMAT(' ABORT -- NUMBER OF RETRIEVAL GASES EXCEEDS',I2)
  606 FORMAT(' ABORT -- NUMBER OF PROFILE RETRIEVALS (NPGAS=',I2,&
         ') EXCEEDS MAXIMUM (MAXPRF=',I2,')')
  610 FORMAT(' READCK1: RETRIEVAL GAS : ', A7, ' NOT IN INPUT LIST *** ABORT')
  611 FORMAT(' COLUMN SCALE FACTOR: ', F10.5)
  612 FORMAT(6F12.4)
  613 FORMAT(' RELATIVE UNCERTAINTIES OF THE A PRIORI PROFILE')

  614 FORMAT(' HALF WIDTH HALF HEIGHT (KM) OF GAUSSIAN INTERLAYER CORRELATION :',ES11.4,/,&
             ' MINIMUM, MAXIMUM ALTITUDE (KM) FOR OFF-DIAGONAL ELEMENTS       : ',2F10.3 )

  615 FORMAT(' HALF WIDTH HALF HEIGHT (KM) OF EXPONENTIAL INTERLAYER CORRELATION : ',ES11.4,/,&
             ' MINIMUM, MAXIMUM ALTITUDE (KM) FOR OFF-DIAGONAL ELEMENTS          : ',2F10.3 )

  617 FORMAT( " READING IN",I3," x",I3," COVARIANCE MATRIX FROM FILE : ", A )
  618 FORMAT( " READING IN",I3," x",I3," INVERSE COVARIANCE MATRIX FROM FILE : ", A )
  619 FORMAT( " ILOGRETRIEVAL FLAG : ", I2)
  620 FORMAT(/,' HALF WIDTH OF INTEGRATION INTERVAL(CM-1)   : ', F10.7 )
 ! 622 FORMAT(  ' LINESHAPE MODEL                          : ', I5, /, &
 !              ' 1-VOIGT, 2-GALATRY, 0-GALATRY IF B0 EXISTS' )
  630 FORMAT(/,'NO GASES BEING RETRIEVED.')
  650 FORMAT(/,' MAXIMUM NUMBER OF ITERATIONS               : ', I5)
  651 FORMAT(' CONVERGENCE VARIABLE MUST BE GREATER THEN 0')
      RETURN

      END SUBROUTINE READCK1


      SUBROUTINE READCK2( CPNAM )

      CHARACTER(LEN=14), DIMENSION(5) :: CPNAM

! --- TEMPERATURE RETRIEVAL
      WRITE(16,120) IFTEMP
      IF( IFTEMP )THEN
         WRITE(16,121)
         WRITE(16,612) TSIGMA(1:NLAYERS)
      ENDIF

! --- FORWARD MODEL PARAMETERS
      IF( .NOT. F_EAPOD ) IEAP = 0
      IF( .NOT. F_EPHASE ) IEPHS = 0
      WRITE(16, 101)
      WRITE(16, 102) IFCO, FPS, F_EAPOD, IEAP, NEAP, F_EPHASE, IEPHS, NEPHS, IEMISSION

      IF( IEMISSION /= 0 )THEN
         WRITE(16,103)
         WRITE(16,104) EMISSION_T_BACK, EMISSION_OBJECT, EMISSION_NORM
      END IF

! --- RETRIEVAL SWITCHES
! --- DEFEAT ISPARM IF F_WSHIFT IS NOT SET
      IF( .NOT. F_WSHIFT ) ISPARM = 0
      WRITE (16, 105)
      WRITE (16, 106) F_RTSOL(4), F_WSHIFT, ISPARM, F_BACKG, NBACK, IFDIFF, IFPHASE, F_RTAPOD, &
                      F_RTPHASE, F_LM

! --- INITIAL SCALES AND VARIANCES FOR FITTED PARAMETERS
      WRITE (16, 109)
      WRITE (16, 110) WSHFT, SWSHFT, BCKSL, SBCKSL, BCKCRV, SBCKCRV, CIPARM(4), &
                        SCPARM(4), PHS, SPHS, SZERO(1), EAPPAR, SEAPPAR, EPHSPAR, SEPHSPAR

      IF( F_LM )THEN
         WRITE(16,107)
         WRITE(16,108) GAMMA_START, GAMMA_DEC, GAMMA_INC, CONVERGENCE
      END IF

! --- SOLAR SPECTRUM PARAMETERS
      IF( IFCO )THEN
! --- DEFINE NAMES OF SOLAR PARAMETERS
         CPNAM(1) = 'Sol - n/a'
         CPNAM(2) = 'Sol - n/a'
         CPNAM(3) = 'Sol - n/a'
         CPNAM(4) = 'SolLnShft'
         CPNAM(5) = 'SolLnStrn'
! --- ADD ONE TO WAVENUMBER SHIFT PARAMETER TO AVOID THE CASE OF ZERO
! --- INITIAL SHIFT
!         CIPARM(:) = CIPARM(:) + 1.D0
!         CPARM(:)  = CIPARM(:)
      ENDIF

! --- PRINT OUT GAS FILES
      IF( F_WRTGASSPC )THEN
         SELECT CASE (GASOUTTYPE)
         CASE (1)
            WRITE(16,130) ' WRITE OUT FINAL GAS FILES ONLY'
         CASE (2)
            WRITE(16,130) ' WRITE OUT GAS FILES FOR ALL ITERATIONS'
         CASE DEFAULT
            WRITE(16,130) ' PARAMETER OUTPUT.WRT_GASFILES.TYPE OUT OF RANGE (1 || 2 ONLY)'
            WRITE(16,130) ' PARAMETER OUTPUT.WRT_GASFILES.TYPE OUT OF RANGE (1 || 2 ONLY)'
            CALL SHUTDOWN
            STOP '2'
         END SELECT
      ENDIF

      RETURN

 101  FORMAT(/,' FORWARD MODEL SWITCHES:')
 102  FORMAT( '  INCLUDE SOLAR LINES                       : ', L5, /, &
              '  INCLUDE PRESSURE SHIFT                    : ', L5, /, &
              '  EFFECTIVE MODULATION FUNCTION TYPE        : ', L5, I5, '   # TERMS : ', I5, /, &
              '  EFFECTIVE PHASE FUNCTION TYPE             : ', L5, I5, '   # TERMS : ', I5, /, &
              '  COMPUTE EMISSION COMPONENT                : ', I5 )

 103  FORMAT(/,' EMISSION PARAMETERS:')
 104  FORMAT( '  BACKGROUND TEMPERATURE                    : ', F12.4, / &
              '  SUN REFLECTED BY                          : ', A5, /, &
              '  NORMALIZATION                             : ', L5)


 105  FORMAT(/,' RETRIEVAL SWITCHES: ')
 106  FORMAT( '  FIT SOLAR SHIFT                           : ', L5, /, &
              '  FIT WAVENUMBER SHIFT                      : ', L5, '   TYPE       : ', I5, /, &
              '  FIT BACKGROUND                            : ', L5, '   TYPE       : ', I5, /, &
              '  FIT DIFFERENT SHIFT BY GAS                : ', L5, /, &
              '  FIT SIMPLE PHASE CORRECTION               : ', L5, /, &
              '  FIT MODULATION FUNCTION                   : ', L5, /, &
              '  FIT PHASE FUNCTION                        : ', L5, /, &
              '  USE LEVENBERG-MARQUARDT                   : ', L5)


 107  FORMAT(/,' LEVENBERG-MARQUARDT PARAMETERS:')
 108  FORMAT( '  GAMMA START VALUE                         : ', F12.3, /, &
              '  DECREASE BY                               : ', F12.6, /, &
              '  INCREASE BY                               : ', F12.6, /, &
              '  CONVERGENCE CRITERION                     : ', F12.6 )

 109  FORMAT(/,' A PRIORI AND UNCERTAINTIES:')
 110  FORMAT( '  WAVENUMBER SHIFT AND VARIANCE             : ', 2F12.7,/, &
              '  BACKGROUND SLOPE AND VARIANCE             : ', 2F12.7,/, &
              '  BACKGROUND CURVATURE AND VARIANCE         : ', 2F12.7,/, &
              '  SOLAR SHIFT AND VARIANCE                  : ', 2F12.7,/, &
              '  SIMPLE PHASE AND VARIANCE                 : ', 2F12.7,/, &
              '  ZERO LEVEL VARIANCE                       : ',  F12.7,/, &
              '  MODULATION FUNCTION SCALE AND VARIANCE    : ', 2F12.7,/, &
              '  PHASE FUNCTION SCALE AND VARIANCE         : ', 2F12.7 )


! 111  FORMAT(/' INITIAL SOLAR WAVENUMBER SHIFT            : ', F12.7)
! 112  FORMAT( ' FIT SOLAR SHIFT FLAG                      : ', L5)

 120  FORMAT(/,' TEMPERATURE RETRIEVAL SWITCH               : ', L5)
 121  FORMAT( ' TEMPERATURE RELATIVE UNCERTAINTIES         :')
 130  FORMAT(/,A)
 612  FORMAT(6F12.4)

      END SUBROUTINE READCK2


      SUBROUTINE READCK3( )

      INTEGER :: I, J, K, L, N, TBCK
      LOGICAL, DIMENSION(NRET) :: INBAND

      DO I=1, NRET
         INBAND(I) = .FALSE.
      ENDDO

      WRITE(16, 100) NBAND

      TBCK = 0
      DO I = 1, NBAND

! --- CONVERT FOV DIAMETER FROM MILLIRADIANS TO SOLID ANGLE SAVE FOVDIA FOR SOLAR
         FOVDIA(I) = OMEGA(I)
         OMEGA(I) = 2.0D0*PI*(1.D0 - COS(1.D-03*OMEGA(I)/2.D0))

         IF( WAVE3(I) .GE. WAVE4(I) )THEN
            WRITE(16,*) 'SFIT4.CTRL: BANDPASS LIMITS OUT OF ORDER'
            WRITE(00,*) 'SFIT4.CTRL: BANDPASS LIMITS OUT OF ORDER'
            CALL SHUTDOWN
            STOP '2'
         ENDIF

         WRITE (16, 101) I
         WRITE (16, 102) WAVE3(I), WAVE4(I), ZSHIFT(I,1), IZERO(I), NRETB(I)
         WRITE (16, 113) OMEGA(I), FOVDIA(I)

! --- CHECK F_ZSHIFT SWITCH AND DEFEAT IZERO IF NECESSARY
         IF( .NOT. F_ZSHIFT(I) ) IZERO(I) = 0
         IF( F_ZSHIFT(I) .AND. IZERO(I) .EQ. 1 ) NKZERO = I

! --- CHECK GASES TO RETRIEVE IN BAND
         K = NRETB(I)
         IF( K .EQ. 0 )then
            WRITE(16, *)' -- NO GASES TO RETRIEVE IN THIS BAND --'
            GOTO 44
         ENDIF
         DO J = 1, K
            IF( J .EQ. 1 )THEN
               WRITE (16, 103) GASB(I,J)
            ELSE
               WRITE (16, 104) GASB(I,J)
            ENDIF
            DO N = 1, NRET
               !print *, i, j, n, k, gasb(i,j)
               IF (GASB(I,J) == NAME(IGAS(N))) GO TO 43
            END DO
            WRITE (16, 105) TRIM(GASB(I,J)), WAVE3(I), WAVE4(I)
            WRITE (*, 105) TRIM(GASB(I,J)), WAVE3(I), WAVE4(I)
            STOP 2 ! IF list is not in the retrieval list, stop
            goto 50
   43       CONTINUE
            IGASB(I,J) = IGAS(N)
            NGASB(I,J) = N
            NGIDX(N,0,I) = 1
            INBAND(N) = .TRUE.
         END DO
         WRITE(16, *)''

! work in progress no vmr or temp retrieval jwh

 50      continue
! --- SWITCH OFF TRETB IN THIS BAND IF IFTEMP IS OFF
         IF( .NOT. IFTEMP ) TRETB(I) = .FALSE.

! --- CHECK TEMPERATURE RETRIEVAL IN THIS BAND
 44      IF( TRETB(I) )THEN
            WRITE(16,110)
            TBCK = 1
         ENDIF

! --- CHECK RETRIEVING ANYTHING IN THIS BAND
         IF(( .NOT. TRETB(I) ) .AND. ( K .EQ. 0 ))THEN
            WRITE(16,*) ' NOT RETRIEVING ANY QUANTITY IN THIS BAND'
            WRITE(00,*) ' NOT RETRIEVING ANY QUANTITY IN THIS BAND'
            CALL SHUTDOWN
            STOP '2'
         ENDIF

! --- CHANNEL PARAMETERS IF EXISTS
         K = NBEAM_OF_BAND(I)
         !print*, k
         IF (K > 0) THEN
            WRITE (16, 107) CHANNEL_MODEL_OF_BAND(I), NBEAM_OF_BAND(I)
            DO L=1, NBEAM_OF_BAND(I)
               WRITE(16,111) CCIPARM(I,L,:) !*CHAN_SCALE(I,K,:)
               WRITE(16,112) SCHAN_SCALE(I,L,:) !CCIPARM(I,K,:)*
            END DO
            !CALL READ_CHANNEL_PARMS (I, K, 16)
         ENDIF

      END DO

      IF( TBCK .EQ. 0 .AND. IFTEMP )THEN
         WRITE(16,*) ' IFTEMP SET BUT NOT IN A BAND'
         WRITE(00,*) ' IFTEMP SET BUT NOT IN A BAND'
         CALL SHUTDOWN
         STOP '2'
      ENDIF

      DO I = 2, NBAND
         IF (WAVE3(I) >= WAVE3(1)) CYCLE
         WRITE(16, *) 'MICRO-WINDOWS MUST BE IN ASCENDING WAVENUMBER ORDER'
         WRITE(00, *) 'MICRO-WINDOWS MUST BE IN ASCENDING WAVENUMBER ORDER'
         CALL SHUTDOWN
         STOP '2'
      END DO

      DO J = 1, NRET
         IF( .NOT. INBAND(J) )THEN
            WRITE(16,114) J, NAME(IGAS(J))
            WRITE(0 ,114) J, NAME(IGAS(J))
            CALL SHUTDOWN
            STOP '2'
         ENDIF
      ENDDO

! --- NUMBER OF SPECTRA, DEFAULT SNR, FIT TOLERANCE
      !WRITE (16, 106) NSPEC, SNR, TOL

      RETURN

 100  FORMAT(/,' NUMBER OF BANDS TO FIT   : ', I5 )

 101  FORMAT(/,' BANDPASS           : ',I5)
 102  FORMAT(  ' WAVENUMBER RANGE                     : ', F12.6, ' - ', F12.6, /, &
               ' ZERO LEVEL SHIFT AND SWITCH          : ', F12.6, ', ', I5, /, &
               ' NUMBER OF RETRIEVAL GASES            : ', I5)


 103  FORMAT(' RETRIEVAL GAS                        : ', A7, $)
 104  FORMAT(' ', A7, $)

 105  FORMAT(/,' ABORT--READCK3 : ', A7, ' IS NOT IN LIST OF RETRIEVAL GASES',&
             ' BANDPASS RANGE =',F10.4,' TO',F10.4,' CM-1')

! 106  FORMAT(/,' NUMBER OF SPECTRA                    : ', I5, /, &
!               ' DEFAULT SIGNAL-TO-NOISE              : ', F10.2, /, &
!               ' TOLERENCE CRITERION                  : ', F10.5, / )

 107  FORMAT(' CHANNEL MODEL AND NUMBER OF BEAMS    : ', A5, I5)
! 108  FORMAT(/,"NOT FITTED GASES INCLUDED IN CALCULATION",/," GAS     SCALE")
! 109  FORMAT(A7,F7.4)
 110  FORMAT(' RETRIEVING TEMPERATURE IN THIS BAND')
 111  FORMAT(' APRIORI                              : ', 4F12.7)
 112  FORMAT(' COVARIANCE                           : ', 4F12.7)
 113  FORMAT(' SOLID ANGLE [STR]                    : ', E12.6, /, &
             ' FOV DIAMETER [MR]                    : ', F12.6 )
 114  FORMAT(/,'RETRIEVAL GAS: ', I3, 2X, A7, ' IS NOT RETRIEVED IN A PARTICULAR BAND.')

      END SUBROUTINE READCK3

      END MODULE READIN
