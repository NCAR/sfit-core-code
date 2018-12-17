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

      MODULE CHANNEL

! PWJ 2004

      USE PARAMS
      USE WRITEOUT
      USE SYNSPEC

      IMPLICIT NONE

      INTEGER, DIMENSION(MAXBND) :: NBEAM_OF_BAND ! NUMBER OF BEAMS IN THIS BAND
      CHARACTER (len=2), DIMENSION(MAXBND) :: CHANNEL_MODEL_OF_BAND
                                            ! TWO VALID VALUES:
                                            ! 'IP' -- INTEFEROGRAM PERTURBATION MODE
                                            ! 'PS' -- PHASE-SHIFTED REFLECTED BEAM
      REAL(DOUBLE), DIMENSION(MAXBND,MAX_NUM_OF_BEAMS,4) :: SCHAN_SCALE, CHAN_SCALE
      REAL(DOUBLE), DIMENSION(MAXBND,MAX_NUM_OF_BEAMS,4) :: CCIPARM
                                            ! CHANNEL SPECTRUM VALUES (pwj)
      REAL(DOUBLE), DIMENSION(MAXBND,MAX_NUM_OF_BEAMS,4) :: CHANNEL_CPARM
                                            ! CHANNEL SPECTRUM VALUES (pwj)
      INTEGER, DIMENSION(MAXBND,MAX_NUM_OF_BEAMS,4) :: CHANNEL_IFIX
                                            ! CHANNEL SPECTRUM VALUES (pwj)
      INTEGER :: FIRST_CHANNEL_PARM_NUM     ! THE POSITION OF FIRST CHANNEL
      ! PARM IN STATE VECTOR  (pwj)
      INTEGER :: NCHAN                      ! MP number of beams retrieved

      CONTAINS

!---------------------------------------------------------------------------
      SUBROUTINE READ_CHANNEL_PARMS(IBAND, NUM_OF_BEAMS, UNIT2)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  AUTHOR:  WUJIAN PENG
!  DATE:  AUG. 4 , 2002
!
!  FUNCTION: THIS SUBROUTINE IS USED TO READ CHANNEL PARAMETERS
!          FOR BANDPASS IBAND FROM FILE UNIT1 AND WRITE THEM
!          INTO FILE UNIT2; THE BANDPASS HAS NUM_OF_BEAMS
!            SETS OF PARAMETERS, EACH SET CONTAINS FOUR PARMS:
!        CHANNEL_CPARM(IBAND,IBEAM,K),K=1,4
!          AND THE CORRESPONDING FLAGS:
!        CHANNEL_IFIX(IBAND,IBEAM,K),K=1,4
!
!  INPUT:  IBAND -- THE CURRENT BANDPASS NUMBER
!        NUM_OF_BEAMS -- NUMBER OF BEAMS IN CURRENT BANDPASS
!        UNIT1 -- INPUT FILE UNIT NUMBER
!        UNIT2 -- OUTPUT FILE UNIT NUMBER
!
!  OUTPUT:  CHANNEL_CPARM -- CHANNEL PARAMETERS
!           CCIPARM       -- INITIAL CHANNEL PARAMETERS
!           CHANNEL_IFIX  -- FLAGS INDICATING WHICH CPARM NEEDS
!           TO BE FIT
!
!  NOTE:   VARIABLES : CHANNEL_CPARM, CHANNEL_CIPARM,
!        CHANNEL_IFIX
!        ARE DEFINED IN 'bands.dat'
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! binput_parse

      INTEGER, INTENT(IN) :: IBAND, NUM_OF_BEAMS, UNIT2

      INTEGER :: J !, K

!-- READ THE PARAMETERS---------------------------------------

      !WRITE (UNIT2, 99) IBAND

      DO J = 1, NUM_OF_BEAMS
         !READ (UNIT1, *) (CCIPARM(IBAND,J,K),K=1,4)
         !WRITE (UNIT2, 100)
         WRITE (UNIT2,  9) J
         WRITE (UNIT2, 11) CCIPARM(IBAND,J,1), CHANNEL_IFIX(IBAND,J,1)
         WRITE (UNIT2, 12) CCIPARM(IBAND,J,2), CHANNEL_IFIX(IBAND,J,2)
         WRITE (UNIT2, 13) CCIPARM(IBAND,J,3), CHANNEL_IFIX(IBAND,J,3)
         WRITE (UNIT2, 14) CCIPARM(IBAND,J,4), CHANNEL_IFIX(IBAND,J,4)

         !READ (UNIT1, *) (CHANNEL_IFIX(IBAND,J,K),K=1,4)
         !WRITE (UNIT2, 200)
         !WRITE (UNIT2, 20) (CHANNEL_IFIX(IBAND,J,K),K=1,4)
      END DO

!-- STORE THE INITIAL PARAMETERS-------------------------------

      CHANNEL_CPARM(IBAND,:NUM_OF_BEAMS,:) = CCIPARM(IBAND,:NUM_OF_BEAMS,:)

      RETURN

   9  FORMAT(' CHANNEL BEAM       : ', I5 )
   11 FORMAT(' INITIAL PEAK AMPLITUDE OF BEAM AND FIT SWITCH   : ', F10.5, I5)
   12 FORMAT(' INITIAL CHANNEL PERIOD                          : ', F10.5, I5)
   13 FORMAT(' INITIAL ZERO PHASE REFERENCE WAVENUMBER         : ', F10.5, I5)
   14 FORMAT(' INITIAL CHANGE IN AMPLITUDE PER WAVENUMBER      : ', F10.5, I5)

!   10 FORMAT(4F10.4)
!   20 FORMAT(4I5)
!   99 FORMAT('INPUT CHANNEL PARAMETER FOR BANDPASS: ',I5)
!  100 FORMAT('INPUT PEAK AMPLITUDE OF REFLECTED BEAM RELATIVE ',' TO PARIMARY',&
!         /,' INPUT SEPARATION(CM-1) BETWEEN CHANNEL PEAKS',/,&
!         ' INPUT ZERO CHANNEL PHASE REFERENCE WAVENUMBER',/,&
!         ' INPUT CHANGE IN CHANNEL AMPLITUDE PER WAVENUMBER')
!  200 FORMAT('INPUT CHANNEL_IFIX FOR EACH PARAMETER(1=YES,0=NO)')

      END SUBROUTINE READ_CHANNEL_PARMS


!---------------------------------------------------------------------------
      SUBROUTINE INSERT_CHANNEL_PARMS(NVAR, PARM, PNAME, SPARM)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  AUTHOR:  WUJIAN PENG
!  DATE:  AUG. 4 , 2002
!
!  FUNCTION: THIS SUBROUTINE IS TO INSERT THE CHANNEL PARMS TO BE FIT
!          INTO THE STATE VECTOR PARM(), AND THE POSITION OF FIRST
!          CHANNEL PARAMETER IN PARM() WILL BE STORED.
!
!  INPUT : NVAR -- THE CURRENT LENGTH OF PARMS WHICH HAVE BEEN FILLED
!              INTO THE STATE VECTOR
!          CHAN_SCALE -- SCALING FACTOR FOR CHANNEL PARAMETERS
!          SCHAN_SCALE -- UNCERTAINTY OF THE SCALING FACTER
!
!  OUTPUT: PARM  -- THE STATE VECTOR TO BE FIT
!          PNAME -- NAMES OF THE STATEVECTOR ELEMENT
!          SPARM -- UNCERTAINTIES OF THE STATEVECTOR INITIAL GUESS
!
!  NOTE  :
!            VARIABLES: NBAND,NBEAM_OF_BAND,
!                       CHANNEL_IFIX,  CHANNEL_CPARM,
!             FIRST_CHANNEL_PARM_NUM
!          ARE DEFINED IN 'bands.dat'
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      INTEGER, INTENT(INOUT)      :: NVAR
      REAL(DOUBLE), INTENT(INOUT) :: PARM(NMAX), SPARM(NMAX)
      CHARACTER, INTENT(INOUT)    :: PNAME(NMAX)*(16)

      INTEGER :: IBAND, IBEAM, K

! --- INSERT THE CHANNEL PARMS TO BE FIT INTO STATE VECTOR PARM()--------
! ---  AND KEEP THE START NUMBER IN THE PARM()---------------------------
! --- USE LACK OF SIGMA VALUES TO SIGNIFY A NON-FIT PARAMETER

      NCHAN = 0
      FIRST_CHANNEL_PARM_NUM = NVAR + 1

      DO IBAND = 1, NBAND
         DO IBEAM = 1, NBEAM_OF_BAND(IBAND)
            DO K = 1, 4
               IF (SCHAN_SCALE(IBAND,IBEAM,K) .GT. TINY(SCHAN_SCALE(IBAND,IBEAM,K))) THEN
                  CHANNEL_IFIX(IBAND,IBEAM,K) = 1
               ELSE
                  CHANNEL_IFIX(IBAND,IBEAM,K) = 0
               END IF
               IF (CHANNEL_IFIX(IBAND,IBEAM,K) .EQ. 0) CYCLE
               NVAR = NVAR + 1
               NCHAN = NCHAN + 1
               SELECT CASE (K)
               CASE (1)
                  !PNAME(NVAR) = 'PEAK_AMP'
                  WRITE(PNAME(NVAR), '(A8,I1,"_",I1)') 'PeakAmp_', IBAND, IBEAM
               CASE (2)
                  !PNAME(NVAR) = 'CHAN_SEP'
                  WRITE(PNAME(NVAR), '(A8,I1,"_",I1)') 'ChanSep_', IBAND, IBEAM
               CASE (3)
                  !PNAME(NVAR) = 'ZERO_PH_REF'
                  WRITE(PNAME(NVAR), '(A8,I1,"_",I1)') 'ZroPhsR_', IBAND, IBEAM
               CASE (4)
                  !PNAME(NVAR) = 'DELTA_PEAK_AMP'
                  WRITE(PNAME(NVAR), '(A8,I1,"_",I1)') 'DlPkAmp_' , IBAND, IBEAM
               END SELECT
               PARM(NVAR)  = 1.0D0
               SPARM(NVAR) = SCHAN_SCALE(IBAND,IBEAM,K)
            END DO
         END DO
      END DO
      RETURN

!--  RETURN TO THE CALLER

      END SUBROUTINE INSERT_CHANNEL_PARMS


!---------------------------------------------------------------------------
      SUBROUTINE RETRIEVE_CHANNEL_PARMS(PARM)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  AUTHOR:  WUJIAN PENG
!  DATE:  AUG. 4 , 2002
!
!  FUNCTION: THIS SUBROUTINE IS TO RETRIEVE THE CHANNEL PARMS TO BE
!          FITTED FROM THE STATE VECTOR PARM(), THE FIRST CHANNEL
!          PARM IN STATE VECTOR PARM() IS INDICATED BY THE VARIABLE
!          'FIRST_CHANNEL_PARM_NUM'
!
!  INPUT:    PARM -- THE STATE VECTOR TO BE FIT
!
!  OUTPUT:   PART OF CHANNEL_CPARM(,,)
!
!  NOTE  :
!            VARIABLES: NBAND,NBEAM_OF_BAND,
!                       CHANNEL_IFIX,  CHANNEL_CPARM, CHANNEL_CIPARM
!       FIRST_CHANNEL_PARM_NUM
!          ARE DEFINED IN 'bands.dat'
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL(DOUBLE) :: PARM(NMAX)

      INTEGER :: IBAND, IBEAM, K, NVAR

!-- RETRIEVE THE CHANNEL PARMS TO BE FIT FROM STATE VECTOR PARM()--------

      NVAR = FIRST_CHANNEL_PARM_NUM - 1

      DO IBAND = 1, NBAND
         DO IBEAM = 1, NBEAM_OF_BAND(IBAND)
            DO K = 1, 4
               IF(CHANNEL_IFIX(IBAND,IBEAM,K) == 0)THEN
                  CHANNEL_CPARM(IBAND,IBEAM,K) = CCIPARM(IBAND,IBEAM,K)
               ELSE
                  NVAR = NVAR + 1
                  CHANNEL_CPARM(IBAND,IBEAM,K) = PARM(NVAR)*CCIPARM(IBAND,IBEAM,K)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      RETURN

      END SUBROUTINE RETRIEVE_CHANNEL_PARMS


!---------------------------------------------------------------------------
      REAL(KIND(0.0D0)) FUNCTION FIT_CHANNEL_PARMS (IBAND, JSCAN, J, BKGND, TCALI, TEMP)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  AUTHOR:  WUJIAN PENG
!  DATE:  AUG. 4 , 2002
!
!  FUNCTION:  THIS FUNCTION IS USED TO GET THE UPDATED SPECTRUM AFTER
!           THE CHANNEL PARMS ARE USED. IT USES TWO CHANNEL MODELS:
!                      IP -- INTERFEROGRAM PERTURBATION MODEL
!                      PS -- PHASE-SHIFTED REFLECTING MODEL
!
!  INPUT:
!      IBAND   -- NUMBER OF CURRENT BANDPASS
!      JSCAN   -- NUMBER OF CURRENT SCAN
!      J       -- INDEX FOR CURRENT MON POINT
!      BKGND   -- THE BACKGROUND FITTING VALUES FOR CURRENT BANDPASS
!      TCALI   -- COMPUTED COMPLEX SPECTRUM BY SUBROUTINE FSPEC2
!
!  OUTPUT:
!              -- THE UPDATED SPECTRUM AFTER USING CHANNEL PARM
!
!  NOTE:  VARIABLES : ZSHIFT, CHANNEL_MODEL_OF_BAND
!         ARE DEFINED IN 'bands.dat'
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      INTEGER, INTENT(IN)               :: IBAND, JSCAN, J
      REAL(DOUBLE), INTENT(IN)          :: BKGND
      REAL(DOUBLE), INTENT(OUT)         :: TEMP
      COMPLEX(DBLE_COMPLEX), INTENT(IN) :: TCALI

      REAL(DOUBLE)          :: THETA, REALPT, COMPX, PSUM, YCTEMP
      CHARACTER(len=2)      :: CH_MOD
      COMPLEX(DBLE_COMPLEX) :: CMOD, CHANNEL_VALUE

      FIT_CHANNEL_PARMS = 0.0D0

!-- GET CHANNEL MODEL OF CURRENT BAND AND
!-- CALL SUBROUTINE TO OBTAIN : REALPT, COMPX, PSUM

      CH_MOD = CHANNEL_MODEL_OF_BAND(IBAND)
      CALL CALCULATE_PARMS (IBAND, J, REALPT, COMPX, PSUM)

!-- FITTING THE CHANNEL SPECTRUM

      IF (CH_MOD=='PS' .OR. CH_MOD=='Ps' .OR. CH_MOD=='pS' .OR. CH_MOD=='ps') &
         THEN                                    !PROCCESS WITH PHASE-SHITED
           !  REFLECTING MODEL
         THETA = ATAN(COMPX/REALPT)
         CMOD = DCMPLX(SQRT(REALPT*REALPT + COMPX*COMPX),0.D0)
         CHANNEL_VALUE = CMOD*DCMPLX(COS(THETA),SIN(THETA))

         YCTEMP = BKGND*(DBLE(TCALI*CHANNEL_VALUE) + ZSHIFT(IBAND,JSCAN))
         FIT_CHANNEL_PARMS = YCTEMP

         IF( F_WRTCHANNEL )THEN
            TEMP = REAL( CHANNEL_VALUE )         !FOR_CHANNEL_VALUE_FILES
         ENDIF

         GO TO 10

      ENDIF

      IF (CH_MOD=='IP' .OR. CH_MOD=='Ip' .OR. CH_MOD=='ip' .OR. CH_MOD=='iP') &
         THEN                                    !PROCESS WITH INTERFEROGRAM-
                !   PERTURBATION MODEL
         YCTEMP = BKGND*(DBLE(TCALI) + ZSHIFT(IBAND,JSCAN)+PSUM)
         FIT_CHANNEL_PARMS = YCTEMP
         IF( F_WRTCHANNEL )TEMP = PSUM   !for_channel_value_files
      ELSE
         WRITE(16, *) ' CHANNEL:FIT_CHANNEL_PARMS: UNKNOWN CHANNEL SPECTRUM MODEL : ', CH_MOD
         WRITE( 0, *) ' CHANNEL:FIT_CHANNEL_PARMS: UNKNOWN CHANNEL SPECTRUM MODEL : ', CH_MOD
         CALL SHUTDOWN
         STOP 1
      ENDIF

   10 CONTINUE

      RETURN

      END FUNCTION FIT_CHANNEL_PARMS


!---------------------------------------------------------------------------
      SUBROUTINE CALCULATE_PARMS(IBAND, J, REALPT, COMPX, PSUM)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                C
!  AUTHOR:  WUJIAN PENG                                          C
!  DATE:  AUG. 4 , 2002                                          C
!                                                                C
! FUNCTION: COMPUTE THREE PARAMETERS :      C
!      REALPT, COMPX, PSUM      C
!         OF BANDPASS IBAND AT POINT J FOR FUNCTION FIT_CHANNEL_PARMS C
!             C
! INPUT:        C
!       IBAND  -- NUMBER OF CURRENT SPECTRA    C
!       J      -- INTEGER INDICATING THE CURRENT FITTING POINT  C
!             C
! OUTPUT:        C
!       COMPX -- A SCALAR  USED FOR NIPLE MODEL: PHASE-SHIFTED  C
!                 REFLECTED BEAMS     C
!       PSUM  -- A SCALAR USED FOR INTERFEROGRAM PERTURBATION  C
!                 MODEL       C
!       REALPT-- A SCALAR USED FOR PHASE_SHIFTED REFLECTED MODEL C
!             C
!  NOTE: VARIABLE ' CHANNEL_CPARM ' IS DEFINED IN  'bands.dat'  C
!                 SPAC IS DEFINED IN 'jmix.dat'   C
!             C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      INTEGER, INTENT(IN)       :: IBAND, J
      REAL(DOUBLE), INTENT(OUT) :: REALPT, COMPX, PSUM

      INTEGER      :: LL, NB
      REAL(DOUBLE) :: WAVE, ABEAM, DELAY, COSPART

!-- INITIALIZE LOCAL VARIABLES----------------------------------------

      WAVE   = WSTART(IBAND) + (J - 1)*SPAC(IBAND)
      REALPT = 1.D0
      COMPX  = 0.D0
      PSUM   = 0.D0
      NB     = NBEAM_OF_BAND(IBAND)

!-- COMPUTE OUTPUT VALUES---------------------------------------------

      DO LL = 1, NB
         ABEAM = CHANNEL_CPARM(IBAND,LL,1)*(1.D0 + CHANNEL_CPARM(IBAND,LL,4)*(WAVE-WSTART(IBAND)))
         DELAY = 6.2831853D0*(WAVE - CHANNEL_CPARM(IBAND,LL,3))/CHANNEL_CPARM(IBAND,LL,2)
         COSPART = ABEAM*COS(DELAY)
         REALPT = REALPT - ABEAM + COSPART
         COMPX = COMPX + ABEAM*SIN(DELAY)
         PSUM = PSUM + COSPART
      END DO

      RETURN

      END SUBROUTINE CALCULATE_PARMS


!---------------------------------------------------------------------------
      SUBROUTINE PRINT_CHANNEL_PARMS(FILE_UNIT)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                C
!  AUTHOR:  WUJIAN PENG                                          C
!  DATE:  AUG. 4 , 2002                                          C
!             C
!  FUNCTION: THIS SUBROUTINE IS TO PRINT ALL CHANNEL PARMS AFTER BEING
!          FIT       C
!  INPUT :   FILE_UNIT -- FILE UNIT  NUMBER OF OUTPUT   C
!             C
!             C
!  OUTPUT:   CHANNEL PARMS --CHANNEL_CPARM(,,)    C
!             C
!  NOTE  :   VARIABLE 'CHANNEL_CPARM' IS DEFINED IN 'bands.dat'  C
!             C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      INTEGER :: FILE_UNIT
      INTEGER :: IBAND, IBEAM, K

!-- PRINT CHANNEL PARMS-------------------------------------------------

      DO IBAND = 1, NBAND
         DO IBEAM = 1, NBEAM_OF_BAND(IBAND)
            WRITE (FILE_UNIT, 99) IBAND, IBEAM
            WRITE (FILE_UNIT, 991)
            WRITE (FILE_UNIT, 100) (CHANNEL_CPARM(IBAND,IBEAM,K),K=1,4)
         END DO
      END DO

      RETURN

   99 FORMAT(/,'RETRIEVED CHANNEL PARAMETERS FOR BANDPASS: ',I2,'    BEAM: ',I5)
  991 FORMAT(' PEAK_AMP      CHAN_SEP      ZERO_PH_REF  ',' CHANGE_IN_CHAN_AMP')
  100 FORMAT(4(f12.5,2X))

      END SUBROUTINE PRINT_CHANNEL_PARMS

      END MODULE CHANNEL
