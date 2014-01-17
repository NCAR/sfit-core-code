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



      MODULE RAYTRACE

      USE PARAMS
      USE MOLCPARAM
      USE DATAFILES
      USE RETVPARAM
      USE ISOTOPE
      USE WRITEOUT

      IMPLICIT NONE

!    Constants from NIST 01/11/2002
      REAL (8), PARAMETER    :: AVOGAD = 6.02214199D+23,        &
                                ALOSMT = SCHMIT,                &
                                GASCON = 8.314472D+07,          &
                                PZERO  = BAR,                   &
                                TZERO  = ZEROC,                 &
                                CLIGHT = 2.99792458D+10

      INTEGER (4), PARAMETER :: MXFSC  = 200,                   & ! NOT USED
                                NFINE  = 207,                   & ! FINE GRID
                                MXLAY  = LAYMAX,                & ! MAX FINISH (100)
                                MXZMD  = 200,                   & ! USR MODEL
                                MXPDIM = LAYMAX+NFINE+MXZMD+2,  & ! MAX PATHS
                                IM2    = MXPDIM-2,              &
                                MXMOL  = MOLTOTAL,              &
                                MXTRAC = 22,                    &
                                NXPBAR = MXLAY * (14+MXMOL)+2,  &
                                NXZOUT = MXLAY * 3 + MXMOL*3


      CHARACTER (LEN=8)     :: HMOLS(MXMOL)
      INTEGER               :: JUNIT(MXMOL)
      REAL (8)              :: USRMIX(MXZMD,MXMOL)
      REAL (8)              :: WMOL(MXMOL)
      REAL (8)              :: AMWT(MXMOL)
      REAL (8)              :: ADCON, ALZERO, AVMWT, AIRMWT

      CHARACTER (LEN=8), DIMENSION(3)      :: HMOD
      CHARACTER (LEN=3), DIMENSION(MXLAY)  :: CTYPE

      INTEGER (4)  :: IRD=-1, IPR=73, IPU=78, IRP=72, IZ
      INTEGER (4)  :: NOPRNT=-1, IMMAX, IMDIM, IBMAX, IBDIM, IOUTMX, IOUTDM, IPMAX
      INTEGER (4)  :: IPHMID, IPDIM, KDIM, KMXNOM, IFINMX
      INTEGER (4)  :: NMRGCALL
      INTEGER (4)  :: IMMAX_B, IMLOW
      INTEGER (4)  :: ICH(4)
      INTEGER (4), DIMENSION(MXLAY) :: IPATH, ITYL

      REAL (8)                     :: DEG, GCAIR, DELTAS, ZMIN, ZMAX, SAMPLE, ADBAR
      REAL (8)                     :: HOBS, RE, REF_LAT

      REAL (8), DIMENSION(MXZMD)   :: ZMDL, PM, TM, TM0, RFNDXM, MWGM, DENW, DRYAIR
      REAL (8), DIMENSION(NFINE)   :: ZFINE
      REAL (8), DIMENSION(MXPDIM)  :: ZPTH, PP, TP, RFNDXP, SP, PPSUM, TPSUM, RHOPSM
      REAL (8), DIMENSION(MXLAY)   :: PBAR, TBAR, WN2L, DVL, WTOTL, ALBL, ZFIN
      REAL (8), DIMENSION(MXLAY)   :: ADBL, AVBL, H2OSL, SECNTA, SOUT, RHOSUM

      REAL (8), DIMENSION(NFINE+MXLAY) :: ZOUT
      REAL (8), DIMENSION(0:MXLAY)     :: ALTZ, PZ, TZ
      REAL (8), DIMENSION(MXMOL)       :: AMTTOT, AMTCUM, ISKIP
      REAL (8), DIMENSION(MXLAY)       :: ZBND, PBND, TBND, ALORNZ, ADOPP, AVOIGT, DRAIRL
      REAL (8)                         :: DENM(MXMOL,MXZMD)
      REAL (8), DIMENSION(MXZMD)       :: RELHUM

      REAL (8), DIMENSION(MXMOL,MXPDIM) :: AMTP, DENP
      REAL (8), DIMENSION(MXMOL,MXLAY)  :: AMOUNT


      REAL (8), ALLOCATABLE :: ZSL(:)

! DEFAULT THAT CAN BE PUT IN BINPUT
      INTEGER (4) :: MODEL_RT   = 7       ! OUR REFMOD + ZPT
      INTEGER (4) :: ITYPE_RT   = 2       ! SLANT PATH
      INTEGER (4) :: NOZERO_RT  = 1
      INTEGER (4) :: NOPRNTF_RT = -1
      INTEGER (4) :: IPUNCH_RT  = -1
      INTEGER (4) :: IFXTYP     = 1       ! IFXTYP IS THE FLAG FOR FIXING THE VALUE OF ITYL
      INTEGER (4) :: MUNITS     = 1
      INTEGER (4) :: IASTRO_RT  = 1

! READ IN FROM STATION.LAYERS
      INTEGER (4) :: IBMAX_B_RT = 1          ! MAX BOUNDARIES INC TOP & BOTTOM
      REAL (8)    :: HSPACE_RT  = 120.D0     ! TOA
      REAL (8)    :: H1F_RT     = 0.D0       ! BOTTOM OF BOTTOM LAYER
      REAL (8)    :: H2F_RT     = 120.D0     ! TOA

! READ IN FROM T15ASC
      REAL (8), DIMENSION(MAXSPE) :: REARTH   = 6390.D0
      REAL (8), DIMENSION(MAXSPE) :: XVB      = 2000.D0
      REAL (8), DIMENSION(MAXSPE) :: REFLAT   = 45.D0
      REAL (8), DIMENSION(MAXSPE) :: APPANG   = 0.D0
      REAL (8), DIMENSION(MAXSPE) :: ASTANG   = 0.D0
      REAL (8), DIMENSION(MAXSPE) :: ASTANG0   = 0.D0

      REAL(DOUBLE), DIMENSION(MAXGAS,LAYMAX) :: FXGAS

! WOULD BE FROM T15ASC FOR HORIZONTAL PATHS...
      REAL (8) :: ANGLEF_RT = 0.0D0
      REAL (8) :: RANGEF_RT = 0.0D0
      REAL (8) :: BETAF_RT  = 0.0D0
      INTEGER (4) :: LENF_RT   = 0
      REAL (8) :: HOBS_RT   = 0.0D0

!     IBDIM IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
!     IOUTDM IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
!     IMDIM IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
!         STORED IN ZMDL (INPUT)
!     IPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH OBTAINED
!         BY MERGING ZMDL AND ZOUT
!     KDIM IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT

      DATA KMXNOM / 7 /

!     DELTAS IS THE NOMINAL SLANT PATH INCREMENT IN KM.
      DATA DELTAS / 5.0D0 /


!     ALZERO IS THE MEAN LORENTZ HALFWIDTH AT PZERO AND 296.0 K.
!     AVMWT IS THE MEAN MOLECULAR WEIGHT USED TO AUTOMATICALLY
!     GENERATE THE LBLRTM BOUNDARIES IN AUTLAY
      DATA ALZERO / 0.04 /, AVMWT / 36.0 /

!     MOLECULAR WEIGHTS

! These get over-written
      DATA AIRMWT / 28.964 /
      DATA   AMWT / 18.015 ,  44.010 , 47.998 , 44.01 ,  &
                    28.011 ,  16.043 , 31.999 , 30.01 ,  &
                    64.06  ,  46.01  , 17.03  , 63.01 ,  &
                    17.00  ,  20.01  , 36.46  , 80.92 ,  &
                   127.91  ,  51.45  , 60.08  , 30.03 ,  &
                    52.46  ,  28.014 , 27.03  , 50.49 ,  &
                    34.01  ,  26.03  , 30.07  , 34.00 ,  &
                    66.01  , 146.05  , 34.08  , 46.03 ,  &
                    33.00  ,  15.99  , 98.    , 30.00 ,  &
                    97.    ,  44.5   , 61*0.0 /

      DATA( ZFINE(IZ), IZ=1, NFINE)/ &
            0.0,       0.2,       0.4,       0.6,       0.8, &
            1.0,       1.2,       1.4,       1.6,       1.8, &
            2.0,       2.2,       2.4,       2.6,       2.8, &
            3.0,       3.2,       3.4,       3.6,       3.8, &
            4.0,       4.2,       4.4,       4.6,       4.8, &
            5.0,       5.2,       5.4,       5.6,       5.8, &
            6.0,       6.2,       6.4,       6.6,       6.8, &
            7.0,       7.2,       7.4,       7.6,       7.8, &
            8.0,       8.2,       8.4,       8.6,       8.8, &
            9.0,       9.2,       9.4,       9.6,       9.8, &
           10.0,      10.2,      10.4,      10.6,      10.8, &
           11.0,      11.2,      11.4,      11.6,      11.8, &
           12.0,      12.2,      12.4,      12.6,      12.8, &
           13.0,      13.2,      13.4,      13.6,      13.8, &
           14.0,      14.2,      14.4,      14.6,      14.8, &
           15.0,      15.2,      15.4,      15.6,      15.8, &
           16.0,      16.2,      16.4,      16.6,      16.8, &
           17.0,      17.2,      17.4,      17.6,      17.8, &
           18.0,      18.2,      18.4,      18.6,      18.8, &
           19.0,      19.2,      19.4,      19.6,      19.8, &
           20.0,      20.2,      20.4,      20.6,      20.8, &
           21.0,      21.2,      21.4,      21.6,      21.8, &
           22.0,      22.2,      22.4,      22.6,      22.8, &
           23.0,      23.2,      23.4,      23.6,      23.8, &
           24.0,      24.2,      24.4,      24.6,      24.8, &
           25.0,      25.2,      25.4,      25.6,      25.8, &
           26.0,      26.2,      26.4,      26.6,      26.8, &
           27.0,      27.2,      27.4,      27.6,      27.8, &
           28.0,      28.2,      28.4,      28.6,      28.8, &
           29.0,      29.2,      29.4,      29.6,      29.8, &
           30.0,      30.5,      31.0,      31.5,      32.0, &
           32.5,      33.0,      33.5,      34.0,      34.5, &
           35.0,      35.5,      36.0,      36.5,      37.0, &
           37.5,      38.0,      38.5,      39.0,      39.5, &
           40.0,      41.0,      42.0,      43.0,      44.0, &
           45.0,      46.0,      47.0,      48.0,      49.0, &
           50.0,      51.0,      52.0,      53.0,      54.0, &
           55.0,      56.0,      57.0,      58.0,      59.0, &
           60.0,      62.0,      64.0,      66.0,      68.0, &
           70.0,      72.0,      74.0,      76.0,      78.0, &
           80.0,      85.0,      90.0,      95.0,     100.0, &
          110.0,     120.0 /


      CONTAINS

!----------------------------------------------------------------------

      SUBROUTINE ATMDRV( )

      IMPLICIT NONE

      INTEGER (4) :: NSCAN, FLAG

!     open default files: read, print, punch
      IRD = 55
      OPEN (IRD,FILE=TFILE(71),STATUS='UNKNOWN')

      IF( F_WRTRAYTC )THEN
         IPR = 66
         OPEN (IPR,FILE='TAPE6',STATUS='UNKNOWN')
         IPU = 7
         IRP = 8
         OPEN (IRP,FILE=TFILE(72),STATUS='UNKNOWN')
      ENDIF
      !WRITE (IPR,900)

! --- FIRST CALL
      FLAG = 1
      READ(IRD,*) NSCAN

! --- NSCAN - NUMBER OF SZA'S TO CALCULATE A RAYTRACE
      !CALL LBLATM( FLAG, NSCAN, ASTANG )

      RETURN

!  900 FORMAT ('1',20X,'*****PROGRAM LBLATM*****v4     ',A10,5X,A10,///)

      END SUBROUTINE ATMDRV


!----------------------------------------------------------------------
      SUBROUTINE READLAYRS( NLAY )

! SEE MKLEV.PRO TO CREATE THE LAYERS THIS CODE READS

! JAN 2011
! READ A 'STATION.LAYERS' FILE FOR THE RETRIEVAL ALT GRIDS

      CHARACTER (LEN=80)    :: BUF
      INTEGER               :: I, NAERR, IDUM, NBND, NLAY
      REAL                  :: RDUM

      CALL FILEOPEN( 71, 3 )

      READ(71,'(A80)') BUF
      WRITE(16,201) BUF
  201 FORMAT(/,' STATION.LAYERS FILE : ', A)
      READ(71,*) NBND


! --- THERE IS 1 LESS LAYER THAN LEVELS - THIS IS NLAY !!!
! --- THIS NEEDS TO MATCH EG SIGMA IN BINPUT
! --- THESE ATM ARRAYS ALL GO HIGH ALT TO LOW
! --- ZBAR IS MIDPOINTS
! --- Z IS THE LOWER BOUNDARIES OF EACH LAYER
      NLAY = NBND -1
      ALLOCATE (ZSL(NBND), ZBAR(NLAY), Z(NLAY), P(NLAY), T(NLAY), PMB(NLAY), STAT=NAERR)
      IF (NAERR /= 0) THEN
         WRITE(16, *) 'RAYTRACE: READLAYERS: COULD NOT ALLOCATE Z AND ZBAR ARRAYS ERROR NUMBER = ', NAERR
         WRITE( 0, *) 'RAYTRACE: READLAYERS: COULD NOT ALLOCATE Z AND ZBAR ARRAYS ERROR NUMBER = ', NAERR
         STOP '3'
      ENDIF

      READ(71,*) BUF

! --- TOA IS AT TOP OF FILE
      DO I=1, NLAY
         READ(71,*) IDUM, ZSL(I), RDUM, RDUM, ZBAR(I)
      ENDDO

      READ(71,*) IDUM, ZSL(NBND)

      WRITE(16,202) '  LAYER BOUNDARIES AND MIDPOINTS'
      DO I=1, NLAY
         WRITE(16,203) I, ZSL(I),  ZBAR(I)
      ENDDO

      WRITE(16,204) ZSL(NBND)
      WRITE(16,205) NLAY

! --- ALLOW INPUT OF TYPE 3 - OCCULTATION PATH TO TANGENT POINT POSSIBLY
!     BELOW OBSERVATION ALTITUDE
      H1F_RT   = -1.
      H2F_RT   = ZSL(1)
      ITYPE_RT = 2
      READ(71,*,END=100) ITYPE_RT
      IF( ITYPE_RT .EQ. 3 )H2F_RT = 0.0
      READ(71,*,END=101) H1F_RT
      IF(( H1F_RT .LT. ZSL(NBND)) .OR. (H1F_RT .GT. ZSL(1)))THEN
         WRITE(16, *) 'RAYTRACE: ITYPE 3: H1 OOR'
         WRITE( 0, *) 'RAYTRACE: ITYPE 3: H1 OOR'
         CALL SHUTDOWN
         STOP '3'
      ENDIF
      GOTO 102

 100 CONTINUE
      H1F_RT     = ZSL(NBND)

 102 CONTINUE
      HSPACE_RT  = ZSL(1)
      IBMAX_B_RT = NBND
      !print *, ITYPE_RT, HSPACE_RT, H1F_RT, H2F_RT, IBMAX_B_RT
      CALL FILECLOSE( 71, 2 )

      RETURN

 101  CONTINUE
      WRITE(16, *) 'RAYTRACE: ITYPE 3 H1 READ ERROR.'
      WRITE( 0, *) 'RAYTRACE: ITYPE 3 H1 READ ERROR.'
      CALL SHUTDOWN
      STOP '3'

 202  FORMAT(A)
 203  FORMAT(I5, 2F12.4)
 204  FORMAT(' BOTTOM OF RETRIEVAL GRID : ', F12.4)
 205  FORMAT(' NUMBER OF LAYERS         : ', I12)

END SUBROUTINE READLAYRS

!----------------------------------------------------------------------
!     path:      %P%
!     revision:  $Revision: 9.8 $
!     created:   $Date: 2006/07/21 17:49:43 $
!     presently: %H%  %T%
!
      SUBROUTINE LBLATM( ITER, NLEV)
!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002 - 2006, Atmospheric & Environmental Research, Inc. (AER).|
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
!**********************************************************************
!
!     LBLATM IS AN ATMOSPHERIC RAY TRACE PROGRAM.
!     IT CREATES AND FORMATS THE ATMOSPHERIC INPUTS FOR THE AFGL
!     LINE-BY-LINE TRANSMITTANCE/RADIANCE PROGRAM LBLRTM.
!
!     SEE THE COMMENTS IN SUBROUTINE ATMPTH FOR DETAILED INSTRUCTIONS O
!     THE USAGE OF THE ATMOSPHERIC INPUTS.
!
!     The geometry was modified for LBLRTM to reflect changes
!     implemented in MODTRAN to solve problems with inconsistent
!     path parameters.
!     These changes include changing some variables and functions to
!     double precision.
!
!**********************************************************************
!-
!-                      STATEMENT FLAGS
!-
!-    LBLATM HAS BEEN STRUCTURED TO HAVE ENHANCED PORTABILITY UNDER
!-    FORTRAN 77.  TWO FLAGS (COLUMN73) HAVE BEEN USED TO FACILITATE
!-    PROGRAM CONVERSION.
!-
!-   &    IDENTIFIES STATEMENTS REQUIRED FOR WORD SIZE LESS THAN 8 CHAR
!-               ALL STATEMENTS FLAGGED WITH & IN COLUMN 73 HAVE
!-               STARTING IN COLUMN 1. THESE TWO CHARACTERS MUST
!-               BE CHANGED TO BLANKS FOR COMPUTERS WITH WORD SIZE
!-               LESS THAN 8 CHARACTERS.
!-
!-   !    IDENTIFIES STATEMENTS REQUIRED TO DOUBLE PRECISION THE
!-               VARIABLES NEEDED FOR CALCULATIONS WHICH NEED MORE
!-               THAN 32 BITS TO OBTAIN SUFFICIENT ACCURACY (I.E.
!-               THE FREQUENCIES). STATEMENTS FLAGGED WITH ! HAVE
!-               STARTING IN COLUMN 1. THESE TWO CHARACTERS SHOULD BE
!-               CHANGED TO BLANKS FOR COMPUTERS HAVING SINGLE
!-               PRECISION LESS THAN 10 SIGNIFICANT DIGITS.
!-
!-   >    IDENTIFIES STATEMENTS THAT MAY BE USEFUL FOR CONVERSION,
!-               TYPICALLY SYSTEM SPECIFIC CALLS (I.E. DATE, TIME,
!-               CPU TIME, RANDOM NUMBER, ETC.).
!-
!----------------------------------------------------------------------
!
!     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
!     MXLAY IS THE MAXIMUM NUMBER OF OUTPUT LAYERS
!     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
!         STORED IN ZMDL (INPUT)
!     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
!         OBTAINED BY MERGING ZMDL AND ZOUT
!     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT

      CHARACTER (LEN=24)   :: HDATE
      CHARACTER (LEN=18)   :: HVRATM

      INTEGER (4) :: IREAD, ITER, ITERMX, JSPEC
      INTEGER     :: NLEV
      REAL (4)    :: TSRT, TSTP

!     ASSIGN CVS VERSION NUMBER TO MODULE
      HVRATM = '$Revision: 9.8 $'
      CALL FDATE(HDATE)
      !PRINT *, HDATE

! --- SETUP OUTPUT VERBOSITY LEVELS USE OLD AND NEW FLAGS
! --- SKIP PUNCH FILE ALTOGETHER
      IF( .NOT. F_WRTRAYTC )THEN
         NOPRNT = -1
         RAYOUTTYPE = 0
      ELSEIF( RAYOUTTYPE .GE. 2 )THEN
         NOPRNT = 2
      ENDIF

      IF( NOPRNT .GE. 2 )WRITE (IPR,900) HDATE, HVRATM

      CALL CPU_TIME(TSRT)

      HMOLS(1:MOLTOTAL) = NAME(1:MOLTOTAL)
      AMWT(1:MOLTOTAL)  = XMASS(1,1:MOLTOTAL)

      KDIM     = MXMOL
      IMDIM    = MXZMD
      IOUTDM   = MXLAY
      IPDIM    = MXPDIM
      IBDIM    = MXZMD
      SAMPLE   = 4.0D0
      ITERMX   = 19


! --- FUNCTIONALITY
! ITER = 0 THIS IS THE FIRST CALL FROM SFIT4  -> READ INPUT
! ITER > 1 FOLLOWS OPT ITER ON TEMP RETRIEVAL

! IREAD = 0 - READS INPUT
! IREAD = 0 - WRITES .MS, .MX, .PT, .SA
! IREAD = 1 - TEMPERATURES HAVE CHANGED - RECALCULATE MODEL WITHOUT READ
! IREAD = 2 - SKIP READ, USE EXISTING MODEL (NEW SZA)
! IREAD = 999 FOR ZENITH CALCULATION AND FINAL WRITE TO .MS FILE

! EACH ITERATION IS WRITTEN TO TAPE6

! --- FIRST TIME THROUGH READ IN TAPE5 INPUT AND TAPE8 ATMOSPHERE (IREAD=0)
!      IF( ITER .LT. 0 .OR. ITER .GT. ITERMX )STOP "LBLATM - ITER GT ITERMX"

! --- NSPEC - NUMBER OF SZA'S TO CALCULATE A RAYTRACE
      DO JSPEC=1, NSPEC

! --- IREAD=0 -> READ INPUT TAPE5 & ATMOS MODEL TAPE8

         IREAD = 2
         IF( JSPEC .EQ. 1 .AND. ITER .EQ. 0 ) IREAD = 0
         IF( JSPEC .EQ. 1 .AND. ITER .GT. 0 ) IREAD = 1

         IF( NOPRNT .GE. 2 )WRITE(IPR,901) "LBLATM - ITER, JSPEC, NSPEC, IREAD, SZA ", ITER, JSPEC, NSPEC, IREAD, ASTANG(JSPEC)

         CALL ATMPTH( JSPEC, IREAD, ASTANG(JSPEC), REARTH(JSPEC), XVB(JSPEC), ITER, NLEV )

      ENDDO

! --- DO ONE MORE THAN NSPEC -> SZA=0.0
      IF( NOPRNT .GE. 2 )WRITE(IPR,901) "LBLATM - ITER, JSPEC, NSPEC, IREAD, SZA ", ITER, NSPEC+1, NSPEC, 999, 0.0D0

      CALL ATMPTH( NSPEC+1, 999,  00.0D0, 6390.D0, 2500.D0, ITER, NLEV )

      CALL CPU_TIME(TSTP)

      WRITE(0,905)  "  RAYTRACE PROCESS TIME : ", TSTP-TSRT
      IF( NOPRNT .GE. 2 )WRITE(IPR,905) "  RAYTRACE PROCESS TIME : ", TSTP-TSRT

      CALL FILECLOSE( IPU, 1 )

      RETURN

 900  FORMAT (2X,'*****PROGRAM LBLATM*****  ',A,5X,A)
 901  FORMAT (/,A,4I4,F10.4,/)
 905  FORMAT (A,2F8.4)

      END SUBROUTINE LBLATM

!     ----------------------------------------------------------------

      SUBROUTINE ATMPTH( JSPEC, IREAD, ASTANG1, RE1, XVBAR1, ITER, NLEV )
!
!**********************************************************************
!
!
!
!
!                  ATMPTH   (ATMOSPHERIC PATH)
!
!
!
!
!                            WILLIAM O. GALLERY
!                          + GAIL   P.  ANDERSON
!                            FRANCIS X. KNEIZYS
!                            JAMES   H. CHETWYND JR.
!                            SHEPARD A. CLOUGH
!
!
!                           +(POINT OF CONTACT FOR THIS PROGRAM)
!
!                                      AIR FORCE GEOPHYSICS LAB
!                                      OPTICAL PHYSICS DIVISION
!                                      HANSCOM AFB
!                                      BEDFORD, MA.  01731
!                                      617-861-4774
!
!
!                                      REVISED:   JULY 1990
!
!**********************************************************************
!
!
!     USER INSTRUCTIONS:
!
!     ATMPTH CALCULATES THE DENSITY WEIGHTED MEAN TEMPERATURE AND
!     PRESSURE AND THE INTEGRATED ABSORBER AMOUNTS (IN MOLECULES
!     CM-2) FOR EACH LAYER ALONG A PATH THROUGH A LAYERED
!     ATMOSPHERE, INCLUDING THE EFFECTS OF REFRACTION AND THE  EARTH'S
!     CURVATURE.  ATMPTH IS DESIGNED TO PREPARE THE ATMOSPHERIC INPUTS
!     TO THE PROGRAM LBLRTM WHICH DOES A LINE-BY-LINE CALCULATION OF
!     ATMOSPHERIC TRANSMITTANCE OR RADIANCE AND IS DESCRIBED IN
!     REFERENCE (1).  THE CONTROL CARDS REQUIRED TO RUN ATMPTH ARE
!     DESCRIBED LATER IN THESE COMMENTS.  A DETAILED DESCRIPTION
!     OF THE ALGORITHM USED HERE AND A DISCUSSION OF THE EFFECTS OF
!     THE EARTH'S CURVATURE AND REFRACTION ARE GIVEN IN REFERENCE (2).
!
!     THE DEFINITIONS AND USES OF THE PATH PARAMETERS ITYPE, H1, H2,
!     ANGLE, RANGE, BETA, AND LEN ARE DESCRIBED IN REFERENCE (2) AND
!     ARE THE SAME AS IN REFERENCE (4).
!
!     THERE ARE SIX BUILT IN ATMOSPHERIC PROFILES WHICH DEFINE THE
!     PRESSURE, TEMPERATURE, AND MIXING RATIOS OF THE 28 MOLECULAR
!     SPECIES INCLUDING H2O, CO2, O3, N2O, CO, CH4, AND O2 ON THE AFGL
!     ATMOSPHERIC LINE PARAMETERS COMPILATION AT 50 STANDARD
!     ALTITUDES.  THESE MODEL ATMOSPHERES ARE DESCRIBED IN
!     REFERENCE (3).  THE USER MAY ALSO INPUT AN ATMOSPHERIC
!     PROFILE AS DESCRIBED LATER (SEE ALSO THE COMMENTS IN
!     THE SUBROUTINE NSMDL). TWENTY-0NE ADDITIONAL MIXING RATIO PROFILE
!     FOR SPECIES CORRESPONDING TO THE MOLECULES ON THE AFGL TRACE GAS
!     COMPILATION ARE INCLUDED.
!
!     THE PRINCIPAL OUTPUT CONSISTS OF THE INTEGRATED ABSORBER AMOUNTS
!     FOR A SET OF LAYERS TO BE INPUT TO THE LINE-BY-LINE CALCULATION.
!     THE NUMBER OF THESE LAYERS REPRESENTS A TRADEOFF BETWEEN ACCURACY
!     AND COMPUTATIONAL SPEED OF THE LINE-BY-LINE CALCULATION.  THE
!     USER HAS THE OPTION OF INPUTTING HIS OWN SET OF LAYER BOUNDARIES
!     OR OF LETTING THE SUBROUTINE AUTLAY GENERATE THESE LAYERS
!     AUTOMATICALLY.  IF THE USER INPUTS BOUNDARY ALTITUDES,  THEY NEED
!     NOT FALL ON THE ATMOSPHERIC PROFILE BOUNDARIES OR INCLUDE THE
!     PATH ENDPOINTS. IF AUTOMATIC LAYERING IS SELECTED, THE USER MAY
!     SPECIFY THE MAXIMUM HALFWIDTH RATIO ACROSS A LAYER AND THE
!     MAXIMUM TEMPERATURE DIFFERENCE ACROSS A LAYER.
!
!     IT IS DIFFICULT TO SPECIFY APRIORI THE RELATIONSHIP BETWEEN
!     THE NUMBER OF LAYERS AND THE ACCURACY:  THE ACCURACY DEPENDS UPON
!     SUCH FACTORS AS THE SPECTRAL REGION, THE DISTRIBUTION OF THE
!     MOLECULES OF INTEREST, THE PARTICULAR PATH TAKEN, AND WHETHER
!     TRANSMITTANCE OR RADIANCE IS CALCULATED. THE LAYERING CREATED
!     BY THE DEFAULT VALUES OF AVTRAT (1.5) AND TDIFF1 (5.0 K) AND
!     TDIFF2 (8.0 K) SHOULD BE CONSIDERED A POINT OF DEPARTURE FOR
!     SUBSEQUENT CALCULATIONS. THE USER SHOULD THEN EXPERIMENT WITH
!     DIFFERENT LAYERING UNTIL THE RESULTS ARE CONSISTENT WITH
!     HIS ACCURACY REQUIREMENTS.
!
!     TO SAVE COMPUTER TIME IN LBLRTM, THE LAYER AMOUNTS ARE ZEROED
!     OUT WHEN
!         1.  THE CUMULATIVE AMOUNT FOR THAT LAYER AND ABOVE IS LESS
!             THAN 0.1 PERCENT OF THE TOTAL,
!         AND
!         2.  A. TRANSMITTANCE IS CALCUALTED (IEMIT = 0)
!             OR
!             B. RADIANCE IS CALCULATED (IEMIT = 1) AND THE PATH IS
!                LOOKING UP ( IPATH = 3)
!     O2 IS  NOT CONSIDERED IN THIS SCHEME.  IF THE ABSORBER
!     FOR A LAYER FOR ALL THE MOLECULES (EXCEPT O2) ARE ZEROED
!     OUT, THEN THAT LAYER AND THOSE ABOVE ARE ELIMINATED
!
!     TO CALCULATE THE AMOUNTS FOR THE TRACE GASES (MOLECULES 8 THROUGH
!     31) THE USER MUST INCREASE NMOL ON CARD 3.1.
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!     OUTPUT :
!
!     THE PRINTED OUTPUT IS ON FILE IPR (DEFAULT=6). SELECTING
!     NOPRNT=1 SUPRESSES THE PRINTING OF THE ATMOSPHERIC PROFILES
!     AND THE LAYER-BY-LAYER RESULTS FOR THE REFRACTED PATH.
!     IF IPUNCH = 1, THEN THE LBLRTM INPUT DATA IS ALSO PUT ON FILE
!     IPU (DEFAULT=7) AND CONSISTS OF A SINGLE CARD IMAGE GIVING THE
!     NUMBER OF LAYERS LMAX AND A 70 CHARACTER FIELD DESCRIBING THE
!     PROFILE AND THE PATH, FOLLOWED BY TWO (OR MORE) CARD IMAGES FOR
!     EACH OF THE LMAX LAYERS
!
!        CARD 2.1    IFORM,LMAX,NMOL,SECNT0,HMOD (1X,I1,I3,I5,F10.6,3A8)
!             IFORM  = COLUMN AMOUNT FORMAT FLAG
!             LMAX   = NUMBER OF LBLRTM LAYERS, MAY DIFFER FROM
!                      IBMAX DEPENDING ON THE PATH.
!             NMOL   = NUMBER OF MOLECULES SELECTED
!             SECNT0 = EFFECTIVE SECANT (SCALE FACTOR) FOR THE AMOUNTS
!             HMOD   = 24 CHARACTER FIELD.
!
!        CARD 2.1.1  PBAR(L),TBAR(L),SECNTK(L),ITYL(L),IPATH(L),
!                    ALTZ(L-1),PZ(L-1),TZ(L-1),ALTZ(L),PZ(L),TZ(L)
!                  (E15.7,2F10.4,A3,I2,1X,F7.2,F8.3,F7.2,F7.2,F8.3,F7.2)
!             PBAR   =  AVERAGE PRESSURE (MB)
!             TBAR   =  AVERAGE TEMPERATURE (K)
!             SECNTK = SCALE FACTOR FOR COLUMN AMOUNT (DEFAULT=0)
!             ITYL  : OVERRIDES THE LBLRTM INTERNAL CALCULATION FOR
!                     ITYPE, NORMALLY LEFT BLANK
!             IPATH : IF THE PATH DOES NOT GO THROUGH A TANGENT HEIGHT,
!                         IF H1.LT.H2   IPATH = 3
!                         IF H1.GT.H2   IPATH = 1
!                      IF THE PATH GOES THROUGH A TANGENT HEIGHT, THEN
!                         FOR THE LAYERS FROM THE TANGENT HEIGHT TO
!                         MIN(H1,H2),   IPATH = 2
!                         FOR THE LAYERS (IF ANY) FROM MIN(H1,H2)
!                         TO H1,  IPATH = 1
!                         FOR THE LAYERS (IF ANY) FROM MIN(H1,H2)
!                         TO H2,  IPATH = 3
!                      FOR A HORIZONTAL PATH,  IPATH = 0
!             ALTZ(L)    UPPER BOUNDARY ALTITUDE (CURRENT LAYER)
!             ALTZ(L-1)  LOWER BOUNDARY ALTITUDE (FOR FIRST LAYER ONLY)
!             PZ(L)      PRESSURE AT ALTZ(L), MB
!             PZ(L-1)    PRESSURE AT ATLZ(L-1),  (FOR FIRST LAYER ONLY)
!             TZ(L)      TEMPERATURE AT ALTZ(L), DEGREES K
!             TZ(L-1)    TEMPERATURE AT ALTZ(L-1),(FOR FIRST LAYER ONLY
!
!        CARD 2.1.2           (AMOUNT(K,L),K=1,7),WBROADL(L)
!                             (1P8E15.7)
!        CARD 2.1.3           (AMOUNT(K,L),K=8,NMOL)
!                             (1P8E15.7)
!             AMOUNT(K)   COLUMN DENSITIES FOR THE K'TH MOLECULAR
!                         SPECIES (MOLECULES CM-2)
!             WBROADL(L)  COLUMN DENSITY FOR BROADENING GASES
!                         (MOLECULES CM-2)
!
!        CARDS 2.1 ARE REPEATED UNITL LMAX LAYERS ARE SPECIFIED.
!
!----------------------------------------------------------------------
!
!  REFERENCES:
!
! (1) LBLRTM - A USERS' GUIDE (AVAILABLE FROM S.A. CLOUGH AT
!                    THE ABOVE ADDRESS)
!        SEE ALSO:
!          FASCODE - FAST ATMOSPHERIC SIGNATURE CODE
!          (SPECTRAL TRANSMITTANCE AND RADIANCE)
!                                                       AFGL-TR-78-0081
!
! (2) AIR MASS COMPUTER PROGRAM FOR ATMOSPHERIC TRANSMITTANCE/RADIANCE:
!     FSCATM
!        W. O. GALLERY, F. X. KNEIZYS, AND S. A. CLOUGH
!                                                       AFGL-TR-83-0065
!
! (3) AFGL ATMOSPHERIC CONSTITUENT PROFILES (0-120 KM)
!        G. P. ANDERSON, S. A. CLOUGH, F.X. KNEIZYS, J. H. CHETWYND
!        AND E. P. SHETTLE
!                                                       AFGL-TR-86-0110
!
! (4) ATMOSPHERIC TRANSMITTANCE/RADIANCE:
!     COMPUTER CODE LOWTRAN 5
!                                                       AFGL-TR-80-0067
!
!**********************************************************************

      CHARACTER (LEN=48)  :: CFORM1,CFORM2
      CHARACTER (LEN=8)   :: COTHER
      CHARACTER (LEN=7)   :: PAFORM(2)
      CHARACTER (LEN=4)   :: HT1HRZ, HT2HRZ, HT1SLT, HT2SLT, PZFORM(5)

      INTEGER (4) :: IPASS, IERROR, LMAX, IAMT, IEMIT, M, IREAD, MODEL
      INTEGER (4) :: ITYPE, NOZERO, IPUNCH, IASTRO, LEN, LENF, IB, ITER
      SAVE           ITYPE, NOZERO, IPUNCH, IASTRO
      INTEGER (4) :: ISTART, IP, LIP, K, IM, IFORM, KASTRO, I2, I
      INTEGER (4) :: L, MOL, LTST, NPTST, IBMAX_B, JSPEC
      INTEGER     :: NLEV

      REAL (8)    :: AVTRAT, TDIFF1, TDIFF2, HSPACE, AVRATS, TDIF1S, TDIF2S
      SAVE           HSPACE, MODEL,IBMAX_B
      REAL (8)    :: SECNT0, H1F, H2F, RANGEF, ZH, ANGLEF, BETAF, TH
      REAL (8)    :: RHOBAR, HIP, ZINT, TIP, AMTAIR, WVIP, SUMAMT, RATP, DRAIR
      REAL (8)    :: PPH2O, DENAIR, DRY_AIR, AIRTOT, FAC !DTEMP, RATIO, ZETA
      REAL (8)    :: AIRMAS, SUMN2, SUMRS, PWTD, TWTD, FACTOR, SUMWK
      REAL (8)    :: FRH2O, ALFCOR, OLDDV, PTST, AST, APP, ASTANG1, RE1, XVBAR1

      REAL (8)    :: PH, A, H1, H2, ANGLE, RANGE, BETA, PHI, ALTD1, ALTD2
      REAL (8)    :: XVBAR, HMIN, HMAX, BENDNG, HMID, WTOT, T296, AIRMS1
      SAVE           H1, XVBAR, HMIN

      REAL (8), DIMENSION(MXZMD)  :: DENSAVE
      REAL (8), DIMENSION(2)      :: TTMP, WVTMP, PTMP, ZTMP
      REAL (8), DIMENSION(MXMOL)  :: WMT

      DATA AVRATS / 1.5D0 /,TDIF1S / 5.0D0 /,TDIF2S / 8.0D0 /

      DATA COTHER / 'OTHER   '/
      DATA HT1HRZ / ' AT '/,HT2HRZ / ' KM '/,HT1SLT / ' TO '/,HT2SLT / ' KM '/
      DATA PZFORM / 'F8.6','F8.5','F8.4','F8.3','F8.2'/
      DATA PAFORM / '1PE15.7','  G15.7'/
      DATA CFORM1 / '(1PE15.7,0PF10.2,10X,A3,I2,1X,2(F7.3,F8.3,F7.2))'/
      !DATA CFORM1 / '(1PE15.7,0PF10.2,10X,A3,I2,1X,2F7.3,F8.3,F7.2)'/
      DATA CFORM2 / '(  G15.7,0PF10.2,10X,A3,I2,23X,(F7.3,F8.3,F7.2))'/

      DATA IERROR / 0 /, IPASS / 0 /

      DATA T296 / STDTEMP /

!     IAMT = 1: CALCULATE AMOUNTS, IAMT = 2: DO NOT CALCULATE AMOUNTS
      DATA IAMT / 1 /

!     AIRMS1 IS ONE AIRMASS OR THE TOTAL AMOUNT FOR A VERTICAL PATH
!     FROM GROUND TO SPACE
      DATA AIRMS1 / 2.153D25 /

! --- IREAD = 0 - FIRST TIME THROUGH AND READ INPUT AND ATMOSPHERE
! ---       > 0 - NEXT PASSES NO READ INPUT AND ATMOSPHERE
! ---       = 999 - LAST TIME, SZA = 0.0, CLOSE MS FILE
! --- NOPRNT = 3 - HIGH PRIORITY PRINT
! ---        = 2 - LOWER PRIORITY
! ---        = 1 - LOWER STILL
! ---        = 0 - LOWER STILL

      IF( IREAD .EQ. 0 )THEN
         WRITE(16,*) ' '
         WRITE(16,890) " MAX LAYERS FOR USER MODEL GRID    : ", MXZMD
         WRITE(16,890) " MAX LAYERS FOR FINAL OUTPUT       : ", LAYMAX
         WRITE(16,890) " NUM LAYERS FOR INTERNAL FINE GRID : ", NFINE
         WRITE(16,890) " MAX LAYERS FOR RAYTRACING GRID    : ", MXLAY
         NMRGCALL = 0
      ENDIF

      WRITE(16,891)  "  ITER, JSPEC, ASTRO SOLAR ANGLE   : ", ITER, JSPEC, ASTANG1

      SECNT0 = 1.0
      IEMIT  = 0

! --- FROM LBLATM
!      IF( IREAD .GT. 0 )THEN
!         NOPRNT = 3
!      ENDIF

      DEG = 180.0D0/PI

!     GCAIR IS THE GAS CONSTANT FOR RHO IN MOL CM(-3), P IN MB, AND T IN K
      GCAIR = 1.0D-3*GASCON/AVOGAD

!     ADCON IS THE CONSTANT FOR THE DOPPLER HALFWIDTH
      ADCON = SQRT(2.0D0 * LOG(2.0D0)*GASCON/CLIGHT**2)

!     ZERO OUT COMMON BLOCKS - EACH CALL
      DRAIRL(:)   = 0.0D0 ! DRY AIR ON ZFIN FOR OUTPUT, CALC IN ATMPTH
      WMT(:)      = 0.0D0 ! TOTAL COLUMN
      ZPTH(:)     = 0.0D0 ! ATMOS ON MODEL + BOUNDRARIES + FINE GRIDS
      PP(:)       = 0.0D0
      TP(:)       = 0.0D0
      SP(:)       = 0.0D0 ! RAYTRACED PATHS
      RFNDXP(:)   = 0.0D0 ! REFRACTIVE INDEX ON PATH
      PPSUM(:)    = 0.0D0 ! SUMS ON PATH
      TPSUM(:)    = 0.0D0
      RHOPSM(:)   = 0.0D0
      DENP(:,:)   = 0.0D0
      AMTP(:,:)   = 0.0D0
      PBAR(:)     = 0.0D0 ! FINAL VALUES FOR LAYERS ZFIN
      TBAR(:)     = 0.0D0
      WN2L(:)     = 0.0D0
      RHOSUM(:)   = 0.0D0
      RFNDXM(:)   = 0.0D0 ! REFRACTIVE INDEX ON MODEL
      DRYAIR(:)   = 0.0D0 ! OF MODEL
      ZOUT(:)     = 0.0D0 ! CALC IN AMERGE
      AMOUNT(:,:) = 0.0D0

      IF (IREAD.LE.0) THEN

! --- READ CONTROL CARD 3.1

!         READ (IRD,900) MODEL,ITYPE,IBMAX_B,NOZERO,NOPRNTF,NMOL,IPUNCH,  &
!                        IFXTYP,MUNITS,RE,HSPACE,XVBAR,dumrd,REF_LAT,    &
!                        IASTRO

         NOZERO  = NOZERO_RT
         !NOPRNT  = NOPRNTF_RT
         IPUNCH  = IPUNCH_RT
         HSPACE  = HSPACE_RT
         ITYPE   = ITYPE_RT
         MODEL   = MODEL_RT
         IASTRO  = IASTRO_RT
         IBMAX_B = IBMAX_B_RT
         IBMAX   = ABS(IBMAX_B)
         RE      = RE1
         XVBAR   = XVBAR1

      ENDIF


!PRINT *, F_WRTRAYTC, RAYOUTTYPE, NOPRNT, IPUNCH
!print *, 'itype ', itype
!print *, 'iter ', iter

      IF (NOPRNT.GE.2) THEN
         WRITE (IPR,902)
         WRITE (IPR,904) MODEL,ITYPE,IBMAX,NOZERO,NOPRNT,NMOL,IPUNCH,   &
                         IFXTYP,MUNITS,RE,HSPACE,XVBAR,REF_LAT,IASTRO
      ENDIF

      M = MODEL
      !IF (NMOL.EQ.0)STOP "ATMPTH - NMOL"
      IF (ITYPE.LE.1.OR.ITYPE.GT.3) GO TO 290
      IF (M.LT.0.OR.M.GT.7) GO TO 290     ! now 7 for old refmod
      IF (IBMAX.GT.IBDIM) GO TO 290
      IF (NMOL.GT.KDIM) GO TO 290

! appears that 78/ipu doe not get opened on T ret wen repeated calls to rayt...
!print *, ipunch, ipass

      IF (IPUNCH.GE.1 .AND. IPASS.EQ.0) THEN
! --- TAPE7 ONLY OPENED IF IPU SELECTED AND FIRST TIME THROUGH
         IPASS = 1
         CALL FILEOPEN ( IPU, 2 )
         WRITE (IPU,905) IPASS
      ENDIF

      IF (RE.NE.0.0) GO TO 60
      RE = 6371.23
      IF (M.EQ.1) RE = 6378.39D0
      IF (M.EQ.4.OR.M.EQ.5) RE = 6356.91D0

   60 CONTINUE

      IF (HSPACE.EQ.0.)  HSPACE  = 120.0D0
      IF (XVBAR.LE.0.)   XVBAR   = 2500.0D0
      IF (REF_LAT .EQ.0) REF_LAT = 45.0D0

      IF (NOPRNT.GE.2) THEN
         WRITE (IPR,906)
         WRITE (IPR,904) MODEL,ITYPE,IBMAX,NOZERO,NOPRNT,NMOL,IPUNCH,   &
                         IFXTYP,MUNITS,RE,HSPACE,XVBAR,REF_LAT,IASTRO
      ENDIF

      IF (ITYPE.EQ.1) THEN

! --- NOT SUPPORTED
! --- HORIZONTAL PATH SELECTED -----------------------------------------

         IF (NOPRNT.GE.0) WRITE (IPR,908)

! --- READ IN CONTROL CARD 3.2

         IF (IREAD.LE.0) READ (IRD,910) H1F,RANGEF
         RANGE  = RANGEF
         ZH     = H1F
         H1     = ZH
         H2     = 0.0D0
         H2F    = H2
         ANGLE  = 0.0D0
         ANGLEF = ANGLE
         BETA   = 0.0D0
         BETAF  = BETA
         LEN    = 0
         LENF   = LEN
         IF (NOPRNT.GE.0)  WRITE (IPR,912) ZH,RANGE

! --- SET UP THE ATMOSPHERIC PROFILE

         !CALL MDLATM (ITYPE,M,IREAD,HSPACE,LMAX)
         CALL MDLATM (ITYPE,M,IREAD,HSPACE)

         IF (IMMAX.EQ.1) THEN
            ZH = ZMDL(1)
            H1F = ZH
            H1 = ZH
            PH = PM(1)
            TH = TM(1)
            RHOBAR = ALOSMT*PH*TZERO/(PZERO*TH)
            DO 70 K = 1, NMOL
               DENP(K,1) = DENM(K,1)
   70       CONTINUE
            DENW(1) = DENM(1,1)
            GO TO 110
         ENDIF

! --- INTERPOLATE ATMOSPHERIC PROFILE DENSITIES TO ZH

         DO 80 IM = 2, IMMAX
            IF (ZH.LT.ZMDL(IM)) GO TO 90
   80    CONTINUE
         IM = IMMAX
   90    CONTINUE
         A = (ZH-ZMDL(IM-1))/(ZMDL(IM)-ZMDL(IM-1))
         CALL EXPINT (PH,PM(IM-1),PM(IM),A)
         TH = TM(IM-1)+(TM(IM)-TM(IM-1))*A
         RHOBAR = ALOSMT*PH*TZERO/(PZERO*TH)
         DO 100 K = 1, NMOL
            CALL EXPINT (DENP(K,1),DENM(K,IM-1),DENM(K,IM),A)
  100    CONTINUE

  110    CONTINUE
         IF (NOPRNT.GE.2) THEN
            WRITE (IPR,914) HMOD,ZH,PH,TH,(HMOLS(K),K=1,NMOL)
            WRITE (IPR,916) RHOBAR,(DENP(K,1),K=1,NMOL)
         ENDIF

!        COMPUTE AMOUNTS FOR A HORIZONTAL PATH

         DO 120 K = 1, NMOL
            AMOUNT(K,1) = DENP(K,1)*RANGE*1.0D+5
  120    CONTINUE
         AMTAIR = RHOBAR*RANGE*1.0D+5
         IF (NOPRNT.GE.2) THEN
            WRITE (IPR,918) HMOD,ZH,PH,TH,RANGE,(HMOLS(K),K=1,NMOL)
            WRITE (IPR,920) AMTAIR,(AMOUNT(K,1),K=1,NMOL)
         ENDIF
         IPATH(1) = 0
         LMAX = 1

         SUMAMT = 0.
         DO 130 K = 1, NMOL
            SUMAMT = SUMAMT+AMOUNT(K,1)
  130    CONTINUE
         WN2L(1) = AMTAIR-SUMAMT

         PBAR(1) = PH
         TBAR(1) = TH
         ALTZ(0) = -RANGE
         ZOUT(1) = ZH
         IOUTMX = 1
         IFINMX = 1
         SECNTA(1) = 1.
         ALTZ(1) = ZH

! --- WRITE ATMOSPHERE TO TAPE7 (IN E15.7 FORMAT)

         IF (IPUNCH.GE.1) THEN
           IFORM = 1
           WRITE (IPU,924) IFORM,LMAX,NMOL,SECNT0,HMOD,RANGE,ZH

!           -------------------------------------
!           > WRITE MOLECULAR INFORMATION IN    <
!           >  - MIXING RATIO IF MUNITS IS 1    <
!           >  - COLUMN DENSITY IF MUNITS IS 0  <
!           -------------------------------------

           IF (MUNITS.EQ.1) THEN
              DRAIR =  WN2L(1)
              DO 135 M = 2,NMOL
                 DRAIR = DRAIR + AMOUNT(M,1)
  135         CONTINUE

!             > IF DRAIR IS ZERO, THEN WRITE OUT AMOUNT ONLY    <
!             > (SINCE AMOUNT ZERO => MIXING RATIO ZERO)        <

              IF (DRAIR.EQ.0 .AND. NOPRNT.GE.0) THEN
                 WRITE (IPU,926) PH,TH,IPATH(1),ZH,ZH,                  &
                                 (AMOUNT(K,1),K=1,7),WN2L(1),           &
                                 (AMOUNT(K,1),K=8,NMOL)
              ELSEIF (NOPRNT.GE.0) THEN
                 WRITE (IPU,926) PH,TH,IPATH(1),ZH,ZH,                  &
                                 (AMOUNT(K,1)/DRAIR,K=1,7),WN2L(1),     &
                                 (AMOUNT(K,1)/DRAIR,K=8,NMOL)
              ENDIF
           ELSE

!             TEST TO MAKE SURE THERE ARE NO FRACTIONAL MOLECULAR
!             AMOUNTS WRITTEN OUT (WILL CAUSE PATH TO ASSUME
!             MIXING RATIO)

              DO 137 K=1,NMOL
                 IF (AMOUNT(K,1).LT.1.) THEN
                    IF (NOPRNT.GE.0) WRITE(IPR,1000) K,AMOUNT(K,1)
                    AMOUNT(K,1) = 0.0
                 ENDIF
  137         CONTINUE

              IF (NOPRNT.GE.0)                                          &
                   WRITE (IPU,926) PH,TH,IPATH(1),ZH,ZH,                &
                                  (AMOUNT(K,1),K=1,7),WN2L(1),          &
                                  (AMOUNT(K,1),K=8,NMOL)
           ENDIF

        ENDIF

      ELSE


! --- SLANT PATH SELECTED-------------------------------------------------------
! --- ITYPE = 2 OR 3: SLANT PATH THROUGH THE ATMOSPHERE

         IF( NOPRNT.GE.0 )WRITE (IPR,930) ITYPE

! --- READ IN CONTROL CARD 3.2 CONTAINING SLANT PATH PARAMETERS

!         IF (IREAD.LE.0) READ (IRD,932) H1F,H2F,ANGLEF,RANGEF,BETAF,LENF,HOBS

! --- SOME ARE DEFAULTS NOW

         H1F    = H1F_RT
         H2F    = H2F_RT
         ANGLEF = ANGLEF_RT
         RANGEF = RANGEF_RT
         BETAF  = BETAF_RT
         LENF   = LENF_RT
         HOBS   = HOBS_RT

         H1     = H1F
         H2     = H2F
         ANGLE  = ANGLEF
! --- NEW TAPE5 - ANGLES ARE IN A LIST AT THE TOP
         ANGLE  = ASTANG1
         RANGE  = RANGEF
         BETA   = BETAF
         LEN    = LENF

! --- ACCOUNT FOR ITYPE=3 AND FINAL RUN FOR VERTICAL MASS PATHS
! --- THIS FORCES ITYPE 2
         IF( IREAD .EQ. 999 )THEN
            H2 = ZMAX
            IF( H1 .GT. HMIN ) H1 = HMIN
         ENDIF

         IF (NOPRNT.GE.0) THEN
            IF (IBMAX_B .LT. 0) THEN
               WRITE (IPR,933) H1,H2,ANGLE,RANGE,BETA,LEN
            ELSE
               WRITE (IPR,934) H1,H2,ANGLE,RANGE,BETA,LEN
            ENDIF
         ENDIF
         !print*,' atmpth top ', h1, h2, angle
! --- GENERATE OR READ IN LBLRTM BOUNDARY LAYERS

         IF (IBMAX.EQ.0) THEN

! --- SELECT AUTOMATIC LAYERING

            IF (IREAD.LE.0) THEN
               READ (IRD,936) AVTRAT,TDIFF1,TDIFF2,ALTD1,ALTD2
               IF (AVTRAT.EQ.0.0) AVTRAT = AVRATS
               IF (TDIFF1.EQ.0.0) TDIFF1 = TDIF1S
               IF (TDIFF2.EQ.0.0) TDIFF2 = TDIF2S
               IF ((ALTD2.LE.0).OR.(ALTD2.LE.ALTD1)) THEN
                  ALTD1 = 0.
                  ALTD2 = 100.
               ENDIF
               IF (NOPRNT.GE.0)                                         &
                    WRITE (IPR,938) AVTRAT,TDIFF1,TDIFF2,ALTD1,ALTD2
               IF (AVTRAT.LE.1.0.OR.TDIFF1.LE.0.0.OR.TDIFF2.LE.0.0)     &
                    GO TO 320
            ENDIF
            GO TO 150

         ENDIF

! --- READ IN LBLRTM BOUNDARY LAYERS

         IF (IREAD.LE.0) THEN
            IF (IBMAX_B .LT. 0) THEN
               READ (IRD,940) (PBND(IB),IB=1,IBMAX)
               IF (NOPRNT.GE.2 )WRITE (IPR,943) (IB,PBND(IB),IB=1,IBMAX)
            ELSE

               !READ (IRD,940) (ZBND(IB),IB=1,IBMAX)

! --- ZBND IS READ IN READLAYRS ALREADY AS Z

!print *, ZSL(1:IBMAX)
!print *, zbnd(1:IBMAX)


               DO I=1, IBMAX
                  ZBND(I) = ZSL(IBMAX-I+1)
               ENDDO


!print *, zbnd(1:IBMAX)

               IF (NOPRNT.GE.2 )WRITE (IPR,942) (IB,ZBND(IB),IB=1,IBMAX)
            ENDIF
         ENDIF

         IF (IBMAX.EQ.0) GO TO 150

         IF (IBMAX_B .GT. 0) THEN
            DO 140 IB = 2, IBMAX
               IF (ZBND(IB).LE.ZBND(IB-1)) GO TO 300
  140       CONTINUE
         ENDIF

         IF (IBMAX_B .LT. 0) THEN
           DO 145 IB = 2,IBMAX
               IF (PBND(IB) .GE. PBND(IB-1)) GO TO 305
  145        CONTINUE
         ENDIF

  150   CONTINUE

! --- THAT WAS END OF TAPE5 READ

! --- SET UP ATMOSPHERIC PROFILE

         IF( IREAD .EQ. 0 .OR. IREAD .EQ. 1)CALL MDLATM (ITYPE,M,IREAD,HSPACE)

! --- INTERPOLATE PBND GRID ONTO ZBND GRID.

! --- TO ENSURE THAT CALCULATED/INPUT ZMDL'S WILL MATCH CALCULATED USER-LEVEL
! --- ALTITUDES, A COMBINATION OF INTERPOLATION AND HYDROSTATICS ARE USED.
! --- ZBND = A * F1(P) + (1 - A) * F2(P), WHERE
! --- F1(P) = INTERPOLATION IN LN(P), F2(P) = HYDROSTATIC CALCULATION

         IF (IBMAX_B .LT. 0) THEN

! --- PRESSURE BOUNDRARIES

            ISTART = 2

            DO 160 IP=1,IBMAX
               PTMP(1)   = 0.0
               TTMP(1)   = 0.0
               WVTMP(1)  = 0.0
               ZTMP(1)   = 0.0
               PTMP(2)   = 0.0
               TTMP(2)   = 0.0
               WVTMP(2)  = 0.0
               ZTMP(2)   = 0.0

               DO 161 LIP=ISTART,IMMAX
                  IF (PBND(IP) .GT. PM(LIP)) GO TO 162
  161             CONTINUE
                  LIP=IMMAX
  162             CONTINUE
                  IF (PBND(IP) .EQ. PM(LIP-1)) THEN
                     ZBND(IP) = ZMDL(LIP-1)
                     TBND(IP) = TM(LIP-1)
                  ELSE
                     IF(PBND(IP) .EQ. PM(LIP)) THEN
                        ZBND(IP) = ZMDL(LIP)
                        TBND(IP) = TM(LIP)
                     ELSE

! PERFORM INTERPOLATION IN LN(PM)
                        HIP =  (ZMDL(LIP)-ZMDL(LIP-1))/ LOG(PM(LIP)/PM(LIP-1))
                        ZINT = ZMDL(LIP-1)+ HIP* LOG(PBND(IP)/PM(LIP-1))

! PERFORM ALTITUDE CALCULATION USING HYDROSTATIC EQUATION
                        PTMP(1) = PM(LIP-1)
                        ZTMP(1) = ZMDL(LIP-1)
                        TTMP(1) = TM(LIP-1)
                        WVTMP(1) = DENW(LIP-1)

                        PTMP(2) = PBND(IP)

                        TIP = (TM(LIP)-TM(LIP-1))/LOG(PM(LIP)/PM(LIP-1))
                        TTMP(2) = TM(LIP-1)+ TIP* LOG(PBND(IP)/PM(LIP-1))

                        WVIP =  (DENW(LIP)-DENW(LIP-1))/ LOG(PM(LIP)/PM(LIP-1))
                        WVTMP(2) =  DENW(LIP-1) + WVIP* LOG(PBND(IP)/PM(LIP-1))
                        CALL CMPALT(2,PTMP,TTMP,WVTMP,ZTMP(1),ZTMP)

! COMBINE THE INTERPOLATION AND THE HYDROSTATIC CALCULATION

                        RATP =  LOG(PBND(IP)/PM(LIP-1))/                &
     &                        LOG(PM(LIP)/PM(LIP-1))

                        A = RATP**3

                        ZBND(IP) = A*ZINT + (1-A)*ZTMP(2)
                        TBND(IP) = TTMP(2)
                     ENDIF
               ENDIF

               ISTART = LIP

  160          CONTINUE

! INTERPOLATE H1, H2 ONTO ALTITUDE GRID
            PTMP(1)   = 0.0
            TTMP(1)   = 0.0
            WVTMP(1)  = 0.0
            ZTMP(1)   = 0.0
            PTMP(2)   = 0.0
            TTMP(2)   = 0.0
            WVTMP(2)  = 0.0
            ZTMP(2)   = 0.0

            DO 166 LIP = 2,IMMAX
               IF (H1 .GT. PM(LIP)) GO TO 167
  166          CONTINUE
               LIP = IMMAX
  167       CONTINUE
            IF (H1 .EQ. PM(LIP-1)) THEN
               H1 = ZMDL(LIP-1)
            ELSE
               IF(H1 .EQ. PM(LIP)) THEN
                  H1 = ZMDL(LIP)
               ELSE

! PERFORM INTERPOLATION IN LN(PM)
                  HIP  = (ZMDL(LIP) - ZMDL(LIP-1))/LOG(PM(LIP)/PM(LIP-1))
                  ZINT = ZMDL(LIP-1) + HIP* LOG(H1/PM(LIP-1))

! PERFORM ALTITUDE CALCULATION USING HYDROSTATIC EQUATION
                  PTMP(1)  = PM(LIP-1)
                  ZTMP(1)  = ZMDL(LIP-1)
                  TTMP(1)  = TM(LIP-1)
                  WVTMP(1) = DENW(LIP-1)

                  PTMP(2) = H1

                  TIP     = (TM(LIP)-TM(LIP-1))/LOG(PM(LIP)/PM(LIP-1))
                  TTMP(2) = TM(LIP-1)+TIP* LOG(H1/PM(LIP-1))

                  WVIP     =  (DENW(LIP)-DENW(LIP-1))/LOG(PM(LIP)/PM(LIP-1))
                  WVTMP(2) =  DENW(LIP-1) + WVIP* LOG(H1/PM(LIP-1))

                  CALL CMPALT(2,PTMP,TTMP, WVTMP,ZTMP(1),ZTMP)

! COMBINE THE INTERPOLATION AND THE HYDROSTATIC CALCULATION

                  RATP =  LOG(H1/PM(LIP-1))/ LOG(PM(LIP)/PM(LIP-1))
                  A    = RATP**3
                  H1   = A*ZINT + (1-A)*ZTMP(2)
               ENDIF
            ENDIF

            IF (H1 .LT. 0.0) THEN
               PRINT 946, H1,ZTMP(1)
               IF (NOPRNT.GE.0) WRITE (IPR,946) H1,ZTMP(1)
               WRITE(16,*)' RAYTRACE : ATMPTH: COMPUTED ALTITUDE VALUE OF H1 IS NEGATIVE'
               WRITE( 0,*)' RAYTRACE : ATMPTH: COMPUTED ALTITUDE VALUE OF H1 IS NEGATIVE'
               CALL SHUTDOWN
               STOP '3'
            ENDIF

            PTMP(1)  = 0.0
            TTMP(1)  = 0.0
            WVTMP(1) = 0.0
            ZTMP(1)  = 0.0

            PTMP(2)  = 0.0
            TTMP(2)  = 0.0
            WVTMP(2) = 0.0
            ZTMP(2)  = 0.0

            DO 168 LIP = 2,IMMAX
               IF (H2 .GT. PM(LIP)) GO TO 169
  168          CONTINUE
               LIP = IMMAX
  169       CONTINUE
            IF (H2 .EQ. PM(LIP-1)) THEN
               H2 = ZMDL(LIP-1)
            ELSE
               IF(H2 .EQ. PM(LIP)) THEN
                  H2 = ZMDL(LIP)
               ELSE
! PERFORM INTERPOLATION IN LN(PM)
                  HIP  = (ZMDL(LIP)-ZMDL(LIP-1))/ LOG(PM(LIP)/PM(LIP-1))
                  ZINT = ZMDL(LIP-1)+HIP* LOG(H2/PM(LIP-1))

! PERFORM ALTITUDE CALCULATION USING HYDROSTATIC EQUATION
                  PTMP(1)  = PM(LIP-1)
                  ZTMP(1)  = ZMDL(LIP-1)
                  TTMP(1)  = TM(LIP-1)
                  WVTMP(1) = DENW(LIP-1)

                  PTMP(2) = H2

                  TIP     = (TM(LIP)-TM(LIP-1))/LOG(PM(LIP)/PM(LIP-1))
                  TTMP(2) = TM(LIP-1)+TIP* LOG(H2/PM(LIP-1))

                  WVIP     =  (DENW(LIP)-DENW(LIP-1))/LOG(PM(LIP)/PM(LIP-1))
                  WVTMP(2) =  DENW(LIP-1) + WVIP* LOG(H2/PM(LIP-1))
                  CALL CMPALT(2,PTMP,TTMP,WVTMP,ZTMP(1),ZTMP)

! COMBINE THE INTERPOLATION AND THE HYDROSTATIC CALCULATION

                  RATP =  LOG(H2/PM(LIP-1)) / LOG(PM(LIP)/PM(LIP-1))
                  A    = RATP**3
                  H2   = A*ZINT + (1-A)*ZTMP(2)

               ENDIF
            ENDIF

            IF (H2 .LT. 0.0) THEN
               PRINT 946, H2,ZTMP(1)
               IF (NOPRNT.GE.0) WRITE (IPR,946) H2,ZTMP(1)
               WRITE(16,*) ' RAYTRACE : ATMPTH: COMPUTED ALTITUDE VALUE OF H2 IS NEGATIVE'
               WRITE( 0,*) ' RAYTRACE : ATMPTH: COMPUTED ALTITUDE VALUE OF H2 IS NEGATIVE'
               CALL SHUTDOWN
               STOP '3'
            ENDIF

! --- END IF IBMAX_B < 0 - PRESSURE BOUDRARIES - INTERPOLATING ONTO Z BOUNDS

         ENDIF

         IF (IBMAX.GE.1) THEN
            IF (ZBND(1).LT.ZMDL(1)) THEN
               IF (NOPRNT.GE.0) WRITE (IPR,944)
               IF (ABS(ZBND(1)-ZMDL(1)).LE.0.0001) THEN
                  ZBND(1) = ZMDL(1)
               ELSE
                  PRINT 946,ZBND(1),ZMDL(1)
                  IF (NOPRNT.GE.0) WRITE (IPR,946) ZBND(1),ZMDL(1)
                  WRITE(16,*) ' RAYTRACE : ATMPTH: BOUNDARIES OUTSIDE OF ATMOS'
                  WRITE( 0,*) ' RAYTRACE : ATMPTH: BOUNDARIES OUTSIDE OF ATMOS'
                  CALL SHUTDOWN
                  STOP '3'
               ENDIF
            ENDIF
         ENDIF

! --- COMPUTE THE REFRACTIVE INDEX PROFILE
! --- RFNDXM IS 1.0-INDEX
! --- EQUATION FOR RFNDXM IS FROM LOWTRAN (REF 3)

         IF (NOPRNT.GE.0) WRITE(IPR,899) 'USING LOWTRAN 6 REFRACTIVE INDEX CALCULATION ON MODEL LEVELS'
         DO 170 IM = 1, IMMAX

            PPH2O = DENM(1,IM)*PZERO*TM(IM)/(TZERO*ALOSMT)

! --- APPROXIMATION TO REFRACTION INDEX (FROM LOWTRAN5)

!           RFNDXM(IM) = ((77.46+0.459E-8*XVBAR**2)*PM(IM)/TM(IM)-
!    *                   (PPH2O/1013.0)*(43.49-0.347E-8*XVBAR**2))*
!    *                   1.0E-6

! --- APPROXIMATION TO REFRACTION INDEX (FROM LOWTRAN6)

            RFNDXM(IM)=((83.42D0+(185.08D0/(1.0D0-(XVBAR/1.14D+5)**2))+       &
                       (4.11D0/(1.0D0-(XVBAR/6.24D+4)**2)))*(PM(IM)*288.15D0)/            &
                       (PZERO*TM(IM))-(43.49D0-(XVBAR/1.7D+4)**2)*(PPH2O/PZERO))*1.0D-06
  170    CONTINUE


! --- PRINT ATMOSPHERIC MODEL PROFILE

! --- DENSITY AT MODEL LEVELS
         IF (NOPRNT.GT.1) THEN
            WRITE (IPR,950) (HMOLS(K),K=1,NMOL)
            !WRITE (IPR,952)
         ENDIF

         DO 180 IM = 1, IMMAX

! --- DENG=DENM(1,IM)*2.989641E-17

            DENAIR = ALOSMT*(PM(IM)/PZERO)*(TZERO/TM(IM))
            DENSAVE(IM) = DENAIR
            IF( NOPRNT .GE. 0 )                                         &
                 WRITE (IPR,954) IM,ZMDL(IM),PM(IM),TM(IM),RFNDXM(IM),  &
                                 DENAIR,(DENM(K,IM),K=1,NMOL)

  180    CONTINUE


! --- MIXING RATIO AT MODEL LEVELS USING DRY AIR
! --- CHANGED FROM PPM TO UNITLESS VMR AT WRITE STMNT

         IF (NOPRNT.GT.1)THEN
            WRITE (IPR,952)
            WRITE (IPR,951) (HMOLS(K),K=1,NMOL)
         ENDIF

         DO 188 IM = 1, IMMAX

            DRY_AIR = DENSAVE(IM) - DENM(1,IM)

            IF( NOPRNT .GT. 1 )                                       &
                WRITE(IPR,954) IM,ZMDL(IM),PM(IM),TM(IM),RFNDXM(IM),  &
                               DENSAVE(IM),(((DENM(K,IM)/DRY_AIR)*1.0D0),K=1,NMOL)

  188    CONTINUE


! --- REDUCE SLANT PATH PARAMETERS TO STANDARD FORM USING FPACK()
! --- CALCULATE APPARENT ZENITH ANGLE ITERATIVELY USING ASTRO()


!--- ITERATION FOR ASTRO TO APPARENT SZA  LABELS 400-4001
         KASTRO = 0
         AST    = ANGLE    ! SAVE ASTRONOMICAL SZA
  400    CONTINUE

         CALL FSCGEO (H1,H2,ANGLE,RANGE,BETA,ITYPE,LEN,HMIN,PHI,IERROR,HOBS)
         IF (IERROR.NE.0) GO TO 310

! --- SET UP LBLRTM LAYER BOUNDARIES

         IF (IBMAX.NE.0) GO TO 200

! --- AUTOMATIC LAYERING SELECTED

         HMAX =   MAX(H1,H2)
         CALL AUTLAY (HMIN,HMAX,XVBAR,AVTRAT,TDIFF1,TDIFF2,ALTD1,ALTD2,IERROR)
         GO TO 220

  200    CONTINUE

! --- USER SUPPLIED LAYERING

!         IF (NOPRNT .GE. 0) WRITE (IPR,956)
!
!! --- CALCULATE HALFWIDTHS - DO WE NEED THIS???
!         IF (IBMAX_B .LT. 0) THEN
!! --- PRESSURE FINAL BOUNDS
!            DO 205 IB = 1, IBMAX
!               CALL HALFWD_P(ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB),ADOPP(IB),AVOIGT(IB))
!  205       CONTINUE
!         ELSE
!! --- ALTITUDE FINAL BOUNDS
!            DO 210 IB = 1, IBMAX
!               CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB),ADOPP(IB),AVOIGT(IB))
!  210       CONTINUE
!         ENDIF
!
  220    CONTINUE
!
!!         IF (NOPRNT .GE. 1) WRITE (IPR,958) ALZERO,AVMWT,XVBAR
!!         DO 230 IB = 1, IBMAX
!!
!!            ZETA  = ALORNZ(IB)/(ALORNZ(IB)+ADOPP(IB))
!!            RATIO = 0.0D0
!!            DTEMP = 0.0D0
!!            IF (IB.NE.IBMAX) RATIO = AVOIGT(IB)/AVOIGT(IB+1)
!!            IF (IB.NE.IBMAX) DTEMP = ABS(TBND(IB)-TBND(IB+1))
!!            IF (NOPRNT .GE. 2) &
!!                WRITE (IPR,960) IB,ZBND(IB),PBND(IB),TBND(IB),ALORNZ(IB), &
!!                                ADOPP(IB),ZETA,AVOIGT(IB),RATIO,DTEMP
!!  230    CONTINUE

         IF (IERROR.NE.0) THEN
            WRITE(16,*) ' RAYTRACE : ATMPTH: IERROR'
            WRITE( 0,*) ' RAYTRACE : ATMPTH: IERROR'
            CALL SHUTDOWN
            STOP '3'
         ENDIF

! --- CALCULATE THE REFRACTED PATH THROUGH THE ATMOSPHERE

         CALL RFPATH (H1,H2,ANGLE,PHI,LEN,HMIN,IAMT,RANGE,BETA,BENDNG)

! --- ITERATION FOR ASTRO TO APPARENT SZA

         IF(KASTRO.EQ.1 .OR. IASTRO.EQ.0 .OR. ABS(ANGLE).LT. 0.00001) GOTO 401
         CALL ASTRO(KASTRO,HMIN,PHI,IERROR,BENDNG,H1,H2,ANGLE,RANGE,BETA,ITYPE,LEN,HOBS)
         GO TO 400

  401    CONTINUE
         !APP = ANGLE


! --- PRINT AMOUNTS BY LAYER AND ACCUMULATE TOTALS

         IF( NOPRNT .GE. 2 )WRITE(IPR,962) (HMOLS(K),K=1,NMOL)
         I2        = IPMAX-1
         AIRTOT    = 0.0D0
         AMTTOT(:) = 0.0D0
         HMID      = MIN(H1,H2)
         DO 260 I = 1, I2
            FAC = 1.0D0
            IF( LEN .EQ. 1 .AND. ZPTH(I+1) .LE. HMID ) FAC = 2.0D0
            AMTAIR = RHOPSM(I) * 1.0D5
            AIRTOT = AIRTOT + FAC * AMTAIR
            DO 250 K = 1, NMOL
               AMTTOT(K) = AMTTOT(K) + FAC * AMTP(K,I)
  250       CONTINUE
            IF( NOPRNT .GE. 2 )THEN
               WRITE(IPR,964) I,ZPTH(I),ZPTH(I+1),AMTAIR,(AMTP(K,I),K=1,NMOL)
               ENDIF
  260    CONTINUE

         IF( NOPRNT .GE. 0 )WRITE(IPR,966) H1,H2,AIRTOT,(AMTTOT(K),K=1,NMOL)

! --- PRINT SUMMARY OF PATH

         AIRMAS = AIRTOT/AIRMS1
         IF (NOPRNT .GE. 0) &
         WRITE (IPR,968) HMOD,H1,H2,ANGLE,RANGE,BETA,PHI,HMIN,BENDNG,LEN,AIRMAS
         IF (ITYPE.EQ.3) ITYPE = 2
         H1F     = H1
         H2F     = H2
         ANGLEF  = ANGLE
         LENF    = LEN

! --- CONDENSE THE AMOUNTS INTO THE LBLRTM OUTPUT LAYERS ZOUT,
! ---  WHICH ARE DEFINED BY THE BOUNDARIES ZBND FROM HMIN TO
! ---  HMAX ALSO, ZERO OUT THE AMOUNT FOR A MOLECULE IF THE
! ---  CUMULATIVE AMOUNT FOR THAT LAYER AND ABOVE IN LESS THAN
! ---  0.1 PERCENT OF THE TOTAL

         CALL FPACK (H1,H2,HMID,LEN,IEMIT,NOZERO)

! --- OUTPUT THE FINAL PROFILE IN COLUMN DENSITY AND MIXING RATIO FROM DRY AIR IN ZFIN GRID

         LMAX = IFINMX-1
         IF( IREAD .EQ. 0 )WRITE(16,890) " NUM LAYERS IN FINAL OUTPUT GRID   : ", LMAX

         IF (NOPRNT.GE.2) THEN
             WRITE (IPR,970) (HMOLS(K),K=1,7),COTHER,(HMOLS(K),K=8,NMOL)
         ENDIF

         IF (IPUNCH.GE.1) THEN
            IFORM = 1
            WRITE(IPU,972)IFORM,LMAX,NMOL,SECNT0,(HMOD(I),I=1,2),H1,H2,ANGLE,LEN
         ENDIF

         IF (IPUNCH.EQ.2) THEN
            WRITE (97) LMAX,NMOL,SECNT0,(HMOD(I),I=1,2),H1,H2,ANGLE,LEN
         ENDIF

         SUMN2 = 0.0D0
         SUMRS = 0.0D0
         PWTD  = 0.0D0
         TWTD  = 0.0D0
         WTOT  = 0.0D0

!        READ FROM/WRITE TO "IFIXTYPE" FILE: IF IFXTYP = -2, USE
!        PRESET ITYL VALUES; IF IFXTYP = 2, CALCULATE AND WRITE OUT
!        ITYL VALUES.

         IF (IFXTYP.EQ.-2) THEN
            OPEN(99,FILE='IFIXTYPE',STATUS='OLD', FORM='FORMATTED')
         ELSEIF (IFXTYP.EQ.2) THEN
            OPEN(99,FILE='IFIXTYPE',STATUS='UNKNOWN', FORM='FORMATTED')
         ENDIF

         DO 280 L = 1, LMAX

            FACTOR = 1.0D0
            IF (IPATH(L).EQ.2) FACTOR = 2.0D0
            SUMWK = 0.0D0

            DO 270 MOL = 1, NMOL
               SUMWK    = SUMWK + AMOUNT(MOL,L)
               WMT(MOL) = WMT(MOL) + AMOUNT(MOL,L)*FACTOR
  270       CONTINUE

            WTOTL(L)  = SUMWK+WN2L(L)
            SUMN2     = SUMN2+WN2L(L)*FACTOR
            SUMRS     = SUMRS+RHOSUM(L)*FACTOR
            WTOT      = WTOT+WTOTL(L)*FACTOR
            PWTD      = PWTD+PBAR(L)*WTOTL(L)*FACTOR
            TWTD      = TWTD+TBAR(L)*WTOTL(L)*FACTOR

!           > DETERMINE ITYL(L), IF DESIRED (WHEN SETTING THE RATIO <
!           > FROM LAYER TO LAYER).  DEFAULT IS ITYL(L) LEFT BLANK, <
!           >                                                       <
!           >                    CTYPE(L) = '   '                   <

            CTYPE(L) = '   '
            IF (IFXTYP.EQ.1) THEN
               FRH2O  = AMOUNT(1,L)/WTOTL(L)
               ALFCOR = (PBAR(L)/PZERO)*SQRT(T296/TBAR(L))
               ADBAR = 3.581155E-07*XVBAR*SQRT(TBAR(L)/AVMWT)
               CALL FIXTYP(IEMIT,FRH2O,ALFCOR,OLDDV,L,CTYPE(L))
               READ(CTYPE(L),1100) ITYL(L)
            ELSEIF (IFXTYP.EQ.2) THEN
               FRH2O  = AMOUNT(1,L)/WTOTL(L)
               ALFCOR = (PBAR(L)/PZERO)*SQRT(T296/TBAR(L))
               ADBAR = 3.581155E-07*XVBAR*SQRT(TBAR(L)/AVMWT)
               CALL FIXTYP(IEMIT,FRH2O,ALFCOR,OLDDV,L,CTYPE(L))
               READ(CTYPE(L),1100) ITYL(L)
               WRITE(99,1100) ITYL(L)
            ELSEIF (IFXTYP.EQ.-2) THEN
               READ(99,1100) ITYL(L)
               WRITE(CTYPE(L),1100) ITYL(L)
            ENDIF

! --- WRITE ATMOSPHERE TO TAPE6 IN COLUMN DENSITY [MOL*CM-2]

            IF (NOPRNT .GE. 2) THEN
               WRITE (IPR,976) L,ZFIN(L),ZFIN(L+1),CTYPE(L),IPATH(L),&
                                  PBAR(L),TBAR(L),RHOSUM(L),            &
                                  (AMOUNT(K,L),K=1,7), WN2L(L),          &
                                  (AMOUNT(K,L),K=8,NMOL)
               ENDIF

! --- WRITE ATMOSPHERE TO TAPE7

            IF (IPUNCH.GE.1) THEN
               LTST = L
               IF (L.EQ.1) LTST = 0
               PTST  =  LOG10(PZ(LTST))
               NPTST = NINT(PTST) + 2
               IF (PTST.LT.0.0) NPTST = 1
               CFORM1(38:41) = PZFORM(NPTST)
               CFORM2(38:41) = PZFORM(NPTST)
               NPTST = 1
               IF (PBAR(L).GE.0.1) NPTST = 2
               CFORM1(2:8) = PAFORM(NPTST)
               CFORM2(2:8) = PAFORM(NPTST)
               IF (L.EQ.1) THEN
                  WRITE (IPU,CFORM1) PBAR(L),TBAR(L),CTYPE(L),IPATH(L), &
     &                               ALTZ(L-1),PZ(L-1),TZ(L-1),         &
     &                               ALTZ(L),  PZ(L),  TZ(L)
               ELSE
                  WRITE (IPU,CFORM2) PBAR(L),TBAR(L),CTYPE(L),IPATH(L), &
     &                               ALTZ(L),  PZ(L),   TZ(L)
               ENDIF

! --- WRITE MOLECULAR INFORMATION TO TAPE7 IN
! ---  MIXING RATIO IF MUNITS IS 1
! ---  COLUMN DENSITY IF MUNITS IS 0

               IF (MUNITS.EQ.1) THEN
                  DRAIRL(L) =  WN2L(L)
                  DO 275 M = 2, NMOL
                     DRAIRL(L) = DRAIRL(L) + AMOUNT(M,L)
  275             CONTINUE

! --- IF DRAIR IS ZERO, THEN WRITE OUT AMOUNT ONLY
                  IF( noprnt .ge.0 )then
                     IF (DRAIRL(L).EQ.0) THEN
                        WRITE (IPU,978) (AMOUNT(K,L),K=1,7),WN2L(L)
                        IF (NMOL.GT.7) WRITE (IPU,978) (AMOUNT(K,L),K=8,NMOL)
                     ELSE
                        WRITE (IPU,978) (AMOUNT(K,L)/DRAIRL(L),K=1,7),WN2L(L)
                        IF (NMOL.GT.7) WRITE (IPU,978) (AMOUNT(K,L)/DRAIRL(L),K=8,NMOL)
                     ENDIF
                  endif
               ELSE

! --- TEST TO MAKE SURE THERE ARE NO FRACTIONAL MOLECULAR
! --- AMOUNTS WRITTEN OUT (WILL CAUSE PATH TO ASSUME MIXING RATIO)

                  DO 277 K=1,NMOL
                     IF (AMOUNT(K,L).LT.1.) THEN
                        IF( noprnt .ge.0 )WRITE(IPR,1000) K,L
                        AMOUNT(K,L) = 0.0
                     ENDIF
  277             CONTINUE

                  IF( noprnt .ge.0 )WRITE (IPU,978) (AMOUNT(K,L),K=1,7),WN2L(L), (AMOUNT(K,L),K=8,NMOL)

               ENDIF

            ENDIF

  280    CONTINUE

! --- CLOSE "IFIXTYPE" FILE, IF USED
         IF (ABS(IFXTYP).EQ.2) THEN
            REWIND(99)
            CLOSE(99)
         ENDIF


! --- WRITE ATMOSPHERE TO TAPE6 IN MIXING RATIO IN DRY AIR

         IF (NOPRNT .GE. 2) &
               WRITE (IPR,973) (HMOLS(K),K=1,7), COTHER, (HMOLS(K),K=8,NMOL)

         DO 285 L = 1, LMAX

            DRAIRL(L) = WN2L(L)
            DO M = 2, NMOL
               DRAIRL(L) = DRAIRL(L) + AMOUNT(M,L)
           ENDDO

! --- IF DRAIR IS ZERO, THEN WRITE OUT AMOUNT ONLY

            IF( DRAIRL(L) .EQ. 0 .AND. NOPRNT.GE.2 )THEN
               WRITE (IPR,976) L,ZFIN(L),ZFIN(L+1),CTYPE(L),IPATH(L),&
                                  PBAR(L),TBAR(L),RHOSUM(L),            &
                                  (AMOUNT(K,L),K=1,7),WN2L(L),          &
                                  (AMOUNT(K,L),K=8,NMOL)
            ELSE IF (NOPRNT.GE.2) THEN
               WRITE (IPR,976) L,ZFIN(L),ZFIN(L+1),CTYPE(L),IPATH(L),&
                                  PBAR(L),TBAR(L),RHOSUM(L),            &
                                  (AMOUNT(K,L)/DRAIRL(L),K=1,7), WN2L(L)/DRAIRL(L),    &
                                  (AMOUNT(K,L)/DRAIRL(L),K=8,NMOL)
            ENDIF

  285   CONTINUE


! TOTAL COLUMNS

         PWTD = PWTD/WTOT
         TWTD = TWTD/WTOT
         L = LMAX
         IF (NOPRNT .GE. 0) THEN
               WRITE (IPR,980) !(HMOLS(K),K=1,7), COTHER, (HMOLS(K),K=8,NMOL)
               WRITE (IPR,984) L,ZFIN(1),ZFIN(L+1),PWTD,TWTD,SUMRS,     &
                               (WMT(K),K=1,7),SUMN2,(WMT(K),K=8,NMOL)
         ENDIF

! WRITE TO ATMOSPHERE FILE (97) INFORMATION NEEDED FOR LEVEL TO LAYER CA

         IF (IPUNCH.EQ.2) THEN
            WRITE (97) IBMAX,(PBAR(L),TBAR(L),L=1,IBMAX-1)
            WRITE (97) (PBND(L),TBND(L),(DENM(K,L),K=1,NMOL),L=1,IBMAX)
            CLOSE (97)
         ENDIF

      ENDIF

      !print*, nlev, lmax
      IF( LMAX .LT. NLEV )NLEV = LMAX
      CALL REFOUT( ITER, JSPEC, IREAD, LMAX, AST, APP, BENDNG )

      RETURN

!     ERROR MESSAGES

  290 WRITE (0,986) MODEL,ITYPE,NMOL,IBMAX
      WRITE(16,*) ' RAYTRACE : ATMPTH: CARD 3.1'
      WRITE( 0,*) ' RAYTRACE : ATMPTH: CARD 3.1'
      CALL SHUTDOWN
      STOP '3'

  300 WRITE (0,988) (ZBND(I),I=1,IBMAX)
      PRINT 988,(ZBND(I),I=1,IBMAX)
      WRITE(16,*) ' RAYTRACE : ATMPTH: ZBND'
      WRITE( 0,*) ' RAYTRACE : ATMPTH: ZBND'
      CALL SHUTDOWN
      STOP '3'

!  301 PRINT 988,(ZTMP(I),I=2,IBMAX)
!      STOP ' USER INPUT LEVELS TOO CLOSE - IBMAX'

  305 WRITE (0,989) (PBND(I),I=1,IBMAX)
      PRINT 989,(PBND(I),I=1,IBMAX)
      WRITE(16,*) ' RAYTRACE : ATMPTH: PBND'
      WRITE( 0,*) ' RAYTRACE : ATMPTH: PBND'
      CALL SHUTDOWN
      STOP '3'

  310 WRITE (0,990)
      WRITE(16,*) ' RAYTRACE : ATMPTH: FSCGEO'
      WRITE( 0,*) ' RAYTRACE : ATMPTH: FSCGEO'
      CALL SHUTDOWN
      STOP '3'

  320 WRITE (0,992) AVTRAT,TDIFF1,TDIFF2
      WRITE(16,*) ' RAYTRACE : AVTRAT,TDIFF'
      WRITE( 0,*) ' RAYTRACE : AVTRAT,TDIFF'
      CALL SHUTDOWN
      STOP '3'

  890 FORMAT (A,I12)
  891 FORMAT( /, A, 2I6, F12.5)

  899 FORMAT (/,A,/)
!  900 FORMAT (7I5,I2,1X,I2,5F10.3,I5)
  902 FORMAT (' CONTROL CARD 3.1: MODEL AND OPTIONS ')
  904 FORMAT (10X,'MODEL   = ',I5,/,10X,'ITYPE   = ',I5,/,10X,        &
     &        'IBMAX   = ',I5,/,10X,'NOZERO  = ',I5,/,10X,'NOPRNT  = ', &
     &        I5,/,10X,'NMOL    = ',I5,/,10X,'IPUNCH  = ',I5,/,10X,     &
     &        'IFXTYP  = ',I5,/,10X,'MUNITS  = ',I5,/,10X,'RE      = ', &
     &        F10.3,' KM',/,10X,'HSPACE  = ',F10.3,' KM',/,10X,         &
     &        'VBAR    = ',F10.3,' CM-1',                               &
     &                     /,10X,'REF_LAT = ',F10.3, ' DEG',            &
     &        /,10X,'IASTRO  = ',I5)
  905 FORMAT('$',I5, 10A8)
  906 FORMAT (/,' CONTROL CARD 3.1 PARAMETERS WITH DEFAULTS:')
  908 FORMAT (//,' HORIZONTAL PATH SELECTED')
  910 FORMAT (F10.3,10X,10X,F10.3)
  912 FORMAT (/,' CONTROL CARD 3.2:',//,10X,'Z     = ',F10.3,' KM',/, &
     &       10X,'RANGE = ',F10.3,' KM')
  914 FORMAT (///' PRESSURE, TEMPERATURE, AND DENSITIES INTERPOLATED',  &
     &        ' FROM THE FOLLOWING ATMOSPHERIC MODEL: ',//,10X,3A8,//,  &
     &        10X, 'Z     = ',F10.3,' KM',/,10X,'P     = ',F10.3,' MB', &
     &        /,10X,'T     = ',F10.3,' K',//,10X,'DENSITIES :',T26,     &
     &        'AIR',(T30,8A10))
  916 FORMAT (T63,'(MOL CM-3)',//,T20,1PE10.3,(T30,8E10.3))
  918 FORMAT ('0SINGLE LAYER INPUT TO LBLRTM',//,10X,'MODEL = ',3A8,/,  &
     &        10X,'Z     = ',F10.3,' KM',/,10X,'P     = ',F10.3,' MB',  &
     &        /,10X,'T     = ',F10.3,' K',/,10X,'RANGE = ',F10.3,' KM', &
     &        //,10X,'AMOUNTS (MOL CM-2):',T36,'AIR',(T32,8A10))
  920 FORMAT (//,T30,1PE10.2,(T30,8E10.2))
!  922 FORMAT (A4)
  924 FORMAT (1X,I1,I3,I5,F10.6,3A8,' * ',F7.3,' KM PATH AT ',F7.3,     &
     &        ' KM ALT')
  926 FORMAT (E15.7,F10.4,10X,I5,1X,F7.3,15X,F7.3,/,(1P8E15.7))
!  928 FORMAT (//,' MULTIPLE SCATTERING TURNED OFF, HMIN = ',F10.6,      &
!     &        ' > HMAXMS = ',F10.6,/)
  930 FORMAT (/,' SLANT PATH SELECTED, ITYPE = ',I5)
!  931 FORMAT (//,' TEST BLOCK ',2I5,A,A,/)

  !932 FORMAT (5F10.4,I5,5X,F10.4)
  933 FORMAT (/' CONTROL CARD 3.2:  SLANT PATH PARAMETERS',/,10X,    &
     &       'H1      = ',F12.6,' MBAR',/,10X,'H2      = ',F12.6,' MBAR'&
     &        /,10X,'ANGLE   = ',F12.6,' DEG',/,10X,'RANGE   = ',F12.6, &
     &        ' KM',/,10X,'BETA    = ',F12.6,' DEG',/,10X,'LEN     = ', &
     &        I10)
  934 FORMAT (/' CONTROL CARD 3.2:  SLANT PATH PARAMETERS',//,10X,    &
     &        'H1      = ',F12.6,' KM',/,10X,'H2      = ',F12.6,' KM',  &
     &        /,10X,'ANGLE   = ',F12.6,' DEG',/,10X,'RANGE   = ',F12.6, &
     &        ' KM',/,10X,'BETA    = ',F12.6,' DEG',/,10X,'LEN     = ', &
     &        I10)
  936 FORMAT (5F10.3)
  938 FORMAT (///,' AUTOLAYERING SELECTED',//,10X,'AVTRAT    = ',F8.2,  &
     &        /,10X,'TDIFF1    = ',F8.2,/,10X,'TDIFF2    = ',F8.2,/,    &
     &        10X,'ALTD1     = ',F8.2,/,10X,'ALTD2     = ',F8.2)
  940 FORMAT (8F10.3)
  942 FORMAT (///,' USER DEFINED BOUNDARIES FOR LBLRTM LAYERS',/,10X,   &
     &        'I',4X,'Z (KM)',//,(10X,I4,F15.8))
  943  FORMAT (///,' USER DEFINED BOUNDARIES FOR LBLRTM LAYERS',/,10X,  &
     &        'I',4X,'P (MB)',//,(10X,I4,F15.8))
  944 FORMAT (' ERROR IN USER INPUT BOUNDARIES ')
  946 FORMAT (' BOUNDARIES ARE OUTSIDE THE RANGE OF THE ATMOSPHERE',/,  &
     &        ' BOUNDARY = ',F10.2,' ATMOSPHERE =',F10.2,/,             &
     &        ' RESET BOUNDARY GT THAN ATMOSPHERE')
!  948 FORMAT ('1ATMOSPHERIC PROFILE SELECTED IS: M = ',I3,5X,3A8)

! --- CHANGED 8 TO 100 FOR WIDE LIST OF MOLS 950, 951, 954, 962, 964, 966, 970, 976, 973, 980, 984
! --- CHANGED 100 TO 200 LINES 976 984 973 978

  950 FORMAT (/,T4,'I',T13,'Z',T24,'P',T38,'T',T46,'REFRACT',T73,       &
              'DENSITY  AT MODEL LEVELS [MOLS CM-3]',/,T46,'INDEX-1',/,T12,'(KM)',T23,  &
              '(MB)',T37,'(K)',T46,'*1.0E6',T63,'AIR',(T68,200(6X,A9)))
  951 FORMAT (/,T4,'I',T13,'Z',T24,'P',T38,'T',T46,'REFRACT',T55,       &
              'DENSITY',T70,'MIXING RATIO ON MODEL LEVELS (BASED UPON DRY AIR)',/&
              T46,'INDEX-1',T56,                                        &
              '(MOL CM-3)'/,T12,                                        &
              '(KM)',T23,                                               &
              '(MB)',T37,'(K)',T46,'*1.0E6',T61,'AIR',(T68,200(6X,A9)))
  952 FORMAT (/)
  954 FORMAT (I4,F11.5,F15.8,F11.5,6P,F11.5,1P,E15.7,(T68,1P,200E15.7))
!  956 FORMAT (///,' HALFWIDTH INFORMATION ON THE USER SUPPLIED ',       &
!     &        'LBLRTM BOUNDARIES',/,' THE FOLLOWING VALUES ARE ',       &
!     &        'ASSUMED:')
!  958 FORMAT (10X,'ALZERO    = ',F9.3,' CM-1 = AVERAGE LORENTZ WIDTH ', &
!     &        'AT STP',/,10X,'AVMWT     = ',F8.2,                       &
!     &        '       = AVERAGE  MOLECULAR WEIGHT',/,10X,               &
!     &        'VBAR      = ',F8.2,'  CM-1 = AVERAGE WAVENUMBER',///,    &
!     &        T5,'I',T12,'Z',T22,'P',T32,'T',T39,'LORENTZ',T49,         &
!     &        'DOPPLER',T61,'ZETA',T70,'VOIGT',T80,'VOIGT',T90,         &
!     &        'TEMP',/,T11,'(KM)',T21,'(MB)',T31,'(K)',T40,'(CM-1)',    &
!     &        T50,'(CM-1)',T70,'(CM-1)',T80,'RATIO',T90,'DIFF (K)',/)
!  960 FORMAT (I5,F10.3,F12.5,F9.2,F9.5,F10.5,F10.3,F10.5,F10.2,F10.1)
  962 FORMAT (/,'INTEGRATED ABSORBER AMOUNTS BY LAYERS IN PATHS - POST APPARENT SZA CALCULATION',//,T5, &
              'I  LAYER BOUNDARIES',T55,'INTEGRATED AMOUNTS ',          &
              '(MOL CM-2)',/,T11,'FROM',T22,'TO',T29,'AIR',T36,         &
              200(1X,A8,1X),/,T11,'(KM)',T21,'(KM)',(T37,200A10))
  964 FORMAT (I5,2F10.3,1P,E10.3,(T36,1P,200E10.3))
  966 FORMAT (/,'TOTAL',F9.3,F10.3,1PE10.3,(T35,1P200E10.3))
  968 FORMAT (/,'SUMMARY OF THE GEOMETRY CALCULATION',//,10X,           &
              'MODEL   = ',4X,3A8,/10X,'H1      = ',F12.6,' KM',/,10X,  &
              'H2      = ',F12.6,' KM',/,10X,'ANGLE   = ',F12.6,' DEG', &
              /,10X,'RANGE   = ',F12.6,' KM',/,10X,'BETA    = ',F12.6,  &
              ' DEG',/,10X,'PHI     = ',F12.6,' DEG',/,10X,             &
              'HMIN    = ',F12.6,' KM',/,10X,'BENDING = ',F12.6,' DEG', &
              /,10X,'LEN     = ',I10,/,10X,'AIRMAS  = ',G12.6,          &
              'RELATIVE TO A VERTICAL PATH , GROUND TO SPACE')
  970 FORMAT (/,'FINAL SET OF LAYERS FOR INPUT TO LBLRTM',/,             &
              ' A LAYER AMOUNT MAY BE SET TO ZERO IF THE CUMULATIVE ',  &
              'AMOUNT FOR THAT LAYER AND ABOVE IS LESS THAN 0.1 ',      &
              'PERCENT',/,' OF THE TOTAL AMOUNT. THIS IS DONE ONLY ',   &
              'FOR THE FOLLOWING CASES',/,5X,'1.  IEMIT = 0 ',          &
              '(TRANSMITTANCE)',/,5X,'2.  IEMIT = 1 (RADIANCE) AND ',   &
              'IPATH = 3 (PATH LOOKING UP)',/,' O2 IS NOT INCLUDED',/,  &
              ' IF THE AMOUNTS FOR ALL THE MOLECULES BUT O2 ARE ',      &
              'ZEROED, THE REMAINING LAYERS ARE ELIMINATED',///,T13,    &
              'LAYER',T23,'I',T25,'I',/,T4,'L',T10,'BOUNDARIES',T23,    &
              'T',T25,'P',T31,'PBAR',T40,'TBAR',T70,                    &
              'INTEGRATED AMOUNTS (MOLS CM-2)',/,T9,'FROM',T18,         &
              'TO',T23,'Y',T25,'T',/,T9,'(KM)',T17,'(KM)',T23,'P',T25,  &
              'H',T31,'(MB)',T41,'(K)',T53,'AIR',(T59,200(6X,A9)))
  972 FORMAT (1X,I1,I3,I5,F10.6,2A8,' H1=',F8.2,' H2=',F8.2,            &
     &        ' ANG=',F8.3,' LEN=',I2)
  973 FORMAT  (/,     &
             T13,'LAYER',T23,'I',T25,'I',/,T4,'L',T10,'BOUNDARIES',    &
             T23,'T',T25,'P',T31,'PBAR',T40,'TBAR',                    &
             T68,'MOLECULAR MIXING RATIOS BY LAYER IN DRY AIR',/,T9,'FROM',       &
             T18,'TO',T23,'Y',T25,'T',/,T9,'(KM)',T17,'(KM)',T23,'P',  &
             T25,'H',T31,'(MB)',T41,'(K)',T53,'AIR',(T59,200(6X,A9)))
!  974 FORMAT ('0',I3,2F8.3,A3,I2,F11.5,F8.2,1X,1P9E15.7)
  976 FORMAT ('0',I3,2F8.3,A3,I2,F11.5,F8.2,1X,1P200E15.7) !,/,(60X,1P100E15.7))
  978 FORMAT (1P200E15.7)
!  980 FORMAT (/,T4,'L  PATH BOUNDARIES',T28,'PBAR',T37,'TBAR',  &
!              T65,'ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH',/,T9,  &
!              'FROM',T18,'TO',/,T9,'(KM)',T17,'(KM)',T28,'(MB)',T38,    &
!              '(K)',T47,'AIR',(T54,100(1X,A9)))
  980 FORMAT (/,T65,'ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH')

!  982 FORMAT ('0',I3,2F8.3,2X,F11.5,F8.2,1X,1P9E10.3)
  984 FORMAT ('0',I3,2F8.3,5X,F11.5,F8.2,1X,1P200E15.4) !,/,(52X,1P8E10.3))
  986 FORMAT (///,' ERROR IN INPUT, CONTROL CARD 3.1: ONE OF THE ',     &
     &        'PARAMETERS MODEL, ITYPE, NMOL IS OUT OF RANGE',//,10X,   &
     &        'MODEL   = ',I5,/,10X,'ITYPE   = ',I5,/,10X,'NMOL    = ', &
     &        I5,10X,' IBMAX =',I5)
  988 FORMAT (///,' ERROR: BOUNDARY ALTITUDES FOR LBLRTM LAYERS ',      &
     &        'ARE NEGATIVE OR NOT IN ASCENDING ORDER',//,5X,' ZBND ',  &
     &        /,(10F10.4))
  989 FORMAT (///,' ERROR: BOUNDARY PRESSURES FOR LBLRTM LAYERS ',      &
     &        'ARE POSITIVE OR NOT IN DESCENDING ORDER',//,5X,' PBND ', &
     &        /,(10F10.4))
  990 FORMAT ('0ERROR FLAG RETURNED FROM FSCGEO:  AN ERROR OCCURED ',   &
     &        'IN PROCESSING THE SLANT PATH PARAMETERS',/,'0PROGRAM ',  &
     &        'STOP')
  992 FORMAT (///,' ERROR: EITHER AVTRAT.LE.1.0 OR TDIFF.LE.0',/,       &
     &        '0PROGRAM STOP  -  AVTRAT = ',E12.6,' TDIFF1 = ',F10.4,   &
     &        ' TDIFF2 = ',F10.4)
 1000 FORMAT ('*** WARNING: ZEROING MOLECULE #',I2.2,' AMOUNT ',        &
     &        'IN LAYER #',I3.3)
 1100 FORMAT (I3)
 !1101 FORMAT (A3)

      END SUBROUTINE ATMPTH

!     ----------------------------------------------------------------

      SUBROUTINE REFOUT( ITER, JSPEC, IREAD, LMAX, AST, APP, BENDNG )

      INTEGER (4)         :: IM, IL, LMAX, ISPEC, IREAD, ITER, JSPEC, NAERR, I, J
      REAL (8)            :: AST, APP, BENDNG, XKONST
      REAL (8), PARAMETER :: CONST = 7.341D+21
      REAL (8)            :: WETAIR(LMAX)
      REAL (8)            :: SCALES(30), XK(100)

      !print*,ITER, JSPEC, IREAD
      DO IL=1, LMAX
         WETAIR(IL) = DRAIRL(IL) + AMOUNT(1,IL)
         !print *, WETAIR(IL), DRAIRL(IL), rhosum(il) !rhosum=wetair
      ENDDO

      APP   = AST - BENDNG
      ISPEC = NINT(AST*100.)
      APPANG(JSPEC) = APP
      !PRINT *, ASTRO, AST, BENDNG, NINT(ASTRO*100.), NINT(AST*100.), ISPEC

! --- SAVE RAYTRACED ATMOSPHERES TO INTERNAL VARIABLES
! --- UPSIDE DOWN!
! --- MASS PATHS ARE TOTAL MASS PATHS IN CM*ATM
      DO IL = 1, LMAX
         CCC(JSPEC,IL)  = TBAR(LMAX-IL+1)*WETAIR(LMAX-IL+1)/CONST
         CORG(JSPEC,IL) = CCC(JSPEC,IL)
      END DO

! --- VERTICAL PATH IS IN MOLEC/CM^2, KODE 999, JSPEC = NSPEC +1
      IF( IREAD .EQ. 999 )THEN
         APPANG(JSPEC) = 0.0
         DO IL = 1, LMAX
            CCC(JSPEC,IL) = WETAIR(LMAX-IL+1)*1.0D0
            CORG(JSPEC,IL) = CCC(JSPEC,IL)
         END DO
      ENDIF
!      write(0,'(a,4i5,2f10.2)') 'refout ',iter, jspec, nspec, lmax, app, astang(jspec)

      IF( ITER .NE. 0 .AND. JSPEC .NE. 1 )RETURN




! --- FOR NOW JSPEC = 1 IS SPECIAL SHOULD HAVE MOST LAYERS - LARGEST SZA (WHEN > 90)
      !IF( JSPEC .EQ. 1 )THEN
      IF( IREAD .EQ. 999 )THEN

!print *, 'update pt '

! --- ZPT
         DO IL = 1, LMAX
            Z(IL)   = ZFIN(LMAX-IL+1)
            PMB(IL) = PBAR(LMAX-IL+1)
            P(IL)   = PMB(IL) / BAR
            T(IL)   = TBAR(LMAX-IL+1)
         ENDDO

!print *, iter
! write(25,'(a,i5,10f10.3)') 'refout ', iter, T(:lmax)
! write(25,*) 'refout ', iter, T(:lmax)

! --- MIXING RATIOS
         DO IM=1, NMOL
            DO IL=1, LMAX
               FXGAS(IM,IL) =  AMOUNT(IM,LMAX-IL+1)/WETAIR(LMAX-IL+1)
               !FXGAS(IM,IL) =  AMOUNT(IM,LMAX-IL+1)/DRAIRL(LMAX-IL+1)
               !WRITE(90,206) (AMOUNT(IM,IL)/WETAIR(IL),IL=LMAX,1,-1)
            ENDDO
         ENDDO

         IF( ITER .NE. 0 )RETURN

! --- SAVE INITIAL WEIGHTED VMR, TEMPERATURE & PRESSURE ARRAYS
         IF( .NOT. ALLOCATED( TORG ))THEN
            ALLOCATE( TORG(LMAX), PORG(LMAX), PMBORG(LMAX), FXORG(NMOL,LMAX), STAT=NAERR )
            IF( NAERR .NE. 0 )THEN
               WRITE(16, *) 'COULD NOT ALLOCATE TORG ARRAY ERROR NUMBER = ', NAERR
               WRITE( 0, *) 'COULD NOT ALLOCATE TORG ARRAY ERROR NUMBER = ', NAERR
               CALL SHUTDOWN
               STOP '3'
            ENDIF
         ENDIF

!print *, 'setting Torg etc in rayt'
         TORG(:LMAX)   = T(:LMAX)
         PORG(:LMAX)   = P(:LMAX)
         PMBORG(:LMAX) = PMB(:LMAX)
         FXORG(:NMOL,:LMAX) = FXGAS(:NMOL,:LMAX)

! --- WRITE OUT STANDARD SFIT INPUT FILES PT, MS, MIX

         IF( .NOT. F_WRTRAYTC )RETURN
! --- PT
         CALL FILEOPEN( 74, 1 )
         DO IL=LMAX, 1, -1
            IF( IL .EQ. LMAX )THEN
               WRITE(74,109) ZFIN(IL), PBAR(IL), TBAR(IL), ZFIN(IL+1)
            ELSE
               WRITE(74,109) ZFIN(IL), PBAR(IL), TBAR(IL)
            ENDIF
         ENDDO

         CALL FILECLOSE( 74, 1 )

! --- MIX
         CALL FILEOPEN( 76, 1 )
         DO IM = 1, NMOL
            IF( IM .EQ. 1 )THEN
                WRITE(76,204) HMOLS(IM), 1, LMAX, NMOL
            ELSE
                WRITE(76,204) HMOLS(IM)
            ENDIF
            WRITE(76,206) (AMOUNT(IM,IL)/DRAIRL(IL),IL=LMAX,1,-1)
            !WRITE(76,206) (AMOUNT(IM,IL)/WETAIR(IL),IL=LMAX,1,-1)
         ENDDO

         CALL FILECLOSE( 76, 1 )

      ENDIF ! jspec=1




      IF( .NOT. F_WRTRAYTC )RETURN

! --- MS & SA
      IF( IREAD .EQ. 0 )THEN
         CALL FILEOPEN( 75, 1 )
         CALL FILEOPEN( 77, 1 )
         WRITE(77,*) TRIM(TAG), ' SELECTION OF SA FOR THIS ALTITUDE GRID'
      ENDIF

      IF( IREAD .NE. 999 )THEN
         WRITE(75,103) ISPEC, LMAX, 1, AST, BENDNG, APP
! --- CONVERT FOR AIR FROM MOLCM-2 TO CM*ATM UNITS
         !WRITE(90,206) (TBAR(IL)*DRAIRL(IL)/CONST,IL=LMAX,1,-1)
         WRITE(75,206) (TBAR(IL)*WETAIR(IL)/CONST,IL=LMAX,1,-1)
      ENDIF

      WRITE(77,103) ISPEC, LMAX, 1, AST, BENDNG, APP

      IF( IREAD .EQ. 999 )THEN

         WRITE(75,103) 0, LMAX, 1, AST, BENDNG, APP
         !WRITE(90,206) (TBAR(IL)*DRAIRL(IL)/CONST,IL=LMAX,1,-1)
         WRITE(75,206) (TBAR(IL)*WETAIR(IL)/CONST,IL=LMAX,1,-1)

         WRITE(75,103) 999, LMAX, 1, AST, BENDNG, APP
         WRITE(75,206) (WETAIR(IL)*1.0D0,IL=LMAX,1,-1)

         CALL FILECLOSE( 75, 1 )

         WRITE(77,207) ZFIN(LMAX+1), (ZFIN(IL),IL=LMAX, 1, -1)

! --- POSSIBLE DIAGONAL VARIANCES
         DATA SCALES / 1., 1.5, 2., 3., 5., 7., 24*0.0 /

         DO J=4, -1, -1
            DO I=1, 6
               XKONST = SCALES(I)/(10.**J)
               !PRINT *, I, J, XKONST
               WRITE(77,209) XKONST*100., ' % / KM'
               WRITE(77,207) (XKONST/((ZFIN(IL+1)-ZFIN(IL))),IL=LMAX, 1, -1)

               WRITE(77,209) XKONST*100., ' % / SQRT(KM)'
               WRITE(77,207) (XKONST/(SQRT(ZFIN(IL+1)-ZFIN(IL))),IL=LMAX, 1, -1)

               WRITE(77,209) XKONST*100., ' % * KM'
               WRITE(77,207) (XKONST*((ZFIN(IL+1)-ZFIN(IL))),IL=LMAX, 1, -1)

               WRITE(77,209) XKONST*100., ' % * SQRT(KM)'
               WRITE(77,207) (XKONST*(SQRT(ZFIN(IL+1)-ZFIN(IL))),IL=LMAX, 1, -1)

               XK(:) = XKONST
               WRITE(77,209) XKONST*100., ' %'
               WRITE(77,207) (XK(IL),IL=1,LMAX)

            ENDDO
         ENDDO

         CALL FILECLOSE( 77, 1 )

      ENDIF

      RETURN

! 900 FORMAT( G12.4, F12.4, F12.6,G12.4)
  103 FORMAT(3I7,3F10.5)
! 105 FORMAT(F10.2,2(1PE10.3,0PF10.3))
  109 FORMAT(F10.2,1PE10.3,0PF10.3,F10.2)
  204 FORMAT(A7, 3X, 5I5 )
! 205 FORMAT(8E10.4)
  206 FORMAT(6(1PE12.4))
  207 FORMAT(6(F12.6))
  209 FORMAT(F10.4,A)

     END SUBROUTINE REFOUT

!     ----------------------------------------------------------------

      !SUBROUTINE MDLATM( ITYPE, MDL, IREAD, HSPACE, LMAX )
      SUBROUTINE MDLATM( ITYPE, MDL, IREAD, HSPACE )

!     *****************************************************************
!     THIS SUBROUTINE LOADS ONE OF THE 6 BUILT IN ATMOSPHERIC PROFILES
!     OR CALLS NSMDL TO READ IN A USER SUPPLIED PROFILE.
!     *****************************************************************

      INTEGER (4) :: I, ISPACE=0, ITYPE, IREAD, MDL !, LMAX
      REAL (8)    :: HSPACE

!     ZMDL BLANK COMMON ALTITUDES FOR LBLRTM BOUNDRIES
!     ZMAX /PARMTR/ HIGHEST LBLRTM ALT
!     ZMIN /PARMTR/ LOWEST LBLRTM ALT
!     ZPTH BLANK COMMON
!     ZST /MLATM/ ORIGINAL LBLRTM ALTITUDES

! MODELS 0 & 7 ONLY!!

! --- WE ONLY READ IN TAPE5 AND TAPE8 ONCE & RUN FOR ANY 3 OF SZA'S
      !IF( IREAD .GT. 0 ) RETURN

      IF( MDL .EQ. 0 .OR. MDL .EQ. 7 )GOTO 40

      IF( MDL .GE. 1 .AND. MDL .LE. 6 )THEN
         WRITE(16, *) "RAYTRACE: MDLATM: NO MODEL 1-6"
         WRITE( 0, *) "RAYTRACE: MDLATM: NO MODEL 1-6"
         CALL SHUTDOWN
         STOP '3'
      ENDIF

      IF( MDL .LT. 0 .OR. MDL .GT. 7 )THEN
         WRITE(16, *) "RAYTRACE: MDLATM: MODEL OUT OF RANGE"
         WRITE( 0, *) "RAYTRACE: MDLATM: MODEL OUT OF RANGE"
         CALL SHUTDOWN
         STOP '3'
      ENDIF

   !40 CALL NSMDL( IREAD, ITYPE, MDL, NOPRNT, LMAX )
   !40 CALL NSMDL( IREAD, MDL, NOPRNT, LMAX )
   40 CALL NSMDL( IREAD, MDL, NOPRNT )

      IF( NOPRNT.GT.0 )WRITE(IPR,*) " MODEL ATM TYPE : ", ITYPE

!   45   IF (IMOLDX.EQ.-99) THEN
!          IF (IMMAX.NE.IBMAX) THEN
!              WRITE(IPR,*) 'ERROR IN ATMOSPHERE SPECIFICATION:'
!              WRITE(IPR,*) '   DESIRED LEVELS MUST MATCH INPUT GRID'
!              WRITE(IPR,*) '   FOR ANALYTIC JACOBIAN CALCULATION'
!              STOP 'ERROR IN LEVEL GRID:  SEE TAPE6'
!          ENDIF
!      ENDIF

      ZMIN = ZMDL(1)
      DO 70 I = 1, IMMAX
         IF (HSPACE+0.001.GT.ZMDL(I)) ISPACE = I
   70 END DO
      IF( ZMDL(ISPACE) .LT. HSPACE )ZMDL(ISPACE) = HSPACE

      IMMAX = ISPACE
      ZMAX  = ZMDL(IMMAX)
      IMLOW = IMMAX
      !print *, 'zmax ',zmax, immax, ispace
      RETURN

!  900 FORMAT (3A8)

      END SUBROUTINE MDLATM

!     ----------------------------------------------------------------

      SUBROUTINE LNGMDL(NMOLK, NLAYK)

!     READS OLDER REFMOD FILE

      LOGICAL               :: FLAG
      INTEGER (4)           :: IUPDN, NLAYK, NMOLK, IM, IL, IGN, J, IERR
      CHARACTER  (LEN=80)   :: BUFFER
      CHARACTER  (LEN=8)    :: THISNAME
      REAL (8), ALLOCATABLE :: RIN(:), ZLNG(:), TPLNG(:,:), GLNG(:,:), YGAS(:)

      !REWIND(IRP)

      ! IUPDN = 0 ASCENDING ORDER
      ! IUPDN = 1 DECENDING ORDER

      CALL FILEOPEN( IRP, 3 )
      READ(IRP,*,ERR=200,END=201) IUPDN, NLAYK, NMOLK
      IF( IUPDN .LT. 0 .OR. IUPDN .GT. 1 )THEN
         WRITE(16,*) 'LNGMDL IUPDN OOR'
         WRITE(00,*) 'LNGMDL IUPDN OOR'
         STOP '3'
      ENDIF
      IF( NMOLK .LT. MOLTOTAL ) THEN
          WRITE(16,*) 'REFERENCE NMOLK (',NMOLK,') LESS THAN MOLTOTAL (',MOLTOTAL,')'
          WRITE(00,*) 'REFERENCE NMOLK (',NMOLK,') LESS THAN MOLTOTAL (',MOLTOTAL,')'
          STOP '3'
      ENDIF

      ALLOCATE( RIN(NLAYK), ZLNG(NLAYK), TPLNG(NLAYK,2), GLNG(NLAYK,NMOLK), YGAS(NMOLK) )
      GLNG(:,:) = 0.0D0

      IF( IUPDN .EQ. 0 )THEN
         READ(IRP,*,ERR=202,END=203) BUFFER
         !PRINT *, BUFFER
         READ(IRP,109,ERR=202,END=203) (ZLNG(IL), IL=1, NLAYK)     !KM
         READ(IRP,*,ERR=202,END=203) BUFFER
         !PRINT *, BUFFER
         READ(IRP,109,ERR=202,END=203) (TPLNG(IL,1), IL=1, NLAYK)  !MBAR
         READ(IRP,*,ERR=202,END=203) BUFFER
         !PRINT *, BUFFER
         READ(IRP,109,ERR=202,END=203) (TPLNG(IL,2), IL=1, NLAYK)  !K
         DO IM=1, NMOLK
            FLAG = .FALSE.
            READ(IRP,108,ERR=204,END=205) IGN, THISNAME, BUFFER
!            PRINT *, THISNAME, BUFFER
            READ(IRP,109,ERR=204,END=205)(RIN(IL), IL=1, NLAYK)    !MIX RATIO
            IF( ADJUSTL(TRIM(THISNAME)) .EQ. NAME(IM) .AND. ADJUSTL(TRIM(THISNAME)) .NE. 'OTHER' )THEN
               GLNG(1:NLAYK,IM) = RIN
               FLAG = .TRUE.
            ELSE
            ! --- CHECK FOR ISO SUBSTITUTION
               ! --- VMR'S ARE IN REFMOD FORM
               IF( USEISO .AND. NISOVMR .NE. NLAYK )THEN
                  WRITE(16,*) ' NISOVMR IS NOT EQUAL TO REFERENCE # LAYERS'
                  WRITE(00,*) ' NISOVMR IS NOT EQUAL TO REFERENCE # LAYERS'
                  CALL SHUTDOWN
                  STOP '3'
               ENDIF
               DO J=1, NISOSEP
                  IF( NEWID(J) .EQ. IM )THEN
                     RIN = NEWVMR(:NLAYK,J)
                     GLNG(1:NLAYK,IM) = RIN
                     !GLNG(1:NLAYK,IM) = GLNG(1:NLAYK,oldid(j))
                     FLAG = .TRUE.
                     CYCLE
                  ENDIF ! NEWID
               ENDDO ! NISO
             ENDIF ! NAME
             ! --- FORMALLY SET TO ZERO VMR'S THAT HAVE NOT BEEN READ IN
             IF( .NOT. FLAG )THEN
                PRINT *, IM, NAME(IM)
                IF( NOPRNT.GT.0 )WRITE(IPR,113) "NO PROFILE FOR GAS :", IM, NAME(IM), "SETTING TO ZEROS."
                GLNG(1:NLAYK,IM) = 0.0D0
             ENDIF
          ENDDO ! NMOLK
      ELSE
         READ(IRP,*,ERR=202,END=203) BUFFER
         !PRINT *, BUFFER
         READ(IRP,109,ERR=202,END=203) (RIN(IL), IL=1, NLAYK)      !KM
         ZLNG = DREV( RIN, NLAYK )
         READ(IRP,*,ERR=202,END=203) BUFFER
         !PRINT *, BUFFER
         READ(IRP,109,ERR=202,END=203) (RIN(IL), IL=1, NLAYK)      !MBAR
         TPLNG(1:NLAYK,1) = DREV( RIN, NLAYK )
         READ(IRP,*,ERR=202,END=203) BUFFER
         !PRINT *, BUFFER
         READ(IRP,109,ERR=202,END=203) (RIN(IL), IL=1, NLAYK)      !K
         TPLNG(1:NLAYK,2) = DREV( RIN, NLAYK )
         DO IM=1, NMOLK
            FLAG = .FALSE.
            READ(IRP,108,ERR=204,END=205) IGN, THISNAME, BUFFER
            !PRINT *, THISNAME, BUFFER
            READ(IRP,109,ERR=204,END=205)(RIN(IL), IL=1, NLAYK)    !MIX RATIO
            IF( ADJUSTL(TRIM(THISNAME)) .EQ. NAME(IM) .AND. ADJUSTL(TRIM(THISNAME)) .NE. 'OTHER' )THEN
               GLNG(1:NLAYK,IM) = DREV( RIN, NLAYK )
               FLAG = .TRUE.
            ELSE
            ! --- CHECK FOR ISO SUBSTITUTION
               ! --- VMR'S ARE IN REFMOD FORM
               IF( USEISO .AND. NISOVMR .NE. NLAYK ) THEN
                  WRITE(16,*) ' NISOVMR = ',NISOVMR,' IS NOT EQUAL TO REFERENCE # LAYERS (=',NLAYK,')'
                  WRITE(00,*) ' NISOVMR = ',NISOVMR,' IS NOT EQUAL TO REFERENCE # LAYERS (=',NLAYK,')'
                  CALL SHUTDOWN
                  STOP '3'
               ENDIF
               DO J=1, NISOSEP
                  IF( NEWID(J) .EQ. IM )THEN
                     RIN = NEWVMR(:NLAYK,J)
                     GLNG(1:NLAYK,IM) = DREV( RIN, NLAYK )
                     !GLNG(1:NLAYK,IM) = GLNG(1:NLAYK,oldid(j))
                     FLAG = .TRUE.
                     CYCLE
                  ENDIF ! NEWID
               ENDDO ! NISO
             ENDIF ! NAME
             ! --- FORMALLY SET TO ZERO VMR'S THAT HAVE NOT BEEN READ IN
             IF( .NOT. FLAG )THEN
                !PRINT *, IM, NAME(IM)
                IF( NOPRNT.GT.0 )WRITE(IPR,113) "NO PROFILE FOR GAS :", IM, NAME(IM), "SETTING TO ZEROS."
                GLNG(1:NLAYK,IM) = 0.0D0
             ENDIF
          ENDDO ! NMOLK
      ENDIF ! UNDN

      CALL FILECLOSE( IRP, 2 )

      WRITE(16,113) " NUM LAYERS FOUND IN USER MODEL    : ", NLAYK

      IF( NOPRNT.GT.0 )WRITE(IPR,114) " LANGLEY FORMATED USER MODEL ATMOSPHERE"
      IF( NOPRNT.GT.0 )WRITE(IPR,111) "ALTITUDE",  "PRESSURE",  "TEMPERAT ", (NAME(IM),IM=1,MOLTOTAL)
      DO IL=1, NLAYK
         IF( NOPRNT.GT.0 )WRITE(IPR, 110) ZLNG(IL), TPLNG(IL,1:2), (GLNG(IL,IM),IM=1,NMOLK)
      ENDDO

! --- SAVE THESE Z P T V'S FROM FILE FOR RAYTRACE
      DO IL=1, NLAYK
         ZMDL(IL) = ZLNG(IL)
         PM(IL)   = TPLNG(IL,1)
         TM(IL)   = TPLNG(IL,2)
         TM0(IL)  = TPLNG(IL,2)
         USRMIX(IL,1:NMOLK) = GLNG(IL,1:NMOLK)
      ENDDO

      IF( ALLOCATED(ZLNG) )DEALLOCATE( RIN, ZLNG, TPLNG, GLNG, YGAS, STAT = IERR )

      RETURN

 200  WRITE(16,*) "200 READ ERROR FIRST LINE : ", TFILE(72)
      WRITE(00,*) "200 READ ERROR FIRST LINE : ", TFILE(72)
      CALL SHUTDOWN
      STOP '3'
 201  WRITE(16,*) "201 EOF ERROR FIRST LINE : ", TFILE(72)
      WRITE(00,*) "201 EOF ERROR FIRST LINE : ", TFILE(72)
      CALL SHUTDOWN
      STOP '3'
 202  WRITE(16,*) "202 READ ERROR Z, P OR T : ", TFILE(72)
      WRITE(00,*) "202 READ ERROR Z, P OR T : ", TFILE(72)
      CALL SHUTDOWN
      STOP '3'
 203  WRITE(16,*) "203 EOF ERROR Z, P OR T : ", TFILE(72)
      WRITE(00,*) "203 EOF ERROR Z, P OR T : ", TFILE(72)
      CALL SHUTDOWN
      STOP '3'
 204  WRITE(16,*) "204 READ ERROR A GAS BLOCK HEADER : ", TFILE(72)
      WRITE(00,*) "204 READ ERROR A GAS BLOCK HEADER  : ", TFILE(72)
      CALL SHUTDOWN
      STOP '3'
 205  WRITE(16,*) "205 EOF ERROR A GAS BLOCK HEADER  : ", TFILE(72)
      WRITE(00,*) "205 EOF ERROR A GAS BLOCK HEADER  : ", TFILE(72)
      CALL SHUTDOWN
      STOP '3'

!101   FORMAT( F10.3,1X,1PE10.3,1X,0PF7.2,8(1PE10.3),8(/,29X,8(1PE10.3)))
!106   FORMAT( 5(E12.5,1X))
!107   FORMAT( 5X,A75)
108   FORMAT( I5,A8,A75)
109   FORMAT( 5(E12.4:,1X))
110   FORMAT( F12.2, 2G12.4, 200E12.4 )
111   FORMAT( 200A12 )
!112   FORMAT( 36X,200(8X,I4))
113   FORMAT( A, I4, 1X, A8, 1X, A)
114   FORMAT( /,A )
!201   FORMAT( 8(1PE10.3))

      END SUBROUTINE LNGMDL

!     ----------------------------------------------------------------

       FUNCTION DREV (X,N) RESULT (Y)
            INTEGER, INTENT (IN) :: N
            REAL (8), DIMENSION(N), INTENT (IN) :: X
            INTEGER              :: I
            REAL (8), DIMENSION(N) :: Y

            DO I=1,N
               Y(I)=X(N-I+1)
            END DO
            RETURN
       END FUNCTION DREV


!     ----------------------------------------------------------------

      !SUBROUTINE NSMDL (IREAD, ITYPE, MDL, NOPRNT, LMAX)
      !SUBROUTINE NSMDL (IREAD, MDL, NOPRNT, LMAX)
      SUBROUTINE NSMDL (IREAD, MDL, NOPRNT )
!     *****************************************************************
!
!
!     NOTES TO USER:
!
!     THIS SUBROUTINE IS FOR READING IN AN ATMOSPHERIC PROFILE
!     CORRESPONDING TO MODEL = 0.  THE PROFILE IS READ IN AFTER
!     CONTROL CARD 3.4
!     MDL = 0 - LBLRTM FORMATTED MODEL ON TAPE8
!     MDL = 7 - REFMOD FORMATTED MODEL ON TAPE8
!
!     CARD 3.4    IMMAX,(HMOD(I),I=1,3)
!                   (I5,3A8)
!
!             IMMAX  NUMBER OF BOUNDARIES FOR THE PROFILE
!
!             HMOD   A 24 CHARACTER HEADER DESCRIBING THE PROFILE
!
!     SEE DETAILS IN RDUNIT ON CARDS 3.5 AND 3.6.1 ... 3.6.N
!
!     *****************************************************************
!
      INTEGER (4)                :: IREAD, MDL, IMTYPE, IM, NOPRNT !, LMAX
      REAL(8), DIMENSION(KMAX+1) :: x, y, y0, b, c, d
      REAL(8), DIMENSION(MXZMD)  :: tm0
!print*, 'lmax rayt :',lmax
      !LMAX = 0
      IF( NOPRNT.GT.0 )WRITE(IPR,900) MDL
!print*, 'nsmdl',IREAD, MDL, NOPRNT, LMAX
      IF( MDL .EQ. 0 )THEN

         IF( NOPRNT .GE. 0 )WRITE (IPR,901)"READING IN LBLRTM FORMAT USER ATMOSPHERE MODEL"
         READ (IRP,905) IMTYPE, IMMAX_B, HMOD
         IMMAX = ABS(IMMAX_B)
         IMLOW = IMMAX

         IF( NOPRNT .GE. 0 )WRITE (IPR,910) IMTYPE, IMMAX, HMOD
         IF( IMMAX .GT. IMDIM )GOTO 30

         DO IM = 1, IMMAX

!     READ IN GENERIC UNITS FOR USER MODEL
            CALL RDUNIT (IM,ZMDL(IM),PM(IM),TM(IM),NMOL)

!     CONVERSION OF GENERIC UNITS TO DENSITIES FOR LBLRTM RUNS
            CALL CONVRT( ZMDL(IM), PM(IM), TM(IM), IM, NMOL, NOPRNT )
            DENW(IM) = DENM(1,IM)

         ENDDO

      ELSE

! MODEL = 7 REFMOD
        IF( IREAD .EQ. 0 )THEN
           IF( NOPRNT .GE. 0 )WRITE (IPR,901)"READING IN LANGLEY FORMAT USER ATMOSPHERE MODEL"
           CALL LNGMDL ( NMOL, IMMAX )
           tm0(:immax) = tm(:immax)

!PRINT *, IREAD, LMAX, KMAX, IMMAX
!PRINT *, ZMDL(:IMMAX)
!PRINT *,''
!PRINT *, TM(:IMMAX)
!PRINT *,''
!PRINT *, PM(:IMMAX)


        ELSE IF( IREAD .EQ. 1 )THEN

!PRINT *, IREAD, LMAX, KMAX, IMMAX
!PRINT *, ZMDL(:IMMAX)
!PRINT *, ''
!PRINT *, TM(:IMMAX)
!PRINT *,''
!PRINT *, PM(:IMMAX)

 x(1:kmax) = DREV(Zbar,KMAX)
 x(kmax+1) = zmdl(immax)

 y(1:kmax) = DREV(T,KMAX)
 y(kmax+1) = tm0(immax)

 y0(1:kmax) = DREV(TORG,KMAX)
 y0(kmax+1) = tm0(immax)

!PRINT *,''
!write(*, '(3f10.3)') (x(im), y(im), y0(im), im=1, kmax+1)

            ! change model T to perturbed T
            CALL spline (KMAX, x, y, b, c, d)

            DO IM = 1, IMMAX

               TM(IM) = seval (KMAX, ZMDL(IM), x, y, b, c, d)

            ENDDO

!PRINT *,''
!write(*, '(3f10.3)') (Zmdl(im), Tm(im), tm0(im), im=1, immax)


        ELSE
           WRITE(16,*) 'NSMDL IREAD OOR'
           WRITE(00,*) 'NSMDL IREAD OOR'
           CALL SHUTDOWN
           STOP '3'
        ENDIF

        IF( NOPRNT .GE. 0 )THEN
           WRITE(IPR,901)"CONVERT UNITS AND CALCULATE AMOUNTS & RH"
           WRITE(IPR,904)
        ENDIF

! --- LOOP OVER LEVELS IN INPUT MODEL
        DO IM = 1, IMMAX

           JUNIT(1:NMOL) = 18
           WMOL(1:NMOL)  = USRMIX(IM,1:NMOL)
           CALL CONVRT( ZMDL(IM), PM(IM), TM(IM), IM, NMOL, NOPRNT )
           DENW(IM)      = DENM(1,IM)

        ENDDO

      ENDIF
!print*,'immaxb ',immax_B
      IF (IMMAX_B .LT. 0) THEN
         CALL CMPALT (IMMAX,PM,TM,DENW,ZMDL(1),ZMDL)
      ENDIF

      DO IM = 2, IMMAX
         IF (ZMDL(IM) .LE. ZMDL(IM-1)) GOTO 35
      ENDDO

      RETURN
!
   30 CONTINUE
      WRITE (0,915) IMMAX,IMDIM
      IF (NOPRNT .GE.0) WRITE (IPR,915) IMMAX,IMDIM

      WRITE(16,*) ' LEVEL ERROR IN NSMDL '
      WRITE(00,*) ' LEVEL ERROR IN NSMDL '
      CALL SHUTDOWN
      STOP '3'

   35 CONTINUE

      IF (NOPRNT.GE.0) WRITE (IPR,920) IM,IM+1,ZMDL(IM),ZMDL(IM+1)

      WRITE(16,*) 'INPUT ALTITUDES NOT IN ASCENDING ORDER'
      WRITE(00,*) 'INPUT ALTITUDES NOT IN ASCENDING ORDER'
      CALL SHUTDOWN
      STOP '3'

!  900 FORMAT (///,' READING IN USER SUPPLIED MODEL ATMOSPHERE',/)
 900  FORMAT(/,' NSMDL : ATMOSPHERE MODEL TYPE : ', I5 )
 901  FORMAT (/,A,/)
 904  FORMAT( /, " ALTITUDE    PRESSURE  TEMPERATURE      RH   WATER_DENS  W_SAT_DENS")
 905  FORMAT (2I5,3A8)
 910  FORMAT (//,10X,'IMTYPE   = ',I5,/,                                &
     &           10X,'IMMAX    = ',I5,/,10X,'PROFILE = ',3A8)
 915  FORMAT (/,' NUMBER OF PROFILE LEVELS IMMAX = ',I5,                &
     &        ' EXCEEDS THE MAXIMUM ALLOWED = ',I5)
 920  FORMAT (///,' ERROR: INPUT ALTITUDES FOR LBLRTM LAYERS ',         &
     &        'ARE NOT IN ASCENDING ORDER',//,5X,                       &
     &        ' ZMDL AT GRID PT I,I+1 = ',I5,I5,/,(2F10.4))
!
      END SUBROUTINE NSMDL
!
!     ----------------------------------------------------------------
!
      SUBROUTINE HEADPR (IPR)

      INTEGER (4) :: IPR, I
!     SUBROUTINE TO WRITE HEADER INFORMATION FOR MODEL  0
!
      WRITE (IPR,900)
      WRITE (IPR,905)
      WRITE (IPR,910) (I,HMOLS(I),I=1,MXMOL)
      WRITE (IPR,915)
!
      RETURN
!
  900 FORMAT (/,'  THE USER HAS ELECTED TO PROVIDE THE REQUIRED',/,     &
     &        '  MODEL ATMOSPHERE SPECIFICATIONS.',/,/,                 &
     &        '  SEE DOCUMENTATION OR "SUBROUTINE RDUNIT" FOR ',/,      &
     &        '  ADDITIONAL INFORMATION.',//)
  905 FORMAT ('  USER OPTIONS FOR PRESSURE AND TEMPERATURE ',//,        &
     &        '               JCHAR   JUNIT ',//,                       &
     &        '    PRESSURE " ",A      10    PRESSURE IN (MB)',/,       &
     &        '                 B      11       "     "  (ATM)',/,      &
     &        '                 C      12       "     "  (TORR)',/,     &
     &        '                1-6    1-6    DEFAULT TO SPECIFIED',     &
     &        ' MODEL ATMOSPHERE',//,                                   &
     &        '    TEMP     " ",A      10    AMBIENT TEMP IN DEG(K)',/, &
     &        '                 B      11       "     "   "   " (C)',/, &
     &        '                1-6    1-6    DEFAULT TO SPECIFIED',     &
     &        ' MODEL ATMOSPHERE',//)
  910 FORMAT (/,' AVAILABLE     ',7('(',I2,')',A8),/,' MOL. SPECIES',   &
     &        (T16,7('(',I2,')',A8)))
  915 FORMAT (/,'  POTENTIAL CHOICE OF UNITS FOR ABOVE SPECIES',/,      &
     &        ' JCHAR = " ",A    - VOLUME MIXING RATIO (PPMV)',/,       &
     &        '       = B        - NUMBER DENSITY (CM-3)',/,            &
     &        '       = C        - MASS MIXING RATIO (GM/KG)',/,        &
     &        '       = D        - MASS DENSITY (GM M-3)',/,            &
     &        '       = E        - PARTIAL PRESSURE (MB)',/,            &
     &        '       = F        - DEW POINT TEMP (K) * H2O ONLY *',/,  &
     &        '       = G        - DEW POINT TEMP (C) * H2O ONLY *',/,  &
     &        '       = H        - RELATIVE HUMIDITY (PERCENT) ',       &
     &                             '*H2O ONLY*',/,                      &
     &        '       = I        - AVAILABLE FOR USER DEFINITION',/,    &
     &        '       = 1-6      - DEFAULT TO SPECIFIED MODEL ',        &
     &        'ATMOSPHERE',/,' JCHAR MUST BE LESS THAN "J"',/)
!
      END SUBROUTINE HEADPR
!
!     ----------------------------------------------------------------
!
      SUBROUTINE RDUNIT (IM,ZMDL1,PM1,TM1,NMOL)
!
!     *******************************************************
!
!       SUBROUTINE DESIGNED TO READ NEW MOLECULAR DATA INPUT
!        PARAMETERS - JCHAR = INPUT KEY (SEE BELOW)
!                     WMOL  = INPUT VALUE FOR LAYER
!
!     ***  ROUTINE ALSO ACCEPTS VARIABLE UNITS ON PRESS AND TEMP
!     ***  THE ASSOCIATED 'JUNIT' DEFINITIONS ARE CONTAINED IN
!               JUNITP, AND JUNITT
!          SEE INPUT KEY BELOW
!
!
!       NMOL = NUMBER OF MOLECULAR SPECIES TO BE CONSIDERED
!               (ORDER IS THAT OF AFGL LINE PARAMETER TAPE)
!
!     FOR MOLECULAR SPECIES ONLY
!
!       JCHAR   JUNIT
!
!     " ",A      10    VOLUME MIXING RATIO (PPMV)
!         B      11    NUMBER DENSITY (CM-3)
!         C      12    MASS MIXING RATIO (GM(K)/KG(AIR))
!         D      13    MASS DENSITY (GM M-3)
!         E      14    PARTIAL PRESSURE (MB)
!         F      15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY
!         G      16     "    "     "  (TD IN T(C)) - H2O ONLY
!         H      17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY
!         I      18    AVAILABLE FOR USER DEFINITION
!        1-6    1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE
!                                                (SEE KEY BELOW)
!
!     ****************************************************************
!     ****************************************************************
!
!     ***** OTHER 'JCHAR' SPECIFICATIONS - JCHARP,JCHART
!
!       JCHAR   JUNIT
!
!      " ",A     10    PRESSURE IN (MB)
!          B     11       "     "  (ATM)
!          C     12       "     "  (TORR)
!         1-6   1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE
!
!      " ",A     10    AMBIENT TEMPERATURE IN DEG(K)
!          B     11       "         "       "  " (C)
!          C     12       "         "       "  " (F)
!         1-6   1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE
!
!     ***** DEFINITION OF "DEFAULT" CHOICES FOR PROFILE SELECTION *****
!
!      FOR THE USER WHO WISHES TO ENTER ONLY SELECTED ORIGINAL
!      VERTICAL PROFILES AND WANTS STANDARD ATMOSPHERE SPECIFICATIONS
!      FOR THE OTHERS, THE FOLLOWING OPTION IS AVAILABLE
!
!     *** JCHAR(P,T OR K) MUST = 1-6 (AS ABOVE)
!
!      FOR MOLECULES 8-35, ONLY US STD PROFILES ARE AVIALABLE
!      THEREFORE, WHEN  'JCHAR(K) = 1-5', JCHAR(K) WILL BE RESET TO 6
!
!     *************************************************************

      CHARACTER (LEN=1) :: JCHAR(MXMOL), JCHARP, JCHART, JLONG

      INTEGER (4)       :: IM, NMOL, JUNITP, JUNITT, K
      INTEGER (4)       :: I
      INTEGER (4), DIMENSION(MXMOL)  ::JOLD

      REAL (8)          :: ZMDL1, PM1, TM1, C1, C2, C3

      DATA JOLD / MXMOL*99 /
      DATA C1 / 18.9766 /,C2 / -14.9595 /,C3 / -2.4388 /

      IF (IM.EQ.0 .AND. NOPRNT.GE.0) CALL HEADPR (IPR)
!
!     *********************************************************
!
!     INPUT READ FOR 'MODEL = 0", I.E. USER-SUPPLIED VERITCAL
!
!     **********************************************************
!
      READ (IRP,900) ZMDL1,PM1,TM1,JCHARP,JCHART,JLONG,(JCHAR(K),K=1,MXMOL)

      JUNITP = JOU(JCHARP)
      JUNITT = JOU(JCHART)

      DO 10 K = 1, NMOL
         JUNIT(K) = JOU(JCHAR(K))
   10 END DO

!     READ IN MOLECLAR INFORMATION AT E15.8 FORMAT FOR FLAG JLONG='L'
      IF (JLONG.EQ.'L') THEN
         READ (IRP,906) (WMOL(K),K=1,NMOL)
      ELSEIF (JLONG.EQ.' ') THEN
         READ (IRP,905) (WMOL(K),K=1,NMOL)
      ELSE
         WRITE(16,*) 'INVALID VALUE FOR JLONG ON RECORD 3.5: ',JLONG
         WRITE(00,*) 'INVALID VALUE FOR JLONG ON RECORD 3.5: ',JLONG
         CALL SHUTDOWN
         STOP '3'
      ENDIF
      IF (IM.EQ.0 .AND. NOPRNT.GE.0) WRITE (IPR,910)

      IF( NOPRNT .GE. 0 )THEN
      IF (JLONG.EQ.'L') THEN
         WRITE (IPR,916) IM,ZMDL,JCHARP,PM,JCHART,TM,                   &
     &                   (K,JCHAR(K),WMOL(K),K=1,NMOL)
      ELSE
         WRITE (IPR,915) IM,ZMDL,JCHARP,PM,JCHART,TM,                   &
     &                   (K,JCHAR(K),WMOL(K),K=1,NMOL)
      ENDIF
      ENDIF

      DO 20 I = 1, NMOL
         JOLD(I) = JUNIT(I)
   20 END DO

      CALL CHECK (PM1,JUNITP,1)
      CALL CHECK (TM1,JUNITT,2)

      RETURN
!
  900 FORMAT (3E10.3,5X,2A1,1X,A1,1X,38A1)
  905 FORMAT (8E10.3)
  906 FORMAT (8E15.8)
  910 FORMAT (//,'  ECHO INPUT PARAMETERS FOR USER PROVIDED MODEL',/,   &
     &        '0   (P : UNIT)=   ',5X,'(T : UNIT)=   ',5X,              &
     &        '(MOLECULE NUMBER : UNIT)=   ')
  915 FORMAT ('0',I4,1X,'(ALT:KM)=',F7.3,4X,'(P:',A1,')=',G11.5,4X,     &
     &        '(T:',A1,')=',F8.3,/,(5X,7(' (',I2,':',A1,')=',1PE10.3)))
  916 FORMAT ('0',I4,1X,'(ALT:KM)=',F7.3,4X,'(P:',A1,')=',G11.5,4X,     &
     &        '(T:',A1,')=',F8.3,/,(5X,7(' (',I2,':',A1,')=',1PE15.8)))
!
      END SUBROUTINE RDUNIT

!
!     ----------------------------------------------------------------
!
      INTEGER (4) FUNCTION JOU (CHAR)

      INTEGER (4)        :: I, INDX, INDX1(22)
      CHARACTER (LEN=1)  :: CHAR, HOLVEC(22)

      DATA (HOLVEC(I),I=1,22) /                                         &
     &                '1','2','3','4','5','6','0','0','0','0',' ','A',  &
     &                'B','C','D','E','F','G','H','I','J','K'/

      DATA (INDX1(I),I=1,22) /                                          &
     &                  1,  2,  3,  4,  5,  6,  0,  0,  0,  0, 10, 10,  &
     &                 11, 12, 13, 14, 15, 16, 17, 18, 19, 20/

      INDX = 0
      DO 10 I = 1, 22
         IF (HOLVEC(I).NE.CHAR) CYCLE !GO TO 10
         INDX = INDX1(I)
         GO TO 20
   10 END DO
   20 IF (INDX.EQ.0) THEN
         WRITE(16,900) CHAR
         WRITE(00,900) CHAR
         CALL SHUTDOWN
         STOP '3'
      ENDIF

      JOU = INDX

      RETURN

  900 FORMAT ('0 INVALID PARAMETER :',2X,A1)

      END FUNCTION JOU

!     ----------------------------------------------------------------

      SUBROUTINE CHECK (A,IA,KEY)

!      UNITS CONVERSION FOR P AND T
!
!     A = P OR T     AND  IA =JUNITP(I.E. MB,ATM,TORR)
!                            =JUNITT(I.E. DEG K OR C)
!                            =JUNITR(I.E. KM,M,OR CM)

      REAL (8)      :: A, PTORR
      INTEGER (4)   :: IA, KEY

      DATA PTORR / 760. /
!
      WRITE(42,*) 'CHECK'
      IF (IA.LE.10) RETURN
!
      GO TO (10,20,30) KEY
!
!     PRESSURE CONVERSIONS
!
   10 IF (IA.EQ.11) THEN
         A = A * PZERO
         RETURN
      ELSEIF (IA.EQ.12) THEN
         A = A*PZERO / PTORR
         RETURN
      ELSE
         WRITE(16,*) ' CHECK(P)'
         WRITE(00,*) ' CHECK(P)'
         CALL SHUTDOWN
         STOP '3'
      ENDIF
!
!     TEMPERATURE COMVERSIONS
!
   20 IF (IA.LE.11) THEN
         A = A + TZERO
         RETURN
      ELSE
         WRITE(16,*) ' CHECK(T)'
         WRITE(00,*) ' CHECK(T)'
         CALL SHUTDOWN
         STOP '3'
      ENDIF
!
!      RANGE CONVERSIONS
!
   30 IF (IA.EQ.11) THEN
         A = A/1.D3
         RETURN
      ELSEIF (IA.EQ.12) THEN
         A = A/1.D5
         RETURN
      ELSE
         WRITE(16,*) ' CHECK(R)'
         WRITE(00,*) ' CHECK(R)'
         CALL SHUTDOWN
         STOP '3'
      ENDIF
!
      END SUBROUTINE CHECK
!
!     ----------------------------------------------------------------
!

      SUBROUTINE CONVRT( Z, P, T, IM, NMOL, NOPRNT )
!
!*************************************************************
!
!        WRITTEN APR, 1985 TO ACCOMMODATE 'JCHAR' DEFINITIONS FOR
!        UNIFORM DATA INPUT -
!
!      JCHAR    JUNIT
!
!    " ",A       10    VOLUME MIXING RATIO (PPMV)
!        B       11    NUMBER DENSITY (CM-3)
!        C       12    MASS MIXING RATIO (GM(K)/KG(AIR))
!        D       13    MASS DENSITY (GM M-3)
!        E       14    PARTIAL PRESSURE (MB)
!        F       15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY
!        G       16     "    "     "  (TD IN T(C)) - H2O ONLY
!        H       17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY
!        I       18    MIXING RATIO
!        J       19    REQUEST DEFAULT TO SPECIFIED MODEL ATMOSPHERE
!
!***************************************************************

      REAL (8)    :: Z, P, T, RHOAIR, C1, C2, C3, B, R
      INTEGER (4) :: IM, NMOL, NOPRNT, K

      DATA C1 / 18.9766 /,C2 / -14.9595 /,C3 / -2.4388 /

      RHOAIR = ALOSMT*(P/PZERO)*(TZERO/T)

!     GET WATER VAPOR DENSITY
      CALL WATVAP ( Z, P, T, IM, JUNIT(1), WMOL(1), DENM(1,IM), NOPRNT )

!     DETERMINE DENSITY OF DRY AIR
      DRYAIR(IM) = RHOAIR - DENM(1,IM)

!     LOOP THROUGH OTHER MOLECULES
      DO 70 K=2,NMOL
         ! --- CHECK FOR 0 AIRMASS IE GASES NOT IN MOLCPARAMS.F90
         IF( AMWT(K) .LT. 0.1 ) THEN
            B = 0.0
            R = 0.0
         ELSE
            B = AVOGAD/AMWT(K)
            R = AIRMWT/AMWT(K)
         ENDIF
         DENM(K,IM) = 0.0

         IF (JUNIT(K).GT.10) GO TO 20
!
!     GIVEN VOL. MIXING RATIO PPM
!
         DENM(K,IM) = WMOL(K)*DRYAIR(IM)*1.0D-6
         CYCLE !GO TO 70
   20    IF (JUNIT(K).NE.11) GO TO 30
!
!     GIVEN NUMBER DENSITY (CM-3)
!
         DENM(K,IM) = WMOL(K)
         CYCLE !GO TO 70
   30    CONTINUE
         IF (JUNIT(K).NE.12) GO TO 40
!
!     GIVEN MASS MIXING RATIO (GM KG-1)
!
         DENM(K,IM) = R*WMOL(K)*1.0D-3*DRYAIR(IM)
         CYCLE !GO TO 70
   40    CONTINUE
         IF (JUNIT(K).NE.13) GO TO 50
!
!     GIVEN MASS DENSITY (GM M-3)
!
         DENM(K,IM) = B*WMOL(K)*1.0D-6
         CYCLE !GO TO 70
   50    CONTINUE
         IF (JUNIT(K).NE.14) GO TO 60
!
!     GIVEN PARTIAL PRESSURE (MB)
!
         DENM(K,IM) = ALOSMT*(WMOL(K)/PZERO)*(TZERO/T)
         CYCLE !GO TO 70
   60    CONTINUE
         IF (JUNIT(K).NE.18) GO TO 61
!
!     JUNIT(18) AVAILABLE FOR USER DEFINITION HERE
!     GIVEN VOL. MIXING RATIO
!
         DENM(K,IM) = WMOL(K)*DRYAIR(IM)

   61    CONTINUE
         IF (JUNIT(K).GT.14 .AND. JUNIT(K) .NE. 18) THEN
            WRITE(16,900) K,JUNIT(K)
            WRITE(00,900) K,JUNIT(K)
            CALL SHUTDOWN
            STOP '3'
         ENDIF

   70 END DO
!
  900 FORMAT (/,'   **** ERROR IN CONVRT ****, JUNIT(',I5,') = ',I5)
!
      RETURN
!
      END SUBROUTINE CONVRT
!
!     ----------------------------------------------------------------
!
      SUBROUTINE WATVAP( Z, P, T, IM, IJUNIT, WMOL1, DENNUM, NOPRNT )
!
!**********************************************************************
!
!        WRITTEN APR, 1985 TO ACCOMMODATE 'JCHAR' DEFINITIONS FOR
!        UNIFORM DATA INPUT -
!
!     JCHAR    JUNIT
!
!    " ",A       10    VOLUME MIXING RATIO (PPMV)
!        B       11    NUMBER DENSITY (CM-3)
!        C       12    MASS MIXING RATIO (GM(K)/KG(AIR))
!        D       13    MASS DENSITY (GM M-3)
!        E       14    PARTIAL PRESSURE (MB)
!        F       15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY
!        G       16     "    "     "  (TD IN T(C)) - H2O ONLY
!        H       17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY
!        I       18    MIXING RATIO
!        J       19    REQUEST DEFAULT TO SPECIFIED MODEL ATMOSPHERE
!
!     THIS SUBROUTINE COMPUTES THE WATERVAPOR NUMBER DENSITY (MOL CM-3)
!     GIVE HUMIDITY  # TD = DEW POINT TEMP(K,C), RH = RELATIVE
!     (PERCENT), PPH2O = WATER VAPOR PARTIAL PRESSURE (MB), DENH2O =
!     WATER VAPOR MASS DENSITY (GM M-3),AMSMIX = MASS MIXING RATIO
!     (GM/KG).
!                     THE FUNCTION DENSAT FOR THE SATURATION
!     WATER VAPOR DENSITY OVER WATER IS ACCURATE TO BETTER THAN 1
!     PERCENT FROM -50 TO +50 DEG C. (SEE THE LOWTRAN3 OR 5 REPORT)
!
!       'JUNIT' GOVERNS CHOICE OF UNITS -
!
!**********************************************************************

      REAL (8)       :: Z, P, T, RHOAIR, A, B, R, WMOL1, DENNUM, C1, C2, C3
      REAL (8)       :: DENSAT, ATEMP, ATD, DENST, RHP
      INTEGER (4)    :: IM, IJUNIT, NOPRNT

      DATA C1 / 18.9766D0 /,C2 / -14.9595D0 /,C3 / -2.4388D0 /

      DENSAT(ATEMP) = ATEMP*B*EXP(C1+C2*ATEMP+C3*ATEMP**2)*1.0D-6

      RHOAIR = ALOSMT*(P/PZERO)*(TZERO/T)
      A = TZERO/T
      B = AVOGAD/AMWT(1)
      R = AIRMWT/AMWT(1)

      IF (IJUNIT.NE.10) GO TO 10

!     GIVEN VOL. MIXING RATIO
!     CONVERT USING DENSITY OF DRY AIR.

      WMOL1 = WMOL1*1.D-06
      DENNUM = (WMOL1/(1.+WMOL1))*RHOAIR
      GO TO 90
   10 IF (IJUNIT.NE.11) GO TO 20
!
!     GIVEN NUMBER DENSITY (CM-3)
!
      DENNUM = WMOL1
      GO TO 90
   20 CONTINUE
      IF (IJUNIT.NE.12) GO TO 30
!
!     GIVEN MASS MIXING RATIO (GM KG-1)

!     CONVERT USING DENSITY OF DRY AIR.  THE FOLLOWING QUADRATIC IS
!
      WMOL1 = WMOL1*R*1.0D-3
      DENNUM = (WMOL1/(1.+WMOL1))*RHOAIR
      GO TO 90
   30 CONTINUE
      IF (IJUNIT.NE.13) GO TO 40
!
!     GIVEN MASS DENSITY (GM M-3)
!
      DENNUM = B*WMOL1*1.0D-6
      GO TO 90
   40 CONTINUE
      IF (IJUNIT.NE.14) GO TO 50
!
!     GIVEN WATER VAPOR PARTIAL PRESSURE (MB)
!
      DENNUM = ALOSMT*(WMOL1/PZERO)*(TZERO/T)
      GO TO 90
   50 CONTINUE
      IF (IJUNIT.NE.15) GO TO 60
!
!     GIVEN DEWPOINT (DEG K)
!
      ATD = TZERO/(WMOL1)
      DENNUM = DENSAT(ATD)*(WMOL1)/T
      GO TO 90
   60 CONTINUE
      IF (IJUNIT.NE.16) GO TO 70
!
!     GIVEN DEWPOINT (DEG C)
!
      ATD = TZERO/(TZERO+WMOL1)
      DENNUM = DENSAT(ATD)*(TZERO+WMOL1)/T
      GO TO 90
   70 CONTINUE
      IF (IJUNIT.NE.17) GO TO 71
!
!     GIVEN RELATIVE HUMIDITY (PERCENT)
!
      DENNUM = DENSAT(A)*(WMOL1/100.0D0)
      GO TO 90

   71 IF (IJUNIT.NE.18) GO TO 80
!
!     JUNIT(18) AVAILABLE FOR USER DEFINITION HERE
!     GIVEN VOL. MIXING RATIO
!
      DENNUM = (WMOL1/(1.+WMOL1))*RHOAIR
      GOTO 90

   80 WRITE (0,900) JUNIT
      WRITE(16,*) 'WATVAP JUNIT OOR'
      WRITE(00,*) 'WATVAP JUNIT OOR'
      CALL SHUTDOWN
      STOP '3'

   90 CONTINUE
      DENST = DENSAT(A)
      RHP = 100.0D0*(DENNUM/DENST)
      RELHUM(IM) = RHP
      IF (NOPRNT .GE. 0) WRITE (IPR,905) Z, P, T, RHP, DENNUM, DENST
      IF (RHP.LE.100.0) GO TO 100
      IF (NOPRNT .GE. 0) WRITE (IPR,910) RHP
  100 CONTINUE

      RETURN

  900 FORMAT (/,'  **** ERROR IN WATVAP ****, JUNIT = ',I5)
! 904  FORMAT( /, "ALTITUDE   PRESSURE  TEMPERATE       RH ")
  905 FORMAT (3G12.4,1X,F6.2,3G12.4)
  910 FORMAT (/,' ****** WARNING (FROM WATVAP) # RELATIVE HUMIDTY = ',  &
     &        G10.3,' IS GREATER THAN 100 PERCENT')
!
      END SUBROUTINE WATVAP
!
!     ----------------------------------------------------------------
!
      SUBROUTINE FSCGEO (H1,H2,ANGLE,RANGE,BETA,ITYPE,LEN,HMIN,PHI,IERROR,HOBS)
!
!     -------------------------------------------------------------
!     THIS ROUTINE WAS MODIFIED FOR LBLRTM TO REFLECT CHANGES
!     IMPLEMENTED IN MODTRAN TO SOLVE PROBLEMS WITH INCONSISTENT
!     PATH PARAMETERS.
!     IT WAS ALSO MODIFIED TO ELIMINATE GOTO STATEMENTS IN ORDER TO
!     MAKE THE PROGRAM EASIER TO UNDERSTAND.
!     THESE CHANGES WERE OBTAINED FROM H. SNELL (MARCH, 1996).
!     -------------------------------------------------------------
!
!     *****************************************************************
!     FSCGEO INTERPRETS THE ALLOWABLE COMBINATIONS OF INPUT PATH
!     PARAMETERS INTO THE STANDARD SET H1,H2,ANGLE,PHI,HMIN, AND LEN.
!     THE ALLOWABLE COMBINATIONS OF INPUT PARAMETERS ARE-
!      FOR ITYPE = 2  (SLANT PATH H1 TO H2)
!       A. H1, H2, AND ANGLE,
!       B. H1, ANGLE, AND RANGE,
!       C. H1, H2, AND RANGE,
!       D. H1, H2, AND BETA -
!      FOR ITYPE = 3 (SLANT PATH H1 TO SPACE, H2 = ZMAX(=100 KM,M=1 TO 6
!       A. H1 AND ANGLE,
!       B. H1 AND HMIN (INPUT AS H2).
!     THE SUBROUTINE ALSO DETECTS BAD INPUT (IMPOSSIBLE GEOMETRY),
!     ITYPE = 2 CASES WHICH INTERSECT THE EARTH, AND RETURNS THESE
!     CASES WITH ERROR FLAGS.
!     THE SUBROUTINE FNDHMN IS CALLED TO CALCULATE HMIN, THE MINIMUM
!     HEIGHT ALONG THE PATH, AND PHI, THE ZENITH ANGLE AT H2, USING THE
!     ATMOSPHERIC PROFILE STORED IN /MDATA/
!     *****************************************************************
!
      INTEGER (4)  :: ITYPE, LEN, IERROR, ITER, ISELCT

      REAL (8)     :: H1, H2, ANGLE, RANGE, BETA, HMIN, PHI, HOBS, H2ST
      REAL (8)     :: ZARG2, ZARG3, ERARG2, ERARG3, RADCONV, SINPHI, SINANGLE
      REAL (8)     :: TOA_ANG, SINTOA_SAT, SINTOA, DIFFANGLE, R1, R2, H_TOA
      REAL (8)     :: TOA_SAT

      ITER = 0

!print*,' fscgeo top ', h1, h2, hmin, itype, angle, zmax

!     CHECK FOR ERROR

      IF ((ITYPE.NE.3).AND.(ITYPE.NE.2)) GOTO 90
!
      IF (ITYPE.EQ.3) THEN
!
!     SLANT PATH TO SPACE
!     NOTE: IF BOTH HMIN AND ANGLE ARE ZERO, THEN ANGLE IS
!           ASSUMED SPECIFIED
!
          IF (H2.EQ.0) THEN
!
!             CASE 3A: H1,SPACE,ANGLE
!
              IF( NOPRNT.GE.0)WRITE (IPR,900)
              H2 = ZMAX
              CALL FNDHMN (H1,ANGLE,H2,LEN,ITER,HMIN,PHI,IERROR)

          ELSE
!
!             CASE 3B: H1,HMIN,SPACE
!
              IF( NOPRNT.GE.0 )WRITE (IPR,905)
              HMIN = H2
              H2 = ZMAX
              IF (H1.LT.HMIN) GO TO 80
              CALL FNDHMN (HMIN,90.0D0,H1,LEN,ITER,HMIN,ANGLE,IERROR)
              CALL FNDHMN (HMIN,90.0D0,H2,LEN,ITER,HMIN,PHI,IERROR)
              IF (HMIN.LT.H1) LEN = 1
          ENDIF
      ENDIF
!
      IF (ITYPE.EQ.2) THEN
!
!       ASSIGN THE VARIABLE ISELCT TO THE FOLLOWING CASES
!       (DEPENDING ON INPUT PARAMETERS):
!
!       -----------------------------------------------
!       H1   H2   ANGLE  RANGE  BETA  =>   CASE  ISELCT
!       -----------------------------------------------
!       X    X      X                       2A     21
!       X           X      X                2B     22
!       X    X             X                2C     23
!       X    X                   X          2D     24
!       -----------------------------------------------
!
         IF (RANGE.GT.0.0) THEN
!
!           MUST BE CASE 2B OR CASE 2C
!
            IF (H2.GT.0.0) THEN
!
!              CASE 2C
!
               ISELCT=23
            ELSEIF (ANGLE.EQ.0.0) THEN
               IF( NOPRNT.GE.0 )WRITE(IPR,1000)
               WRITE(*,1000)
               ISELCT=23
            ELSE
!
!              CASE 2B
!
               ISELCT=22
            ENDIF
         ELSEIF (BETA.GT.0.0) THEN
!
!           CASE 2D (BETA CANNOT BE ZERO)
!
            ISELCT=24
         ELSE
!
!           CASE 2A, SINCE RANGE AND BETA ARE BOTH ZERO
!
            ISELCT=21
         ENDIF
!
         IF (ISELCT.EQ.21) THEN
!
!           CASE 2A: H1, H2, ANGLE
!
            IF (NOPRNT .GE.0) WRITE (IPR,910)
            IF (H1.GE.H2.AND.ANGLE.LE.90.0) GO TO 110
            IF (H1.EQ.0.0.AND.ANGLE.GT.90.0) GO TO 120
            IF (H2.LT.H1.AND.ANGLE.GT.90.0.AND.NOPRNT.GE.0) WRITE (IPR,915) LEN
            H2ST = H2
            CALL FNDHMN (H1,ANGLE,H2,LEN,ITER,HMIN,PHI,IERROR) !here

            IF (H2.NE.H2ST) GO TO 120
         ENDIF
!
         IF (ISELCT.EQ.22) THEN
!
!           CASE 2B: H1, ANGLE, RANGE
!           ASSUME REFRACTION
!
            IF (NOPRNT .GE.0) WRITE (IPR,920)
            CALL NEWH2(H1,H2,ANGLE,RANGE,BETA,LEN,HMIN,PHI)
         ENDIF
!
         IF (ISELCT.EQ.23) THEN
!
!           CASE 2C: H1, H2, RANGE
!
            IF (NOPRNT .GE.0) WRITE (IPR,930)
            IF (ABS(H1-H2).GT.RANGE) GO TO 100
            R1 = H1+RE
            R2 = H2+RE
!
            ZARG2 = (H1**2-H2**2+RANGE**2+2.0D0*RE*(H1-H2))/(2.0D0*R1*RANGE)
            ERARG2 = ABS(ZARG2)-1.0
            IF ((ERARG2.LE.1.0E-6).AND.(ERARG2.GE.0.0)) THEN
               IF (ZARG2.LT.0.0) THEN
                  ZARG2 = -1.0
               ELSE
                  ZARG2 = 1.0
               ENDIF
            ENDIF
            ANGLE = 180.0-ACOS(ZARG2)*DEG
            ZARG3 = (H2**2-H1**2+RANGE**2+2.0D0*RE*(H2-H1))/(2.0D0*R2*RANGE)
            ERARG3 = ABS(ZARG3)-1.0
            IF ((ERARG3.LE.1.0E-6).AND.(ERARG3.GE.0.0)) THEN
               IF (ZARG3.LT.0.0) THEN
                 ZARG3 = -1.0
               ELSE
                 ZARG3 = 1.0
               ENDIF
            ENDIF
            PHI = 180.0-ACOS(ZARG3)*DEG
            BETA = PHI+ANGLE-180.
!
            IF (RANGE.GT.2.0.AND.BETA.GT.0) THEN
               CALL FDBETA (H1,H2,BETA,ANGLE,PHI,LEN,HMIN,IERROR)
            ELSE
               LEN = 0
               IF (ANGLE.GT.90.0.AND.PHI.GT.90.0) LEN = 1
               CALL FNDHMN (H1,ANGLE,H2,LEN,ITER,HMIN,PHI,IERROR)
            ENDIF
         ENDIF
!
         IF (ISELCT.EQ.24) THEN
!
!        CASE 2D: H1, H2, BETA
!
            CALL FDBETA (H1,H2,BETA,ANGLE,PHI,LEN,HMIN,IERROR)
         ENDIF
      ENDIF
!

!     END OF ALLOWED CASES
!
!     TEST IERROR AND RECHECK LEN
!
      IF (IERROR.NE.0) RETURN
      LEN = 0
      IF (HMIN .LT. MIN(H1,H2)) LEN = 1
!
!     REDUCE PATH ENDPOINTS ABOVE ZMAX TO ZMAX
!
      IF (HMIN.GE.ZMAX) GO TO 130
      IF (H1.GT.ZMAX.OR.H2.GT.ZMAX) CALL REDUCE (H1,H2,ANGLE,PHI,ITER)

!
!     AT THIS POINT THE FOLLOWING PARAMETERS ARE DEFINED-
!         H1,H2,ANGLE,PHI,HMIN,LEN
!
!     CALCULATE SIN(PHI) AND SIN(ANGLE) AND OUTPUT
!
      RADCONV  = 2.0D0 * PI / 360.0D0
      SINPHI   = SIN(RADCONV*PHI)
      SINANGLE = SIN(RADCONV*ANGLE)

      IF (NOPRNT .GE. 0)WRITE (IPR,935) H1,H2,ANGLE,SINANGLE,PHI,SINPHI,HMIN,LEN

!     CALCULATE AND OUTPUT GEOMETRY FROM SATELLITE ABOVE 120KM.
!     SUBTRACT FROM 180 DEGREES TO CORRECTLY PLACE ANGLE IN THE
!     3RD QUADRANT.
!
      IF (HOBS.GT.0.) THEN
         IF (H2.GT.H1) THEN
            H_TOA = H2
            SINTOA = SINPHI
            TOA_ANG = PHI
         ELSE
            H_TOA = H1
            SINTOA = SINANGLE
            TOA_ANG = ANGLE
         ENDIF
         SINTOA_SAT = ((RE+H_TOA)/(RE+HOBS))*SINTOA
         TOA_SAT = 180.0D0 - ASIN(SINTOA_SAT)/RADCONV
         SINTOA_SAT = SIN(RADCONV*TOA_SAT)
         DIFFANGLE = TOA_SAT - TOA_ANG
         IF (NOPRNT .GE. 0)WRITE (IPR,937) HOBS,TOA_SAT,SINTOA_SAT,DIFFANGLE
      ENDIF

      RETURN
!
!     ERROR MESSAGES
!
   80 CONTINUE
      IF (NOPRNT .GE. 0)WRITE (IPR,940) H1,HMIN
      GO TO 140
                            !,ITYPE
   90 WRITE (0,945) ITYPE
      GO TO 140
  100 WRITE (0,950) H1,H2,RANGE
      GO TO 140
  110 CONTINUE
      WRITE (0,955) H1,H2,ANGLE
      GO TO 140
  120 WRITE (0,960)
      GO TO 140
  130 WRITE (0,965) ZMAX,H1,H2,HMIN
  140 IERROR = 1
!
      RETURN
!
  900 FORMAT (/,' CASE 3A: GIVEN H1,H2=SPACE,ANGLE')
  905 FORMAT (/,' CASE 3B: GIVEN H1, HMIN, H2=SPACE')
  910 FORMAT (/,' CASE 2A: GIVEN H1, H2, ANGLE')
  915 FORMAT (/,' EITHER A SHORT PATH (LEN=0) OR A LONG PATH ',         &
     &        'THROUGH A TANGENT HEIGHT (LEN=1) IS POSSIBLE: LEN = ',   &
     &        I3)
  920 FORMAT (/,' CASE 2B:, GIVEN H1, ANGLE, RANGE',//,10X,             &
     &        'NOTE: H2 IS COMPUTED FROM H1, ANGLE, AND RANGE ',        &
     &        'ASSUMING REFRACTION')
!  925 FORMAT (/,10X,'CALCULATED H2 IS LESS THAN ZERO:',/,10X,           &
!     &        'RESET H2 = 0.0 AND RANGE = ',F10.3)
  930 FORMAT (/,' CASE 2C: GIVEN H1, H2, RANGE',//,10X,                 &
     &        'NOTE: ANGLE IS COMPUTED FROM H1, H2, AND RANGE ',        &
     &        'ASSUMING NO REFRACTION')
  935 FORMAT (' SLANT PATH PARAMETERS IN STANDARD FORM',/               &
     &        /,10X,'H1         = ',F12.6,' KM',                        &
     &        /,10X,'H2         = ',F12.6,' KM',                        &
     &        /,10X,'ANGLE      = ',F12.6,' DEG',                       &
     &        /,10X,'SIN(ANGLE) = ',F12.6,                              &
     &        /,10X,'PHI        = ',F12.6,' DEG',                       &
     &        /,10X,'SIN(PHI)   = ',F12.6,                              &
     &        /,10X,'HMIN       = ',F12.6,' KM',                        &
     &        /,10X,'LEN        = ',I10)
  937 FORMAT (///,' SLANT PATH PARAMETERS AT SATELLITE',/               &
     &        /,10X,'H_SAT        = ',F12.6,' KM',                      &
     &        /,10X,'PHI_SAT      = ',F12.6,' DEG'                      &
     &        /,10X,'SIN(PHI_SAT) = ',F12.6,                            &
     &        /,10X,'PHI_SAT-PHI  = ',F12.6,' DEG')
  940 FORMAT ('0FSCGEO: CASE 3B (H1,HMIN,SPACE): ERROR IN INPUT DATA',  &
     &        //,10X,'H1 = ',F12.6,'    IS LESS THAN HMIN = ',F12.6)
  945 FORMAT ('0FSCGEO: ERROR IN INPUT DATA, ITYPE NOT EQUAL TO ',      &
     &        ' 2, OR 3.   ITYPE = ',I10,E23.14)
  950 FORMAT ('0FSCGEO: CASE 2C (H1,H2,RANGE): ERROR IN INPUT DATA',    &
     &        //,10X,'ABS(H1-H2) GT RANGE;  H1 = ',F12.6,'    H2 = ',   &
     &        F12.6,'    RANGE = ',F12.6)
  955 FORMAT ('0FSCGEO: CASE 2A (H1,H2,ANGLE): ERROR IN INPUT DATA',    &
     &        //,10X,'H1 = ',F12.6,'    IS GREATER THAN OR EQUAL TO',   &
     &        ' H2 = ',F12.6,/,10X,'AND ANGLE = ',F12.6,'    IS LESS',  &
     &        ' THAN OR EQUAL TO 90.0')
  960 FORMAT ('0FSCGEO: ITYPE = 2: SLANT PATH INTERSECTS THE EARTH',    &
     &        ' AND CANNOT REACH H2')
  965 FORMAT (' FSCGEO:  THE ENTIRE PATH LIES ABOVE THE TOP ZMAX ',     &
     &        'OF THE ATMOSPHERIC PROFILE',//,10X,'ZMAX = ',G12.6,5X,   &
     &        '  H1 = ',G12.6,5X,'  H2 = ',G12.6,'  HMIN = ',G12.6)
!
 1000 FORMAT (/3X, 'AMBIGUOUS INPUTS:',/3X,'H1 AND RANGE ARE BOTH > 0', &
     &    /3X,'BUT H2 AND ANGLE = 0',//3X,'PATH COULD BE 2B OR 2C',     &
     &    //5X,'WILL ASSUME 2C',//3X,                                   &
     &    'CHANGE IN FSCGEO IF 2B IS DESIRED')

      END SUBROUTINE FSCGEO
!
!     ----------------------------------------------------------------
!
      SUBROUTINE REDUCE (H1,H2,ANGLE,PHI,ITER)
!
!     *****************************************************************
!     ZMAX IS THE HIGHEST LEVEL IN THE ATMOSPHERIC PROFILE STORED IN
!     COMMON /MDATA/.  IF H1 AND/OR H2 ARE GREATER THAN ZMAX, THIS
!     SUBROUTINE REDUCES THEM TO ZMAX AND RESETS ANGLE AND/OR PHI
!     AS NECESSARY. THIS REDUCTION IS NECESSARY,FOR EXAMPLE FOR
!     SATELLITE ALTITUDES, BECAUSE (1) THE DENSITY PROFILES ARE
!     POORLY DEFINED ABOVE ZMAX AND (2) THE CALCULATION TIME FOR
!     PATHS ABOVE ZMAX CAN BE EXCESSIVE ( EG. FOR GEOSYNCRONOUS
!     ALTITUDES)
!     *****************************************************************
!
      REAL (8)     :: H1, H2, ANGLE, PHI, SH, GAMMA, CPATH, CZMAX, ZMAX, ANGMAX
      INTEGER (4)  :: ITER
      REAL (8), EXTERNAL :: ANDEX

      IF (H1.LE.ZMAX.AND.H2.LE.ZMAX) RETURN

      CALL FINDSH (H1,SH,GAMMA)
      CPATH = ANDEX(H1,SH,GAMMA)*(RE+H1)*SIN(ANGLE/DEG)

      CALL FINDSH (ZMAX,SH,GAMMA)
      CZMAX = ANDEX(ZMAX,SH,GAMMA)*(RE+ZMAX)
      ANGMAX = 180.0D0-ASIN(CPATH/CZMAX)*DEG

      IF (H1.LE.ZMAX) GO TO 10
      H1 = ZMAX
      ANGLE = ANGMAX
   10 CONTINUE
      IF (H2.LE.ZMAX) GO TO 20
      H2 = ZMAX
      PHI = ANGMAX
   20 CONTINUE
      IF (ITER.EQ.0.AND.NOPRNT.GE.0) WRITE (IPR,900) ZMAX,ANGMAX
!
      RETURN
!
  900 FORMAT (///,' FROM SUBROUTINE REDUCE : ',/,10X,'ONE OR BOTH OF',  &
     &        ' H1 AND H2 ARE ABOVE THE TOP OF THE ATMOSPHERIC ',       &
     &        'PROFILE ZMAX = ',F10.3,'  AND HAVE BEEN RESET TO ZMAX.', &
     &        /,10X,'ANGLE AND/OR PHI HAVE ALSO BEEN RESET TO THE ',    &
     &        'ZENITH ANGLE AT ZMAX = ',F10.3,' DEG')
!
      END SUBROUTINE REDUCE
!
!     ----------------------------------------------------------------
!
      SUBROUTINE FDBETA (H1,H2,BETAS,ANGLE,PHI,LEN,HMIN,IERROR)
!
!     *****************************************************************
!     GIVEN H1,H2,AND BETA (THE EARTH CENTERED ANGLE) THIS SUBROUTINE
!     CALCULATES THE INITIAL ZENITH ANGLE AT H1 THROUGH AN ITERATIVE
!     PROCEDURE
!     *****************************************************************

      INTEGER (4)  :: IERROR, LEN, ITERMX, IBMSAV, IFLAG, IAMTB, IORDER, ITER
      REAL (8)     :: H1, H2, BETAS, ANGLE, PHI, HMIN, BETA1, BETA2, BETAP
      REAL (8)     :: RA, RB, SG, ANGLE1, ANGLE2, BETA, DBETA, ANGLS1, ANGLEP
      REAL (8)     :: TOLRNC, BETD, ZER, HA, HB, RANGE, ANGLS2, BENDNG, HMING
      REAL (8)     :: DANG, TEMP !, DERIV, DC

      DATA TOLRNC / 5.0D-3 /,ITERMX / 10 /,BETD / 0.04D0 /
      DATA ZER / 0. /

      BETA = BETAS
      IFLAG = 0
      IF (H1.LE.H2) THEN
          IORDER = 1
          HA = H1
          HB = H2
      ELSE
          IORDER = -1
          HA = H2
          HB = H1
      ENDIF
!
!     IF AUTOLAYERING SELECTED(IBMAX = 0) THEN SET UP DUMMY
!     LBLRTM OUTPUT LAYERS
!
      IBMSAV = IBMAX
      IF (IBMAX.EQ.0) THEN
          IBMAX = 2
          ZBND(1) = ZMIN
          ZBND(2) = ZMAX
      ENDIF
!
!     SET PARAMETER TO SUPRESS CALCULATION OF AMOUNTS
!
      IAMTB = 2
!
!     GUESS AT ANGLE, INTEGRATE TO FIND BETA, TEST FOR
!     CONVERGENCE, AND ITERATE
!     FIRST GUESS AT ANGLE: USE THE GEOMETRIC SOLUTION (NO REFRACTION)
!
      IF (NOPRNT .GE. 0)WRITE (IPR,900)
      ITER = 0
      RA = RE+HA
      RB = RE+HB
      SG = SQRT((HA-HB)**2+4.0D0*RA*RB*(SIN(BETA/(2.0D0*DEG)))**2)
      ANGLE1 = 180.0D0-ACOS((HA**2-HB**2+2.0D0*RE*(HA-HB)+SG**2)            &
     &         /(2.0D0*RA*SG))*DEG
      HMIN = HA
      IF (ANGLE1.GT.90.0D0) HMIN = RA*SIN(ANGLE1/DEG)-RE
      HMING = HMIN
      ANGLS1 = ANGLE1
      CALL FNDHMN (HA,ANGLS1,HB,LEN,ITER,HMIN,PHI,IERROR)
      LEN = 0
      IF (HMIN.LT.HA) LEN = 1
      CALL RFPATH (HA,HB,ANGLS1,PHI,LEN,HMIN,IAMTB,RANGE,BETA1,BENDNG)
      IF (NOPRNT .GE. 0)WRITE (IPR,905) ITER,ANGLS1,BETA,ZER,SG,HMING,ZER,ZER
!
!     OBTAIN DERIVATIVE
!
      SG = SQRT((HA-HB)**2+4.0D0*RA*RB*(SIN((BETA+BETD)/(2.0D0*DEG)))**2)
      ANGLEP = 180.0D0-ACOS((HA**2-HB**2+2.0D0*RE*(HA-HB)+SG**2)            &
     &         /(2.0D0*RA*SG))*DEG
      DANG = ANGLE1-ANGLEP
      IF (HMIN.LT.0.0) THEN
          IFLAG = 1
          HMIN = 0.0
          CALL FNDHMN (HMIN,90.0D0,HA,LEN,ITER,HMIN,ANGLS1,IERROR)
      ENDIF
      ITER = 1
      LEN = 0
      IF (ANGLE1.GT.90.0) LEN = 1
      CALL FNDHMN (HA,ANGLS1,HB,LEN,ITER,HMIN,PHI,IERROR)
      LEN = 0
      IF (HMIN.LT.HA) LEN = 1
      CALL RFPATH (HA,HB,ANGLS1,PHI,LEN,HMIN,IAMTB,RANGE,BETA1,BENDNG)
      DBETA = BETA-BETA1
      IF (NOPRNT .GE. 0)WRITE (IPR,905) ITER,ANGLS1,BETA1,DBETA,RANGE,HMIN,PHI,BENDNG
      IF (IFLAG.EQ.1.AND.BETA1.LT.BETA) GO TO 90
   50 CONTINUE
      ANGLEP = ANGLE1-DANG
      LEN = 0
      IF (ANGLEP.GT.90.0) LEN = 1
      CALL FNDHMN (HA,ANGLEP,HB,LEN,ITER,HMIN,PHI,IERROR)
      LEN = 0
      IF (HMIN.LT.HA) LEN = 1
      CALL RFPATH (HA,HB,ANGLEP,PHI,LEN,HMIN,IAMTB,RANGE,BETAP,BENDNG)
      IF (ABS(BETA1-BETAP).LT.TOLRNC) GO TO 60
      ITER = ITER+1
      !DC = BETAP-BETA1
      !DERIV = -DC/BETD
      ANGLE2 = ANGLE1+(ANGLE1-ANGLEP)*(BETA-BETA1)/(BETA1-BETAP)
      ANGLS2 = ANGLE2
      LEN = 0
      IF (ANGLE2.GT.90.0) LEN = 1
      CALL FNDHMN (HA,ANGLS2,HB,LEN,ITER,HMIN,PHI,IERROR)
      LEN = 0
      IF (HMIN.LT.HA) LEN = 1
      CALL RFPATH (HA,HB,ANGLS2,PHI,LEN,HMIN,IAMTB,RANGE,BETA2,BENDNG)
      DBETA = BETA-BETA2
      IF (NOPRNT .GE. 0)WRITE (IPR,905) ITER,ANGLS2,BETA2,DBETA,RANGE,HMIN,PHI,BENDNG
      IF (BETA2.LT.BETA.AND.HMIN.LT.0.0) GO TO 90
      ANGLE1 = ANGLE2
      ANGLS1 = ANGLE1
      BETA1 = BETA2
      IF (ABS(BETA-BETA2).LT.TOLRNC) GO TO 70
      IF (ITER.GT.ITERMX) GO TO 100
      GO TO 50
   60 ANGLE2 = ANGLEP
      ANGLS2 = ANGLE2
      BETA = BETAP
   70 CONTINUE
      IF (HMIN.LT.0.0) GO TO 90
!
!     CONVERGED TO A SOLUTION
!
      ANGLE = ANGLE2
      BETA = BETA2
!
!     ASSIGN ANGLE AND PHI TO PROPER H1 AND H2
!
      IF (IORDER.NE.1) THEN
          TEMP = PHI
          PHI = ANGLE
          ANGLE = TEMP
      ENDIF
      IBMAX = IBMSAV
      BETAS = BETA
!
      RETURN
!
!     ERROR MESSAGES
!
   90 CONTINUE
      WRITE (0,910)
      GO TO 110
  100 CONTINUE
      WRITE (0,915) H1,H2,BETA,ITER,ANGLE1,BETA1,ANGLE2,BETA2
!
  110 IERROR = 1
!
      RETURN
!
  900 FORMAT (///,' CASE 2D: GIVEN H1, H2,  BETA:',//,                  &
     &        ' ITERATE AROUND ANGLE UNTIL BETA CONVERGES',//,          &
     &        ' ITER    ANGLE',T21,'BETA',T30,'DBETA',T40,'RANGE',      &
     &        T51,'HMIN',T61,'PHI',T70,'BENDING',/,T10,'(DEG)',T21,     &
     &        '(DEG)',T30,'(DEG)',T41,'(KM)',T51,'(KM)',T60,'(DEG)',    &
     &        T71,'(DEG)',/)
  905 FORMAT (I5,3F10.4,2F10.3,2F10.4)
  910 FORMAT ('0FDBETA, CASE 2D(H1,H2,BETA): REFRACTED TANGENT ',       &
     &        'HEIGHT IS LESS THAN ZERO-PATH INTERSECTS THE EARTH',     &
     &        //,10X,'BETA IS TOO LARGE FOR THIS H1 AND H2')
  915 FORMAT ('0FDBETA, CASE 2D (H1,H2,BETA): SOLUTION DID NOT ',       &
     &        ' CONVERGE',//,10X,'H1 = ',F12.6,'    H2 = ',F12.6,       &
     &        '    BETA = ',F12.6,'    ITERATIONS = ',I4,//,10X,        &
     &        'LAST THREE ITERATIONS ',//,(10X,'ANGLE = ',F15.9,        &
     &        '    BETA = ',F15.9))
!
      END SUBROUTINE FDBETA
!
!     ----------------------------------------------------------------
!
      SUBROUTINE FNDHMN (H1,ANGLE,H2,LEN,ITER,HMIN,PHI,IERROR)
!
!     *****************************************************************
!     THIS SUBROUTINE CALCULATES THE MINIMUM ALTITUDE HMIN ALONG
!     THE REFRACTED PATH AND THE FINAL ZENITH ANGLE PHI.
!     THE PARAMETER LEN INDICATES WHETHER THE PATH GOES THROUGH
!     A TANGENT HEIGHT (LEN=1) OR NOT (LEN=0).  IF ANGLE > 90 AND
!     H1 > H2, THEN LEN CAN EITHER BE 1 OR 0, AND THE CHOICE IS
!     LEFT TO THE USER.
!     THE (INDEX OF REFRACTION - 1.0) IS MODELED AS AN EXPONENTIAL
!     BETWEEN THE LAYER BOUNDARIES, WITH A SCALE HEIGHT SH AND AN
!     AMOUNT AT THE GROUND GAMMA.
!     CPATH IS THE REFRACTIVE CONSTANT FOR THIS PATH AND
!     EQUALS  INDEX(H1)*(RE+H1)*SIN(ANGLE).
!     *****************************************************************

      REAL (8)     :: H1, ANGLE, H2, HMIN, PHI, DERIV, DC
      REAL (8)     :: CPATH, CRFRCT, SH, GAMMA, CT1, CTP, CH2, CMIN
      REAL (8)     :: DH, ETA, H, HTP, HT1
      INTEGER (4)  :: N, LEN, ITER, IERROR
      REAL (8), EXTERNAL :: ANDEX
!      REAL (8) :: ANDEX

      DATA DH / 0.2D0 /,ETA / 5.0D-7 /

! --- ETA MAY BE TOO SMALL FOR SOME COMPUTERS. TRY 1.0E-7 FOR 32 BIT WORD MACHINES

      CRFRCT(H) = ( RE + H )*ANDEX( H, SH, GAMMA )

      N = 0
      CALL FINDSH( H1, SH, GAMMA )
      CPATH = CRFRCT(H1)*SIN(ANGLE/DEG)
      CALL FINDSH( H2, SH, GAMMA )
      CH2 = CRFRCT(H2)
      IF (ABS(CPATH/CH2).GT.1.0) GO TO 70
      IF (ANGLE.LE.90.0) THEN
          LEN = 0
          HMIN = H1
          GO TO 60
      ENDIF
      IF (H1.LE.H2) LEN = 1
      IF (LEN.NE.1) THEN
          LEN = 0
          HMIN = H2
          GO TO 60
      ENDIF

!     LONG PATH THROUGH A TANGENT HEIGHT.
!     SOLVE ITERATIVELY FOR THE TANGENT HEIGHT HT.
!     HT IS THE HEIGHT FOR WHICH  INDEX(HT)*(RE+HT) = CPATH.

      CALL FINDSH( 0.0D0, SH, GAMMA )
      CMIN = CRFRCT(0.0D0)
!
!     FOR BETA CASES (ITER>0), ALLOW FOR HT < 0.0
!
      IF (ITER.EQ.0.AND.CPATH.LT.CMIN) GO TO 50
      HT1 = H1*SIN(ANGLE/DEG)+(SIN(ANGLE/DEG)-1.0)*RE
!
!     ITERATE TO FIND HT
!
   30 CONTINUE
      N = N+1
      CALL FINDSH (HT1,SH,GAMMA)
      CT1 = CRFRCT(HT1)
      IF (ABS((CPATH-CT1)/CPATH).LT.ETA) GO TO 40
      IF (N.GT.15) GO TO 80
      HTP = HT1-DH
      CALL FINDSH (HTP,SH,GAMMA)
      CTP = CRFRCT(HTP)
      DERIV=(CT1-CTP)/DH
      HT1=HT1+(CPATH-CT1)/DERIV
      GO TO 30
   40 CONTINUE
      HMIN=HT1
      GO TO 60
   50 CONTINUE
!
!     TANGENT PATH INTERSECTS EARTH
!
      H2 = 0.0
      HMIN = 0.0
      LEN = 0
      CH2 = CMIN
      WRITE (0,900) H1,ANGLE
   60 CONTINUE
!
!     CALCULATE THE ZENITH ANGLE PHI AT H2
!
      PHI = ASIN(CPATH/CH2)*DEG
      IF (ANGLE.LE.90.0.OR.LEN.EQ.1) PHI = 180.0-PHI
!
      RETURN
!
!     H2 LT TANGENT HEIGHT FOR THIS H1 AND ANGLE
!
   70 CONTINUE
      WRITE (0,905)
      IERROR = 2
!
      RETURN
!
   80 CONTINUE
      DC = CPATH-CT1
      WRITE (0,910) N,CPATH,CT1,DC,HT1
!
      WRITE(16,*) ' FNDHMN '
      WRITE(00,*) ' FNDHMN '
      CALL SHUTDOWN
      STOP '3'

  900 FORMAT (///,' TANGENT PATH WITH H1 = ',F10.3,' AND ANGLE = ',     &
     &        F10.3,' INTERSECTS THE EARTH',//,10X,                     &
     &        'H2 HAS BEEN RESET TO 0.0 AND LEN TO 0')
  905 FORMAT ('0H2 IS LESS THAN THE TANGENT HEIGHT FOR THIS PATH ',     &
     &        'AND CANNOT BE REACHED')
  910 FORMAT (///,'0FROM SUBROUTINE FNDHMN :',//,10X,                   &
     &        'THE PROCEEDURE TO FIND THE TANGENT HEIGHT DID NOT ',     &
     &        'CONVERG AFTER ',I3,'  ITERATIONS',//,10X,'CPATH   = ',   &
     &        F12.5,' KM',//,10X,'CT1     = ',F12.5,' KM',//,10X,       &
     &        'DC      = ',E12.3,' KM',//,10X,'HT1     = ',F12.5,' KM')
!
      END SUBROUTINE FNDHMN
!
!     ----------------------------------------------------------------
!
      SUBROUTINE FINDSH (H,SH,GAMMA)
!
!     *****************************************************************
!     GIVEN AN ALTITUDE H, THIS SUBROUTINE FINDS THE LAYER BOUNDARIES
!     Z(I1) AND Z(I2) WHICH CONTAIN H,  THEN CALCULATES THE SCALE
!     HEIGHT (SH) AND THE VALUE AT THE GROUND (GAMMA+1) FOR THE
!     REFRACTIVITY (INDEX OF REFRACTION -1)
!     *****************************************************************

      REAL (8)      :: H, SH, GAMMA
      INTEGER (4)   :: IM, I1, I2

      DO 10 IM = 2, IMMAX
         I2 = IM
         IF (ZMDL(IM).GE.H) GO TO 20
   10 END DO
      I2 = IMMAX
   20 CONTINUE
      I1 = I2-1
      CALL SCALHT (ZMDL(I1),ZMDL(I2),RFNDXM(I1),RFNDXM(I2),SH,GAMMA)
!
      RETURN
!
      END SUBROUTINE FINDSH
!
!     ----------------------------------------------------------------
!
      SUBROUTINE SCALHT( Z1, Z2, RFNDX1, RFNDX2, SH, GAMMA )
!
!     *****************************************************************
!     THIS SUBROUTINE CALCULATES THE SCALE HEIGHT SH OF THE (INDEX OF
!     REFRACTION-1.0) FROM THE VALUES OF THE INDEX AT THE ALTITUDES Z1
!     AND Z2 ( Z1 < Z2). IT ALSO CALCULATES THE EXTRAPOLATED VALUE
!     GAMMA OF THE (INDEX-1.0) AT Z = 0.0
!     *****************************************************************
!

      REAL (8) :: Z1, Z2, RFNDX1, RFNDX2, SH, GAMMA
      REAL (8) :: RF1, RF2, RATIO

      RF1 = RFNDX1+1.0D-20
      RF2 = RFNDX2+1.0D-20
      RATIO = RF1/RF2
      IF (ABS(RATIO-1.0D0).LT.1.0D-05) GO TO 10
      SH = (Z2-Z1)/ LOG(RATIO)
      GAMMA = RF1*(RF2/RF1)**(-Z1/(Z2-Z1))
      GO TO 20
   10 CONTINUE
!
!     THE VARIATION IN THE INDEX OF REFRACTION WITH HEIGHT IS
!     INSIGNIFICANT OR ZERO
!
      SH = 0.0
      GAMMA = RFNDX1
   20 CONTINUE
!
      RETURN
!
      END SUBROUTINE SCALHT


! ----------------------------------------------------------------
!
 !     REAL (8) FUNCTION ANDEX (H,SH,GAMMA)
!
!     DOUBLE PRECISION VERSION OF ANDEX - NEEDED FOR IMPROVED GEOMETRY
!
!     *****************************************************************
!     COMPUTES THE INDEX OF REFRACTION AT HEIGHT H, SH IS THE
!     SCALE HEIGHT, GAMMA IS THE VALUE AT H=0 OF THE REFRACTIVITY =
!     INDEX-1
!     *****************************************************************
!
 !     REAL (8) :: H, SH, GAMMA
 !
 !     IF (SH.EQ.0.0) THEN
 !        ANDEX = 1.0D0 + GAMMA
 !     ELSE
 !        ANDEX = 1.0D0 + GAMMA*EXP(-H/SH)
 !     ENDIF
!
 !     RETURN
!
 !     END FUNCTION ANDEX
!
! ----------------------------------------------------------------
!
      REAL (8) FUNCTION RADRF (H,SH,GAMMA)
!
!     DOUBLE PRECISION VERSION OF RADREF - NEEDED FOR IMPROVED GEOMETRY
!
!     *****************************************************************
!     COMPUTES THE RADIUS OF CURVATURE OF THE REFRACTED RAY FOR
!     A HORIZONTAL PATH.  RADREF = ANDEX/ D(ANDEX)/D(RADIUS)
!     *****************************************************************
!
      REAL (8) :: H, SH, GAMMA, BIGNUM
      DATA BIGNUM / 1.0D36 /

      IF (SH.EQ.0.0) GO TO 10
      RADRF = SH*(1.0D0 + EXP(H/SH)/GAMMA)
!
      RETURN
!
   10 RADRF = BIGNUM
!
      RETURN
!
      END FUNCTION RADRF

!     ----------------------------------------------------------------
!
      SUBROUTINE RFPATH (H1,H2,ANGLE,PHI,LEN,HMIN,IAMT,RANGE,BETA, BENDNG)
!
!     -------------------------------------------------------------
!     THIS ROUTINE WAS MODIFIED FOR LBLRTM TO REFLECT CHANGES
!     IMPLEMENTED IN MODTRAN TO SOLVE PROBLEMS WITH INCONSISTENT
!     PATH PARAMETERS.
!     IT WAS ALSO MODIFIED TO ELIMINATE GOTO STATEMENTS IN ORDER TO
!     MAKE THE PROGRAM EASIER TO UNDERSTAND.
!     THESE CHANGES WERE OBTAINED FROM H. SNELL (MARCH, 1996).
!     -------------------------------------------------------------
!
!     *****************************************************************
!     THIS SUBROUTINE TRACES THE REFRACTED RAY FROM H1 WITH AN
!     INITIAL ZENITH ANGLE ANGLE TO H2 WHERE THE ZENITH ANGLE IS PHI,
!     AND CALCULATES THE ABSORBER AMOUNTS (IF IAMT.EQ.1) ALONG
!     THE PATH.  IT STARTS FROM THE LOWEST POINT ALONG THE PATH
!     (THE TANGENT HEIGHT HMIN IF LEN = 1 OR HA = MIN(H1,H2) IF LEN = 0
!     AND PROCEEDS TO THE HIGHEST POINT.  BETA AND RANGE ARE THE
!     EARTH CENTERED ANGLE AND THE TOTAL DISTANCE RESPECTIVELY
!     FOR THE REFRACTED PATH FROM H1 TO H2, AND BENDNG IS THE TOTAL
!     BENDING ALONG THE PATH
!     *****************************************************************
!
!
      REAL (8)     :: DS, DBEND, S, SINAI, COSAI, CPATH, SH, GAMMA, HA
      REAL (8)     :: ANGLEA, RHOBAR, THETA, DBETA, PBAR1, TBAR1
      REAL (8)     :: H1, H2, ANGLE, PHI, HMIN, RANGE, BETA, BENDNG
      INTEGER (4)  :: LEN, IAMT, I_2, IORDER, J2, J, IHLOW, IHIGH
      REAL (8), EXTERNAL :: ANDEX
      CHARACTER (LEN=2) :: HLOW(2)

      DATA HLOW / 'H1','H2'/

      DATA I_2/2/

!     REORDER H1 AND H2 TO HA AND HB (HA .LE. HB)

      IF (H1.LE.H2) THEN
          IORDER = 1
          HA = H1
          !HB = H2
          ANGLEA = ANGLE
      ELSE
          IORDER = -1
          HA = H2
          !HB = H1
          ANGLEA = PHI
      ENDIF
!
!     MERGE THE ATMOSPHERIC PROFILE STORED IN ZMDL WITH H1,H2,(HMIN) AN
!     THE BOUNDARIES ZBND
!

!print*, 'rfpath ', h1, h2, hmin

      CALL AMERGE (H1,H2,HMIN,LEN)
      IF (IAMT.EQ.1.AND.NOPRNT.GE.0) WRITE (IPR,900)
!
!     CALCULATE CPATH SEPERATELY FOR LEN = 0,1
!
      IF (LEN.EQ.0) THEN
          CALL FINDSH ( HA, SH, GAMMA)
          CPATH = ( RE + HA )*ANDEX( HA, SH, GAMMA )*SIN( ANGLEA/DEG )
      ELSE
          CALL FINDSH ( HMIN, SH, GAMMA)
          CPATH = ( RE + HMIN )*ANDEX( HMIN, SH, GAMMA )
      ENDIF
!
      BETA = 0.0
      S = 0.0
      BENDNG = 0.0
      IF (LEN.EQ.1) THEN
!
!     TANGENT PATH
!
          IF (IORDER.EQ.-1) THEN
              IHLOW = 2
          ELSE
              IHLOW = 1
          ENDIF
          IF (IAMT.EQ.1.AND.NOPRNT.GE.0) WRITE (IPR,905) HLOW(IHLOW)
          SINAI = 1.0D0
          COSAI = 0.0D0
          THETA = 90.0D0
      ELSE
!
!     SHORT PATH
!
!     ANGLEA IS THE ZENITH ANGLE AT HA IN DEG
!     SINAI IS SIN OF THE INCIDENCE ANGLE
!     COSAI IS CARRIED SEPERATELY TO AVOID A PRECISION PROBLEM
!     WHEN SINAI IS CLOSE TO 1.0
!
          THETA = ANGLEA
          IF (ANGLEA.LE.45.0) THEN
              SINAI = SIN(ANGLEA/DEG)
              COSAI = -COS(ANGLEA/DEG)
          ELSE
              SINAI = COS((90.0-ANGLEA)/DEG)
              COSAI = -SIN((90.0-ANGLEA)/DEG)
          ENDIF
          IF (IORDER.EQ.-1) THEN
              IHLOW = 2
          ELSE
              IHLOW = 1
          ENDIF
          IHIGH = MOD(IHLOW,I_2)+1
          IF (IAMT.EQ.1.AND.NOPRNT.GE.0)                                &
     &        WRITE (IPR,910) HLOW(IHLOW),HLOW(IHIGH)
      ENDIF
!
!     LOOP OVER THE LAYERS
!
      J2 = IPMAX-1
      DO 100 J = 1, J2
         CALL SCALHT (ZPTH(J),ZPTH(J+1),RFNDXP(J),RFNDXP(J+1),SH,GAMMA)
         CALL ALAYER (J,SINAI,COSAI,CPATH,SH,GAMMA,IAMT,DS,DBEND)
         DBEND = DBEND*DEG
         PHI = ASIN(SINAI)*DEG
         DBETA = THETA-PHI+DBEND
         PHI = 180.0D0-PHI
         S = S+DS
         BENDNG = BENDNG+DBEND
         BETA = BETA+DBETA
         IF (IAMT.EQ.1) THEN
             PBAR1 = PPSUM(J)/RHOPSM(J)
             TBAR1 = TPSUM(J)/RHOPSM(J)
             RHOBAR = RHOPSM(J)/DS

             IF (NOPRNT.GE.0) WRITE (IPR,915) J,ZPTH(J),ZPTH(J+1),      &
     &            THETA,DS,S,DBETA,BETA,PHI,DBEND,BENDNG,PBAR1,          &
     &            TBAR1,RHOBAR

         ENDIF
         THETA = 180.0-PHI
!
         IF (LEN.EQ.1) THEN
!
!            FOR TANGENT PATHS, DOUBLE THE QUANTITIES BENDNG,BETA,
!            AND S FOR THE SYMMETRIC PART OF THE PATH
!
             IF ((J+1).EQ.IPHMID) THEN
                 BENDNG = 2.0*BENDNG
                 BETA = 2.0*BETA
                 S = 2.0*S
                 IF (IAMT.EQ.1.AND.NOPRNT.GE.0)                         &
     &                         WRITE (IPR,920) S,BETA,BENDNG
                 IF (IPHMID.NE.IPMAX) THEN
                     IF (IORDER.EQ.-1) THEN
                         IHLOW = 2
                     ELSE
                         IHLOW = 1
                     ENDIF
                     IHIGH = MOD(IHLOW,I_2)+1
                     IF (IAMT.EQ.1.AND.NOPRNT.GE.0)                     &
     &                   WRITE (IPR,910) HLOW(IHLOW),HLOW(IHIGH)
                 ENDIF
             ENDIF
         ENDIF
  100 END DO
      IF (IORDER.EQ.-1) PHI = ANGLEA
      RANGE = S
!
      RETURN
!
  900 FORMAT (/,'CALCULATION OF THE REFRACTED PATH THROUGH THE ',        &
     &        'ATMOSPHERE',/,T5,'I',T14,'ALTITUDE',T30,'THETA',T38,   &
     &        'DRANGE',T47,'RANGE',T57,'DBETA',T65,'BETA',T76,'PHI',    &
     &        T84,'DBEND',T91,'BENDING',T102,'PBAR',T111,'TBAR',T119,   &
     &        'RHOBAR',/,T11,'FROM',T22,'TO',/,T11,'(KM)',T21,'(KM)',   &
     &        T30,'(DEG)',T39,'(KM)',T48,'(KM)',T57,'(DEG)',T65,        &
     &        '(DEG)',T75,'(DEG)',T84,'(DEG)',T92,'(DEG)',T102,'(MB)',  &
     &        T112,'(K)',T117,'(MOL CM-3)',/)
  905 FORMAT (' ',T10,'TANGENT',T20,A2,/,T10,'HEIGHT',/)
  910 FORMAT (' ',T14,A2,' TO ',A2,/)
  915 FORMAT (' ',I4,2F10.3,10F9.3,1PE9.2)
  920 FORMAT ('0',T10,'DOUBLE RANGE, BETA, BENDING',/,T10,              &
     &        'FOR SYMMETRIC PART OF PATH',T44,F9.3,T62,F9.3,T89,       &
     &        F9.3,/)
!
      END SUBROUTINE RFPATH
!
!     ----------------------------------------------------------------
!
      SUBROUTINE AMERGE (H1,H2,HMIN,LEN)
!
!     *****************************************************************
!     AMERGE CREATES A SET OF LAYER BOUNDARIES ZOUT WHICH INCLUDES
!     HMIN, (HMID), HMAX AND ALL OF ZBND BETWEEN HMIN AND HMAX.
!     ZOUT DEFINES THE LAYERS FOR THE LBLRTM CALCULATION.
!     ZOUT IS THEN MERGED WITH THE ATMOSPHERIC PROFILE IN ZMDL INTO ZPTH
!     INTERPOLATING TO THE LEVELS ZOUT WHEN NECESSARY.  THE RAY
!     TRACE IS CALCULATED USING THE PROFILE IN ZPTH.

! ZBND - OUTPUT GRID SPECIFIED IN TAPE5
! ZMDL - INPUT MODEL GRID WITH TM & PM
! ZFIN - ZBND WITH TOP AND BOTTOM - FINISHED OUTPUT GRID
! ZFINE - BUILT IN FINE GRID
! ZPTH - ZMDL + ZFINE + ZFIN FOR RAYTRACING
! ZOUT - ZFINE + ZFIN (INTERMEDIATE)

!     *****************************************************************

      INTEGER (4) :: LEN, I_2, IHMAX, I1, IB, IH, IOUT, IM, K, IP, JM

      REAL (8)    ::  H1, H2, HMIN, ZH(3), TOL, HMID, HMAX, A

      ZFIN = 0.0D0

      DATA TOL / 5.0D-4 /, I_2 / 2 /


!print *, 'amerge zbnd ', zbnd
!print *, 'amerge h1 h2 ', h1, h2

      NMRGCALL = NMRGCALL +1
      !PRINT *, 'NMRGCALL ', NMRGCALL

!     HMID .EQ. MINIMUM OF H1, H2

      HMID =   MIN(H1,H2)
      HMAX =   MAX(H1,H2)
      IHMAX = 2
      ZH(1) = HMIN
      IF (LEN.EQ.0) THEN
          ZH(2) = HMAX
      ELSE
          ZH(2) = HMID
          IF (ABS(H1-H2).LT.TOL) H1 = H2
          IF (H1.NE.H2) THEN
              IHMAX = 3
              ZH(3) = HMAX
          ENDIF
      ENDIF

!     MERGE ZH AND ZBND BETWEEN ZH(1) AND ZH(IHMAX) TO CREATE ZFIN
! --- THIS IS OUR FINAL OUTPUT INCLUDING A TOP & INTERNAL BOUNDS & BOTTOM

      ZFIN(1) = ZH(1)
      DO 30 I1 = 1, IBMAX
         IF (ABS(ZBND(I1)-ZH(1)).LT.TOL) ZH(1) = ZBND(I1)
         IF (ZBND(I1).GT.ZH(1)) GO TO 40
   30 END DO
      I1 = IBMAX
   40 CONTINUE

!     ZBND(I1) IS SMALLEST ZBND .GT. ZH(1)

      IOUT = 1
      IB = I1
      IH = 2
   50 CONTINUE
      IOUT = IOUT+1
      IF (IB.GT.IBMAX) GO TO 60
      IF (ABS(ZBND(IB)-ZH(IH)).LT.TOL) ZH(IH) = ZBND(IB)
      IF (ZBND(IB).LT.ZH(IH)) GO TO 70
      IF (ZBND(IB).EQ.ZH(IH)) IB = IB+1

!     INSERT ZH(IH)

   60 CONTINUE
      ZFIN(IOUT) = ZH(IH)
      IH = IH+1
      IF (IH.GT.IHMAX) GO TO 80
      GO TO 50

!     INSERT ZBND(IB)

   70 CONTINUE
      ZFIN(IOUT) = ZBND(IB)
      IB = IB+1
      GO TO 50
!
   80 CONTINUE

      IFINMX = IOUT

      !WRITE(*,901) "ZFIN", IFINMX, ZFIN(1:IFINMX)


! REALLY WE WANT A FINER GRID THAN THE ZBND & ZMDL MERGE
! BUT STILL WANT IT TO START & STOP ON HMIN & HMAX
! AND WE WANT ZBND POINTS IN THE GRID AT THE END TO EXTRACT THEM IN THE STD WAY
! HERE WE INSERT ZFINE INTO ZFIN

!     MERGE ZFINE AND ZFIN BETWEEN ZH(1) AND ZH(IHMAX) TO CREATE ZOUT

      !WRITE(*,901) "ZH ", 3, ZH

      ZOUT(1) = ZFIN(1)
      IB   = 1
      IOUT = 1
      I1   = 0
 23   I1   = I1 +1
         IF( I1 .GT. NFINE ) GOTO 22
         IF( ZFINE(I1) .LT. ZOUT(1) )GOTO 23
 21      IF( ABS(ZFINE(I1) - ZFIN(IB)) .GT. TOL .AND. ZFINE(I1) .LT. ZFIN(IB) )THEN
            ZOUT(IOUT) = ZFINE(I1)
            IOUT = IOUT +1
            GOTO 23
         ENDIF
         IF( ABS(ZFINE(I1) - ZFIN(IB)) .LT. TOL )THEN
            ZOUT(IOUT) = ZFINE(I1)
            IF( IB .GE. IFINMX )GOTO 22
            IOUT = IOUT +1
            IB   = IB   +1
            GOTO 23
         ENDIF
         IF( ZFINE(I1) .GT. ZFIN(IB) .AND. IB .LE. IFINMX )THEN
            ZOUT(IOUT) = ZFIN(IB)
            IF( IB .GE. IFINMX )GOTO 22
            IOUT = IOUT +1
            IB   = IB   +1
            GOTO 21
         ENDIF

 22   CONTINUE
      IOUTMX = IOUT

      !WRITE(*,901) "ZOUT", IOUTMX , ZOUT(1:IOUTMX)

!     NOW MERGE ZOUT AND ZMDL INTO ZPTH (FROM ZOUT(1) TO ZOUT(IOUTMX))
!     AND INTERPOLATE PRESSURE, TEMPERATURE, AND DENSITY WHEN
!     NECESSARY

!     FIND SMALLEST ZMDL .GT. HMIN

      DO 90 IM = 1, IMMAX
         IF (ZMDL(IM).GE.HMIN) GO TO 100
   90 END DO
      WRITE(16,900) HMIN
      WRITE(00,900) HMIN
      CALL SHUTDOWN
      STOP '3'
  100 CONTINUE
      IPHMID = 0
      IP = 0
      IOUT = 1
  110 CONTINUE
      IP = IP+1
      IF (IP.GT.IPDIM) THEN
         WRITE(16,905) IPDIM
         WRITE( 0,905) IPDIM
         CALL SHUTDOWN
         STOP '3'
      ENDIF
      !print *, 'zpth ', zpth
      !print*, 'im ', im, zmdl(im), iout, zout(iout), ip, zpth(ip), ioutmx

      IF (IM.GT.IMMAX) GO TO 130
      IF (ABS(ZOUT(IOUT)-ZMDL(IM)).LT.TOL) ZMDL(IM) = ZOUT(IOUT)
      IF (ZOUT(IOUT).LT.ZMDL(IM)) GO TO 130
      IF (ZOUT(IOUT).EQ.ZMDL(IM)) IOUT = IOUT+1
!
!     INSERT ZMDL(IM)
!
      ZPTH(IP)   = ZMDL(IM)
      PP(IP)     = PM(IM)
      TP(IP)     = TM(IM)
      RFNDXP(IP) = RFNDXM(IM)
      DO 120 K = 1, NMOL
         DENP(K,IP) = DENM(K,IM)
  120 END DO
      IM = IM+1
      IF (ABS(ZPTH(IP)-HMID).LT.TOL) HMID = ZPTH(IP)
      IF (ZPTH(IP).EQ.HMID) IPHMID = IP
      IF (ABS(ZPTH(IP)-ZOUT(IOUTMX)).LT.TOL) ZOUT(IOUTMX) = ZPTH(IP)
      IF (ZPTH(IP).EQ.ZOUT(IOUTMX)) GO TO 150
      GO TO 110
!
!     INSERT LEVEL FROM ZOUT(IOUT) AND INTERPOLATE
!
  130 CONTINUE
      ZPTH(IP) = ZOUT(IOUT)
      JM = IM
      JM = MAX(JM,I_2)
      A = (ZOUT(IOUT)-ZMDL(JM-1))/(ZMDL(JM)-ZMDL(JM-1))
      CALL EXPINT (PP(IP),PM(JM-1),PM(JM),A)
      TP(IP) = TM(JM-1)+(TM(JM)-TM(JM-1))*A
      CALL EXPINT (RFNDXP(IP),RFNDXM(JM-1),RFNDXM(JM),A)

      DO 140 K = 1, NMOL
         CALL EXPINT (DENP(K,IP),DENM(K,JM-1),DENM(K,JM),A)
  140 END DO
      IF (ABS(ZPTH(IP)-HMID).LT.TOL) ZPTH(IP) = HMID
      IF (ZPTH(IP).EQ.HMID) IPHMID = IP
      IOUT = IOUT+1
      IF (ABS(ZPTH(IP)-ZOUT(IOUTMX)).LT.TOL) ZPTH(IP) = ZOUT(IOUTMX)
      IF (ZPTH(IP).EQ.ZOUT(IOUTMX)) GO TO 150
      GO TO 110

  150 CONTINUE


      !PRINT *, IP, ZPTH(IP),IOUTMX,  ZOUT(IOUTMX)

      IPMAX = IP
      IF( NMRGCALL .LE. 1 .AND. NOPRNT.GE.0)THEN
         WRITE(IPR,906)  " BOUNDRARIES DETERMINED FOR OUTPUT : ", IFINMX
         WRITE(IPR,906)  " NUM BOUNDRARIES WITH FINE GRID    : ", IOUTMX
         WRITE(IPR,906)  " NUM BOUNDRARIES FOR PATHS CALC    : ", IPMAX
      ENDIF

      !WRITE(*,901) "ZPTH", IPMAX, ZPTH(1:IPMAX)

      RETURN
!
  900 FORMAT ('0FROM AMERGE- ATMOSPHERIC PROFILE IN ZMDL DOES NOT',     &
     &        ' EXTEND UP TO HMIN = ',E12.5)
!  901 FORMAT(A,2X,I5,500G12.4)
  905 FORMAT ('0FROM AMERGE- MERGING THE ATMOSPHERIC PROFILE AND THE ', &
     &        'LBLRTM BOUNDARIES INTO ZPTH(IPDIM) EXCEEDS THE ',        &
     &        'DIMENSION IPDIM = ',I5)
 906  FORMAT( A, I6 )

      END SUBROUTINE AMERGE

!     ----------------------------------------------------------------

      SUBROUTINE ALAYER (J,SINAI,COSAI,CPATH,SH,GAMMA,IAMT,S,BEND)

!     -------------------------------------------------------------
!     THIS ROUTINE WAS MODIFIED FOR LBLRTM TO REFLECT CHANGES
!     IMPLEMENTED IN MODTRAN TO SOLVE PROBLEMS WITH INCONSISTENT
!     PATH PARAMETERS.
!     IT WAS ALSO MODIFIED TO ELIMINATE GOTO STATEMENTS IN ORDER TO
!     MAKE THE PROGRAM EASIER TO UNDERSTAND.
!     THESE CHANGES WERE OBTAINED FROM H. SNELL (MARCH, 1996).
!     -------------------------------------------------------------
!
!     *****************************************************************
!     THIS SUBROUTINE TRACES THE OPTICAL RAY THROUGH ONE LAYER FROM
!     Z1 TO Z2 AND IF IAMT.NE.2 CALCULATES THE INTEGRATED ABSORBER
!     AMOUNTS FOR THE LAYER. SINAI IS THE SIN OF THE INITIAL INCIDENCE
!     ANGLE (= 180 - ZENITH ANGLE). COSAI IS CARRIED SEPERATELY TO
!     AVOID A PRECISION PROBLEM NEAR SINAI = 1. CPATH IS THE CONSTANT
!     OF REFRACTION FOR THE PATH = INDEX*RADIUS*SINAI, SH AND GAMMA ARE
!     THE SCALE HEIGHT AND THE AMOUNT AT THE GROUND FOR THE REFRACTIVIT
!     (= 1-INDEX OF REFRACTION), S IS THE REFRACTED PATH LENGTH THROUGH
!     THE LAYER, BETA IS THE EARTH CENTERED ANGLE, AND BEND IS THE
!     BENDING THROUGH THE LAYER. IAMT CONTROLS WHETHER AMOUNTS ARE
!     CALCULATED OR NOT.
!     *****************************************************************

      INTEGER (4) :: J, IAMT, N, K

      REAL (8)    :: S, BEND, DS, DBEND, W1, W2, W3, DSDX1, DSDX2, DSDX3
      REAL (8)    :: DBNDX1, DBNDX2, DBNDX3, R1, R2, R3, X1, X2, X3, RATIO1, RATIO2
      REAL (8)    :: RATIO3, SINAI1, SINAI2, SINAI3, COSAI1, COSAI2, COSAI3
      REAL (8)    :: CPATH, DX, DH, SINAI, COSAI, D31, D32, D21, DHMIN, GAMMA
      REAL (8)    :: SH, EPSILN, Z1, Z2, H1, H2, H3, Y1=0.0, Y3, PA=0.0, PB=0.0, TA, TB
      REAL (8)    :: RHOA=0.0, RHOB, DZ=0.0, HP=0.0, HRHO=0.0, DSDZ
      REAL (8), EXTERNAL :: ANDEX
      REAL (8), DIMENSION(MXMOL) :: HDEN, DENA, DENB

      DATA EPSILN / 1.0D-5 /

!     INITIALIZE VARIABLES FOR THE CALCULATION OF THE PATH

      N = 0
      Z1 = ZPTH(J)
      Z2 = ZPTH(J+1)
      H1 = Z1
      R1 = RE+H1
      DHMIN = DELTAS**2/(2.0*R1)
      SINAI1 = SINAI
      COSAI1 = COSAI
      IF((1.0D0-SINAI).LT.EPSILN ) Y1 = COSAI1**2/2.0D0+COSAI1**4/8.0D0+COSAI1**6*3.0D0/48.0D0
      Y3 = 0.0
      X1 = -R1*COSAI1
      RATIO1 = R1/RADRF(H1,SH,GAMMA)
      DSDX1 = 1.0D0/(1.0D0-RATIO1*SINAI1**2)
      DBNDX1 = DSDX1*SINAI1*RATIO1/R1
      S = 0.0D0
      BEND = 0.0D0
      IF (IAMT.NE.2) THEN

!         INITIALIZE THE VARIABLES FOR THE CALCULATION OF THE ABSORBER AMOUNTS

          PA = PP(J)
          PB = PP(J+1)
          IF (PB.EQ.PA) THEN
             WRITE(16,*) 'ALAYER: PRESSURES IN ADJOINING LAYERS MUST DIFFER', PB
             WRITE(00,*) 'ALAYER: PRESSURES IN ADJOINING LAYERS MUST DIFFER', PB
             CALL SHUTDOWN
             STOP '3'
          ENDIF
          TA = TP(J)
          TB = TP(J+1)
          RHOA = PA/(GCAIR*TA)
          RHOB = PB/(GCAIR*TB)
          DZ = ZPTH(J+1)-ZPTH(J)
          HP = -DZ/ LOG(PB/PA)
          IF (ABS( RHOB/RHOA - 1.0D0 ).GE.EPSILN) THEN
              HRHO = -DZ/ LOG(RHOB/RHOA)
          ELSE
              HRHO = 1.0D30
          ENDIF
          DO 40 K = 1, NMOL
              DENA(K) = DENP(K,J)
              DENB(K) = DENP(K,J+1)
              IF ((DENA(K).EQ.0.0D0.OR.DENB(K).EQ.0.0D0).OR.                &
     &            (ABS(1.0-DENA(K)/DENB(K)).LE.EPSILN)) THEN
!
!                 USE LINEAR INTERPOLATION
!
                  HDEN(K) = 0.0D0
              ELSE
!
!                 USE EXPONENTIAL INTERPOLATION
!
                  HDEN(K) = -DZ / LOG(DENB(K)/DENA(K))
              ENDIF
   40     CONTINUE
      ENDIF
!
!     LOOP THROUGH PATH
!     INTEGRATE PATH QUANTITIES USING QUADRATIC INTEGRATION WITH
!     UNEQUALLY SPACED POINTS
!
   60 CONTINUE
      N = N+1
      DH = -DELTAS*COSAI1
      DH = MAX(DH,DHMIN)
      H3 = H1+DH
      IF (H3.GT.Z2) H3 = Z2
      DH = H3-H1
      R3 = RE+H3
      H2 = H1+DH/2.0
      R2 = RE+H2
      SINAI2 = CPATH/(ANDEX( H2, SH, GAMMA )*R2)
      SINAI3 = CPATH/(ANDEX( H3, SH, GAMMA )*R3)
      RATIO2 = R2/RADRF( H2, SH, GAMMA )
      RATIO3 = R3/RADRF( H3, SH, GAMMA )
      IF ((1.0D0-SINAI2).LE.EPSILN) THEN

!        NEAR A TANGENT HEIGHT, COSAI = -SQRT(1-SINAI**2) LOSES
!        PRECISION. USE THE FOLLOWING ALGORITHM TO GET COSAI.

         Y3 = Y1+(SINAI1*(1.0D0-RATIO1)/R1+4.0D0*SINAI2*(1.0D0-RATIO2)/R2+    &
     &        SINAI3*(1.0D0-RATIO3)/R3)*DH/6.0D0
         COSAI3 = -SQRT(2.0D0*Y3-Y3**2)
         X3 = -R3*COSAI3
         DX = X3-X1
         W1 = 0.5D0*DX
         W2 = 0.0D0
         W3 = 0.5D0*DX
      ELSE
         COSAI2 = -SQRT(1.0D0-SINAI2**2)
         COSAI3 = -SQRT(1.0D0-SINAI3**2)
         X2 = -R2*COSAI2
         X3 = -R3*COSAI3
!
!        CALCULATE WEIGHTS
!
         D31 = X3-X1
         D32 = X3-X2
         D21 = X2-X1
         IF (D32.EQ.0.0.OR.D21.EQ.0.0) THEN
            W1 = 0.5D0*D31
            W2 = 0.0D0
            W3 = 0.5D0*D31
         ELSE
            W1 = (2.0D0-D32/D21)*D31/6.0D0
            W2 = D31**3/(D32*D21*6.0D0)
            W3 = (2.0D0-D21/D32)*D31/6.0D0
         ENDIF
      ENDIF
      DSDX2 = 1.0D0/(1.0D0-RATIO2*SINAI2**2)
      DSDX3 = 1.0D0/(1.0D0-RATIO3*SINAI3**2)
      DBNDX2 = DSDX2*SINAI2*RATIO2/R2
      DBNDX3 = DSDX3*SINAI3*RATIO3/R3
!
!     INTEGRATE
!
      DS = W1*DSDX1+W2*DSDX2+W3*DSDX3
      S = S+DS
      DBEND = W1*DBNDX1+W2*DBNDX2+W3*DBNDX3
      BEND = BEND+DBEND
      IF (IAMT.NE.2) THEN
!
!         CALCULATE AMOUNTS
!
         DSDZ = DS/DH
         PB = PA*EXP(-DH/HP)
         RHOB = RHOA*EXP(-DH/HRHO)
         IF ((DH/HRHO).GE.EPSILN) THEN
            PPSUM(J)  = PPSUM(J)+DSDZ*(HP/(1.0D0+HP/HRHO))*(PA*RHOA-PB*RHOB)
            TPSUM(J)  = TPSUM(J)+DSDZ*HP*(PA-PB)/GCAIR
            RHOPSM(J) = RHOPSM(J)+DSDZ*HRHO*(RHOA-RHOB)
         ELSE
            PPSUM(J)  = PPSUM(J)+0.5D0*DS*(PA*RHOA+PB*RHOB)
            TPSUM(J)  = TPSUM(J)+0.5D0*DS*(PA+PB)/GCAIR
            RHOPSM(J) = RHOPSM(J)+0.5D0*DS*(RHOA+RHOB)
         ENDIF
         DO 130 K = 1, NMOL
            IF ((HDEN(K).EQ.0.0).OR.(ABS(DH/HDEN(K)).LT.EPSILN)) THEN
!
!                 LINEAR INTERPOLATION
!                 1.0E05 FACTOR CONVERTS UNITS KM TO CM
!
               DENB(K)=DENP(K,J)+(DENP(K,J+1)-DENP(K,J))*(H3-Z1)/DZ
               AMTP(K,J) = AMTP(K,J)+0.5D0*(DENA(K)+DENB(K))*DS*1.0D5
            ELSE
!
!                 EXPONENTIAL INTERPOLATION
!
               DENB(K) = DENP(K,J)*EXP(-(H3-Z1)/HDEN(K))
               AMTP(K,J) = AMTP(K,J)+DSDZ*HDEN(K)*(DENA(K)-DENB(K))*1.0D5
            ENDIF
  130    CONTINUE
         PA = PB
         RHOA = RHOB
         DO 140 K = 1, NMOL
            DENA(K) = DENB(K)
  140    CONTINUE
      ENDIF
!
      IF (H3.LT.Z2) THEN
         H1 = H3
         R1 = R3
         SINAI1 = SINAI3
         RATIO1 = RATIO3
         Y1 = Y3
         COSAI1 = COSAI3
         X1 = X3
         DSDX1 = DSDX3
         DBNDX1 = DBNDX3
      ELSE
         SINAI = SINAI3
         COSAI = COSAI3
         SP(J) = S
         RETURN
      ENDIF
!
      GO TO 60
!
      END SUBROUTINE ALAYER
!
!     ----------------------------------------------------------------
!
      SUBROUTINE AUTLAY( HMIN, HMAX, XVBAR, AVTRAT, TDIFF1, TDIFF2, ALTD1,    &
     &                   ALTD2, IERROR)
!
!     *****************************************************************
!     THIS SUBROUTINE AUTOMATICALLY SELECTS A SET OF LBLRTM BOUNDARY
!     LEVELS WHICH SATISFY THE FOLLOWING TWO TESTS:
!          1. THE RATIO OF THE VOIGT HALFWIDTHS BETWEEN BOUNDARIES
!             IS LESS THAN OR EQUAL TO AVTRAT, AND
!          2. THE TEMPERATURE DIFFERENCE BETWEEN BOUNDARIES IS
!             LESS THAN OR EQUAL TO TDIFF
!     TDIFF VARIES FROM TDIFF1 AT HMIN TO TDIFF2 AT HMAX,
!     WITH EXPONENTIAL INTERPOLATION BETWEEN
!     THESE BOUNDARIES ARE ROUNDED DOWN TO THE NEAREST TENTH KM
!     NOTE THAT THESE TESTS APPLY TO THE LAYER BOUNDARIES
!     NOT TO THE AVERAGE VALUES FROM ONE LAYER TO THE NEXT.
!     *****************************************************************

      INTEGER (4) :: IERROR, IM, IB, IBM1, IPASS, IHMIN=0, IND

      REAL (8)    :: XVBAR, HMIN, HMAX, AVTRAT, TDIFF1, TDIFF2, ALTD1, ALTD2
      REAL (8)    :: HTOP, TMIN, TMAX, X, Y, ALOGY, TDIFF, ZZ, FAC, ZX, ZROUND
      REAL (8)    :: P, T, AL, AD, ZBNDTI, ALOGX

      REAL (8), DIMENSION(MXZMD)  :: AVTM

!     FUNCTION ZROUND ROUNDS THE ALTITUDE Z DOWN TO THE
!     NEAREST TENTH KM
      ZROUND(ZX) = 0.1D0* REAL( INT(10.0*ZX), 8)

      HMIN = MAX(HMIN,ZMDL(1))
!
      DO 10 IM = 2, IMMAX
         IHMIN = IM
         IF (ZMDL(IM).GT.HMIN) GO TO 20
   10 END DO
   20 CONTINUE
      HTOP = HMAX
      HTOP = MIN(HTOP,ZMAX)
      IM = IHMIN-1
      ZZ = ZMDL(IM)
      CALL HALFWD (ZZ,XVBAR,P,T,AL,AD,AVTM(IM))
      IB = 1
      ZBND(IB) = HMIN
      IM = IHMIN
      CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB),         &
     &             ADOPP(IB),AVOIGT(IB))
!
!     BEGIN IM LOOP
!
   30 CONTINUE
      IB = IB+1
      IF (IB.GT.IBDIM) GO TO 90
      IBM1 = IB-1
      TMIN = TBND(IBM1)
      TMAX = TBND(IBM1)
      IND = 0
!
!     BEGIN IB LOOP
!
   40 CONTINUE
      IPASS = 0
      ZBND(IB) = ZMDL(IM)
      ZBNDTI = ZMDL(IM)
      IF (ZBND(IB).GE.HTOP) ZBND(IB) = HTOP
      CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB), ADOPP(IB),AVOIGT(IB))
      AVTM(IM) = AVOIGT(IB)
!
!     TEST THE RATIO OF THE VOIGT WIDTHS AGAINST AVTRAT
!
      IF ((AVOIGT(IB-1)/AVOIGT(IB)).LT.AVTRAT) GO TO 50
!
!     ZMDL(IM) FAILS THE HALFWIDTH RATIO TEST
!
      IPASS = 1
      AVOIGT(IB) = AVOIGT(IB-1)/AVTRAT
      X = AVTM(IM)/AVTM(IM-1)
      ALOGX = 1.0D0 - X
      IF (ABS(ALOGX).LT.0.001) THEN
         ZBND(IB) = (ZMDL(IM)+ZMDL(IM-1))/2.
         GO TO 50
      ELSE
         ALOGX =  LOG(X)
      ENDIF
      Y = AVOIGT(IB)/AVTM(IM-1)
      ALOGY = 1.-Y
      IF (ABS(ALOGY).GT.0.001) ALOGY =  LOG(Y)
      ZBND(IB) = ZMDL(IM-1)+(ZMDL(IM)-ZMDL(IM-1))*ALOGY/ALOGX
   50 CONTINUE
!
!     TEST THE TEMPERATURE DIFFERENCE AGAINST TDIFF
!
      FAC = (ZBND(IB-1)-ALTD1)/(ALTD2-ALTD1)
      CALL EXPINT (TDIFF,TDIFF1,TDIFF2,FAC)
      IF (TM(IM).GT.TMAX) THEN
         IND = 1
         TMAX = TM(IM)
      ENDIF
      IF (TM(IM).LT.TMIN) THEN
         IND = 2
         TMIN = TM(IM)
      ENDIF
      IF (TMAX-TMIN.LE.TDIFF) GO TO 60
      IF (IND.EQ.1) TBND(IB) = TMIN+TDIFF
      IF (IND.EQ.2) TBND(IB) = TMAX-TDIFF
!
!     ZBND(IB) FAILS THE TEMPERATURE DIFFERENCE TEST
!
      IPASS = 2
      IF (ABS(TM(IM)-TM(IM-1)).LT.0.0001) THEN
         ZBNDTI = (ZMDL(IM)+ZMDL(IM-1))/2.
      ELSE
         ZBNDTI = ZMDL(IM-1)+(ZMDL(IM)-ZMDL(IM-1))*                     &
     &            (TBND(IB)-TM(IM-1))/(TM(IM)-TM(IM-1))
      ENDIF
   60 CONTINUE
      IF (ZBNDTI.LT.ZBND(IB)) ZBND(IB) = ZBNDTI
      IF (ZBND(IB).GE.HTOP) THEN
         ZBND(IB) = HTOP
         IF (ZBND(IB)-ZBND(IB-1).LE.0.1) THEN
            IB = IB-1
            ZBND(IB) = HTOP
            CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB),   &
     &                   ADOPP(IB),AVOIGT(IB))
         ENDIF
         GO TO 80
      ENDIF
      IF (IPASS.NE.0) GO TO 70
!
!     BOTH HALFWIDTH AND TEMPERATURE TEST PASS FOR ZBND(IB) = ZMDL(IM),
!     NOW TRY ZBND(IB) = ZMDL(IM+1)
!
      IM = IM+1
      GO TO 40
   70 CONTINUE
!
!     ONE OF THE TESTS FAILED AND A NEW BOUNDRY ZBND WAS PRODUCED
!
      ZBND(IB) = ZROUND(ZBND(IB))
      CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB),         &
     &             ADOPP(IB),AVOIGT(IB))
      GO TO 30
   80 CONTINUE
      IBMAX = IB
      IF (NOPRNT .GE. 0)WRITE (IPR,900) AVTRAT,TDIFF1,HMIN,TDIFF2,HMAX
!
      RETURN
!
   90 CONTINUE
      WRITE (0,905) IBDIM
      IBMAX = IBDIM
      IERROR = 5
!
      RETURN
!
  900 FORMAT (///,                                                      &
     &        ' LBLRTM LAYER BOUNDARIES PRODUCED BY THE AUTOMATIC ',    &
     &        'LAYERING ROUTINE AUTLAY',/,' THE USER SHOULD EXAMINE ',  &
     &        'THESE BOUNDARIES AND MODIFY THEM IF APPROPRIATE',/,      &
     &        ' THE FOLLOWING PARAMETERS ARE USED:',//,10X,             &
     &        'AVTRAT    = ',F8.2,'       = MAX RATIO OF VOIGT WIDTHS', &
     &        /,10X,'TDIFF1    = ',F8.2,'       = MAX TEMP DIFF AT ',   &
     &        F4.0,' KM',/10X,'TDIFF2    = ',F8.2,                      &
     &        '       = MAX TEMP DIFF AT ',F4.0,' KM')
  905 FORMAT (///,' ERROR IN AUTLAY:',/,5X,'THE NUMBER OF ',            &
     &        'GENERATED LAYER BOUNDARIES EXCEEDS THE DIMENSION IBDIM', &
     &        ' OF THE ARRAY ZBND.  IBDIM = ',I5,/,5X,'PROBABLE CAUSE', &
     &        ': EITHER AVTRAT AND/OF TDIFF ARE TOO SMALL',/,5X,        &
     &        'THE GENERATED LAYERS FOLLOW')
!
      END SUBROUTINE AUTLAY
!
!     ----------------------------------------------------------------
!
      SUBROUTINE HALFWD_P( XVBAR, P, T, ALORNZ1, ADOPP1, AVOIGT1 )
!
!     *****************************************************************
!     GIVEN AN PRESSURE AND TEMP. AND AVERAGE WAVENUMBER VBAR, THIS
!     SUBROUTINE
!     CALCULATES THE LORENTZ, THE DOPPLER, AND THE VOIGT HALFWIDTHS
!     (AT HALFHEIGHT) ALORNZ, ADOPP, AND AVOIGT RESPECTIVELY FOR
!     THE ALTITUDE Z
!     AN AVERAGE LORENTZ WIDTH ALZERO AND AN AVERAGE MOLECULAR
!     WEIGHT AVMWT ARE ASSUMED
!     *****************************************************************!
!
      REAL(8)     :: XVBAR, P, T, ALPHAL, ALPHAD, ALPHAV, V, AL, AD
      REAL(8)     :: ALORNZ1, ADOPP1, AVOIGT1

!     FUNCTIONS
!     ALZERO IS AT 1013.25 MB AND 296.0 K
!
      ALPHAL(P,T) = ALZERO*(P/PZERO)*SQRT(296.0D0/T)
      ALPHAD(T,V) = ADCON*V*SQRT(T/AVMWT)
      ALPHAV(AL,AD) = 0.5D0*(AL+SQRT(AL**2+4.0D0*AD**2))

      ALORNZ1  = ALPHAL(P,T)
      ADOPP1   = ALPHAD(T,XVBAR)
      AVOIGT1  = ALPHAV(ALORNZ1,ADOPP1)
!
      RETURN
!
      END SUBROUTINE HALFWD_P
!
!     ----------------------------------------------------------------
!
!
!     ----------------------------------------------------------------
!
      SUBROUTINE HALFWD (Z,XVBAR,P,T,ALORNZ1,ADOPP1,AVOIGT1)
!
!     *****************************************************************
!     GIVEN AN ALTITUDE Z AND AN AVERAGE WAVENUMBER VBAR, THIS
!     SUBROUTINE INTERPOLATES P AND T FROM THE PROFILE IN ZMDL  AND
!     CALCULATES THE LORENTZ, THE DOPPLER, AND THE VOIGT HALFWIDTHS
!     (AT HALFHEIGHT) ALORNZ, ADOPP, AND AVOIGT RESPECTIVELY FOR
!     THE ALTITUDE Z
!     AN AVERAGE LORENTZ WIDTH ALZERO AND AN AVERAGE MOLECULAR
!     WEIGHT AVMWT ARE ASSUMED
!     *****************************************************************

      REAL(8)     :: Z, XVBAR, P, T, ALPHAL, ALPHAD, ALPHAV, V, AL, AD
      REAL(8)     :: ALORNZ1, ADOPP1, AVOIGT1, FAC
      INTEGER (4) :: I2, IM

!     FUNCTIONS
!     ALZERO IS AT 1013.25 MB AND 296.0 K
!
      ALPHAL(P,T)    = ALZERO*(P/PZERO)*SQRT(296.0D0/T)
      ALPHAD(T,V)    = ADCON*V*SQRT(T/AVMWT)
      ALPHAV(AL,AD)  = 0.5D0*(AL+SQRT(AL**2+4.0D0*AD**2))

      DO 10 I2 = 2, IMMAX
         IM = I2
         IF (ZMDL(IM).GE.Z) GO TO 20
   10 END DO
      IM = IMMAX
   20 CONTINUE
      FAC = (Z-ZMDL(IM-1))/(ZMDL(IM)-ZMDL(IM-1))
      CALL EXPINT (P,PM(IM-1),PM(IM),FAC)
      T = TM(IM-1)+(TM(IM)-TM(IM-1))*FAC
      ALORNZ1  = ALPHAL(P,T)
      ADOPP1   = ALPHAD(T,XVBAR)
      AVOIGT1  = ALPHAV(ALORNZ1,ADOPP1)
!
      RETURN
!
      END SUBROUTINE HALFWD
!
!     ----------------------------------------------------------------
!
      SUBROUTINE FPACK (H1,H2,HMID,LEN,IEMIT,NOZERO)
!
!     *****************************************************************
!     FPACK TAKES THE AMOUNTS STORED IN THE LAYERS DEFINED BY ZPTH AND
!     PACKS THEM INTO THE LAYERS DEFINED BY ZFIN.  IT ALSO ZEROS OUT
!     LAYER AMOUNTS IF THE AMOUNT FOR THAT LAYER AND ABOVE IS LESS
!     THAN 0.1 PERCENT OF THE TOTAL FOR THAT MOLECULE, UNLESS THE
!     NOZERO OPTION IS SELECTED.
!     *****************************************************************

      INTEGER (4) :: LEN, IEMIT, NOZERO, I2, IOUT, IP, K, L, L2, LMAX, ISKPT
      INTEGER (4) :: NMOL_MAX
      REAL (8)    :: H1, H2, HMID, SUMAMT, FAC, TOL !,tbound
      DATA TOL / 5.0D-4 /

      I2 = IPMAX-1
      IOUT = 1
      PZ(0) = PP(1)
      TZ(0) = TP(1)

!     IF ENTRY IN TAPE5 FOR TBOUND < 0, USE TZ(O) AS BOUNDARY
!     TEMPERATURE

!      IF (TBOUND.LT.0.) TBOUND = TZ(0)
!
      DO 20 IP = 1, I2

         PBAR(IOUT) = PBAR(IOUT)+PPSUM(IP)
         TBAR(IOUT) = TBAR(IOUT)+TPSUM(IP)
         RHOSUM(IOUT) = RHOSUM(IOUT)+RHOPSM(IP)
         SOUT(IOUT) = SOUT(IOUT)+SP(IP)
         DO 10 K = 1, NMOL
            AMOUNT(K,IOUT) = AMOUNT(K,IOUT)+AMTP(K,IP)
   10    CONTINUE

         IF (ABS(ZPTH(IP+1)-ZFIN(IOUT+1)).LT.TOL) THEN
            PZ(IOUT) = PP(IP+1)
            TZ(IOUT) = TP(IP+1)
            IOUT = IOUT+1
         ENDIF

   20 END DO

      IF (IOUT.NE.IFINMX) GO TO 110

!     CALCULATE THE DENSITY WEIGHTED PRESSURE AND TEMPERATURE AND
!     ZERO OUT LAYER AMOUNTS AFTER 99.9 PERCENT OF THE TOTAL

      ISKIP(7) = 0

      DO 30 K = 1, NMOL
         AMTCUM(K) = 0.0
         ISKIP(K) = 0
         IF (AMTTOT(K).EQ.0.0) ISKIP(K) = 1
   30 END DO
      L2 = IFINMX-1
      LMAX = L2
      DO 90 L = 1, L2
         PBAR(L) = PBAR(L)/RHOSUM(L)
         TBAR(L) = TBAR(L)/RHOSUM(L)
!
!     ADJUST RHOSUM FOR THE PATH LENGTH IN CM NOT KM
!
         RHOSUM(L) = RHOSUM(L)*1.0D+5
!
         SUMAMT = 0.0D0
         DO 40 K = 1, NMOL
            SUMAMT = SUMAMT + AMOUNT(K,L)
   40    CONTINUE
         WN2L(L) = RHOSUM(L)-SUMAMT

! --- CALCULATE 'EFFECTIVE SECANT' SECNTA

         SECNTA(L) = SOUT(L)/(ZFIN(L+1)-ZFIN(L))
         IF (L.EQ.1) ALTZ(0) = ZFIN(1)
         ALTZ(L) = ZFIN(L+1)

!     SET  IPATH

         IF (LEN.EQ.1) GO TO 50
         IF (H1.LT.H2) IPATH(L) = 3
         IF (H1.GT.H2) IPATH(L) = 1
         GO TO 60
   50    CONTINUE
         IF (ZFIN(L).LT.HMID) IPATH(L) = 2
         IF (ZFIN(L).GE.HMID.AND.H1.GT.H2) IPATH(L) = 1
         IF (ZFIN(L).GE.HMID.AND.H1.LT.H2) IPATH(L) = 3
   60    CONTINUE
!
!     TEST FOR ZEROING OF AMOUNTS
!
         ISKPT = 0
         NMOL_MAX = NMOL
         IF (ISKIP(7).EQ.1) NMOL_MAX = NMOL - 1
         FAC = 1.0D0
         IF (IPATH(L).EQ.2) FAC = 2.0D0
!
         DO 80 K = 1, NMOL

            IF (NOZERO.EQ.1) GO TO 70

            IF (ISKIP(K).NE.1) THEN
               IF (K.EQ.7 .OR. (IEMIT.EQ.1.AND.IPATH(L).NE.3)) GO TO 70
               IF (((AMTTOT(K)-AMTCUM(K))/AMTTOT(K)).GT.0.001) GO TO 70
            ENDIF
!
!     ZERO OUT THIS AMOUNT
!
            ISKIP(K) = 1
            AMOUNT(K,L) = 0.0
            ISKPT = ISKPT+1
!
!     IF ALL BUT O2 ARE ZEROED, ELIMINATE ALL HIGHER LAYERS
!
            IF (ISKPT.GE.(NMOL_MAX))  GO TO 100
   70       CONTINUE
            AMTCUM(K) = AMTCUM(K)+FAC*AMOUNT(K,L)
   80    CONTINUE
         LMAX = L
   90 END DO
  100 CONTINUE
      IFINMX = LMAX+1
!
      RETURN
!
  110 WRITE (0,900) IOUT,IFINMX
      WRITE(16,*) ' ERROR FPACK '
      WRITE( 0,*) ' ERROR FPACK '
      CALL SHUTDOWN
      STOP '3'

  900 FORMAT ('0FROM FPACK-  ERROR, IOUT = ',I5,'  DOES NOT MATCH ',    &
     &        'IFINMX = ',I5)
!
      END SUBROUTINE FPACK


!     ----------------------------------------------------------------


      SUBROUTINE FIXTYP(IEMIT,FRH2O,ALFCOR,OLDDV,L,CINP)

!     *****************************************************************
!     THIS SUBROUTINE CALCULATES ITYL, THE ITYPE (RATIO OF DV FROM
!     ONE LAYER TO THE NEXT) FOR EACH LAYER FOR OUTPUT TO TAPE7, IF
!     DESIRED (IFXTYP = 1).
!     *****************************************************************

      CHARACTER (LEN=3) :: CINP

      INTEGER (4) :: I_2, IEMIT, L, ISCAL, IDV, ITYPE

      REAL (8)    :: DV, AVBAR, FRH2O, ALFCOR, OLDDV, H2OSLF, ALBAR, SCAL
      REAL (8)    :: TYPMAX, TTYPE

      DATA I_2/2/

      DV = 0.

!     CORRECT FOR WATER SELF BROADENING

      H2OSLF = (1.0D0 - FRH2O + 5.0D0 * FRH2O )
      ALBAR = ALZERO*ALFCOR*H2OSLF

      AVBAR = 0.5D0*(ALBAR+SQRT(ALBAR*ALBAR+4.0D0*ADBAR*ADBAR))

      DV = AVBAR/SAMPLE
!
      TTYPE = 0.0D0
      ITYPE = 99
!
!     DV IS ASSUMED TO BE LESS THAN 1
!     SET DV TO 3 SIGNIFICANT FIGURES
!
      IF (L.EQ.1) THEN
         ISCAL =  NINT(LOG10(DV)-3.0D0)
         SCAL  = 10.**ISCAL
         IDV   = NINT((DV/SCAL)+0.5D0)
!
!        SET IDV TO BE EVEN
!
         IF (MOD(IDV,I_2).GT.0) IDV = IDV+1
         DV = SCAL* REAL(IDV)
      ELSE
         TTYPE = OLDDV/DV
         TYPMAX = 2.5
         TYPMAX = 15
!print *, 'ttype ',ttype, typmax
         IF (TTYPE.GT.TYPMAX) THEN
!print *, l, iscal, scal, dv, avbar
!print *, 'ttype ',ttype, typmax
            WRITE(16,*) "FIXTYPE - TTPYE"
            WRITE(00,*) "FIXTYPE - TTPYE"
            CALL SHUTDOWN
            STOP '3'

         ELSEIF (TTYPE.GE.1.2) THEN
!
!           TTYPE IS BETWEEN 1.2 AND TYPMAX
!
            DV = OLDDV
            ITYPE = NINT(1.0D0/(TTYPE-1.0D0)+0.5D0)
            IF (ITYPE.EQ.3) ITYPE = 2
            DV = OLDDV* REAL(ITYPE)/ REAL(ITYPE+1)
         ELSEIF (TTYPE.GE.0.8) THEN
!
!           TTYPE IS BETWEEN 0.8 AND 1.2 (SET TO 1.0)
!
            DV = OLDDV
            ITYPE = 0
         ELSE
!
!           TTYPE IS LESS THAN 0.8
!
            DV = OLDDV
            ITYPE = 0
            IF (IEMIT.NE.1) THEN
               ITYPE = NINT(TTYPE/(1.0D0-TTYPE)+0.5D0)
               DV = DV* REAL(ITYPE+1)/ REAL(ITYPE)
               ITYPE = -ITYPE
            ENDIF
         ENDIF
      ENDIF
!
      OLDDV = DV
!
      WRITE(CINP,900) ITYPE
!
      RETURN
!
  900 FORMAT(I3)
!
      END SUBROUTINE FIXTYP
!
!     ----------------------------------------------------------------
!
!
! -------------------------------------------------------------------
!
      SUBROUTINE NEWH2(H1,H2,ANGLE,RANGE,BETA,LEN,HTAN,PHI)
!
!     CHANGED FOR LBLRTM TO CORRECT GEOMETRY PROBLEMS
!
!     THIS ROUTINE DETERMINES H2,BETA, TANGENT HEIGHT AND LEN.
!     ADOPTED FROM THE MODTRAN2 GEOMETRY PACKAGE
!
!     INPUTS ARE: H1, ZENTIH ANGLE (ANGLE) AND RANGE.
!     LEN = 1 IF THE PATH GOES THROUGH HTAN.
!
!     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO FASE01
!     MXLAY IS THE MAXIMUM NUMBER OF OUTPUT LAYERS
!     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
!         STORED IN ZMDL (INPUT)
!     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
!         OBTAINED BY MERGING ZMDL AND ZOUT
!     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
!
      INTEGER (4)  :: LEN, J, JMAX=0

      REAL (8)     :: H1, H2, ANGLE, RANGE, BETA, HTAN, PHI
      REAL (8)     :: CPATH, CPJ, CPJ1, SH, GAMMA, RE2
      REAL (8)     :: ZJ1, ZJ !, CRFRCT, H
      REAL (8), EXTERNAL :: ANDEX
      !CRFRCT(H) = ( RE2 + H )*ANDEX( H, SH, GAMMA )

      RE2=RE
!     COMPUTE CPATH OR PATH CONSTANT
      CALL FINDSH( H1, SH, GAMMA )
!      CPATH = CRFRCT(H1)*SIN(ANGLE/DEG)
      CPATH = ( RE2 + H1 )*ANDEX( H1, SH, GAMMA )*SIN(ANGLE/DEG)

!     ANGLE = 90 AT H1 IMPLIES THAT H1 = TANGENT HEIGHT
      IF (ANGLE.EQ.90.0) THEN
          HTAN=H1
      ELSE
          DO 100 J=1,IMMAX
              IF (H1.GE.ZMDL(J)) JMAX=J
  100     CONTINUE
          JMAX=JMAX+1
          ZJ1=ZMDL(JMAX)
          !CPJ1=CRFRCT(ZJ1)
          CPJ1 = ( RE2 + ZJ1 )*ANDEX( ZJ1, SH, GAMMA )
          HTAN=-1.0
          DO 200 J=JMAX,1,-1
              IF (HTAN.LT.0.0) THEN
                  IF (J.EQ.1) THEN
                      HTAN=0.0
                  ELSE
                      CPJ=CPJ1
                      ZJ=ZJ1
                      ZJ1=ZMDL(J-1)
                      !CPJ1=CRFRCT(ZJ1)
                      CPJ1 = ( RE2 + ZJ1 )*ANDEX( ZJ1, SH, GAMMA )
                      IF ((CPATH.LE.CPJ).AND.(CPATH.GE.CPJ1)) THEN
                          HTAN=RTBIS(ZJ1,CPJ1,ZJ,CPJ,CPATH)
                      ENDIF
                  ENDIF
              ENDIF
  200     CONTINUE
      ENDIF
!
!     FIND H2, BETA AND LEN
!
      CALL FNDPTH(CPATH,H1,HTAN,H2,RANGE,BETA,LEN,ANGLE,PHI)
!
!     ENSURE LEN IS NOT RESET IN FSCGEO IF DIRECT PATH
      IF (LEN.EQ.0) HTAN=H2
!
!     IF (ANGLE .LE. 90.0) HTAN CARRIES HMIN NOT HTAN
      IF (ANGLE .LE. 90.0) HTAN = MIN(H1,H2)
!
      RETURN
!
      END SUBROUTINE NEWH2
!
! ----------------------------------------------------------------
!
      REAL (8) FUNCTION RTBIS(X1,CX1,X2,CX2,CPATH)
!
!     THIS FUNCTION FINDS THE ROOT OF
!            FUNC(X) = X*REFRACTIVE INDEX - CPA
!
!     THE ROOT IS ACTUALLY THE TANGENT HEIGHT, BETWEEN X1 AND X2.
!     THIS ROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS, ET AL.
!
      INTEGER (4)   :: J, JMAX

      REAL (8)      :: CX1, CX2, CPATH, F, FMID, SH, GAMMA
      REAL (8)      :: X1, X2, XACC, DX, XMID
      REAL (8), EXTERNAL :: ANDEX
      DATA XACC / 1.0D-5 /

      PARAMETER (JMAX=40)

      FMID=CX2-CPATH
      F=CX1-CPATH
      IF(F*FMID.GE.0.)THEN
         WRITE(16,*) 'ROOT MUST BE BRACKETED FOR BISECTION.'
         WRITE(00,*) 'ROOT MUST BE BRACKETED FOR BISECTION.'
         CALL SHUTDOWN
         STOP '3'
      ENDIF
      IF(F.LT.0.)THEN
         RTBIS = X1
         DX    = X2 - X1
      ELSE
         RTBIS = X2
         DX    = X1 - X2
      ENDIF
      DO 11 J=1,JMAX
         DX   = DX * 0.5D0
         XMID = RTBIS + DX
         CALL FINDSH( XMID, SH, GAMMA )
         FMID = ANDEX( XMID, SH, GAMMA )*( XMID + RE ) - CPATH
         IF(FMID.LE.0.0D0)RTBIS = XMID
         IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
   11 END DO
!
!     COMES HERE IF UNABLE TO SOLVE.
!
      IF (ABS(CX2) .LT. ABS(CX1)) THEN
         RTBIS = X2
      ELSE
         RTBIS = X1
      ENDIF
      RETURN

      END FUNCTION RTBIS
!
! ----------------------------------------------------------------
!
      SUBROUTINE FNDPTH( CPATH, H1, HTAN, H2, RANGEI, BETA, LEN, ANGLE, PHI )
!
!     THIS ROUTINE DETERMINES H2, BETA AND LEN.
!     INPUTS ARE H1, HTAN (TANGENT HEIGHT), RANGE (RANGEI) AND
!     THE PATH CONSTANT, CPATH.
!     RANGEO IS THE OUTPUT RANGE WHICH SHOULD EQUAL THE INPUT RANGE.
!
!
      REAL (8)      :: SSAVE, STHETA, CAPRJ, PNTGRN, CTHETA=0.0, CTHET1, DX, H1, H2
      REAL (8)      :: HTAN, RANGEI, BETA, ANGLE, PHI, DR, RANGEO, R1, R2, DZ, Z
      REAL (8)      :: DRNG, DBETA, R, DIFF, CPATH, SH, GAMMA, RX, RATIO, RPLDR
      REAL (8)      :: PERP, BASE, Z2
      REAL (8), EXTERNAL :: ANDEX
      INTEGER (4)   :: I, LEN

      DATA DR / 0.005D0 /

      IF (RANGEI .LT. DR)THEN
         WRITE(16,*) 'STOPPED IN FNDPTH'
         WRITE(00,*) 'STOPPED IN FNDPTH'
         CALL SHUTDOWN
         STOP '3'
      ENDIF
!     (RANGEI .LT. DR) SHOULD NOT HAPPEN; SO THIS CHECK IS REDUNDANT.

      RANGEO = 0
      BETA = 0
      DO 200 I = 1, 2
!
!         IF (ANGLE .LE. 90.0000 .AND. I .EQ. 1) GO TO 200
         IF (ANGLE .LE. 90.0000 .AND. I .EQ. 1) CYCLE

!        IF (ANGLE .LE. 90.0000) THE PATH DOES NOT GO THROUGH HTAN.
!        IF (ANGLE .LE. 90.0000) THE I = 1 CALCULATION SHOULD NOT BE DON
!        IF (ANGLE .LE. 90.0000) FOR I = 2, R1 = H1
!
         IF (I .EQ. 1) THEN
            R1 = H1
            R2 = HTAN
         ELSEIF (I .EQ. 2) THEN
!            IF (HTAN .LT. 0.001 .AND. ANGLE .GT. 90) GO TO 200
            IF (HTAN .LT. 0.001 .AND. ANGLE .GT. 90) CYCLE
!
!           IF (HTAN APPROXIMATELY 0) THEN YOU ARE ABOUT TO HIT THE EART
!
            R2 = ZMAX
            IF (ANGLE .LE. 90.0000) THEN
               R1 = H1
            ELSE
               R1 = HTAN
            ENDIF
         ENDIF
         IF (R2 .LT. R1) THEN
            DZ = -DR
         ELSE
            DZ = DR
         ENDIF

         Z = R1
         DO 100 WHILE (Z.LT.R2)
            Z2=Z
            R=Z+RE
            CALL FINDSH( Z2, SH, GAMMA )
            RX = ANDEX( Z2, SH, GAMMA)
            STHETA = CPATH/(RX*R)
            IF (STHETA .GT. 1.0D0) STHETA = 1.0D0
            IF (STHETA .LT.-1.0D0) STHETA =-1.0D0
            SSAVE = STHETA
            CTHETA = SQRT(1.0D0-STHETA**2)
            IF (R1 .GT. R2) CTHETA = -CTHETA
!
!           IF (R1 .GT. R2) THEN CTHETA IS NEGATIVE BECAUSE THETA .GT. 9
!
            RATIO=-(RX*SH)/(RX-1.0D0)
            CAPRJ = -R/RATIO
            PNTGRN = 1.0D0/(1.0D0-CAPRJ*STHETA*STHETA)
            RPLDR = R+DZ
            Z2 = Z+DZ
            CALL FINDSH( Z2, SH, GAMMA )
            RX = ANDEX( Z2, SH, GAMMA )
            STHETA = CPATH/(RX*RPLDR)
            CTHET1 = CTHETA
            CTHETA = SQRT(1.0D0-STHETA**2)
            IF (R1 .GT. R2) CTHETA = -CTHETA
            DX=CTHETA*DZ+(CTHETA-CTHET1)*R
            DRNG = PNTGRN*DX
            RANGEO = RANGEO + DRNG
!
            DBETA = (((SSAVE+STHETA)*0.5D0) * (PNTGRN*DX)) /               &
     &              (Z-0.5D0*DZ+RE)
            BETA = BETA+DBETA
            IF (RANGEO .GE. RANGEI) THEN
               DIFF = (RANGEI-(RANGEO-DRNG))
               H2 = Z + (DZ/DRNG)*DIFF
               BETA = BETA*DEG
               IF (I .EQ. 2) THEN
                  LEN = 1
                  IF (ANGLE .LE. 90.0000D0) LEN = 0
                  IF (H2 .LT. HTAN) THEN
!
!                    THIS WILL BE THE CASE IF I = 2, AND YOU HAVE
!                    GONE THROUGH THE R-LOOP BARELY (ONLY) ONCE.
!
                     H2 = HTAN
                     LEN = 0
                  ENDIF
               ELSE
                  LEN = 0
               ENDIF
!
!              CORRECTION FOR VERY SHORT PATHS; HERE IT IS ABOUT 5 KM
!
               IF (RANGEI .LT. 5.0D0 .AND. RANGEO/RANGEI .GT. 1.05D0) THEN
!
!                 CALCULATE BETA BY STARIGHT LINE GEOMETRY.
!
                  PERP  = SIN(ANGLE/DEG)*RANGEI
                  BASE = COS(ANGLE/DEG)*RANGEI + RE+H1
                  BETA = ATAN(PERP/BASE)*DEG
                  RANGEO = RANGEI
!
!                 H2 = BASE - RE
!
                  H2 = COS(ANGLE/DEG)*RANGEI+H1
               ENDIF
               PHI = 180.0D0 - ACOS(CTHETA)*DEG
               RETURN
            ENDIF
            Z=Z+DZ
  100    CONTINUE
  200 END DO
!
!     COMES HERE IF YOU HAVE REACHED ZMAX, BUT YOUR RANGEI IS STILL
!     NOT EQUAL TO OUTPUT VALUE.
!     IN THIS CASE DO THE FOLLOWING.
!
      RANGEI = RANGEO
      H2 = ZMAX
      IF (ANGLE .LE. 90) THEN
         LEN = 0
      ELSE
         LEN = 1
      ENDIF
      IF (HTAN .LT. 0.001 .AND. ANGLE .GT. 90) THEN
!
!        YOU HAVE HIT THE EARTH IF YOU ARE AT THIS POINT OF THE CODE
!
         LEN = 0
         H2  = 0.0D0
      ENDIF
      BETA = BETA*DEG
      PHI = 180.0D0 - ACOS(CTHETA)*DEG
!
      RETURN
      END SUBROUTINE FNDPTH
!
!     ----------------------------------------------------------------
!
!
!     ----------------------------------------------------------------

      SUBROUTINE CMPALT( ILVL, PM, TM, DENW, REF_Z, ZMDL )

!**************************************************************
!     AUTHOR: TONY CLOUGH, JENNIFER DELAMERE, JOHN WARDEN
!             JANUARY 2001
!     PROGRAM TO CALCULATE ALTITUDE LEVEL (ZMDL) GIVEN
!     PRESSURE (PM), TEMPERATURE (TM) AND THE NUMBER DENSITY
!     OF WATER VAPOR (DENW) USING THE HYDROSTATIC EQUATION
!
!     INPUT:
!      A) PRESSURE (MBAR)
!      B) TEMPERATURE (KELVIN)
!      C) NUMBER DENSITY OF WATER VAPOR
!
!     OUTPUT:
!      A) ALTITUDE (KM)
!      IDEAL GAS LAW: CRIDOR (1996)
!**************************************************************
!
!      REAL PM(MXZMD),TM(MXZMD),DENW(MXZMD),ZMDL(MXZMD)
!      REAL H2O_MIXRAT(MXZMD),COMP_FACTOR(MXZMD),ZTEMP(MXZMD)

      INTEGER (4) :: ILVL, J, I

      REAL (8) :: Y, CHI0, T0, DT, C1, C2, C3, A, B, ALPHA, BTZ, XINT_TOT, REF_Z
      REAL (8) :: CA0, CA1, CA2, CB0, CB1, CC0, CC1, CD, CE, G0, CHIM, GAVE, DCHI
      REAL (8) :: XMASS_H2O, XMASS_DRY, XMASS_RATIO, RGAS, TOTAL_AIR, DRY_AIR

      REAL (8), DIMENSION(ILVL)  :: PM, TM, DENW, ZMDL, H2O_MIXRAT, COMP_FACTOR, ZTEMP

      DATA CA0/1.58123E-6/, CA1/-2.9331E-8/, CA2/1.1043E-10/
      DATA CB0/5.707E-6/, CB1/-2.051E-8/
      DATA CC0/1.9898E-4/, CC1/-2.376E-6/
      DATA CD/1.83E-11/, CE/-0.0765E-8/

      DATA XMASS_H2O/0.018015/,XMASS_DRY/0.0289654/
      DATA RGAS/ 8.31441 /, BTZ/ 1.380662E-23 /

! CALCULATE GRAVITY AT REFERENCE LATITUDE AT SURFACE
      G0 = 9.80612D0 - 0.02586D0*COS(2.0D0*PI*REF_LAT/180.0D0)

! CALCULATE THE NUMBER DENSITY OF TOTAL AIR MOLECULES [MOLEC/CM^3]
! CALCULATE THE COMPRESSIBILITY FACTOR (COMP_FAC) FOR THE
! IDEAL GAS LAW
      XMASS_RATIO = XMASS_H2O/XMASS_DRY
      DO 10 J=1,ILVL
         DT             = TM(J) - TZERO
         TOTAL_AIR      = PM(J)*1.0D-4 / (BTZ*TM(J))
         DRY_AIR        = TOTAL_AIR - DENW(J)
         H2O_MIXRAT(J)  = DENW(J) / DRY_AIR
         CHIM           = XMASS_RATIO*H2O_MIXRAT(J)
         COMP_FACTOR(J) = 1. - (PM(J)*100.0D0/TM(J)) *                  &
     &        (CA0 + CA1*DT + CA2*DT**2 +                               &
     &        (CB0 + CB1*DT)*CHIM + (CC0 + CC1*DT)*CHIM**2) +           &
     &        (CD + CE*CHIM**2)*(PM(J)*100.D0/TM(J))**2

   10 END DO

! CONVERT REFERENCE ALTITUDE TO METERS

      ZTEMP(1) = REF_Z * 1000.0D0
      ZMDL(1)  = REF_Z

      DO 20 I=1, ILVL - 1
         GAVE = G0 * (RE/(RE+ZTEMP(I)/1000.0D0))**2
         Y    =  LOG(PM(I+1)/PM(I))

         IF (Y .NE. 0.0) THEN
            CHI0 = H2O_MIXRAT(I)
            DCHI = (H2O_MIXRAT(I+1)-H2O_MIXRAT(I))/Y

            T0 = TM(I)
            DT = (TM(I+1) - TM(I))/Y

            C1 = T0 + T0*CHI0
            C2 = T0*DCHI + DT*CHI0 + DT
            C3 = DT*DCHI

            B = 1 + XMASS_RATIO*CHI0
            A = XMASS_RATIO*DCHI
            ALPHA = A/B

            IF ( ABS(ALPHA*Y) .GE. 0.01) THEN
               WRITE(16,*)'LAYER TOO THICK'
               WRITE(00,*)'LAYER TOO THICK'
               CALL SHUTDOWN
               STOP '3'
            ENDIF

            XINT_TOT = C1*Y + 0.5*(C2-C1*ALPHA)*Y**2 +                  &
     &           0.3333D0*(C3-C2*ALPHA+C1*ALPHA**2)*Y**3
            XINT_TOT =  -XINT_TOT*RGAS/(XMASS_DRY*GAVE*B)

            ZTEMP(I+1) = ZTEMP(I) + XINT_TOT*COMP_FACTOR(I)
            ZMDL(I+1) = ZTEMP(I+1)/1000.0D0
        ELSE
           ZTEMP(I+1) = ZMDL(I)*1000.0D0
           ZMDL(I+1) = ZMDL(I)
        ENDIF
   20    CONTINUE

         RETURN

      END SUBROUTINE CMPALT

!     -------------------------------------------------------------

      SUBROUTINE EXPINT (X,X1,X2,A)
!
!**********************************************************************
!     THIS SUBROUTINE EXPONENTIALLY INTERPOLATES X1 AND X2 TO X BY
!     THE FACTOR A
!**********************************************************************
!
      REAL (8)    :: X, X1, X2, A

      IF (X1.EQ.0.0.OR.X2.EQ.0.0) GO TO 10
      X = X1*(X2/X1)**A
!
      RETURN
!
   10 X = X1+(X2-X1)*A
!
      RETURN
!
      END SUBROUTINE EXPINT


! ----------------------------------------------------------------------
      SUBROUTINE ASTRO(KASTRO,HMIN,PHI,IERROR,BENDNG,H1,H2,ANGLE,RANGE,BETA,ITYPE,LEN,HOBS)

!   THIS SUBROUTINE TREATS ANGLE AS AN ASTRONOMICAL ANGLE ON INPUT AND
!   ITERATIVELY CALCULATES THE CORRESPONDING APPARENT ZENITH ANGLE.
!   THE APPARENT ZENITH ANGLE IS THEN STORED IN ANGLE
!   EPS-CONVERGENCE PARAMETER
!   JACK C LARSEN SASC

      INTEGER (4) :: IERROR, ITYPE, LEN, NOSP, K1, I, J, K, KASTRO, IEQ

      REAL (8)    :: HMIN, PHI, BENDNG, H1, H2, ANGLE, RANGE, BETA, HOBS
      REAL (8)    :: EPS, AST, APP
      REAL (8), DIMENSION(3)  :: EQ, EQ2, RO

      IF (NOPRNT .GE. 0)WRITE(IPR,103)
      EPS    = 1.D-8
      AST    = ANGLE
      KASTRO = 1
      NOSP   = NOPRNT
      !NOPRNT = 1
!   SET FIRST GUESS FOR APPARENT ANGLE
      APP    = ANGLE-BENDNG
      EQ2(1) = 0.0D0
      EQ2(2) = 0.0D0
      EQ2(3) = 0.0D0
      RO(1)  = APP
      K1     = 1
! THE *DO 30* LOOP ALLOWS 10 ITERATIONS FOR CONVERGENCE
      DO 30 J=1,10
         DO 20 K=K1,3
            IF(J.GT.1 .AND. K.EQ.2) GO TO 15
            BETA  = 0.D0
            RANGE = 0.D0
            CALL FSCGEO (H1,H2,APP,RANGE,BETA,ITYPE,LEN,HMIN,PHI,IERROR, HOBS)
            CALL RFPATH(H1,H2,APP,PHI,LEN,HMIN,0,RANGE,BETA,BENDNG)
            EQ(K)=RO(K)+BENDNG
!   RELATIVE ERROR
            EQ2(K)=(AST-EQ(K))/AST
            IF(ABS(EQ2(K)) .LE. EPS)GO TO 35
            IF(K.EQ.3)GO TO 25
            IF(K.EQ.1)RO(2) = RO(1) + 0.2D0*BENDNG
!   SECANT METHOD
   15       IF((K.EQ.2).AND.(ABS(EQ2(2)-EQ2(1)).LT. 1.D-12)) PRINT 101, J,K
            IF(K.EQ.2) RO(3)=RO(2)-EQ2(2)*(RO(2)-RO(1))/(EQ2(2)-EQ2(1))
            APP=RO(K+1)
   20    CONTINUE
! EQ2 IS CHECKED TO SEE IF IT MEETS THE CONVERGENCE CRITERIA.
   25    IF(ABS(EQ2(3)) .LE. EPS)GO TO 35
!   MOVE CLOSEST VALUES INTO EQ2(1),RO(1)
!   SHIFT LAST VALUE INTO EQ2(2),RO(2)
         IEQ = 2
         IF( ABS(EQ2(1)) .LT. ABS(EQ2(2)) )IEQ = 1
         EQ2(1) = EQ2(IEQ)
         EQ2(2) = EQ2(3)
         RO(1)  = RO(IEQ)
         RO(2)  = RO(3)
         K1     = 2
   30 END DO
      PRINT 107
   35 CONTINUE
      RANGE = 0.0D0
      BETA  = 0.0D0
      IF (NOPRNT .GE. 0)WRITE(IPR,102)AST,APP,BENDNG,J
      WRITE( 16,102)AST,APP,BENDNG,J
      IF (NOPRNT .GE. 0)WRITE(IPR,104)
      NOPRNT = NOSP
      ANGLE  = APP
      DO 50 I=1,IM2
         PPSUM(I)    = 0.0D0
         TPSUM(I)    = 0.0D0
         !RHOPSMT(I) = 0.0D0
         RHOPSM(I)   = 0.0D0
         DO 60 K=1,KDIM
            AMTP(K,I) = 0.0D0
   60    CONTINUE
   50 END DO

      RETURN

  101 FORMAT(2X,'DIVISION BY ZERO IN ANGLE PREDICTOR',/                 &
     &2X,4(I4,2X))
  102 FORMAT(/,1X,'ASTRONOMICAL ZENITH ANGLE =',F12.8,/,                  &
     &         1X,'    APPARENT ZENITH ANGLE =',F12.8,/,                  &
     &         1X,'             BENDNG ANGLE =',F12.8,/,                  &
     &         1X,'                    AFTER  ',I5,' ITERATIONS IN ASTRO')
  103 FORMAT(//,'BEGIN APPARENT ZENITH ANGLE CALCULATIONS...')
  104 FORMAT(1X,'APPARENT ZENITH ANGLE CALCULATIONS FINISHED.',/,/,/)
  107 FORMAT(1X,'EQUATION 2 DID NOT CONVERGE IN 10 ITERATIONS--HALT')

      END SUBROUTINE ASTRO

      SUBROUTINE spline (n, x, y, b, c, d)

      INTEGER n
      REAL (8) x (n), y (n), b (n), c (n), d (n)
!
!  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!
!    for  x(i) .le. x .le. x(i+1)
!
!  input..
!
!    n = the number of data points or knots (n.ge.2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!
!  output..
!
!    b, c, d  = arrays of spline coefficients as defined above.
!
!  using  p  to denote differentiation,
!
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!
!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!
!
      INTEGER nm1, ib, i
      REAL (8) t
!
      nm1 = n - 1
      IF (n.lt.2) return
      IF (n.lt.3) goto 50
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
      d (1) = x (2) - x (1)
      c (2) = (y (2) - y (1) ) / d (1)
      DO 10 i = 2, nm1
         d (i) = x (i + 1) - x (i)
         b (i) = 2. * (d (i - 1) + d (i) )
         c (i + 1) = (y (i + 1) - y (i) ) / d (i)
         c (i) = c (i + 1) - c (i)
   10 END DO
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
      b (1) = - d (1)
      b (n) = - d (n - 1)
      c (1) = 0.
      c (n) = 0.
      IF (n.eq.3) goto 15
      c (1) = c (3) / (x (4) - x (2) ) - c (2) / (x (3) - x (1) )
      c (n) = c (n - 1) / (x (n) - x (n - 2) ) - c (n - 2) / (x (n - 1) &
      - x (n - 3) )
      c (1) = c (1) * d (1) **2 / (x (4) - x (1) )
      c (n) = - c (n) * d (n - 1) **2 / (x (n) - x (n - 3) )
!
!  forward elimination
!
   15 DO 20 i = 2, n
         t = d (i - 1) / b (i - 1)
         b (i) = b (i) - t * d (i - 1)
         c (i) = c (i) - t * c (i - 1)
   20 END DO
!
!  back substitution
!
      c (n) = c (n) / b (n)
      DO 30 ib = 1, nm1
         i = n - ib
         c (i) = (c (i) - d (i) * c (i + 1) ) / b (i)
   30 END DO
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
      b (n) = (y (n) - y (nm1) ) / d (nm1) + d (nm1) * (c (nm1) + 2. *  &
      c (n) )
      DO 40 i = 1, nm1
         b (i) = (y (i + 1) - y (i) ) / d (i) - d (i) * (c (i + 1)      &
         + 2. * c (i) )
         d (i) = (c (i + 1) - c (i) ) / d (i)
         c (i) = 3. * c (i)
   40 END DO
      c (n) = 3. * c (n)
      d (n) = d (n - 1)
      RETURN
!
   50 b (1) = (y (2) - y (1) ) / (x (2) - x (1) )
      c (1) = 0.
      d (1) = 0.
      b (2) = b (1)
      c (2) = 0.
      d (2) = 0.
      RETURN
      END SUBROUTINE spline

      REAL (8) function seval (n, u, x, y, b, c, d)

      INTEGER n
      REAL (8) u, x (n), y (n), b (n), c (n), d (n)
!
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
      INTEGER i, j, k
      REAL (8) dx
      DATA i / 1 /
      IF (i.ge.n) i = 1
      IF (u.lt.x (i) ) goto 10
      IF (u.le.x (i + 1) ) goto 30
!
!  binary search
!
   10 i = 1
      j = n + 1
   20 k = (i + j) / 2
      IF (u.lt.x (k) ) j = k
      IF (u.ge.x (k) ) i = k
      IF (j.gt.i + 1) goto 20
!
!  evaluate spline
!
   30 dx = u - x (i)
      seval = y (i) + dx * (b (i) + dx * (c (i) + dx * d (i) ) )
      RETURN
      END FUNCTION seval





      END MODULE RAYTRACE
