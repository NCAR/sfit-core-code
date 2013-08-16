      MODULE BANDPARAM

      USE PARAMS

      IMPLICIT NONE

!      LOGICAL, DIMENSION(MAXBND) :: EMISSION_NORM

      CHARACTER (LEN=7), DIMENSION(MAXBND,MOLMAX) :: GASB = ''

      CHARACTER(80), DIMENSION(MAXSPE) :: STITLE
      INTEGER      :: NATMOS = 0                ! SUM OF ALL OBSERVED SPECTRA DATA POINTS
      INTEGER      :: NBAND  = 0                ! NUMBER OF UNIQUE BANDPASSES
      INTEGER      :: NZERO  = 0
      INTEGER      :: NFITS  = 0                ! TOTAL NUMBER OF SPECTRAL SEGMENTS FITTED
      INTEGER      :: NKZERO = 0

      LOGICAL, DIMENSION(MAXBND)      :: TRETB
      INTEGER, DIMENSION(MAXBND)      :: NBANDS = -1

      REAL(DOUBLE) :: DLINES = 4.0D0            ! INTERVAL OUTSIDE REGION FOR INPUT OF LINE DATA [cm-1]
      REAL(DOUBLE) :: DELNU                     ! HALF-WIDTH OF INTERVAL FOR CALCULATING EACH LINE SHAPE
      REAL(DOUBLE) :: SNR = 1.0                 ! DEFAULT YIKES!

      REAL(DOUBLE), DIMENSION(MAXBND) :: WAVE3  ! BEGINNING WAVENUMBER FOR INTEGRATION
      REAL(DOUBLE), DIMENSION(MAXBND) :: WAVE4  ! ENDING WAVENUMBER FOR INTEGRATION
      REAL(DOUBLE), DIMENSION(MAXBND) :: WAVE5  ! BEGINNING WAVENUMBER FOR LINE LIST
      REAL(DOUBLE), DIMENSION(MAXBND) :: WAVE6  ! ENDING WAVENUMBER FOR LINE LIST
      REAL(DOUBLE), DIMENSION(MAXBND) :: WSTART ! FIRST W# CALCULATED READING T15 IN BAND
      REAL(DOUBLE), DIMENSION(MAXBND) :: WSTOP  ! LAST  W# CALCULATED READING T15 IN BAND
      REAL(DOUBLE), DIMENSION(MAXBND) :: WMON   ! STARTING wave# of band ???? WAVE1
      REAL(DOUBLE), DIMENSION(MAXBND) :: SPAC   ! POINT SPACING IN BAND
      REAL(DOUBLE), DIMENSION(MAXBND) :: DN     ! POINT SPACING FOR MONOCHROMATIC CALCULATIONS
      REAL(DOUBLE), DIMENSION(MAXBND) :: PMAX   ! MAXIMUM PATH DIFFERENCE [CM]
      REAL(DOUBLE), DIMENSION(MAXBND) :: PMAX0   ! MAXIMUM PATH DIFFERENCE [CM]
      REAL(DOUBLE), DIMENSION(MAXBND) :: OMEGA  ! FIELD OF VIEW IN [MRAD]
      REAL(DOUBLE), DIMENSION(MAXBND) :: OMEGA0 ! FIELD OF VIEW IN [MRAD]
      REAL(DOUBLE), DIMENSION(MAXBND) :: WAVFAC ! FACTOR OF WAVENAMUBER FOR EACH BAND

      REAL(DOUBLE), DIMENSION(MAXBND,MAXSPE)  :: SCNSNR  ! SNR BY SCAN AND BAND

      INTEGER, DIMENSION(MAXBND,MOLMAX)       :: IGASB ! ID OF GASES RETRIEVED IN BAND I
      INTEGER, DIMENSION(MAXBND,MOLMAX)       :: NGASB ! NRET INDEX OF GAS RETRIEVED IN BAND I
      INTEGER, DIMENSION(MOLMAX,0:2,0:MAXBND) :: NGIDX = 0
																	! INDEX REGION IN PARM OF GAS I AND FLAG BAND
																	! I,0,K IS FLAG THAT GAS I IS IN BAND K
																	! I,1,K IS LOWER INDEX IN PARM
																	! I,2,K IS UPPER INDEX IN PARM

      INTEGER, DIMENSION(MAXBND,MAXSPE) :: ISCAN ! NOT!!!! SPECTRUM INDEX FOR THE I BAND J SCAN (1000*SZA)
      INTEGER, DIMENSION(MAXSPE) :: ISPEC        ! SPECTRUM CODE FOR THE J SCAN (1000*SZA)
      INTEGER, DIMENSION(2,MAXBND,MAXSPE) :: ISCNDX ! LOW, HIGH INDEX OF SPECTRA IN K BY BAND AND SCAN

      INTEGER, DIMENSION(MAXBND) :: NSPAC ! RATIO OF OBS TO MONO POINT SPACING
      INTEGER, DIMENSION(MAXBND) :: NM    ! NUMBER OF POINTS IN TCALC BY BAND
      INTEGER, DIMENSION(MAXBND) :: NPRIM ! NUMBER OF POINTS IN OBS SPEC BY BAND
      INTEGER, DIMENSION(MAXBND) :: NRETB ! NUMBER OF MOLEC RETRIEVED BY BAND
      INTEGER, DIMENSION(MAXSPE) :: NSCAN ! NUMBER OF SCANS TO BE FIT FOR A GIVEN BANDPASS (IBAND)
      INTEGER, DIMENSION(MOLMAX) :: MINZ  ! MINIMUM ALTITUDE RETRIEVED FOR
                                          ! RETRIEVAL GAS I
      INTEGER, DIMENSION(MOLMAX) :: MAXZ  ! MAXIMUM ALTITUDE RETRIEVED FOR
                                          ! RETRIEVAL GAS I
      INTEGER, DIMENSION(MOLMAX) :: NUMZ  ! NUMBER OF ALTITUDES TO BE FIT FOR RETRIEVAL GAS I
!        THESE NEXT 2 MAY NOT BE USED
      INTEGER, DIMENSION(MOLMAX) :: IRETS ! SPECTRUM NUMBER OF FIRST SCAN RETRIEVED FOR GAS I
      INTEGER, DIMENSION(MOLMAX) :: IRETF ! SPECTRUM NUMBER OF FINAL SCAN RETRIEVED FOR GAS I


      SAVE

      END MODULE BANDPARAM
