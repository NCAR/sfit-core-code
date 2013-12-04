      MODULE PARAMS

      INTEGER, PARAMETER :: BYTE_LOG = SELECTED_INT_KIND(2)
      INTEGER, PARAMETER :: SHORT_LOG = SELECTED_INT_KIND(4)
      INTEGER, PARAMETER :: LONG_LOG = SELECTED_INT_KIND(18)
      INTEGER, PARAMETER :: BYTE = SELECTED_INT_KIND(2)
      INTEGER, PARAMETER :: SHORT = SELECTED_INT_KIND(4)
      INTEGER, PARAMETER :: LONG = SELECTED_INT_KIND(18)
      INTEGER, PARAMETER :: DOUBLE = SELECTED_REAL_KIND(13)
      INTEGER, PARAMETER :: EXTENDED = SELECTED_REAL_KIND(30)
      INTEGER, PARAMETER :: DOUBLE_EXT = SELECTED_REAL_KIND(50)
      INTEGER, PARAMETER :: DBLE_COMPLEX = SELECTED_REAL_KIND(13)
      INTEGER, PARAMETER :: EXT_COMPLEX = SELECTED_REAL_KIND(30)

!  --- VERSION STRING - NO SPACES: USE AS FIRST WORD IN OUTPUT FILE TAG
!      CHARACTER(len=10) :: BUILDDATE = BDATE
      CHARACTER (LEN=78),  PARAMETER :: VERSION = 'SFIT4:V0.9.5.0.PRE_RELEASE:20131202'
      CHARACTER (LEN=100)            :: TAG              ! OUTPUT FILE TAG IS VERSION + RUNTIME

!      LOGICAL, PARAMETER :: BUG = .TRUE.
      LOGICAL, PARAMETER :: BUG = .FALSE.
      LOGICAL, PARAMETER :: ANALYTIC_K = .TRUE.
!      LOGICAL, PARAMETER :: ANALYTIC_K = .FALSE.

      ! mp: For emission some of the predefault numbers seem to small (LNMAX, MAXGAS,LAYMAX)

      INTEGER, PARAMETER :: MMAX = 2*131072       ! MAXIMUM NUMBER OF SPECTRAL DATA POINTS
      INTEGER, PARAMETER :: NMAX = 255            ! MAXIMUM NUMBER OF FITTING PARAMETERS
      INTEGER, PARAMETER :: MOLMAX = 10           ! MAXIMUM NUMBER OF RETRIEVAL GASES
      INTEGER, PARAMETER :: LAYMAX = 100          ! MAXIMUM NUMBER OF ATMOSPHERIC LAYERS
      INTEGER, PARAMETER :: MAXSPE = 40           ! MAXIMUM NUMBER OF SPECTRA (SCANS)
      INTEGER, PARAMETER :: MAXGAS = 100          ! MAXIMUM NUMBER OF GASES IN BANDPASS
      INTEGER, PARAMETER :: LNMAX = 8*131072      ! MAXIMUM NUMBER OF LINES
      INTEGER, PARAMETER :: MAXPRF = 10            ! MAXIMUM NUMBER OF PROFILE RETRIEVALS
      INTEGER, PARAMETER :: MAXBND = 125          ! MAXIMUM NUMBER OF BANDPASSES
!      INTEGER MAXCROSS                           ! MAXIMUM NUMBER OF MONOCHROMATIC
!                                                 !  POINTS/BANDPASS (NOW SET AT RUNTIME)
      INTEGER, PARAMETER :: MAXSNR = 20           ! MAXIMUM NUMBER OF ALTERNATE S/N VALUES
      INTEGER, PARAMETER :: LNMXCO = 131072       ! MAXIMUM NUMBER OF SOLAR CO LINES
      INTEGER, PARAMETER :: MAXEAP = 2000         ! EMPIRICAL APODIZTION POINTS
      INTEGER, PARAMETER :: MOLTOTAL = 99         ! TOTAL NUMBER OF POSSIBLE MOLECULES
      INTEGER, PARAMETER :: MAX_NUM_OF_BEAMS = 20 ! MAX NUMBER OF BEAMS FOR A BANDPASS  !pwj
                                                  ! WUJIAN PENG ADDED ON JULY, 2002     !PWJ
                                                  ! INCREASED TO 20 BY MATHIAS PALM !MP
      INTEGER, PARAMETER :: ISOMAX = 15           ! MAXIMUM NUMBER OF SEPARATED ISOTOPE SPECIES
      INTEGER, PARAMETER :: IFLNMSZ = 255         ! MAX LENGTH OF DATAFILE NAMES


! Numerical Constants
! Phys. Today Aug 2001


      REAL(DOUBLE), PARAMETER :: PI = 3.141592653589793D0
      REAL(DOUBLE), PARAMETER :: SCHMIT = 2.6867775D+19
      REAL(DOUBLE), PARAMETER :: ALOGSQ = 0.8325546111576978D0
      REAL(DOUBLE), PARAMETER :: PISQ = 1.7724538509083950D0
      REAL(DOUBLE), PARAMETER :: RFACTOR = 3.581165292D-7
      REAL(DOUBLE), PARAMETER :: RCONST2 = 1.4387752D0
      REAL(DOUBLE), PARAMETER :: STDTEMP = 296.0D0
      REAL(DOUBLE), PARAMETER :: ZEROC = 273.150D0
      REAL(DOUBLE), PARAMETER :: BAR = 1013.250D0
! Bronstein 1995 - german edition
      REAL(DOUBLE), PARAMETER :: c_boltz = 1.380648D-23   ! Boltzmann's constant
      REAL(DOUBLE), PARAMETER :: v_light = 299792458.0D0  ! speed of light
      REAL(DOUBLE), PARAMETER :: c_planck = 6.58211928D-34! Planck's constant





!  ********************************************************************
!   CMPEPSILON is the epsilon used in floating point comparisons
!              instead of checking for equality between floating point
!              numbers.
!              instead of if (x .eq. y) we use
!              if (ABS(x-y) .LE. (ABS(x+y)*CMPEPSILON))
      REAL(DOUBLE), PARAMETER :: CMPEPSILON = 0.0000000000000005

      REAL(DOUBLE):: PLANCK_C1,PLANCK_C2

      END MODULE PARAMS
