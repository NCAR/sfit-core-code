      MODULE DATAFILES

      USE PARAMS

      IMPLICIT NONE

      CHARACTER*(IFLNMSZ) :: TFILE(100), LINDIR
!
!      								 TFIL07, &      ! - OBSOLETE- PT INPUT
!                              TFIL08, &      ! PBPFILE - OUTPUT SPECTRA, FITS AND DIFFERENCES
!                              TFIL09, &      ! ISOTOPE INPUT FILE - OBSOLETE MIX
!                              TFIL10, &      ! SOLAR LINE DATA
!                              TFIL11, &      ! FINAL SOLAR SPECTRUM
!                              TFIL12, &      ! - OBSOLETE - MASS PATH
!                              TFIL13, &      !
!                              TFIL14, &      ! INTERNAL TO LINEPARAM - CFGL FILE
!                              TFIL15, &      ! T15ASC SPECTRA INPUT
!                              TFIL16, &      ! DETAIL - THIS NAME CANNOT BE CHANGES IN SFT4.CTL FILE
!                              TFIL17, &      ! STATEVEC.MXF FINAL MIXING RATIOS
!                              TFIL18, &      ! STATEVEC
!                              TFIL19, &      ! SYNTHETIC SPECTRUM OUTPUT
!                              TFIL20, &      ! SUMMARY
!                              TFIL21, &      !
!                              TFIL22, &      ! STATEVEC.PRC FINAL PARTIAL COLUMNS
!                              TFIL23, &      ! EMPIRICAL MODULATION FUNCTION
!                              TFIL24, &      ! EMPIRICAL PHASE FUNCTION
!                              TFIL30, &      ! CHNSPEC.OUT1 - COMPUTE CHANNEL SPECTRUM
!                              TFIL31, &      ! - OBSOLETE - STATEVEC.ST - SHORT OUTPUT
!                              TFIL40, &      ! CHNSPEC.OUT2 - COMPUTE CHANNEL SPECTRUM
!                              TFIL61, &      ! SA OUTPUT
!                              TFIL62, &      ! SAINV INPUT
!                              TFIL63, &      ! COMPLETE SA OUTPUT
!                              TFIL64, &      ! FULL SA INPUT
!                              TFIL66, &      ! K.OUT JACOBIAN, SA, SE, SAINV
!                              TFIL67, &      ! SE.OUT
!                              TFIL68, &      ! KT.OUT INVERSE JACOBIAN OUTPUT
!                              TFIL70, &      ! DETAIL.OPT - OE DETAIL OUTPUT
!                              TFIL71, &      ! STATION.LAYERS INPUT LEVELS AND MIDPOINTS
!                              TFIL72, &      ! REFERENCE.PRF - ZPT + REFMOD INPUT VMR PROFILES
!                              TFIL73, &      ! RAYTRACE.OUT AKA TAPE6
!                              TFIL74, &      ! PRESSURE - TEMPERATURE
!                              TFIL75, &      ! MASS PATHS
!                              TFIL76, &      ! MIXING RATIOS
!                              TFIL77, &      ! LAYER BASED SA
!                              TFIL78, &      ! RAYTRACE PUNCH OUTPUT - NOT USED
!                              TFIL81, &      ! AK OUTPUT
!                              TFIL82, &      ! MEASUREMENT ERROR OUTPUT
!                              TFIL83, &      ! SMOOTHING ERROR OUTPUT
!                              TFIL84, &      ! AK SMOOTHING ERROR OUTPUT
!                              TFIL85, &      ! AK EIGENVECTORS ERROR OUTPUT
!                              TFIL87, &      ! TABLE OF APRIORI PROFILES AFTER RAYTRACE
!                              TFIL88, &      ! TABLE OF RETRIEVED PROFILES AFTER RAYTRACE
!                              TFIL89, &      ! STATEVECTOR PARAMETER / FACTORS BY ITERATION OUTPUT
!                              TFIL90, &      ! KB MATRIX
!                              TFIL91, &      ! SPECTRA BY ITERATION
!                              TFIL92, &      ! AB MATRIX
!                              TFIL93, &      ! GAIN MATRIX
!                              LINDIR, &

      INTEGER :: NCHAR

      CONTAINS

!----------------------------------------------------------------------
      SUBROUTINE FILESETUP

!      INTEGER :: I

! --- SET DEFAULT FILENAMES

! --- ASCII SPECTRA OUTPUT FITTED, OBSERVED & DIFFERENCE
      TFILE(08) = 'pbpfile'

! --- ISOTOPE SEPARATION FILE
      TFILE(09) = 'isotope.input'

! --- SOLAR INPUT DATA FILE
      TFILE(10) = 'solar_data.input'

! --- SOLAR OUTPUT FINAL SPECTRUMFILE
      TFILE(11) = 'solspec.dat'

! --- ASCII SPECTRA INPUT
      TFILE(15) = 't15asc.4'

      TFILE(16) = 'sfit4.dtl'

! --- STATEVEC
      TFILE(18) = 'statevec'

! --- SET MIXOUTPUT FILENAME AS STATE VECTOR FILENAME + .MXF
      TFILE(17) = TRIM(TFILE(18))//'.mxf'

! --- SAVED SYNTHETIC SPECTRUM
      TFILE(19) = 'synspec.out'

! --- SUMMARY
      TFILE(20) = 'summary'

! --- SET PARTICAL COLUMN OUTPUT FILENAME AS STATE VECTOR FILENAME + .PRC
      TFILE(22) = TRIM(TFILE(18))//'.prc'

! --- EMPIRICAL MODULATION FUNCTION
      TFILE(23) = 'ils.dat'

! --- EMPIRICAL PHASE FUNCTION
      TFILE(24) = 'ils.dat'

! --- CHANNEL OUTPUT
      TFILE(30) = 'chnspec1.output'

      TFILE(40) = 'chnspec2.output'

! --- SET SHORT FILENAME AS SUMMARY FILENAME + .ST
!      TFILE(31) = TRIM(TFILE(20))//'.st'

! --- DEFAULT FIRST RETRIEVAL GAS SA OUTPUT FILENAME
! --- OUTPUT DEPENDS ON NAME SA .NORM, .VMR, .PCOL
      TFILE(61) = 'sa.out'

! --- DEFAULT SAINV INPUT FILENAME
      TFILE(62) = 'sainv.input'

! --- COMPLETE AS IMPLEMENTED SA OUTPUT FILENAME
      TFILE(63) = 'sa.complete'

! --- COMPLETE A POSTERIORI COVARIANCE SHAT
      TFILE(64) = 'shat.complete'

! --- DEFAULT K,SE,SA MATRICES OUTPUT FILENAME
      TFILE(66) = 'k.out'

! --- SE OUTPUT
      TFILE(67) = 'seinv.out'

! --- DEFAULT K TRANSPOSED FILENAME
      TFILE(68) = 'kt.out'

! --- DEFAULT SAINV MATRIX FILENAME
      TFILE(69) = 'sainv.out'

! ---  LEVENBERG-MARQUARDT DETAILS
      TFILE(70) = 'detail.opt'

! --- LAYERING SCHEME - OUTPUT FROM WACCM PROFILES
      TFILE(71) = 'station.layers'

! --- REFERENCE VMR PROFILES
      TFILE(72) = 'reference.prf'

! --- RAYTRACE DETAILED OUTPUT AKA TAPE6
      TFILE(73) = 'raytrace.out'

! --- RAYTRACE OUTPUT THAT MATCHES OLD SFIT2 INPUT - DEPRICATED
      TFILE(74) = 'raytrace.pt'
      TFILE(75) = 'raytrace.ms'
      TFILE(76) = 'raytrace.mix'
      TFILE(77) = 'raytrace.sa'

! --- MORE RAYTRACE OUTPUT
      TFILE(78) = 'raytrace.pnch'

! --- RESERVED FOR GASOUT NAME CHANGES - SEE FRWDMDL.F90
      !TFILE(80)

! --- AVERAGING KERNELS OF THE TARGET GAS
      TFILE(81) = 'ak.out'

! --- MEASUREMENT ERROR GIVEN SNR AS IMPLEMENTED
      TFILE(82) = 'smeas.target'

! --- SMOOTHING ERROR GIVEN SA AS IMPLEMENTED
      TFILE(83) = 'ssmooth.target'

! --- TABLE OF APRIORI PROFILES
      TFILE(87) = 'aprfs.table'

! --- TABLE OF APRIORI AND RETRIEVED PROFILES
      TFILE(88) = 'rprfs.table'

! --- FITTED PARAMETER AKA STATEVECTOR BY ITERATION
      TFILE(89) = 'parm.vectors'

! --- DERIVATIVE OF PARAMETER B'S
      TFILE(90) = 'kb.out'

      TFILE(91) = 'spec.out'

! --- DY#KB
      TFILE(92) = 'ab.out'

! --- DY
      TFILE(93) = 'd.complete'

      RETURN

      END SUBROUTINE FILESETUP


      SUBROUTINE FILEOPEN (INDEX, SW )

      INTEGER, INTENT(IN) :: INDEX, SW
      INTEGER             :: IOS

! --- OPEN NORMAL FILES 31, 18, 20, 73

      IF(( INDEX .LT. 1 ) .OR. ( INDEX .GT. 100 ))THEN
         PRINT *, 'DATAFILES: FILEOPEN: LUN INDEX : ', INDEX, '  OUT OF RANGE.'
         STOP 'ERROR'
      ENDIF

      SELECT CASE (SW)
         CASE (1)
! --- NEW FILE FOR OUTPUT
            OPEN(UNIT=INDEX, FILE=TFILE(INDEX), STATUS='UNKNOWN', IOSTAT=IOS, POSITION='ASIS')
            IF( IOS .NE. 0 )THEN
               WRITE(16,106) INDEX, TRIM(TFILE(INDEX)), IOS
               WRITE( 0,106) INDEX, TRIM(TFILE(INDEX)), IOS
               STOP 'FILEOPEN ERROR...SEE DETAIL FILE'
            ENDIF

         CASE( 2 )
! --- NEW FILE FOR OUTPUT - REPLACE
            OPEN(UNIT=INDEX, FILE=TFILE(INDEX), STATUS='REPLACE', IOSTAT=IOS, POSITION='ASIS')
            IF( IOS .NE. 0 )THEN
               WRITE(16,106) INDEX, TRIM(TFILE(INDEX)), IOS
               WRITE( 0,106) INDEX, TRIM(TFILE(INDEX)), IOS
               STOP 'FILEOPEN ERROR...SEE DETAIL FILE'
            ENDIF

         CASE( 3 )
! --- OLD FILE FOR INPUT
            OPEN(UNIT=INDEX, FILE=TFILE(INDEX), STATUS='OLD', IOSTAT=IOS, POSITION='ASIS')
            IF( IOS .NE. 0 )THEN
               WRITE(16,105) INDEX, TRIM(TFILE(INDEX)), IOS
               WRITE( 0,105) INDEX, TRIM(TFILE(INDEX)), IOS
               STOP 'FILEOPEN ERROR...SEE DETAIL FILE'
            ENDIF

         CASE( 4 )
! --- FILE FOR OUTPUT - APPEND
            OPEN(UNIT=INDEX, FILE=TFILE(INDEX), STATUS='OLD', IOSTAT=IOS, POSITION='APPEND')
            IF( IOS .NE. 0 )THEN
               WRITE(16,106) INDEX, TRIM(TFILE(INDEX)), IOS
               WRITE( 0,106) INDEX, TRIM(TFILE(INDEX)), IOS
               STOP 'FILEOPEN ERROR...SEE DETAIL FILE'
            ENDIF

         CASE DEFAULT
             PRINT *, 'DATAFILES: FILEOPEN: SWITCH, ', SW, '  OUT OF RANGE.'
             STOP 'ERROR'
       END SELECT


 105  FORMAT(/,' FILEOPEN: INPUT FILE OPEN ERROR-UNIT : ',I5, ' FILENAME: "',A,'"', ' IOSTAT: ', I5)
 106  FORMAT(/,' FILEOPEN: OUTPUT FILE OPEN ERROR-UNIT : ',I5, ' FILENAME: "',A,'"', ' IOSTAT: ', I5)


      END SUBROUTINE FILEOPEN


      SUBROUTINE FILECLOSE( INDEX, SW )

      INTEGER, INTENT(IN) :: INDEX, SW
      LOGICAL             :: FILOPEN = .FALSE.
      INTEGER             :: IOS = 0

      IF(( INDEX .LT. 1 ) .OR. ( INDEX .GT. 100 ))THEN
         PRINT *, 'DATAFILES: FILECLOSE: LUN INDEX : ', INDEX, '  OUT OF RANGE.'
         STOP 'ERROR'
      ENDIF


      SELECT CASE( SW )
         CASE( 1 )
! --- OUTPUT FILE
            INQUIRE( UNIT=INDEX, OPENED=FILOPEN )
            IF( FILOPEN )CLOSE(UNIT=INDEX, IOSTAT=IOS)
            IF( IOS .NE. 0 )THEN
               WRITE(16,106) INDEX, TRIM(TFILE(INDEX)), IOS
               WRITE( 0,106) INDEX, TRIM(TFILE(INDEX)), IOS
               STOP 'FILECLOSE ERROR...SEE DETAIL FILE'
            ENDIF

         CASE( 2 )
! --- INPUT FILE
            INQUIRE( UNIT=INDEX, OPENED=FILOPEN )
            IF( FILOPEN )CLOSE(UNIT=INDEX, IOSTAT=IOS)
            IF( IOS .NE. 0 )THEN
               WRITE(16,105) INDEX, TRIM(TFILE(INDEX)), IOS
               WRITE( 0,105) INDEX, TRIM(TFILE(INDEX)), IOS
               STOP 'FILECLOSE ERROR...SEE DETAIL FILE'
            ENDIF


         CASE DEFAULT
             PRINT *, 'DATAFILES: FILECLOSE: SWITCH, ', SW, '  OUT OF RANGE.'
             STOP 'ERROR'
       END SELECT

 105  FORMAT(/,' FILECLOSE: INPUT FILE CLOSE ERROR-UNIT : ',I5, ' FILENAME: "',A,'"', ' IOSTAT: ', I5)
 106  FORMAT(/,' FILECLOSE: OUTPUT FILE CLOSE ERROR-UNIT : ',I5, ' FILENAME: "',A,'"', ' IOSTAT: ', I5)


      END SUBROUTINE FILECLOSE

      END MODULE DATAFILES
