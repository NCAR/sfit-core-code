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

      MODULE RETVPARAM

      USE PARAMS

      CHARACTER (LEN=1)               :: EMISSION_OBJECT  ! 'M' FOR MOON, OTHERS NOT YET INCLUDED

      CHARACTER (LEN=7), DIMENSION(MOLTOTAL) :: GAS

      LOGICAL                         :: F_RTPHASE = .FALSE.
      LOGICAL                         :: F_RTAPOD = .FALSE.
      INTEGER(4)                      :: NMOL  ! MAX MOLECULES READ IN SEE RAYTRACE:LNGMDL
      INTEGER                         :: IFPS
      LOGICAL                         :: EMISSION, EMISSION_NORM, IFTEMP=.FALSE., FPS, RAYTONLY
      LOGICAL                         :: CONVERGE=.FALSE., DIVWARN=.FALSE.

      LOGICAL, DIMENSION(MOLMAX)      :: CORRELATE = .FALSE.  ! OFF-AXIS CORRELATION FOR PROFILE RETRIEVAL
      CHARACTER (len=3), DIMENSION(MOLMAX)      :: REGMETHOD = 'OEM'  ! METHOD OF CONSTRAINT (OEM AND TP SUPPORTED)
      INTEGER                         :: NRPRFGAS             ! NR OF GASES RETRIEVED AS PROFILES
      LOGICAL, DIMENSION(MOLMAX)      :: IFPRF = .FALSE.      ! PROFILE RETRIEVAL
      LOGICAL, DIMENSION(MOLMAX)      :: IFPRF_KB = .FALSE.   ! PROFILE FOR KB MATRIX CALCULATION
      LOGICAL                         :: IFDIFF = .FALSE.     ! DIFFERENTIAL WAVENUMBER SHIFT
      LOGICAL                         :: IFCALCSE = .FALSE.   ! USE RMS AS 1/SE ON PER SCAN
      LOGICAL                         :: F_WSHIFT             ! RETRIEVAL OF WAVENUMBER SHIFT
      LOGICAL                         :: F_BACKG              ! RETRIEVAL OF BACKGROUND SLOPE OR CURVE
      LOGICAL, DIMENSION(MAXBND)      :: F_ZSHIFT
      LOGICAL                         :: F_CONTABS = .FALSE. ! RETRIEVE CONTINUUM
      INTEGER, DIMENSION(MAXBND)      :: IZERO
      REAL(DOUBLE), DIMENSION(MAXBND,MAXSPE) :: ZSHIFT ! ZERO SHIFTS BY BAND AND SCAN
      REAL(DOUBLE), DIMENSION(MAXBND) :: SZERO  !

      INTEGER                         :: IEMISSION
      REAL(DOUBLE)                    :: EMISSION_T_BACK ! CONTAINS THE BACKGROUND TEMEPRATURE IF EMISSION

      INTEGER, DIMENSION(MAXSPE)      :: KZTAN
      INTEGER, DIMENSION(MOLMAX)      :: ISHIFT
      INTEGER, DIMENSION(MAXGAS)      :: ICODE, IFMIX, ISCODE
      INTEGER, DIMENSION(MOLMAX)      :: IRET
      INTEGER, DIMENSION(MOLMAX)      :: IGAS
      INTEGER, DIMENSION(MOLMAX)      :: IFOFF = 0
      INTEGER :: IFLINE, IFSZA, IFFOV


      REAL(DOUBLE), DIMENSION(MOLMAX) :: COLSF  = 0.0D0  ! SCALE FACTOR FOR APRIORI VMR COLUMN RETRIEVAL
      REAL(DOUBLE), DIMENSION(MOLMAX) :: SCOLSF = 0.0D0  ! SIGMA FOR APRIORI VMR COLUMN RETRIEVAL
      REAL(DOUBLE), DIMENSION(LAYMAX) :: TSIGMA      ! SIGMA FOR TEMPERATURE RETRIEVAL

      INTEGER :: ICOUNT, NRET=0, NGAS, ISPARM, NBACK=1, NBKFIT, &
                 NSHIFT, NSPEC, NDIFF, NPHASE
      INTEGER :: NTEMP=0, NTEMP1=0, NILINE=0,NPLINE=0,NTLINE=0,NRLGAS=0
      INTEGER :: NSOLAR=0, NSOLAR1=0

      LOGICAL :: IFPHASE = .FALSE.
      LOGICAL :: LOG_STATEV(MAXGAS) = .FALSE.
      INTEGER :: ILOGRETRIEVAL(MAXGAS)
      INTEGER :: ITRMAX = 0
      INTEGER :: NLAYERS = 0                             ! NUMBER OF LAYERS READ IN FROM SFIT4.CTL FILE

      REAL(DOUBLE), DIMENSION(LAYMAX)        :: PMASMX ! MAXIMUM MASS
      REAL(DOUBLE), DIMENSION(MOLMAX,LAYMAX) :: XORG   ! INITIAL VALUES MIXING RATIOS - FITTED
      REAL(DOUBLE), DIMENSION(MOLMAX,LAYMAX) :: X      ! WORKING VALUES MIXING RATIOS - FITTED
      REAL(DOUBLE), DIMENSION(MAXGAS,LAYMAX) :: XGAS   ! INITIAL VALUES MIXING RATIOS - ALL
      INTEGER                                :: KMAX

!      REAL(DOUBLE) DIMENSION( LAYMAX, MAXPRF ) :: ALTFIT
!      REAL(DOUBLE) DIMENSION( LAYMAX, MAXPRF ) :: XINIT
      REAL(DOUBLE), DIMENSION(LAYMAX,MAXPRF) :: XFIT
      REAL(DOUBLE), DIMENSION(LAYMAX,MAXPRF) :: SIG
      REAL(DOUBLE), DIMENSION(MOLMAX)        :: ZWID
      REAL(DOUBLE), DIMENSION(MOLMAX)        :: ZGMIN
      REAL(DOUBLE), DIMENSION(MOLMAX)        :: ZGMAX

      REAL(DOUBLE), ALLOCATABLE :: CCC(:,:), CORG(:,:)    ! MASS PATH
      REAL(DOUBLE), ALLOCATABLE :: P(:)         ! WEIGHTED PRESSURE (ATM)
      REAL(DOUBLE), ALLOCATABLE :: PORG(:)      ! WEIGHTED PRESSURE (ATM) ORIGINAL
      REAL(DOUBLE), ALLOCATABLE :: T(:)         ! WEIGHTED TEMPERATURE (K)
      REAL(DOUBLE), ALLOCATABLE :: TORG(:)      ! WEIGHTED TEMPERATURE ARRAY (K) ORIGINAL
      REAL(DOUBLE), ALLOCATABLE :: PMB(:)       ! WEIGHTED PRESSURE (MB)
      REAL(DOUBLE), ALLOCATABLE :: PMBORG(:)    ! WEIGHTED PRESSURE (MB) ORIGINAL
      REAL(DOUBLE), ALLOCATABLE :: FXORG(:,:)   ! WEIGHTED MIXING RATIOS ALL ORIGINAL

      REAL(DOUBLE), ALLOCATABLE :: Z(:)         ! RETRIEVAL BOUNDRARIES INCLUDE TOP AND BOTTOM
      REAL(DOUBLE), ALLOCATABLE :: ZBAR(:)      ! RETRIEVAL GRID MIDPOINTS

      REAL(DOUBLE), DIMENSION (:,:), ALLOCATABLE :: SNR_CLC, SNR_THE


! CELL - Only one cell for now

      INTEGER :: NCELL = 0                      ! NUMBER OF CELL OPTICAL PATHS, READ FROM SFIT4.CTL
      INTEGER :: NPATH = 0                      ! NUMBER OF OPTICAL PATHS CELL & ATMOSPHERE = NCELL + NLAYERS
      INTEGER, DIMENSION(2) :: NCELLS           ! CELL PATH ID - not implemented!!!!!!!
      LOGICAL, DIMENSION(MOLMAX)      :: IFCELL = .FALSE.      ! THIS IS A CELL GAS IN ARRAY 1 TO NRET
      CHARACTER (LEN=7), DIMENSION(2) :: CGAS
      INTEGER, DIMENSION(2)           :: CGASID
      REAL(DOUBLE),  DIMENSION(2)     :: CTEMP
      REAL(DOUBLE),  DIMENSION(2)     :: CPRES
      REAL(DOUBLE),  DIMENSION(2)     :: CVMR
      REAL(DOUBLE),  DIMENSION(2)     :: CPATH
      REAL(DOUBLE), DIMENSION(MOLMAX) :: L1LAMBDA = -1.0D0     ! Regularisation strength
                                                               ! for L1 (TP) method

      LOGICAL :: F_MEAS_TRANSMIS
      INTEGER :: NUM_FILTER
      REAL (DOUBLE), ALLOCATABLE :: FILTERTRANS (:,:)


    CONTAINS

      real(double) function FTRANS(WNUM)
        ! INTERPOLATES THE CURVE IN FILTERTRANS TO THE WAVENUMBER WNUM

        implicit none

        integer :: nr, ind_min, ind_max
        real(double) :: wnum
        logical :: flag

        flag = .false.
        do nr = 1,num_filter
           if (filtertrans(1,nr).le.wnum &
                .and.filtertrans(1,nr+1).gt.wnum) then
              ind_min = nr
              ind_max = nr+1
              flag = .true.
              exit
           end if
        end do

        if (.not.flag) then
           WRITE(16,*) "FILTER MEASUREMENT DOES NOT COVER MEASUREMENT, SET TO 1.0"
           WRITE(0,*) "FILTER MEASUREMENT DOES NOT COVER MEASUREMENT, SET TO 1.0"
           ftrans = 1.0d0
           return
        end if


        ftrans = ((wnum - filtertrans(1,ind_min)) * filtertrans(2,ind_max) &
             + ((filtertrans(1,ind_max) - wnum) * filtertrans(2,ind_min))) &
             / (filtertrans(1,ind_max) - filtertrans(1,ind_min))

        return
      end function FTRANS


      SUBROUTINE RELEASE_MEM_RTP

      IF( ALLOCATED( P ))DEALLOCATE( P )
      IF( ALLOCATED( PORG ))DEALLOCATE( PORG )
      IF( ALLOCATED( T ))DEALLOCATE( T )
      IF( ALLOCATED( TORG ))DEALLOCATE( TORG )
      IF( ALLOCATED( PMB ))DEALLOCATE( PMB )
      IF( ALLOCATED( PMBORG ))DEALLOCATE( PMBORG )
      IF( ALLOCATED( Z ))DEALLOCATE( Z )
      IF( ALLOCATED( ZBAR ))DEALLOCATE( ZBAR )
      IF( ALLOCATED( SNR_CLC ))DEALLOCATE( SNR_CLC )
      IF( ALLOCATED( SNR_THE ))DEALLOCATE( SNR_THE )
      IF( ALLOCATED( CCC ))DEALLOCATE( CCC )
      IF( ALLOCATED( CORG ))DEALLOCATE( CORG )
      IF( ALLOCATED( FILTERTRANS ))DEALLOCATE( FILTERTRANS )

      END SUBROUTINE RELEASE_MEM_RTP

      END MODULE RETVPARAM
