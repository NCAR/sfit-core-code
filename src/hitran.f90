module hitran

   use params
   use datafiles
   use bandparam
   use lineparam
   use isotope

   implicit none

   integer, parameter     :: nhit=99, ngal=2, ncia=2, nglines=100000, flagoff=280,nlmx=5, nsdv=5 
   integer, parameter     :: nlmlines=100000, nsdlines=100000
   real(8), parameter     :: weps = 1.0d-6
   integer                :: hnml, gnml, lnml, snml
   integer                :: stlun, map(nhit)
   integer :: nhit_files, nlm_files, ngal_files, nsdv_files
   character (len=255), dimension(nhit) :: hitran_files
   character (len=255), dimension(ngal) :: gal_files
   character (len=255), dimension(nsdv) :: sdv_files
   character (len=255), dimension(nlmx) :: lm_files
   
   type, public :: hitranfile
      integer             :: mo        ! molecule id from subdir name for this file --- replaces map()
      integer             :: flag      ! 0-hitran, 1-cia f&s
      integer             :: lun       ! unit # for this record
      character (len=300) :: buf       ! read buffer
   end type hitranfile

   type, public :: galatrydata
      integer            :: n          ! number of data lines / file(gas)
      integer            :: lun        ! unit for this file
      character (len=255) :: buf        ! read buffer
      integer            :: mo(nglines)    ! mol id
      integer            :: is(nglines)    ! isotope id #
!      real(8)            :: nu(nglines)    ! wavenumber
      real(8)            :: bt(nglines)    ! intensity [cm-1/(molec/cm-2)]
      character(len=60)  :: qa(nglines)    ! quanta data
      ! The format of the quanta fields depends on the species. Refer to
      ! Rothman, L. S. et.al.
      ! The HITRAN 2004 molecular spectroscopic database
      ! Journal Of Quantitative Spectroscopy & Radiative Transfer, 2005, 96, 139-204
      real(8)            :: g0_air(nglines)
      real(8)            :: td_g0_air(nglines)
      real(8)            :: beta(nglines)
      real(4)            :: g0_self(nglines)
      real(4)            :: s_air(nglines)
      real(4)            :: g2_air(nglines)
      real(4)            :: ts_air(nglines)
      real(4)            :: lm_air(nglines) ! line mixing coefficients
      real(4)            :: lm_t1(nglines)  ! extra parameters of F. Hase to model
      real(4)            :: lm_t2(nglines)  ! temperature dependency
   end type galatrydata

   type, public :: linemixfile
      integer            :: n          ! number of data lines / file(gas)
      integer            :: ist        ! record start
      integer            :: lun        ! unit for this file
      character (len=64) :: buf        ! read buffer
      integer            :: mo         ! mol id
   end type linemixfile

   type :: linemixdata
      integer            :: mo         ! mol id
      real(8)            :: nu         ! wavenumber
      real(4)            :: dt(3)      ! intensity [cm-1/(molec/cm-2)]
   end type linemixdata

! this map is obsolete
! the linelist directory structure is the key to the gas names and the molecule id numbers
! those id's and names must be the same in the reference.prf file
! eg a files containing hitran lines is read from one subdir in linelist then the molid will be changed
! to the 2digit integer 0NN of the subdir name and assumed to be for gas 0NN_abcdef
! --- CODE NUMBERS TO CONVERT FROM HITRAN TO ATMOS / SFIT
!     INDEX IS HITRAN ID NUMBER, VALUE IS SFIT
!     SFIT     HITRAN         SFIT        PSEUDOLINES FILE MOLID/ISO
!     #        # NAME         # NAME
      DATA MAP / &
      1,    &! 1 H2O
      2,    &! 2 CO2
      3,    &! 3 O3
      4,    &! 4 N2O
      5,    &! 5 CO
      6,    &! 6 CH4
      7,    &! 7 O2
      8,    &! 8 NO
      9,    &! 9 SO2
      10,   &! 10 NO2
      11,   &! 11 NH3
      12,   &! 12 HNO3
      13,   &! 13 OH
      14,   &! 14 HF
      15,   &! 15 HCL
      16,   &! 16 HBR
      17,   &! 17 HI
      18,   &! 18 CLO
      19,   &! 19 OCS
      20,   &! 20 H2CO
      21,   &! 21 HOCL
      41,   &! 22 N2          HO2
      28,   &! 23 HCN         H202
      30,   &! 24 CH3CL       HONO
      23,   &! 25 H2O2        HO2NO2
      40,   &! 26 C2H2        N2O5
      38,   &! 27 C2H6        CLONO2
      0,    &! 28 PH3         HCN
      36,   &! 29 COF2        CH3F
      50,   &! 30 SF6         CH3CL
      47,   &! 31 H2S         CF4
      46,   &! 32 HCOOH       CCL2F2
      22,   &! 33 HO2         CCL3F3
      0,    &! 34 O           CH3CCL3
      35,   &! 35 CLONO2      CCL4          PS 35/1
      0,    &! 36 NO+         COF2
      0,    &! 37 HOBR        COCLF
      39,   &! 38 C2H4        C2H6
      64,   &! 39 CH3OH*      C2H4
      44,   &! 40 CH3BR*      C2H2
      69,   &! 41 CH3CN*      N2
      31,   &! 42 CF4*        CHF2CL
      43,   &! 43             COCL2
      44,   &! 44             CH3BR
      45,   &! 45             CH3I
      46,   &! 46             HCOOH
      47,   &! 47             H2S
      48,   &! 48             CHCL2F
      49,   &! 49             O2CIA
      50,   &! 50             SF6
      51,   &! 51             NF3
      52,   &! 52             OTHER
      53,   &! 53             OTHER
      54,   &! 54             OTHER
      55,   &! 55             OTHER
      56,   &! 56             OTHER
      57,   &! 57             OTHER
      58,   &! 58             OCLO
      59,   &! 59             F134A
      60,   &! 60             C3H8
      61,   &! 61             F142B
      62,   &! 62             CFC113
      63,   &! 63             F141B
      64,   &! 64             CH3OH
      65,   &! 65             CH3CNPL
      66,   &! 66             C2H6PL
      67,   &! 67             PAN
      68,   &! 68             CH3CHO
      69,   &! 69             CH3CN
      70,   &! 70             OTHER
      71,   &! 71             OTHER
      72,   &! 72             OTHER
      73,   &! 73             OTHER
      74,   &! 74             OTHER
      75,   &! 75             OTHER
      76,   &! 76             OTHER
      77,   &! 77             OTHER
      78,   &! 78             OTHER
      79,   &! 79             OTHER
      80,   &! 80             OTHER
      81,   &! 81             OTHER
      82,   &! 82             OTHER
      83,   &! 83             OTHER
      84,   &! 84             OTHER
      85,   &! 85             OTHER
      86,   &! 86             OTHER
      87,   &! 87             OTHER
      88,   &! 88             OTHER
      89,   &! 89             OTHER
      90,   &! 90             OTHER
      91,   &! 91             OTHER
      92,   &! 92             OTHER
      93,   &! 93             OTHER
      94,   &! 94             OTHER
      95,   &! 95             OTHER
      96,   &! 96             OTHER
      97,   &! 97             OTHER
      98,   &! 98             OTHER
      99    &! 99             OTHER
      /

! map to 1 record in the ascii linelist file made by hbin
! byte range  quantity
! 1   - 160   hitran - see HITRAN docs
! 161 - 172   galatry beta
! 173 - 184   sdv gamma0
! 185 - 196   sdv gamma2
! 197 - 208   sdv eta2
! 209 - 220   lmx ltk1
! 221 - 232   lmx ltk2
! 233 - 244   lmx ylm
! 245 - 256   ---
! 257 - 268   ---
! 269 - 280   ---
! 281 - 288   logical flags for: galatry, fcia, scia, sdv, lmix, 6, 7, 8
! 289 - 300   <nothing>

! from lineparam.f90: GALATRY_FLAG=1, FCIA_FLAG=2, SCIA_FLAG=3, SDV_FLAG=4, LM_FLAG=5

end module hitran

