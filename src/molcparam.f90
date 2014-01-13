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

      MODULE molcparam

! June 2013
! add isotopes to N2O # 6 7 8 / from Toth via ggg linelist

! sfit4 v003.92
! remove pre-separated isotopes
! 53 ch3d
! 49 hdo
! 54 o3668:
! 55 o3686:
! 56 o3667:
! 57 o3676:

! moltotal 70->99 May 2010(sfit2) then June 2012(sfit4 v003.9)
! added PSL:

!id      name        mass           nmode    symbol      # lines     result
!71      CH3COOH     60.05          18                   14957       insert.3
!72      C5H8        68.12          33                   3001        insert.2
!73      MVK         70.09          27       C4H6O       9986        insert.5
!74      MACR        70.09          27       C4H6O       8968        insert.4
!75      C3H6        42.08          21                   9990        insert.0
!76      C4H8        56.11          30                   64965       insert.1

! November 2009
!     Updated for sfit v394 & HITRAN 2008

! 01/09/07
!     Added CH3CHO

! 12/09/05
!     Added completed PAN

! 10/31/05
!     Added CH3OH, CH3CN, C2H6PL, and PAN (PAN not complete)

! 7/18/05
!     Changed mass of 6th isotope of CO from 20 to 30
!     Changed to match raytrace 2.05
!           CFCLO to COCLF
!           CCL2O to COCL2
!           CHFCL2 to CHCL2F

! Changed made 6/22/04 at DU
! Changed COCLF to CFCLO and COCL2 to CCL2O

      USE params

      IMPLICIT NONE

      INTEGER, PARAMETER                      :: NI=9
      REAL(DOUBLE), DIMENSION(MOLTOTAL)       :: TDEP    ! DEFAULT COEFFICIENT TEMP.
                                                         ! DEP. AIR H-W
      REAL(DOUBLE), DIMENSION(NI,MOLTOTAL)    :: XMASS   ! MOLECULAR MASS FOR ALL ISOTOPES
      INTEGER, DIMENSION(MOLTOTAL)            :: NHIISO  ! NUMBER OF ISOTOPES
      CHARACTER (LEN=7), DIMENSION(MOLTOTAL)  :: NAME    ! NAMES OF GASES (CHAR*7)
      REAL(DOUBLE)                            :: TAUMIN  ! CROSS-SECTION CUTOFF
      INTEGER, PRIVATE                        :: J,JISO,IMOLM

      DATA TAUMIN/ 1.D-6 /

!     REVISION DATE: 21 Mar 2005
!      Updated TDEP values for NO, OH, CLO, H2O2, and F142B.
!
!     REVISION DATE: 11 Feb 2002
!      copied dblocks from sfit109g to this file
!
!     REVISION DATE:  8 NOV 2001:
!      added pseudo lines for 11 molecules so:
!      updated NMODE, IVIB for molecs: 26,27,31,32,33,35,42,43,50,61,62
!      and 39 though no new lines
!     REVISION DATE:  2 NOV 2001:
!      update IVIB for CCL2F2 from Goeff's
!      isotopomer.dat file
!      All molecular mass values in this file should
!      be actual, dummy values for pseudolines are set
!      in binput
!       REVISION DATE:  30 OCT 2001, moved to 70 molecules
!       REVISION DATE:  7/10/01, added isotope 6 to CO
!      added isotopes 2,3 to H2S
!      added isotopes 3 to HDO
!      copied data from sfit2:data.f
!       REVISION DATE:  3/8/01, 3/26/01 Added two H2O isotopes,
!               one HDO isotope and one OCS isotope.
!       REVISION DATE:  FEB 19, 1992

!  --- MOLECULES IN ATMOS ORDER

      DATA (NAME(J),J=1,MOLTOTAL) / &
        'H2O',     'CO2',     'O3',     'N2O',     'CO',      &
        'CH4',     'O2',      'NO',     'SO2',     'NO2',     & ! 10
        'NH3',     'HNO3',    'OH',     'HF',      'HCL',     &
        'HBR',     'HI',      'CLO',    'OCS',     'H2CO',    & ! 20
        'HOCL',    'HO2',     'H2O2',   'HONO',    'HO2NO2',  &
        'N2O5',    'CLONO2',  'HCN',    'CH3F',    'CH3CL',   & ! 30
        'CF4',     'CCL2F2',  'CCL3F',  'CH3CCL3', 'CCL4',    &
        'COF2',    'COCLF',   'C2H6',   'C2H4',    'C2H2',    & ! 40
        'N2',      'CHF2CL',  'COCL2',  'CH3BR',   'CH3I',    &
        'HCOOH',   'H2S',     'CHCL2F', 'O2CIA',   'SF6',     & ! 50
        'NF3',     'OTHER',   'OTHER',  'OTHER',   'OTHER',   &
        'OTHER',   'OTHER',   'OCLO',   'F134A',   'C3H8',    & ! 60
        'F142B',   'CFC113',  'F141B',  'CH3OH',   'CH3CNPL', &
        'C2H6PL',  'PAN',     'CH3CHO ','CH3CN',   'OTHER',   & ! 70
        'CH3COOH', 'C5H8',    'MVK',    'MACR',    'C3H6',    &
        'C4H8',    'OTHER',   'OTHER',  'OTHER',   'OTHER',   & ! 80
        'OTHER',   'OTHER',   'OTHER',  'OTHER',   'OTHER',   &
        'OTHER',   'OTHER',   'OTHER',  'OTHER',   'OTHER',   & ! 90
        'OTHER',   'OTHER',   'OTHER',  'OTHER',   'OTHER',   &
        'OTHER',   'OTHER',   'OTHER',  'OTHER'               / ! 99

      DATA (TDEP(J),J=1,MOLTOTAL)/ &
            1.5,       1.0,       1.5,       1.0,       1.0, &
            1.5,       1.0,      1.23,       1.5,       1.5, & ! 10
            1.5,       1.5,      1.05,       1.0,       1.0, &
            1.0,       1.0,      1.23,       1.0,       1.5, & ! 20
            1.5,       1.5,       2.0,       1.5,       2.0, &
            2.0,      2.23,       1.0,       1.5,       1.5, & ! 30
            1.5,       1.5,       1.5,       1.5,       1.5, &
            1.5,       1.5,       1.9,       1.5,       1.0, & ! 40
            1.0,       1.5,       1.5,       1.5,       1.5, &
            1.5,       1.5,       1.5,       1.0,       1.5, & ! 50
            1.5,       0.0,       0.0,       0.0,       0.0, &
            0.0,       0.0,       1.5,       1.5,       1.5, & ! 60
            2.0,       2.0,       1.5,       1.5,       1.5, &
            1.9,       1.5,       1.5,       1.5,       0.0, & ! 70
            1.5,       1.5,       1.5,       1.5,       1.5, &
            1.5,       0.0,       0.0,       0.0,       0.0, & ! 80
            0.0,       0.0,       0.0,       0.0,       0.0, &
            0.0,       0.0,       0.0,       0.0,       0.0, & ! 90
            0.0,       0.0,       0.0,       0.0,       0.0, &
            0.0,       0.0,       0.0,       0.0             / ! 99

      DATA ((XMASS(JISO,IMOLM),JISO=1,NI),IMOLM=1,10)/ &
       18.D0,  20.D0,  19.D0,  19.D0,  21.D0,  20.D0,  3*0.D0,       & !H2O
       44.D0,  45.D0,  46.D0,  45.D0,  47.D0,  46.D0,  48.D0,  47.D0, 49.D0, & !CO2
       48.D0,  50.D0,  50.D0,  49.D0,  49.D0, 4*0.D0,                & !O3
!       48.D0,  8*0.D0,                                               & !O3
       44.D0,  45.D0,  45.D0,  46.D0,  45.D0,  46.D0,  47.D0,  47.D0, 0.D0,  & !N2O
       28.D0,  29.D0,  30.D0,  29.D0,  31.D0,  30.D0, 3*0.D0,        & !CO
       16.D0,  17.D0,  17.D0,  18.D0,  5*0.D0,                       & !CH4
       32.D0,  34.D0,  33.D0, 6*0.D0,                                & !O2
       30.D0,  31.D0,  32.D0, 6*0.D0,                                & !NO
       64.D0,  66.D0, 7*0.D0,                                        & !SO2
       46.D0, 8*0.D0/                                                  !NO2

      DATA ((XMASS(JISO,IMOLM),JISO=1,NI),IMOLM=11,20)/ &
       17.D0,  18.D0, 7*0.D0,                                        & !NH3
       63.D0, 8*0.D0,                                                & !HNO3
       17.D0,  19.D0,  18.D0, 6*0.D0,                                & !OH
       20.D0,  22.D0, 7*0.D0,                                        & !HF
       36.D0,  38.D0, 38.D0,  40.D0, 5*0.D0,                         & !HCL
       80.D0,  82.D0, 7*0.D0,                                        & !HBR
      128.D0, 8*0.D0,                                                & !HI
       51.D0,  53.D0, 7*0.D0,                                        & !CLO
       60.D0,  62.D0,  61.D0,  61.D0,  62.D0,  4*0.D0,               & !OCS
       30.D0,  31.D0,  32.D0, 6*0.D0/                                  !H2CO

      DATA ((XMASS(JISO,IMOLM),JISO=1,NI),IMOLM=21,30)/ &
       52.D0,  54.D0, 7*0.D0,                                        & !HOCL
       33.D0,  8*0.D0,                                               & !HO2
       34.D0,  8*0.D0,                                               & !H2O2
       47.D0,  8*0.D0,                                               & !HONO
       79.D0,  8*0.D0,                                               & !HO2NO2
      108.D0,  8*0.D0,                                               & !N2O5      PL
       97.D0,  99.D0, 7*0.D0,                                        & !CLONO2    PL
       27.D0,  28.D0,  28.D0, 6*0.D0,                                & !HCN
       34.D0,  8*0.D0,                                               & !CH3F
       50.D0,  52.D0, 7*0.D0/                                          !CH3CL

      DATA ((XMASS(JISO,IMOLM),JISO=1,NI),IMOLM=31,40)/ &
       88.D0,  8*0.D0,                                               & !CF4       PL
      1.21D2,  8*0.D0,                                               & !CCL2F2    PL
      1.36D2,  8*0.D0,                                               & !CCL3F     PL
      1.44D2,  8*0.D0,                                               & !CH3CCL3
      1.52D2,  8*0.D0,                                               & !CCL4      PL
       66.D0,  8*0.D0,                                               & !COF2
       82.D0,  8*0.D0,                                               & !COCLF     PL
       30.D0,  31.D0, 7*0.D0,                                        & !C2H6
       28.D0,  29.D0, 7*0.D0,                                        & !C2H4
       26.D0,  27.D0, 7*0.D0/                                          !C2H2

      DATA ((XMASS(JISO,IMOLM),JISO=1,NI),IMOLM=41,50)/ &
       28.D0,  8*0.D0,                                               & !N2
       86.D0,  8*0.D0,                                               & !CHF2CL
       98.D0,  8*0.D0,                                               & !COCL2
       95.D0,  97.0D0, 7*0.D0,                                       & !CH3BR
      1.40D2,  8*0.D0,                                               & !CH3I
       46.D0,  8*0.D0,                                               & !HCOOH
       34.D0,  36.D0, 35.D0, 6*0.D0,                                 & !H2S
      1.03D2,  8*0.D0,                                               & !CHCL2F
       32.D0,  32.0D0, 7*0.D0,                                       & !O2CIA     PL
      1.46D2,  8*0.D0/                                                 !SF6       PL

      DATA ((XMASS(JISO,IMOLM),JISO=1,NI),IMOLM=51,60)/ &
       71.D0,  8*0.D0,                                               & !NF3       PL
       9*0.D0,                                                       & !OTHER
       9*0.D0,                                                       & !OTHER
       9*0.D0,                                                       & !OTHER
       9*0.D0,                                                       & !OTHER
       9*0.D0,                                                       & !OTHER
       9*0.D0,                                                       & !OTHER
       67.D0,  8*0.D0,                                               & !OCLO
       83.D0,  8*0.D0,                                               & !F134A
       44.D0,  8*0.D0/                                                 !C3H8      PL

      DATA ((XMASS(JISO,IMOLM),JISO=1,NI),IMOLM=61,70)/ &
      100.D0,  8*0.D0,                                               & !F142B     PL
      187.D0,  8*0.D0,                                               & !CFC113    PL
      117.D0,  8*0.D0,                                               & !F141B     PL
     32.04D0,  8*0.D0,                                               & !CH3OH
     41.05D0,  8*0.D0,                                               & !CH3CNPL   PL
       30.D0,  8*0.D0,                                               & !C2H6PL    PL
      121.D0,  8*0.D0,                                               & !PAN       PL
       44.D0,  8*0.D0,                                               & !CH3CHO    PL
     41.05D0,  8*0.D0,                                               & !CH3CN     PL
        0.D0,  8*0.D0/                                                 !OTHER

      DATA ((XMASS(JISO,IMOLM),JISO=1,NI),IMOLM=71,80)/ &
       60.D0,  8*0.D0,                                               & !CH3COOH   PL
       68.D0,  8*0.D0,                                               & !C5H8      PL
       70.D0,  8*0.D0,                                               & !MVK       PL
       70.D0,  8*0.D0,                                               & !MACR      PL
       42.D0,  8*0.D0,                                               & !C3H6      PL
       56.D0,  8*0.D0,                                               & !C4H8      PL
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0/                                                 !OTHER

      DATA ((XMASS(JISO,IMOLM),JISO=1,NI),IMOLM=81,90)/ &
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0/                                                 !OTHER

      DATA ((XMASS(JISO,IMOLM),JISO=1,NI),IMOLM=91,MOLTOTAL)/ &
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0,                                               & !OTHER
        0.D0,  8*0.D0/                                                 !OTHER 99


      DATA NHIISO/ &
        6,  9,  5,  8,  6,  4,  3,  3,  2,  1,                       & !1-10
        2,  1,  3,  2,  4,  2,  1,  2,  5,  3,                       & !11-20
        2,  1,  1,  1,  1,  1,  2,  3,  1,  2,                       & !21-30
        1,  1,  1,  1,  1,  1,  1,  2,  2,  2,                       & !31-40
        1,  1,  1,  2,  1,  1,  3,  1,  2,  1,                       & !41-50
        1,  1,  0,  0,  0,  0,  0,  1,  1,  1,                       & !51-60
        1,  1,  1,  1,  1,  1,  1,  1,  1,  0,                       & !61-70
        1,  1,  1,  1,  1,  1,  0,  0,  0,  0,                       & !71-80
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,                       & !81-90
        0,  0,  0,  0,  0,  0,  0,  0,  0/                             !91-99

      END MODULE MOLCPARAM
