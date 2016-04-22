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

module tips

! Nov 2009 updated from HITRAN release version TIPS_2009

  USE params, ONLY:  DOUBLE

  IMPLICIT NONE

  INTEGER, PARAMETER :: NMOL = 38
  INTEGER, PARAMETER :: MAX_ISO = 20
  INTEGER, DIMENSION(NMOL,MAX_ISO) :: ISO82
  INTEGER, DIMENSION(NMOL) :: ISONM
  INTEGER :: ITIPS, JTIPS
!-----------------------------------------------
!
!    The number of isotopes for a particular molecule:
  DATA (ISONM(ITIPS),ITIPS=1,NMOL)/ &
!     H2O, CO2, O3, N2O, CO, CH4, O2,
       6,   9,  18,   5,  6,   3,  3, &
!      NO, SO2, NO2, NH3, HNO3, OH, HF, HCl, HBr, HI,
       3,   2,   1,   2,    1,  3,  1,   2,   2,  1, &
!     ClO, OCS, H2CO, HOCl, N2, HCN, CH3Cl, H2O2, C2H2, C2H6, PH3 !4/24/97
       2,   5,    3,    2,  1,   3,     2,    1,    2,    1,   1, &
!     COF2, SF6, H2S, HCOOH, HO2, O, ClONO2,  NO+, HOBr, C2H4
        1,   1,   3,     1,   1, 1,      2,    1,    2,    2/
!


!       H2O
  DATA (ISO82(1,JTIPS),JTIPS=1,6)/ 161, 181, 171, 162, 182, 172/

!       CO2
  DATA (ISO82(2,JTIPS),JTIPS=1,9)/ 626, 636, 628, 627, 638, 637, 828, 728, 727/

!       O3
  DATA (ISO82(3,JTIPS),JTIPS=1,18)/ 666, 668, 686, 667, 676, 886, 868, 678, 768, &
       786, 776, 767, 888, 887, 878, 778, 787, 777/

!       N2O
  DATA (ISO82(4,JTIPS),JTIPS=1,5)/ 446, 456, 546, 448, 447/

!       CO
  DATA (ISO82(5,JTIPS),JTIPS=1,6)/ 26, 36, 28, 27, 38, 37/

!       CH4
  DATA (ISO82(6,JTIPS),JTIPS=1,3)/ 211, 311, 212/

!       O2
  DATA (ISO82(7,JTIPS),JTIPS=1,3)/ 66, 68, 67/

!       NO
  DATA (ISO82(8,JTIPS),JTIPS=1,3)/ 46, 56, 48/

!       SO2
  DATA (ISO82(9,JTIPS),JTIPS=1,2)/ 626, 646/

!       NO2
  DATA (ISO82(10,JTIPS),JTIPS=1,1)/ 646/

!       NH3
  DATA (ISO82(11,JTIPS),JTIPS=1,2)/ 4111, 5111/

!       HNO3
  DATA (ISO82(12,JTIPS),JTIPS=1,1)/ 146/

!       OH
  DATA (ISO82(13,JTIPS),JTIPS=1,3)/ 61, 81, 62/

!       HF
  DATA (ISO82(14,JTIPS),JTIPS=1,1)/ 19/

!       HCl
  DATA (ISO82(15,JTIPS),JTIPS=1,2)/ 15, 17/

!       HBr
  DATA (ISO82(16,JTIPS),JTIPS=1,2)/ 19, 11/

!       HI
  DATA (ISO82(17,JTIPS),JTIPS=1,1)/ 17/

!       ClO
  DATA (ISO82(18,JTIPS),JTIPS=1,2)/ 56, 76/

!       OCS
  DATA (ISO82(19,JTIPS),JTIPS=1,5)/ 622, 624, 632, 623, 822/

!       H2CO
  DATA (ISO82(20,JTIPS),JTIPS=1,3)/ 126, 136, 128/

!       HOCl
  DATA (ISO82(21,JTIPS),JTIPS=1,2)/ 165, 167/

!       N2
  DATA (ISO82(22,JTIPS),JTIPS=1,1)/ 44/

!       HCN
  DATA (ISO82(23,JTIPS),JTIPS=1,3)/ 124, 134, 125/

!       GH3Cl
  DATA (ISO82(24,JTIPS),JTIPS=1,2)/ 215, 217/

!       H2O2
  DATA (ISO82(25,JTIPS),JTIPS=1,1)/ 1661/

!       C2H2
  DATA (ISO82(26,JTIPS),JTIPS=1,2)/ 1221, 1231/

!       C2H6
  DATA (ISO82(27,JTIPS),JTIPS=1,1)/ 1221/

!       PH3
  DATA (ISO82(28,JTIPS),JTIPS=1,1)/ 1111/

!       COF2
  DATA (ISO82(29,JTIPS),JTIPS=1,1)/ 269/

!       SF6
  DATA (ISO82(30,JTIPS),JTIPS=1,1)/ 29/

!       H2S
  DATA (ISO82(31,JTIPS),JTIPS=1,3)/ 121, 141, 131/

!       HCOOH
  DATA (ISO82(32,JTIPS),JTIPS=1,1)/ 126/

!       HO2
  DATA (ISO82(33,JTIPS),JTIPS=1,1)/ 166/

!       O
  DATA (ISO82(34,JTIPS),JTIPS=1,1)/ 6/

!       ClONO2
  DATA (ISO82(35,JTIPS),JTIPS=1,2)/ 5646, 7646/

!       NO+
  DATA (ISO82(36,JTIPS),JTIPS=1,1)/ 46/

!       HOBr
  DATA (ISO82(37,JTIPS),JTIPS=1,2)/ 169, 161/

!       C2H4
  DATA (ISO82(38,JTIPS),JTIPS=1,2)/ 221, 231/



  CHARACTER (LEN=6), DIMENSION(0:NMOL) :: MOLID
!-----------------------------------------------
!
  DATA (MOLID(ITIPS),ITIPS=0,NMOL)/ '   All', '   H2O', '   CO2', '    O3', &
       '   N2O', '    CO', '   CH4', '    O2', '    NO', '   SO2', '   NO2', &
       '   NH3', '  HNO3', '    OH', '    HF', '   HCl', '   HBr', '    HI', &
       '   ClO', '   OCS', '  H2CO', '  HOCl', '    N2', '   HCN', ' CH3Cl', &
       '  H2O2', '  C2H2', '  C2H6', '   PH3', '  COF2', '   SF6', '   H2S', &
       ' HCOOH', '   HO2', '     O', 'ClONO2', '   NO+', '  HOBr', '  C2H4'/
!

  INTEGER, PARAMETER :: NT = 119
  REAL(DOUBLE), DIMENSION(NT) :: TDAT
  DATA TDAT/ 60., 85., 110., 135., 160., 185., 210., 235., 260., 285., 310.&
       , 335., 360., 385., 410., 435., 460., 485., 510., 535., 560., 585., &
       610., 635., 660., 685., 710., 735., 760., 785., 810., 835., 860., 885.&
       , 910., 935., 960., 985., 1010., 1035., 1060., 1085., 1110., 1135., &
       1160., 1185., 1210., 1235., 1260., 1285., 1310., 1335., 1360., 1385., &
       1410., 1435., 1460., 1485., 1510., 1535., 1560., 1585., 1610., 1635., &
       1660., 1685., 1710., 1735., 1760., 1785., 1810., 1835., 1860., 1885., &
       1910., 1935., 1960., 1985., 2010., 2035., 2060., 2085., 2110., 2135., &
       2160., 2185., 2210., 2235., 2260., 2285., 2310., 2335., 2360., 2385., &
       2410., 2435., 2460., 2485., 2510., 2535., 2560., 2585., 2610., 2635., &
       2660., 2685., 2710., 2735., 2760., 2785., 2810., 2835., 2860., 2885., &
       2910., 2935., 2960., 2985., 3010./


contains

!***********************
!
      SUBROUTINE BD_TIPS_2008(MOL, TEMP, ISO, GI, QT)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!***********************
!    Subroutine BD_TIPS_2003 written by R.R. Gamache
!
!
!...date last changed 25        March, 2005
!  --  UPDATES  --
!    18 December 2003 better vibrational fundamentals for PH3
!
!    This program calculates the total internal
!    partition sum (TIPS) for a given molecule, isotopic species, and
!    temperature.  Current limitations are the molecular species on the
!    HITRAN molecular database and the temperature range 70 - 3000 K.
!
!...This program calculates the TIPS by 4-point LaaGrange interpolation
!
!..  JQSRT - HITRAN Special Issue, 2003.
!..  J. Fischer(a) R.R. Gamache(a&), A. Goldman(b),L.S. Rothman(c), and
!..
!..  (a)  Department of Environmental, Earth, and Atmospheric Sciences,
!..       University of Massachusetts Lowell, Lowell, MA 01854, U.S.A.
!..
!..  (b)  Department of Physics, University of Denver, Denver, CO 80208, U.S.A.
!..
!..  (c)  Harvard-Smithsonian Center for Astrophysics, 60 Garden St, Cambridge, MA 02138 USA
!..
!..  (d)  Laboratoire de Photophysique Mol�culaire, Universit� Paris Sud, 91405 Orsay, FRANCE
!..
!..  &  Corresponding author. Email address: Robert_Gamache@uml.edu
!..  Abstract
!..   Total internal partition sums (TIPS) are calculated for all molecular species in
!..  the 2000 HITRAN database.  In addition, the TIPS for 13 other isotopomers/isotopologues
!..  of ozone and carbon dioxide are presented.  The calculations address the corrections
!..  suggested by Goldman et al. (JQSRT 2000;66:55-86).  The calculations consider the
!..  temperature range 70-3000 K to be applicable to a variety of remote sensing needs.
!..  The method of calculation for each molecular species is stated and
!..  data from the literature are discussed.  A new method of recall for the partition sums,
!..  Lagrange 4-point interpolation, is developed.  This method, unlike
!..  the TIPS code, allows all molecular species to be considered.
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  15:55:42   3/17/05
!...Switches: -yf12
!++
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: MOL ! HITRAN molecule number
      INTEGER  :: ISO             ! isotopomer index
      REAL(DOUBLE)  :: TEMP       ! temperature in K
      REAL(DOUBLE)  :: GI         ! state independent degeneracy factor
      REAL(DOUBLE)  :: QT         ! total internal partition sum
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      !CHARACTER :: STOPNGO*30
!-----------------------------------------------
!
!
!     GO TO PARTICULAR MOLECULE:

      SELECT CASE (MOL)
      CASE (1)
         CALL QT_H2O (TEMP, ISO, GI, QT)
      CASE (2)
         CALL QT_CO2 (TEMP, ISO, GI, QT)
      CASE (3)
         CALL QT_O3 (TEMP, ISO, GI, QT)
      CASE (4)
         CALL QT_N2O (TEMP, ISO, GI, QT)
      CASE (5)
         CALL QT_CO (TEMP, ISO, GI, QT)
      CASE (6)
         CALL QT_CH4 (TEMP, ISO, GI, QT)
      CASE (7)
         CALL QT_O2 (TEMP, ISO, GI, QT)
      CASE (8)
         CALL QT_NO (TEMP, ISO, GI, QT)
      CASE (9)
         CALL QT_SO2 (TEMP, ISO, GI, QT)
      CASE (10)
         CALL QT_NO2 (TEMP, ISO, GI, QT)
      CASE (11)
         CALL QT_NH3 (TEMP, ISO, GI, QT)
      CASE (12)
         CALL QT_HNO3 (TEMP, ISO, GI, QT)
      CASE (13)
         CALL QT_OH (TEMP, ISO, GI, QT)
      CASE (14)
         CALL QT_HF (TEMP, ISO, GI, QT)
      CASE (15)
         CALL QT_HCL (TEMP, ISO, GI, QT)
      CASE (16)
         CALL QT_HBR (TEMP, ISO, GI, QT)
      CASE (17)
         CALL QT_HI (TEMP, ISO, GI, QT)
      CASE (18)
         CALL QT_CLO (TEMP, ISO, GI, QT)
      CASE (19)
         CALL QT_OCS (TEMP, ISO, GI, QT)
      CASE (20)
         CALL QT_H2CO (TEMP, ISO, GI, QT)
      CASE (21)
         CALL QT_HOCL (TEMP, ISO, GI, QT)
      CASE (22)
         CALL QT_N2 (TEMP, ISO, GI, QT)
      CASE (23)
         CALL QT_HCN (TEMP, ISO, GI, QT)
      CASE (24)
         CALL QT_CH3CL (TEMP, ISO, GI, QT)
      CASE (25)
         CALL QT_H2O2 (TEMP, ISO, GI, QT)
      CASE (26)
         CALL QT_C2H2 (TEMP, ISO, GI, QT)
      CASE (27)
         CALL QT_C2H6 (TEMP, ISO, GI, QT)
      CASE (28)
         CALL QT_PH3 (TEMP, ISO, GI, QT)
      CASE (29)
         CALL QT_COF2 (TEMP, ISO, GI, QT)
      CASE (30)
         CALL QT_SF6 (TEMP, ISO, GI, QT)
      CASE (31)
         CALL QT_H2S (TEMP, ISO, GI, QT)
      CASE (32)
         CALL QT_HCOOH (TEMP, ISO, GI, QT)
      CASE (33)
         CALL QT_HO2 (TEMP, ISO, GI, QT)
      CASE (34)
!        ...not applicable to O
         GI = 0
         QT = 0.
      CASE (35)
         CALL QT_CLONO2 (TEMP, ISO, GI, QT)
      CASE (36)
         CALL QT_NOP (TEMP, ISO, GI, QT)
      CASE (37)
         CALL QT_HOBR (TEMP, ISO, GI, QT)
      CASE (38)
         CALL QT_C2H4 (TEMP, ISO, GI, QT)
      CASE DEFAULT
!        ... return 0s if any other molecule
         GI = 0
         QT = 0.
      END SELECT
      RETURN
      END SUBROUTINE BD_TIPS_2008


!***********************
!
      SUBROUTINE BD_TIPS_2008_SFIT2(MOL, TEMP, ISO, GI, QT)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!***********************
!    Subroutine BD_TIPS_2008 written by R.R. Gamache
!
!
! 2008 TIPS from HITRAN 08 merged into sfit2 394 version preserving molecule
! mapping Nov 2009
!
! using HITRAN 08
!
!...date last changed 15 May, 2008
!  --  UPDATES  --
!     5 May, 2008 better high temperature CH4 partition sums, see the paper
!     Ch. Wenger, J.P. Champion, V. Boudon, "The partition sum of methane
!     at high temperature," J Quant Spectros Radiat Transfer 109,
!     Issue 16, November 2008, Pages 2697-2706
!
!
!  sfit2 version uses sfit2 mapping of molecules instead of standard
!  HITRAN mapping  4/27/2005
!
!
!
!...date last changed 25        March, 2005
!  --  UPDATES  --
!    18 December 2003 better vibrational fundamentals for PH3
!
!    This program calculates the total internal
!    partition sum (TIPS) for a given molecule, isotopic species, and
!    temperature.  Current limitations are the molecular species on the
!    HITRAN molecular database and the temperature range 70 - 3000 K.
!
!...This program calculates the TIPS by 4-point LaaGrange interpolation
!
!..  JQSRT - HITRAN Special Issue, 2003.
!..  J. Fischer(a) R.R. Gamache(a&), A. Goldman(b),L.S. Rothman(c), and
!..
!..  (a)  Department of Environmental, Earth, and Atmospheric Sciences,
!..       University of Massachusetts Lowell, Lowell, MA 01854, U.S.A.
!..
!..  (b)  Department of Physics, University of Denver, Denver, CO 80208, U.S.A.
!..
!..  (c)  Harvard-Smithsonian Center for Astrophysics, 60 Garden St, Cambridge, MA 02138 USA
!..
!..  (d)  Laboratoire de Photophysique Mol�culaire, Universit� Paris Sud, 91405 Orsay, FRANCE
!..
!..  &  Corresponding author. Email address: Robert_Gamache@uml.edu
!..  Abstract
!..   Total internal partition sums (TIPS) are calculated for all molecular species in
!..  the 2000 HITRAN database.  In addition, the TIPS for 13 other isotopomers/isotopologues
!..  of ozone and carbon dioxide are presented.  The calculations address the corrections
!..  suggested by Goldman et al. (JQSRT 2000;66:55-86).  The calculations consider the
!..  temperature range 70-3000 K to be applicable to a variety of remote sensing needs.
!..  The method of calculation for each molecular species is stated and
!..  data from the literature are discussed.  A new method of recall for the partition sums,
!..  Lagrange 4-point interpolation, is developed.  This method, unlike
!..  the TIPS code, allows all molecular species to be considered.
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  15:55:42   3/17/05
!...Switches: -yf12
!++
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: MOL ! HITRAN molecule number
      INTEGER  :: ISO             ! isotopomer index
      REAL(DOUBLE)  :: TEMP       ! temperature in K
      REAL(DOUBLE)  :: GI         ! state independent degeneracy factor
      REAL(DOUBLE)  :: QT         ! total internal partition sum
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      !CHARACTER :: STOPNGO*30
!-----------------------------------------------
!
!
!     GO TO PARTICULAR MOLECULE:

      SELECT CASE (MOL)
      CASE (1)
         CALL QT_H2O (TEMP, ISO, GI, QT)
      CASE (2)
         CALL QT_CO2 (TEMP, ISO, GI, QT)
      CASE (3)
         CALL QT_O3 (TEMP, ISO, GI, QT)
      CASE (4)
         CALL QT_N2O (TEMP, ISO, GI, QT)
      CASE (5)
         CALL QT_CO (TEMP, ISO, GI, QT)
      CASE (6)
         CALL QT_CH4 (TEMP, ISO, GI, QT)
      CASE (7)
         CALL QT_O2 (TEMP, ISO, GI, QT)
      CASE (8)
         CALL QT_NO (TEMP, ISO, GI, QT)
      CASE (9)
         CALL QT_SO2 (TEMP, ISO, GI, QT)
      CASE (10)
         CALL QT_NO2 (TEMP, ISO, GI, QT)
      CASE (11)
         CALL QT_NH3 (TEMP, ISO, GI, QT)
      CASE (12)
         CALL QT_HNO3 (TEMP, ISO, GI, QT)
      CASE (13)
         CALL QT_OH (TEMP, ISO, GI, QT)
      CASE (14)
         CALL QT_HF (TEMP, ISO, GI, QT)
      CASE (15)
         CALL QT_HCL (TEMP, ISO, GI, QT)
      CASE (16)
         CALL QT_HBR (TEMP, ISO, GI, QT)
      CASE (17)
         CALL QT_HI (TEMP, ISO, GI, QT)
      CASE (18)
         CALL QT_CLO (TEMP, ISO, GI, QT)
      CASE (19)
         CALL QT_OCS (TEMP, ISO, GI, QT)
      CASE (20)
         CALL QT_H2CO (TEMP, ISO, GI, QT)
      CASE (21)
         CALL QT_HOCL (TEMP, ISO, GI, QT)
      CASE (22)
         CALL QT_HO2 (TEMP, ISO, GI, QT)
      CASE (23)
         CALL QT_H2O2 (TEMP, ISO, GI, QT)
      CASE (28)
         CALL QT_HCN (TEMP, ISO, GI, QT)
      CASE (30)
         CALL QT_CH3CL (TEMP, ISO, GI, QT)
      CASE (36)
         CALL QT_COF2 (TEMP, ISO, GI, QT)
      CASE (38)
         CALL QT_C2H6 (TEMP, ISO, GI, QT)
      CASE (39)
         CALL QT_C2H4 (TEMP, ISO, GI, QT)
      CASE (40)
         CALL QT_C2H2 (TEMP, ISO, GI, QT)
      CASE (41)
         CALL QT_N2 (TEMP, ISO, GI, QT)
      CASE (46)
         CALL QT_HCOOH (TEMP, ISO, GI, QT)
      CASE (47)
         CALL QT_H2S (TEMP, ISO, GI, QT)
      CASE (49)
!        ... HDO maps to H2O isotope 4, 5, 6
         CALL QT_H2O (TEMP, ISO+3, GI, QT)
      CASE (53)
!        ... CH3D maps to CH4 isotope 3, 4
         CALL QT_CH4 (TEMP, ISO+2, GI, QT)
      CASE (54)
!        ... O3668 maps to O3 isotope 2
         CALL QT_O3 (TEMP, 2, GI, QT)
      CASE (55)
!        ... O3686 maps to O3 isotope 3
         CALL QT_O3 (TEMP, 3, GI, QT)
      CASE (56)
!        ... O3667 maps to O3 isotope 4
         CALL QT_O3 (TEMP, 4, GI, QT)
      CASE (57)
!        ... O3676 maps to O3 isotope 5
         CALL QT_O3 (TEMP, 5, GI, QT)

!  CLONO2 is cross section
!      CASE (27)
!         CALL QT_CLONO2 (TEMP, ISO, GI, QT)

!  SF6 is cross section
!      CASE (50) THEN
!         CALL QT_SF6 (TEMP, ISO, GI, QT)

!  PH3 not in sfit2
!      CASE (28)
!         CALL QT_PH3 (TEMP, ISO, GI, QT)

!  O not in sfit2
!      CASE (34)
!        ...not applicable to O
!         GI = 0
!         QT = 0.

! NOP not in sfit2
!      CASE (36)
!         CALL QT_NOP (TEMP, ISO, GI, QT)

! HOBR not in sfit2
!      CASE (37)
!         CALL QT_HOBR (TEMP, ISO, GI, QT)

      CASE DEFAULT
!        ... return 0s if any other molecule since they are not
!            handled in tips
         GI = 0
         QT = 0.
      END SELECT
      RETURN
      END SUBROUTINE BD_TIPS_2008_SFIT2



!**************************************
      SUBROUTINE CLEAR
!**************************************
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
!-----------------------------------------------
!
      DO I = 1, 24
         WRITE (*, '(1X)')
      END DO
      RETURN
      END SUBROUTINE CLEAR
!
!     *****************
      SUBROUTINE QT_H2O(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(6) :: XGJ
      REAL(DOUBLE), DIMENSION(6,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1., 1., 6., 6., 6., 36/
!...      H2O
!...        --       161
      DATA (QOFT(1,J),J=1,119)/ 0.16824E+02, 0.27771E+02, 0.40408E+02, &
         0.54549E+02, 0.70054E+02, 0.86817E+02, 0.10475E+03, 0.12380E+03, &
         0.14391E+03, 0.16503E+03, 0.18714E+03, 0.21021E+03, 0.23425E+03, &
         0.25924E+03, 0.28518E+03, 0.31209E+03, 0.33997E+03, 0.36883E+03, &
         0.39870E+03, 0.42959E+03, 0.46152E+03, 0.49452E+03, 0.52860E+03, &
         0.56380E+03, 0.60015E+03, 0.63766E+03, 0.67637E+03, 0.71631E+03, &
         0.75750E+03, 0.79999E+03, 0.84380E+03, 0.88897E+03, 0.93553E+03, &
         0.98353E+03, 0.10330E+04, 0.10840E+04, 0.11365E+04, 0.11906E+04, &
         0.12463E+04, 0.13037E+04, 0.13628E+04, 0.14237E+04, 0.14863E+04, &
         0.15509E+04, 0.16173E+04, 0.16856E+04, 0.17559E+04, 0.18283E+04, &
         0.19028E+04, 0.19793E+04, 0.20581E+04, 0.21391E+04, 0.22224E+04, &
         0.23080E+04, 0.24067E+04, 0.24975E+04, 0.25908E+04, 0.26867E+04, &
         0.27853E+04, 0.28865E+04, 0.29904E+04, 0.30972E+04, 0.32068E+04, &
         0.33194E+04, 0.34349E+04, 0.35535E+04, 0.36752E+04, 0.38001E+04, &
         0.39282E+04, 0.40597E+04, 0.41945E+04, 0.43327E+04, 0.44745E+04, &
         0.46199E+04, 0.47688E+04, 0.49215E+04, 0.50780E+04, 0.52384E+04, &
         0.54027E+04, 0.55710E+04, 0.57434E+04, 0.59200E+04, 0.61008E+04, &
         0.62859E+04, 0.64754E+04, 0.66693E+04, 0.68679E+04, 0.70710E+04, &
         0.72788E+04, 0.74915E+04, 0.77090E+04, 0.79315E+04, 0.81590E+04, &
         0.83917E+04, 0.86296E+04, 0.88728E+04, 0.91214E+04, 0.93755E+04, &
         0.96351E+04, 0.99005E+04, 0.10171E+05, 0.10448E+05, 0.10731E+05, &
         0.11020E+05, 0.11315E+05, 0.11617E+05, 0.11924E+05, 0.12238E+05, &
         0.12559E+05, 0.12886E+05, 0.13220E+05, 0.13561E+05, 0.13909E+05, &
         0.14263E+05, 0.14625E+05, 0.14995E+05, 0.15371E+05, 0.15755E+05, &
         0.16147E+05/
!...        --       181
      DATA (QOFT(2,J),J=1,119)/ 0.15960E+02, 0.26999E+02, 0.39743E+02, &
         0.54003E+02, 0.69639E+02, 0.86543E+02, 0.10463E+03, 0.12384E+03, &
         0.14412E+03, 0.16542E+03, 0.18773E+03, 0.21103E+03, 0.23531E+03, &
         0.26057E+03, 0.28681E+03, 0.31406E+03, 0.34226E+03, 0.37130E+03, &
         0.40135E+03, 0.43243E+03, 0.46456E+03, 0.49777E+03, 0.53206E+03, &
         0.56748E+03, 0.60405E+03, 0.64179E+03, 0.68074E+03, 0.72093E+03, &
         0.76238E+03, 0.80513E+03, 0.84922E+03, 0.89467E+03, 0.94152E+03, &
         0.98982E+03, 0.10396E+04, 0.10909E+04, 0.11437E+04, 0.11982E+04, &
         0.12543E+04, 0.13120E+04, 0.13715E+04, 0.14328E+04, 0.14959E+04, &
         0.15608E+04, 0.16276E+04, 0.16964E+04, 0.17672E+04, 0.18401E+04, &
         0.19151E+04, 0.19922E+04, 0.20715E+04, 0.21531E+04, 0.22370E+04, &
         0.23232E+04, 0.24118E+04, 0.25030E+04, 0.25967E+04, 0.26929E+04, &
         0.27918E+04, 0.28934E+04, 0.29978E+04, 0.31050E+04, 0.32151E+04, &
         0.33281E+04, 0.34441E+04, 0.35632E+04, 0.36854E+04, 0.38108E+04, &
         0.39395E+04, 0.40715E+04, 0.42070E+04, 0.43459E+04, 0.44883E+04, &
         0.46343E+04, 0.47840E+04, 0.49374E+04, 0.50946E+04, 0.52558E+04, &
         0.54209E+04, 0.55900E+04, 0.57632E+04, 0.59407E+04, 0.61224E+04, &
         0.63084E+04, 0.64988E+04, 0.66938E+04, 0.68933E+04, 0.70975E+04, &
         0.73064E+04, 0.75202E+04, 0.77389E+04, 0.79625E+04, 0.81913E+04, &
         0.84252E+04, 0.86644E+04, 0.89089E+04, 0.91588E+04, 0.94143E+04, &
         0.96754E+04, 0.99422E+04, 0.10215E+05, 0.10493E+05, 0.10778E+05, &
         0.11068E+05, 0.11365E+05, 0.11668E+05, 0.11977E+05, 0.12293E+05, &
         0.12616E+05, 0.12945E+05, 0.13281E+05, 0.13624E+05, 0.13973E+05, &
         0.14330E+05, 0.14694E+05, 0.15066E+05, 0.15445E+05, 0.15831E+05, &
         0.16225E+05/
!...        --       171
      DATA (QOFT(3,J),J=1,119)/ 0.95371E+02, 0.16134E+03, 0.23750E+03, &
         0.32273E+03, 0.41617E+03, 0.51722E+03, 0.62540E+03, 0.74036E+03, &
         0.86185E+03, 0.98970E+03, 0.11238E+04, 0.12642E+04, 0.14097E+04, &
         0.15599E+04, 0.17159E+04, 0.18777E+04, 0.20453E+04, 0.22188E+04, &
         0.23983E+04, 0.25840E+04, 0.27760E+04, 0.29743E+04, 0.31792E+04, &
         0.33907E+04, 0.36091E+04, 0.38346E+04, 0.40672E+04, 0.43072E+04, &
         0.45547E+04, 0.48100E+04, 0.50732E+04, 0.53446E+04, 0.56244E+04, &
         0.59128E+04, 0.62100E+04, 0.65162E+04, 0.68317E+04, 0.71567E+04, &
         0.74915E+04, 0.78363E+04, 0.81914E+04, 0.85571E+04, 0.89335E+04, &
         0.93211E+04, 0.97200E+04, 0.10131E+05, 0.10553E+05, 0.10988E+05, &
         0.11435E+05, 0.11895E+05, 0.12368E+05, 0.12855E+05, 0.13356E+05, &
         0.13870E+05, 0.14399E+05, 0.14943E+05, 0.15502E+05, 0.16076E+05, &
         0.16666E+05, 0.17272E+05, 0.17895E+05, 0.18534E+05, 0.19191E+05, &
         0.19865E+05, 0.20557E+05, 0.21267E+05, 0.21996E+05, 0.22744E+05, &
         0.23512E+05, 0.24299E+05, 0.25106E+05, 0.25935E+05, 0.26784E+05, &
         0.27655E+05, 0.28547E+05, 0.29462E+05, 0.30400E+05, 0.31361E+05, &
         0.32345E+05, 0.33353E+05, 0.34386E+05, 0.35444E+05, 0.36527E+05, &
         0.37637E+05, 0.38772E+05, 0.39934E+05, 0.41124E+05, 0.42341E+05, &
         0.43587E+05, 0.44861E+05, 0.46165E+05, 0.47498E+05, 0.48862E+05, &
         0.50256E+05, 0.51682E+05, 0.53139E+05, 0.54629E+05, 0.56152E+05, &
         0.57708E+05, 0.59299E+05, 0.60923E+05, 0.62583E+05, 0.64279E+05, &
         0.66011E+05, 0.67779E+05, 0.69585E+05, 0.71429E+05, 0.73312E+05, &
         0.75234E+05, 0.77195E+05, 0.79197E+05, 0.81240E+05, 0.83325E+05, &
         0.85452E+05, 0.87622E+05, 0.89835E+05, 0.92093E+05, 0.94395E+05, &
         0.96743E+05/
!...        --       162
      DATA (QOFT(4,J),J=1,119)/ 0.75792E+02, 0.12986E+03, 0.19244E+03, &
         0.26253E+03, 0.33942E+03, 0.42259E+03, 0.51161E+03, 0.60619E+03, &
         0.70609E+03, 0.81117E+03, 0.92132E+03, 0.10365E+04, 0.11567E+04, &
         0.12820E+04, 0.14124E+04, 0.15481E+04, 0.16891E+04, 0.18355E+04, &
         0.19876E+04, 0.21455E+04, 0.23092E+04, 0.24791E+04, 0.26551E+04, &
         0.28376E+04, 0.30268E+04, 0.32258E+04, 0.34288E+04, 0.36392E+04, &
         0.38571E+04, 0.40828E+04, 0.43165E+04, 0.45584E+04, 0.48089E+04, &
         0.50681E+04, 0.53363E+04, 0.56139E+04, 0.59009E+04, 0.61979E+04, &
         0.65049E+04, 0.68224E+04, 0.71506E+04, 0.74898E+04, 0.78403E+04, &
         0.82024E+04, 0.85765E+04, 0.89628E+04, 0.93618E+04, 0.97736E+04, &
         0.10199E+05, 0.10637E+05, 0.11090E+05, 0.11557E+05, 0.12039E+05, &
         0.12535E+05, 0.13047E+05, 0.13575E+05, 0.14119E+05, 0.14679E+05, &
         0.15257E+05, 0.15851E+05, 0.16464E+05, 0.17094E+05, 0.17743E+05, &
         0.18411E+05, 0.19098E+05, 0.19805E+05, 0.20532E+05, 0.21280E+05, &
         0.22049E+05, 0.22840E+05, 0.23652E+05, 0.24487E+05, 0.25345E+05, &
         0.26227E+05, 0.27132E+05, 0.28062E+05, 0.29016E+05, 0.29997E+05, &
         0.31002E+05, 0.32035E+05, 0.33094E+05, 0.34180E+05, 0.35295E+05, &
         0.36438E+05, 0.37610E+05, 0.38812E+05, 0.40044E+05, 0.41306E+05, &
         0.42600E+05, 0.43926E+05, 0.45284E+05, 0.46675E+05, 0.48100E+05, &
         0.49559E+05, 0.51053E+05, 0.52583E+05, 0.54148E+05, 0.55750E+05, &
         0.57390E+05, 0.59067E+05, 0.60783E+05, 0.62539E+05, 0.64334E+05, &
         0.66170E+05, 0.68047E+05, 0.69967E+05, 0.71929E+05, 0.73934E+05, &
         0.75983E+05, 0.78078E+05, 0.80217E+05, 0.82403E+05, 0.84636E+05, &
         0.86917E+05, 0.89246E+05, 0.91625E+05, 0.94053E+05, 0.96533E+05, &
         0.99064E+05/
!...        --       182
      DATA (QOFT(5,J),J=1,119)/ 0.82770E+02, 0.13749E+03, 0.20083E+03, &
         0.27176E+03, 0.34955E+03, 0.43370E+03, 0.52376E+03, 0.61944E+03, &
         0.72050E+03, 0.82679E+03, 0.93821E+03, 0.10547E+04, 0.11763E+04, &
         0.13031E+04, 0.14350E+04, 0.15723E+04, 0.17150E+04, 0.18633E+04, &
         0.20172E+04, 0.21770E+04, 0.23429E+04, 0.25149E+04, 0.26934E+04, &
         0.28784E+04, 0.30702E+04, 0.32690E+04, 0.34750E+04, 0.36885E+04, &
         0.39096E+04, 0.41386E+04, 0.43758E+04, 0.46213E+04, 0.48755E+04, &
         0.51386E+04, 0.54109E+04, 0.56927E+04, 0.59841E+04, 0.62856E+04, &
         0.65973E+04, 0.69197E+04, 0.72529E+04, 0.75973E+04, 0.79533E+04, &
         0.83210E+04, 0.87009E+04, 0.90933E+04, 0.94985E+04, 0.99168E+04, &
         0.10348E+05, 0.10794E+05, 0.11254E+05, 0.11728E+05, 0.12217E+05, &
         0.12722E+05, 0.13242E+05, 0.13778E+05, 0.14331E+05, 0.14900E+05, &
         0.15486E+05, 0.16091E+05, 0.16713E+05, 0.17353E+05, 0.18012E+05, &
         0.18691E+05, 0.19389E+05, 0.20108E+05, 0.20847E+05, 0.21607E+05, &
         0.22388E+05, 0.23191E+05, 0.24017E+05, 0.24866E+05, 0.25738E+05, &
         0.26633E+05, 0.27553E+05, 0.28498E+05, 0.29468E+05, 0.30464E+05, &
         0.31486E+05, 0.32536E+05, 0.33612E+05, 0.34716E+05, 0.35849E+05, &
         0.37011E+05, 0.38202E+05, 0.39424E+05, 0.40676E+05, 0.41959E+05, &
         0.43274E+05, 0.44622E+05, 0.46002E+05, 0.47416E+05, 0.48864E+05, &
         0.50348E+05, 0.51866E+05, 0.53421E+05, 0.55012E+05, 0.56640E+05, &
         0.58307E+05, 0.60012E+05, 0.61757E+05, 0.63541E+05, 0.65366E+05, &
         0.67233E+05, 0.69141E+05, 0.71092E+05, 0.73087E+05, 0.75125E+05, &
         0.77209E+05, 0.79338E+05, 0.81513E+05, 0.83736E+05, 0.86006E+05, &
         0.88324E+05, 0.90693E+05, 0.93111E+05, 0.95580E+05, 0.98100E+05, &
         0.10067E+06/
!...        --       172
      DATA (QOFT(6,J),J=1,119)/ 0.49379E+03, 0.82021E+03, 0.11980E+04, &
         0.16211E+04, 0.20851E+04, 0.25870E+04, 0.31242E+04, 0.36949E+04, &
         0.42977E+04, 0.49317E+04, 0.55963E+04, 0.62911E+04, 0.70164E+04, &
         0.77722E+04, 0.85591E+04, 0.93777E+04, 0.10228E+05, 0.11112E+05, &
         0.12030E+05, 0.12983E+05, 0.13971E+05, 0.14997E+05, 0.16061E+05, &
         0.17163E+05, 0.18306E+05, 0.19491E+05, 0.20719E+05, 0.21991E+05, &
         0.23309E+05, 0.24673E+05, 0.26086E+05, 0.27549E+05, 0.29064E+05, &
         0.30631E+05, 0.32254E+05, 0.33932E+05, 0.35669E+05, 0.37464E+05, &
         0.39321E+05, 0.41242E+05, 0.43227E+05, 0.45279E+05, 0.47399E+05, &
         0.49589E+05, 0.51852E+05, 0.54189E+05, 0.56602E+05, 0.59094E+05, &
         0.61666E+05, 0.64320E+05, 0.67058E+05, 0.69883E+05, 0.72796E+05, &
         0.75801E+05, 0.78899E+05, 0.82092E+05, 0.85382E+05, 0.88773E+05, &
         0.92266E+05, 0.95863E+05, 0.99568E+05, 0.10338E+06, 0.10731E+06, &
         0.11135E+06, 0.11551E+06, 0.11979E+06, 0.12419E+06, 0.12871E+06, &
         0.13337E+06, 0.13815E+06, 0.14307E+06, 0.14812E+06, 0.15331E+06, &
         0.15865E+06, 0.16412E+06, 0.16975E+06, 0.17553E+06, 0.18146E+06, &
         0.18754E+06, 0.19379E+06, 0.20020E+06, 0.20678E+06, 0.21352E+06, &
         0.22044E+06, 0.22753E+06, 0.23480E+06, 0.24226E+06, 0.24990E+06, &
         0.25773E+06, 0.26575E+06, 0.27397E+06, 0.28239E+06, 0.29102E+06, &
         0.29985E+06, 0.30889E+06, 0.31814E+06, 0.32762E+06, 0.33731E+06, &
         0.34724E+06, 0.35739E+06, 0.36777E+06, 0.37840E+06, 0.38926E+06, &
         0.40038E+06, 0.41174E+06, 0.42335E+06, 0.43523E+06, 0.44737E+06, &
         0.45977E+06, 0.47245E+06, 0.48540E+06, 0.49863E+06, 0.51214E+06, &
         0.52595E+06, 0.54005E+06, 0.55444E+06, 0.56914E+06, 0.58415E+06, &
         0.59947E+06/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_H2O


!
!     *****************
      SUBROUTINE QT_CO2(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(9) :: XGJ
      REAL(DOUBLE), DIMENSION(9,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1., 2., 1., 6., 2., 12., 1., 6., 1/
!...      CO2
!...        --       626
      DATA (QOFT(1,J),J=1,119)/ 0.53642E+02, 0.75947E+02, 0.98292E+02, &
         0.12078E+03, 0.14364E+03, 0.16714E+03, 0.19160E+03, 0.21731E+03, &
         0.24454E+03, 0.27355E+03, 0.30456E+03, 0.33778E+03, 0.37343E+03, &
         0.41170E+03, 0.45280E+03, 0.49692E+03, 0.54427E+03, 0.59505E+03, &
         0.64948E+03, 0.70779E+03, 0.77019E+03, 0.83693E+03, 0.90825E+03, &
         0.98440E+03, 0.10656E+04, 0.11522E+04, 0.12445E+04, 0.13427E+04, &
         0.14471E+04, 0.15580E+04, 0.16759E+04, 0.18009E+04, 0.19334E+04, &
         0.20739E+04, 0.22225E+04, 0.23798E+04, 0.25462E+04, 0.27219E+04, &
         0.29074E+04, 0.31032E+04, 0.33097E+04, 0.35272E+04, 0.37564E+04, &
         0.39976E+04, 0.42514E+04, 0.45181E+04, 0.47985E+04, 0.50929E+04, &
         0.54019E+04, 0.57260E+04, 0.60659E+04, 0.64221E+04, 0.67952E+04, &
         0.71859E+04, 0.75946E+04, 0.80222E+04, 0.84691E+04, 0.89362E+04, &
         0.94241E+04, 0.99335E+04, 0.10465E+05, 0.11020E+05, 0.11598E+05, &
         0.12201E+05, 0.12828E+05, 0.13482E+05, 0.14163E+05, 0.14872E+05, &
         0.15609E+05, 0.16376E+05, 0.17173E+05, 0.18001E+05, 0.18861E+05, &
         0.19754E+05, 0.20682E+05, 0.21644E+05, 0.22643E+05, 0.23678E+05, &
         0.24752E+05, 0.25865E+05, 0.27018E+05, 0.28212E+05, 0.29449E+05, &
         0.30730E+05, 0.32055E+05, 0.33426E+05, 0.34845E+05, 0.36312E+05, &
         0.37828E+05, 0.39395E+05, 0.41015E+05, 0.42688E+05, 0.44416E+05, &
         0.46199E+05, 0.48041E+05, 0.49942E+05, 0.51902E+05, 0.53925E+05, &
         0.56011E+05, 0.58162E+05, 0.60379E+05, 0.62664E+05, 0.65019E+05, &
         0.67444E+05, 0.69942E+05, 0.72515E+05, 0.75163E+05, 0.77890E+05, &
         0.80695E+05, 0.83582E+05, 0.86551E+05, 0.89605E+05, 0.92746E+05, &
         0.95975E+05, 0.99294E+05, 0.10271E+06, 0.10621E+06, 0.10981E+06, &
         0.11351E+06/
!...        --       636
      DATA (QOFT(2,J),J=1,119)/ 0.10728E+03, 0.15189E+03, 0.19659E+03, &
         0.24164E+03, 0.28753E+03, 0.33486E+03, 0.38429E+03, 0.43643E+03, &
         0.49184E+03, 0.55104E+03, 0.61449E+03, 0.68263E+03, 0.75589E+03, &
         0.83468E+03, 0.91943E+03, 0.10106E+04, 0.11085E+04, 0.12137E+04, &
         0.13266E+04, 0.14477E+04, 0.15774E+04, 0.17163E+04, 0.18649E+04, &
         0.20237E+04, 0.21933E+04, 0.23743E+04, 0.25673E+04, 0.27729E+04, &
         0.29917E+04, 0.32245E+04, 0.34718E+04, 0.37345E+04, 0.40132E+04, &
         0.43087E+04, 0.46218E+04, 0.49533E+04, 0.53041E+04, 0.56749E+04, &
         0.60668E+04, 0.64805E+04, 0.69171E+04, 0.73774E+04, 0.78626E+04, &
         0.83736E+04, 0.89114E+04, 0.94772E+04, 0.10072E+05, 0.10697E+05, &
         0.11353E+05, 0.12042E+05, 0.12765E+05, 0.13523E+05, 0.14317E+05, &
         0.15148E+05, 0.16019E+05, 0.16930E+05, 0.17883E+05, 0.18879E+05, &
         0.19920E+05, 0.21008E+05, 0.22143E+05, 0.23328E+05, 0.24563E+05, &
         0.25852E+05, 0.27195E+05, 0.28594E+05, 0.30051E+05, 0.31568E+05, &
         0.33146E+05, 0.34788E+05, 0.36496E+05, 0.38271E+05, 0.40115E+05, &
         0.42031E+05, 0.44021E+05, 0.46086E+05, 0.48230E+05, 0.50453E+05, &
         0.52759E+05, 0.55150E+05, 0.57628E+05, 0.60195E+05, 0.62854E+05, &
         0.65608E+05, 0.68459E+05, 0.71409E+05, 0.74461E+05, 0.77618E+05, &
         0.80883E+05, 0.84258E+05, 0.87746E+05, 0.91350E+05, 0.95073E+05, &
         0.98918E+05, 0.10289E+06, 0.10698E+06, 0.11121E+06, 0.11558E+06, &
         0.12008E+06, 0.12472E+06, 0.12950E+06, 0.13443E+06, 0.13952E+06, &
         0.14475E+06, 0.15015E+06, 0.15571E+06, 0.16143E+06, 0.16732E+06, &
         0.17338E+06, 0.17962E+06, 0.18604E+06, 0.19264E+06, 0.19943E+06, &
         0.20642E+06, 0.21360E+06, 0.22098E+06, 0.22856E+06, 0.23636E+06, &
         0.24436E+06/
!...        --       628
      DATA (QOFT(3,J),J=1,119)/ 0.11368E+03, 0.16096E+03, 0.20833E+03, &
         0.25603E+03, 0.30452E+03, 0.35442E+03, 0.40640E+03, 0.46110E+03, &
         0.51910E+03, 0.58093E+03, 0.64709E+03, 0.71804E+03, 0.79422E+03, &
         0.87607E+03, 0.96402E+03, 0.10585E+04, 0.11600E+04, 0.12689E+04, &
         0.13857E+04, 0.15108E+04, 0.16449E+04, 0.17883E+04, 0.19416E+04, &
         0.21054E+04, 0.22803E+04, 0.24668E+04, 0.26655E+04, 0.28770E+04, &
         0.31021E+04, 0.33414E+04, 0.35956E+04, 0.38654E+04, 0.41516E+04, &
         0.44549E+04, 0.47761E+04, 0.51160E+04, 0.54755E+04, 0.58555E+04, &
         0.62568E+04, 0.66804E+04, 0.71273E+04, 0.75982E+04, 0.80944E+04, &
         0.86169E+04, 0.91666E+04, 0.97446E+04, 0.10352E+05, 0.10990E+05, &
         0.11660E+05, 0.12363E+05, 0.13101E+05, 0.13874E+05, 0.14683E+05, &
         0.15531E+05, 0.16418E+05, 0.17347E+05, 0.18317E+05, 0.19332E+05, &
         0.20392E+05, 0.21499E+05, 0.22654E+05, 0.23859E+05, 0.25116E+05, &
         0.26426E+05, 0.27792E+05, 0.29214E+05, 0.30695E+05, 0.32236E+05, &
         0.33840E+05, 0.35508E+05, 0.37242E+05, 0.39045E+05, 0.40917E+05, &
         0.42862E+05, 0.44881E+05, 0.46977E+05, 0.49152E+05, 0.51407E+05, &
         0.53746E+05, 0.56171E+05, 0.58683E+05, 0.61286E+05, 0.63981E+05, &
         0.66772E+05, 0.69661E+05, 0.72650E+05, 0.75742E+05, 0.78940E+05, &
         0.82246E+05, 0.85664E+05, 0.89196E+05, 0.92845E+05, 0.96613E+05, &
         0.10050E+06, 0.10452E+06, 0.10867E+06, 0.11295E+06, 0.11736E+06, &
         0.12191E+06, 0.12661E+06, 0.13145E+06, 0.13643E+06, 0.14157E+06, &
         0.14687E+06, 0.15232E+06, 0.15794E+06, 0.16372E+06, 0.16968E+06, &
         0.17580E+06, 0.18211E+06, 0.18859E+06, 0.19526E+06, 0.20213E+06, &
         0.20918E+06, 0.21643E+06, 0.22388E+06, 0.23154E+06, 0.23941E+06, &
         0.24750E+06/
!...        --       627
      DATA (QOFT(4,J),J=1,119)/ 0.66338E+03, 0.93923E+03, 0.12156E+04, &
         0.14938E+04, 0.17766E+04, 0.20676E+04, 0.23705E+04, 0.26891E+04, &
         0.30267E+04, 0.33866E+04, 0.37714E+04, 0.41839E+04, 0.46267E+04, &
         0.51023E+04, 0.56132E+04, 0.61618E+04, 0.67508E+04, 0.73827E+04, &
         0.80603E+04, 0.87863E+04, 0.95636E+04, 0.10395E+05, 0.11284E+05, &
         0.12233E+05, 0.13246E+05, 0.14326E+05, 0.15477E+05, 0.16702E+05, &
         0.18005E+05, 0.19390E+05, 0.20861E+05, 0.22422E+05, 0.24077E+05, &
         0.25832E+05, 0.27689E+05, 0.29655E+05, 0.31734E+05, 0.33931E+05, &
         0.36250E+05, 0.38698E+05, 0.41280E+05, 0.44002E+05, 0.46869E+05, &
         0.49886E+05, 0.53062E+05, 0.56400E+05, 0.59909E+05, 0.63594E+05, &
         0.67462E+05, 0.71521E+05, 0.75777E+05, 0.80238E+05, 0.84911E+05, &
         0.89804E+05, 0.94925E+05, 0.10028E+06, 0.10588E+06, 0.11173E+06, &
         0.11785E+06, 0.12423E+06, 0.13090E+06, 0.13785E+06, 0.14510E+06, &
         0.15265E+06, 0.16053E+06, 0.16873E+06, 0.17727E+06, 0.18615E+06, &
         0.19540E+06, 0.20501E+06, 0.21501E+06, 0.22540E+06, 0.23619E+06, &
         0.24740E+06, 0.25904E+06, 0.27112E+06, 0.28365E+06, 0.29664E+06, &
         0.31012E+06, 0.32409E+06, 0.33856E+06, 0.35356E+06, 0.36908E+06, &
         0.38516E+06, 0.40180E+06, 0.41902E+06, 0.43683E+06, 0.45525E+06, &
         0.47429E+06, 0.49397E+06, 0.51431E+06, 0.53532E+06, 0.55702E+06, &
         0.57943E+06, 0.60256E+06, 0.62644E+06, 0.65107E+06, 0.67648E+06, &
         0.70269E+06, 0.72972E+06, 0.75758E+06, 0.78629E+06, 0.81588E+06, &
         0.84636E+06, 0.87775E+06, 0.91008E+06, 0.94337E+06, 0.97763E+06, &
         0.10129E+07, 0.10492E+07, 0.10865E+07, 0.11249E+07, 0.11644E+07, &
         0.12050E+07, 0.12467E+07, 0.12896E+07, 0.13337E+07, 0.13789E+07, &
         0.14255E+07/
!...        --       638
      DATA (QOFT(5,J),J=1,119)/ 0.22737E+03, 0.32194E+03, 0.41671E+03, &
         0.51226E+03, 0.60963E+03, 0.71017E+03, 0.81528E+03, 0.92628E+03, &
         0.10444E+04, 0.11707E+04, 0.13061E+04, 0.14518E+04, 0.16085E+04, &
         0.17772E+04, 0.19588E+04, 0.21542E+04, 0.23644E+04, 0.25903E+04, &
         0.28330E+04, 0.30934E+04, 0.33726E+04, 0.36717E+04, 0.39918E+04, &
         0.43342E+04, 0.47001E+04, 0.50907E+04, 0.55074E+04, 0.59515E+04, &
         0.64244E+04, 0.69276E+04, 0.74626E+04, 0.80310E+04, 0.86344E+04, &
         0.92744E+04, 0.99528E+04, 0.10671E+05, 0.11432E+05, 0.12236E+05, &
         0.13086E+05, 0.13984E+05, 0.14932E+05, 0.15932E+05, 0.16985E+05, &
         0.18096E+05, 0.19265E+05, 0.20495E+05, 0.21788E+05, 0.23148E+05, &
         0.24576E+05, 0.26075E+05, 0.27648E+05, 0.29298E+05, 0.31027E+05, &
         0.32839E+05, 0.34736E+05, 0.36721E+05, 0.38798E+05, 0.40970E+05, &
         0.43240E+05, 0.45611E+05, 0.48087E+05, 0.50671E+05, 0.53368E+05, &
         0.56180E+05, 0.59111E+05, 0.62165E+05, 0.65347E+05, 0.68659E+05, &
         0.72107E+05, 0.75694E+05, 0.79425E+05, 0.83303E+05, 0.87334E+05, &
         0.91522E+05, 0.95872E+05, 0.10039E+06, 0.10507E+06, 0.10994E+06, &
         0.11498E+06, 0.12021E+06, 0.12563E+06, 0.13125E+06, 0.13707E+06, &
         0.14309E+06, 0.14933E+06, 0.15579E+06, 0.16247E+06, 0.16938E+06, &
         0.17653E+06, 0.18392E+06, 0.19156E+06, 0.19946E+06, 0.20761E+06, &
         0.21604E+06, 0.22473E+06, 0.23371E+06, 0.24298E+06, 0.25254E+06, &
         0.26240E+06, 0.27258E+06, 0.28307E+06, 0.29388E+06, 0.30502E+06, &
         0.31651E+06, 0.32834E+06, 0.34052E+06, 0.35307E+06, 0.36599E+06, &
         0.37929E+06, 0.39298E+06, 0.40706E+06, 0.42155E+06, 0.43645E+06, &
         0.45178E+06, 0.46753E+06, 0.48373E+06, 0.50038E+06, 0.51748E+06, &
         0.53506E+06/
!...        --       637
      DATA (QOFT(6,J),J=1,119)/ 0.13267E+04, 0.18785E+04, 0.24314E+04, &
         0.29888E+04, 0.35566E+04, 0.41426E+04, 0.47550E+04, 0.54013E+04, &
         0.60886E+04, 0.68232E+04, 0.76109E+04, 0.84574E+04, 0.93678E+04, &
         0.10348E+05, 0.11402E+05, 0.12536E+05, 0.13755E+05, 0.15065E+05, &
         0.16471E+05, 0.17980E+05, 0.19598E+05, 0.21330E+05, 0.23184E+05, &
         0.25166E+05, 0.27283E+05, 0.29543E+05, 0.31953E+05, 0.34521E+05, &
         0.37256E+05, 0.40164E+05, 0.43256E+05, 0.46541E+05, 0.50026E+05, &
         0.53723E+05, 0.57641E+05, 0.61790E+05, 0.66180E+05, 0.70823E+05, &
         0.75729E+05, 0.80910E+05, 0.86378E+05, 0.92145E+05, 0.98224E+05, &
         0.10463E+06, 0.11137E+06, 0.11846E+06, 0.12592E+06, 0.13375E+06, &
         0.14198E+06, 0.15062E+06, 0.15969E+06, 0.16920E+06, 0.17916E+06, &
         0.18959E+06, 0.20052E+06, 0.21196E+06, 0.22392E+06, 0.23642E+06, &
         0.24949E+06, 0.26314E+06, 0.27740E+06, 0.29227E+06, 0.30779E+06, &
         0.32398E+06, 0.34085E+06, 0.35842E+06, 0.37673E+06, 0.39579E+06, &
         0.41563E+06, 0.43626E+06, 0.45772E+06, 0.48003E+06, 0.50322E+06, &
         0.52730E+06, 0.55232E+06, 0.57829E+06, 0.60524E+06, 0.63320E+06, &
         0.66219E+06, 0.69226E+06, 0.72342E+06, 0.75571E+06, 0.78916E+06, &
         0.82380E+06, 0.85966E+06, 0.89678E+06, 0.93518E+06, 0.97490E+06, &
         0.10160E+07, 0.10585E+07, 0.11023E+07, 0.11477E+07, 0.11946E+07, &
         0.12430E+07, 0.12929E+07, 0.13445E+07, 0.13977E+07, 0.14526E+07, &
         0.15093E+07, 0.15677E+07, 0.16280E+07, 0.16901E+07, 0.17541E+07, &
         0.18200E+07, 0.18880E+07, 0.19579E+07, 0.20300E+07, 0.21042E+07, &
         0.21805E+07, 0.22591E+07, 0.23400E+07, 0.24232E+07, 0.25087E+07, &
         0.25967E+07, 0.26871E+07, 0.27801E+07, 0.28757E+07, 0.29739E+07, &
         0.30747E+07/
!...        --       828
      DATA (QOFT(7,J),J=1,119)/ 0.60334E+02, 0.85430E+02, 0.11058E+03, &
         0.13590E+03, 0.16167E+03, 0.18821E+03, 0.21588E+03, 0.24502E+03, &
         0.27595E+03, 0.30896E+03, 0.34431E+03, 0.38225E+03, 0.42301E+03, &
         0.46684E+03, 0.51397E+03, 0.56464E+03, 0.61907E+03, 0.67753E+03, &
         0.74027E+03, 0.80753E+03, 0.87961E+03, 0.95676E+03, 0.10393E+04, &
         0.11275E+04, 0.12217E+04, 0.13222E+04, 0.14293E+04, 0.15434E+04, &
         0.16648E+04, 0.17940E+04, 0.19312E+04, 0.20769E+04, 0.22315E+04, &
         0.23954E+04, 0.25691E+04, 0.27529E+04, 0.29474E+04, 0.31530E+04, &
         0.33702E+04, 0.35995E+04, 0.38414E+04, 0.40965E+04, 0.43654E+04, &
         0.46484E+04, 0.49464E+04, 0.52598E+04, 0.55892E+04, 0.59353E+04, &
         0.62988E+04, 0.66803E+04, 0.70804E+04, 0.74998E+04, 0.79394E+04, &
         0.83998E+04, 0.88817E+04, 0.93859E+04, 0.99132E+04, 0.10464E+05, &
         0.11040E+05, 0.11642E+05, 0.12270E+05, 0.12925E+05, 0.13609E+05, &
         0.14321E+05, 0.15064E+05, 0.15838E+05, 0.16643E+05, 0.17482E+05, &
         0.18355E+05, 0.19263E+05, 0.20207E+05, 0.21188E+05, 0.22208E+05, &
         0.23267E+05, 0.24366E+05, 0.25508E+05, 0.26692E+05, 0.27921E+05, &
         0.29195E+05, 0.30516E+05, 0.31886E+05, 0.33304E+05, 0.34773E+05, &
         0.36294E+05, 0.37869E+05, 0.39499E+05, 0.41185E+05, 0.42929E+05, &
         0.44732E+05, 0.46596E+05, 0.48522E+05, 0.50513E+05, 0.52569E+05, &
         0.54692E+05, 0.56884E+05, 0.59146E+05, 0.61481E+05, 0.63890E+05, &
         0.66375E+05, 0.68937E+05, 0.71578E+05, 0.74301E+05, 0.77107E+05, &
         0.79998E+05, 0.82976E+05, 0.86043E+05, 0.89201E+05, 0.92452E+05, &
         0.95799E+05, 0.99242E+05, 0.10278E+06, 0.10643E+06, 0.11018E+06, &
         0.11403E+06, 0.11799E+06, 0.12206E+06, 0.12625E+06, 0.13055E+06, &
         0.13497E+06/
!...        --       728
      DATA (QOFT(8,J),J=1,119)/ 0.70354E+03, 0.99615E+03, 0.12893E+04, &
         0.15846E+04, 0.18848E+04, 0.21940E+04, 0.25162E+04, 0.28554E+04, &
         0.32152E+04, 0.35991E+04, 0.40099E+04, 0.44507E+04, 0.49242E+04, &
         0.54332E+04, 0.59802E+04, 0.65681E+04, 0.71996E+04, 0.78776E+04, &
         0.86050E+04, 0.93847E+04, 0.10220E+05, 0.11114E+05, 0.12070E+05, &
         0.13091E+05, 0.14182E+05, 0.15345E+05, 0.16585E+05, 0.17906E+05, &
         0.19311E+05, 0.20805E+05, 0.22393E+05, 0.24078E+05, 0.25865E+05, &
         0.27760E+05, 0.29768E+05, 0.31893E+05, 0.34140E+05, 0.36516E+05, &
         0.39025E+05, 0.41674E+05, 0.44469E+05, 0.47416E+05, 0.50520E+05, &
         0.53789E+05, 0.57229E+05, 0.60847E+05, 0.64650E+05, 0.68645E+05, &
         0.72840E+05, 0.77242E+05, 0.81859E+05, 0.86699E+05, 0.91770E+05, &
         0.97081E+05, 0.10264E+06, 0.10846E+06, 0.11454E+06, 0.12090E+06, &
         0.12754E+06, 0.13447E+06, 0.14171E+06, 0.14927E+06, 0.15715E+06, &
         0.16536E+06, 0.17392E+06, 0.18284E+06, 0.19213E+06, 0.20179E+06, &
         0.21185E+06, 0.22231E+06, 0.23319E+06, 0.24450E+06, 0.25625E+06, &
         0.26845E+06, 0.28112E+06, 0.29427E+06, 0.30791E+06, 0.32206E+06, &
         0.33674E+06, 0.35196E+06, 0.36772E+06, 0.38406E+06, 0.40098E+06, &
         0.41850E+06, 0.43663E+06, 0.45539E+06, 0.47480E+06, 0.49488E+06, &
         0.51564E+06, 0.53710E+06, 0.55928E+06, 0.58219E+06, 0.60586E+06, &
         0.63029E+06, 0.65553E+06, 0.68157E+06, 0.70844E+06, 0.73616E+06, &
         0.76476E+06, 0.79424E+06, 0.82464E+06, 0.85597E+06, 0.88826E+06, &
         0.92153E+06, 0.95580E+06, 0.99108E+06, 0.10274E+07, 0.10648E+07, &
         0.11033E+07, 0.11429E+07, 0.11837E+07, 0.12256E+07, 0.12687E+07, &
         0.13131E+07, 0.13586E+07, 0.14055E+07, 0.14536E+07, 0.15031E+07, &
         0.15539E+07/
!...        --       727
      DATA (QOFT(9,J),J=1,119)/ 0.20518E+04, 0.29051E+04, 0.37601E+04, &
         0.46209E+04, 0.54961E+04, 0.63969E+04, 0.73353E+04, 0.83227E+04, &
         0.93698E+04, 0.10486E+05, 0.11681E+05, 0.12962E+05, 0.14337E+05, &
         0.15815E+05, 0.17403E+05, 0.19110E+05, 0.20942E+05, 0.22909E+05, &
         0.25018E+05, 0.27278E+05, 0.29699E+05, 0.32290E+05, 0.35060E+05, &
         0.38019E+05, 0.41177E+05, 0.44545E+05, 0.48135E+05, 0.51957E+05, &
         0.56023E+05, 0.60346E+05, 0.64938E+05, 0.69812E+05, 0.74981E+05, &
         0.80461E+05, 0.86264E+05, 0.92406E+05, 0.98902E+05, 0.10577E+06, &
         0.11302E+06, 0.12067E+06, 0.12875E+06, 0.13726E+06, 0.14622E+06, &
         0.15566E+06, 0.16559E+06, 0.17604E+06, 0.18702E+06, 0.19855E+06, &
         0.21066E+06, 0.22336E+06, 0.23669E+06, 0.25065E+06, 0.26528E+06, &
         0.28061E+06, 0.29664E+06, 0.31342E+06, 0.33096E+06, 0.34930E+06, &
         0.36845E+06, 0.38845E+06, 0.40933E+06, 0.43111E+06, 0.45383E+06, &
         0.47751E+06, 0.50219E+06, 0.52790E+06, 0.55466E+06, 0.58252E+06, &
         0.61151E+06, 0.64166E+06, 0.67300E+06, 0.70558E+06, 0.73943E+06, &
         0.77458E+06, 0.81108E+06, 0.84896E+06, 0.88827E+06, 0.92904E+06, &
         0.97131E+06, 0.10151E+07, 0.10605E+07, 0.11076E+07, 0.11563E+07, &
         0.12068E+07, 0.12590E+07, 0.13130E+07, 0.13689E+07, 0.14267E+07, &
         0.14865E+07, 0.15483E+07, 0.16121E+07, 0.16781E+07, 0.17462E+07, &
         0.18165E+07, 0.18892E+07, 0.19641E+07, 0.20415E+07, 0.21213E+07, &
         0.22036E+07, 0.22884E+07, 0.23759E+07, 0.24661E+07, 0.25590E+07, &
         0.26547E+07, 0.27533E+07, 0.28549E+07, 0.29594E+07, 0.30670E+07, &
         0.31778E+07, 0.32918E+07, 0.34090E+07, 0.35296E+07, 0.36536E+07, &
         0.37812E+07, 0.39123E+07, 0.40470E+07, 0.41855E+07, 0.43278E+07, &
         0.44739E+07/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_CO2


!
!     *****************
      SUBROUTINE QT_O3(T, ISO, GSI, QT)          ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(18) :: XGJ
      REAL(DOUBLE), DIMENSION(18,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1., 1., 1., 6., 6., 1., 1., 6., 6., 6., 36., 1., 1., 6., 6., &
         36., 1., 6./
!...       O3
!...        --       666
      DATA (QOFT(1,J),J=1,119)/ 0.30333E+03, 0.51126E+03, 0.75274E+03, &
         0.10241E+04, 0.13236E+04, 0.16508E+04, 0.20068E+04, 0.23935E+04, &
         0.28136E+04, 0.32703E+04, 0.37672E+04, 0.43082E+04, 0.48975E+04, &
         0.55395E+04, 0.62386E+04, 0.69996E+04, 0.78272E+04, 0.87264E+04, &
         0.97026E+04, 0.10761E+05, 0.11907E+05, 0.13146E+05, 0.14485E+05, &
         0.15929E+05, 0.17484E+05, 0.19158E+05, 0.20957E+05, 0.22887E+05, &
         0.24956E+05, 0.27172E+05, 0.29541E+05, 0.32072E+05, 0.34773E+05, &
         0.37652E+05, 0.40718E+05, 0.43979E+05, 0.47444E+05, 0.51123E+05, &
         0.55026E+05, 0.59161E+05, 0.63540E+05, 0.68172E+05, 0.73069E+05, &
         0.78240E+05, 0.83698E+05, 0.89453E+05, 0.95517E+05, 0.10190E+06, &
         0.10862E+06, 0.11569E+06, 0.12311E+06, 0.13091E+06, 0.13909E+06, &
         0.14767E+06, 0.15666E+06, 0.16608E+06, 0.17594E+06, 0.18626E+06, &
         0.19706E+06, 0.20834E+06, 0.22012E+06, 0.23242E+06, 0.24526E+06, &
         0.25866E+06, 0.27262E+06, 0.28717E+06, 0.30233E+06, 0.31811E+06, &
         0.33453E+06, 0.35161E+06, 0.36937E+06, 0.38784E+06, 0.40702E+06, &
         0.42694E+06, 0.44762E+06, 0.46909E+06, 0.49135E+06, 0.51444E+06, &
         0.53838E+06, 0.56318E+06, 0.58887E+06, 0.61548E+06, 0.64303E+06, &
         0.67153E+06, 0.70102E+06, 0.73153E+06, 0.76306E+06, 0.79566E+06, &
         0.82934E+06, 0.86413E+06, 0.90006E+06, 0.93716E+06, 0.97545E+06, &
         0.10150E+07, 0.10557E+07, 0.10977E+07, 0.11411E+07, 0.11858E+07, &
         0.12318E+07, 0.12792E+07, 0.13281E+07, 0.13784E+07, 0.14302E+07, &
         0.14835E+07, 0.15384E+07, 0.15948E+07, 0.16529E+07, 0.17126E+07, &
         0.17740E+07, 0.18371E+07, 0.19020E+07, 0.19686E+07, 0.20371E+07, &
         0.21074E+07, 0.21797E+07, 0.22538E+07, 0.23300E+07, 0.24081E+07, &
         0.24883E+07/
!...        --       668
      DATA (QOFT(2,J),J=1,119)/ 0.64763E+03, 0.10916E+04, 0.16073E+04, &
         0.21870E+04, 0.28271E+04, 0.35272E+04, 0.42900E+04, 0.51197E+04, &
         0.60225E+04, 0.70057E+04, 0.80771E+04, 0.92455E+04, 0.10520E+05, &
         0.11911E+05, 0.13427E+05, 0.15079E+05, 0.16878E+05, 0.18834E+05, &
         0.20960E+05, 0.23267E+05, 0.25767E+05, 0.28472E+05, 0.31397E+05, &
         0.34553E+05, 0.37957E+05, 0.41620E+05, 0.45559E+05, 0.49790E+05, &
         0.54327E+05, 0.59187E+05, 0.64387E+05, 0.69944E+05, 0.75877E+05, &
         0.82203E+05, 0.88943E+05, 0.96114E+05, 0.10374E+06, 0.11184E+06, &
         0.12043E+06, 0.12954E+06, 0.13918E+06, 0.14939E+06, 0.16018E+06, &
         0.17159E+06, 0.18362E+06, 0.19632E+06, 0.20970E+06, 0.22380E+06, &
         0.23863E+06, 0.25423E+06, 0.27063E+06, 0.28786E+06, 0.30594E+06, &
         0.32490E+06, 0.34478E+06, 0.36561E+06, 0.38743E+06, 0.41026E+06, &
         0.43413E+06, 0.45909E+06, 0.48517E+06, 0.51241E+06, 0.54084E+06, &
         0.57049E+06, 0.60141E+06, 0.63365E+06, 0.66722E+06, 0.70219E+06, &
         0.73858E+06, 0.77644E+06, 0.81581E+06, 0.85674E+06, 0.89927E+06, &
         0.94345E+06, 0.98932E+06, 0.10369E+07, 0.10863E+07, 0.11375E+07, &
         0.11906E+07, 0.12457E+07, 0.13027E+07, 0.13618E+07, 0.14229E+07, &
         0.14862E+07, 0.15517E+07, 0.16194E+07, 0.16894E+07, 0.17618E+07, &
         0.18366E+07, 0.19139E+07, 0.19937E+07, 0.20761E+07, 0.21612E+07, &
         0.22490E+07, 0.23395E+07, 0.24330E+07, 0.25293E+07, 0.26286E+07, &
         0.27309E+07, 0.28363E+07, 0.29449E+07, 0.30568E+07, 0.31720E+07, &
         0.32905E+07, 0.34125E+07, 0.35381E+07, 0.36672E+07, 0.38000E+07, &
         0.39366E+07, 0.40770E+07, 0.42213E+07, 0.43696E+07, 0.45220E+07, &
         0.46785E+07, 0.48392E+07, 0.50043E+07, 0.51737E+07, 0.53476E+07, &
         0.55261E+07/
!...        --       686
      DATA (QOFT(3,J),J=1,119)/ 0.31656E+03, 0.53355E+03, 0.78557E+03, &
         0.10688E+04, 0.13815E+04, 0.17235E+04, 0.20960E+04, 0.25011E+04, &
         0.29420E+04, 0.34223E+04, 0.39459E+04, 0.45172E+04, 0.51408E+04, &
         0.58213E+04, 0.65639E+04, 0.73735E+04, 0.82555E+04, 0.92152E+04, &
         0.10259E+05, 0.11391E+05, 0.12619E+05, 0.13949E+05, 0.15387E+05, &
         0.16940E+05, 0.18614E+05, 0.20417E+05, 0.22357E+05, 0.24440E+05, &
         0.26675E+05, 0.29070E+05, 0.31633E+05, 0.34374E+05, 0.37299E+05, &
         0.40420E+05, 0.43746E+05, 0.47285E+05, 0.51049E+05, 0.55047E+05, &
         0.59289E+05, 0.63788E+05, 0.68554E+05, 0.73598E+05, 0.78932E+05, &
         0.84568E+05, 0.90519E+05, 0.96796E+05, 0.10341E+06, 0.11039E+06, &
         0.11772E+06, 0.12544E+06, 0.13356E+06, 0.14208E+06, 0.15103E+06, &
         0.16041E+06, 0.17026E+06, 0.18057E+06, 0.19137E+06, 0.20268E+06, &
         0.21450E+06, 0.22687E+06, 0.23979E+06, 0.25328E+06, 0.26736E+06, &
         0.28206E+06, 0.29738E+06, 0.31336E+06, 0.33000E+06, 0.34733E+06, &
         0.36537E+06, 0.38414E+06, 0.40366E+06, 0.42396E+06, 0.44505E+06, &
         0.46696E+06, 0.48971E+06, 0.51332E+06, 0.53782E+06, 0.56323E+06, &
         0.58958E+06, 0.61689E+06, 0.64518E+06, 0.67448E+06, 0.70482E+06, &
         0.73623E+06, 0.76872E+06, 0.80234E+06, 0.83710E+06, 0.87303E+06, &
         0.91017E+06, 0.94853E+06, 0.98816E+06, 0.10291E+07, 0.10713E+07, &
         0.11149E+07, 0.11599E+07, 0.12063E+07, 0.12541E+07, 0.13034E+07, &
         0.13542E+07, 0.14066E+07, 0.14606E+07, 0.15161E+07, 0.15733E+07, &
         0.16322E+07, 0.16928E+07, 0.17552E+07, 0.18194E+07, 0.18854E+07, &
         0.19532E+07, 0.20230E+07, 0.20947E+07, 0.21684E+07, 0.22441E+07, &
         0.23219E+07, 0.24018E+07, 0.24838E+07, 0.25680E+07, 0.26545E+07, &
         0.27432E+07/
!...        --       667
      DATA (QOFT(4,J),J=1,119)/ 0.37657E+04, 0.63472E+04, 0.93454E+04, &
         0.12715E+05, 0.16435E+05, 0.20502E+05, 0.24929E+05, 0.29742E+05, &
         0.34975E+05, 0.40668E+05, 0.46868E+05, 0.53624E+05, 0.60990E+05, &
         0.69018E+05, 0.77768E+05, 0.87296E+05, 0.97666E+05, 0.10894E+06, &
         0.12118E+06, 0.13446E+06, 0.14885E+06, 0.16441E+06, 0.18123E+06, &
         0.19938E+06, 0.21894E+06, 0.23998E+06, 0.26261E+06, 0.28690E+06, &
         0.31295E+06, 0.34084E+06, 0.37068E+06, 0.40256E+06, 0.43659E+06, &
         0.47287E+06, 0.51151E+06, 0.55262E+06, 0.59632E+06, 0.64272E+06, &
         0.69194E+06, 0.74412E+06, 0.79937E+06, 0.85783E+06, 0.91963E+06, &
         0.98492E+06, 0.10538E+07, 0.11265E+07, 0.12031E+07, 0.12837E+07, &
         0.13686E+07, 0.14579E+07, 0.15517E+07, 0.16502E+07, 0.17536E+07, &
         0.18621E+07, 0.19758E+07, 0.20949E+07, 0.22196E+07, 0.23501E+07, &
         0.24866E+07, 0.26292E+07, 0.27783E+07, 0.29339E+07, 0.30963E+07, &
         0.32658E+07, 0.34425E+07, 0.36266E+07, 0.38184E+07, 0.40181E+07, &
         0.42260E+07, 0.44422E+07, 0.46671E+07, 0.49008E+07, 0.51437E+07, &
         0.53959E+07, 0.56578E+07, 0.59296E+07, 0.62116E+07, 0.65040E+07, &
         0.68071E+07, 0.71213E+07, 0.74468E+07, 0.77838E+07, 0.81328E+07, &
         0.84939E+07, 0.88676E+07, 0.92541E+07, 0.96536E+07, 0.10067E+08, &
         0.10493E+08, 0.10934E+08, 0.11390E+08, 0.11860E+08, 0.12345E+08, &
         0.12846E+08, 0.13363E+08, 0.13895E+08, 0.14445E+08, 0.15011E+08, &
         0.15595E+08, 0.16196E+08, 0.16815E+08, 0.17453E+08, 0.18110E+08, &
         0.18786E+08, 0.19482E+08, 0.20198E+08, 0.20934E+08, 0.21691E+08, &
         0.22470E+08, 0.23270E+08, 0.24093E+08, 0.24939E+08, 0.25807E+08, &
         0.26699E+08, 0.27616E+08, 0.28556E+08, 0.29522E+08, 0.30514E+08, &
         0.31531E+08/
!...        --       676
      DATA (QOFT(5,J),J=1,119)/ 0.18608E+04, 0.31363E+04, 0.46177E+04, &
         0.62826E+04, 0.81202E+04, 0.10129E+05, 0.12316E+05, 0.14693E+05, &
         0.17277E+05, 0.20089E+05, 0.23153E+05, 0.26492E+05, 0.30133E+05, &
         0.34103E+05, 0.38430E+05, 0.43145E+05, 0.48277E+05, 0.53858E+05, &
         0.59920E+05, 0.66497E+05, 0.73624E+05, 0.81336E+05, 0.89671E+05, &
         0.98668E+05, 0.10836E+06, 0.11880E+06, 0.13002E+06, 0.14207E+06, &
         0.15500E+06, 0.16884E+06, 0.18365E+06, 0.19947E+06, 0.21636E+06, &
         0.23438E+06, 0.25356E+06, 0.27398E+06, 0.29568E+06, 0.31873E+06, &
         0.34318E+06, 0.36911E+06, 0.39656E+06, 0.42561E+06, 0.45632E+06, &
         0.48877E+06, 0.52302E+06, 0.55914E+06, 0.59722E+06, 0.63732E+06, &
         0.67952E+06, 0.72390E+06, 0.77055E+06, 0.81954E+06, 0.87097E+06, &
         0.92491E+06, 0.98146E+06, 0.10407E+07, 0.11027E+07, 0.11677E+07, &
         0.12356E+07, 0.13066E+07, 0.13807E+07, 0.14582E+07, 0.15390E+07, &
         0.16233E+07, 0.17113E+07, 0.18029E+07, 0.18984E+07, 0.19978E+07, &
         0.21012E+07, 0.22089E+07, 0.23208E+07, 0.24372E+07, 0.25581E+07, &
         0.26837E+07, 0.28141E+07, 0.29494E+07, 0.30898E+07, 0.32354E+07, &
         0.33864E+07, 0.35428E+07, 0.37049E+07, 0.38728E+07, 0.40466E+07, &
         0.42264E+07, 0.44125E+07, 0.46050E+07, 0.48040E+07, 0.50098E+07, &
         0.52224E+07, 0.54420E+07, 0.56689E+07, 0.59031E+07, 0.61449E+07, &
         0.63943E+07, 0.66517E+07, 0.69172E+07, 0.71909E+07, 0.74731E+07, &
         0.77639E+07, 0.80635E+07, 0.83721E+07, 0.86900E+07, 0.90172E+07, &
         0.93541E+07, 0.97008E+07, 0.10058E+08, 0.10424E+08, 0.10802E+08, &
         0.11190E+08, 0.11589E+08, 0.11999E+08, 0.12420E+08, 0.12853E+08, &
         0.13298E+08, 0.13755E+08, 0.14223E+08, 0.14705E+08, 0.15199E+08, &
         0.15706E+08/
!...        --       886
      DATA (QOFT(6,J),J=1,119)/ 0.67639E+03, 0.11401E+04, 0.16787E+04, &
         0.22843E+04, 0.29532E+04, 0.36856E+04, 0.44842E+04, 0.53545E+04, &
         0.63030E+04, 0.73381E+04, 0.84686E+04, 0.97040E+04, 0.11054E+05, &
         0.12530E+05, 0.14143E+05, 0.15903E+05, 0.17823E+05, 0.19915E+05, &
         0.22190E+05, 0.24663E+05, 0.27346E+05, 0.30254E+05, 0.33400E+05, &
         0.36800E+05, 0.40469E+05, 0.44423E+05, 0.48678E+05, 0.53251E+05, &
         0.58160E+05, 0.63423E+05, 0.69058E+05, 0.75085E+05, 0.81524E+05, &
         0.88395E+05, 0.95719E+05, 0.10352E+06, 0.11181E+06, 0.12063E+06, &
         0.12999E+06, 0.13991E+06, 0.15043E+06, 0.16157E+06, 0.17335E+06, &
         0.18580E+06, 0.19895E+06, 0.21283E+06, 0.22746E+06, 0.24288E+06, &
         0.25911E+06, 0.27619E+06, 0.29415E+06, 0.31301E+06, 0.33283E+06, &
         0.35362E+06, 0.37542E+06, 0.39827E+06, 0.42221E+06, 0.44726E+06, &
         0.47348E+06, 0.50089E+06, 0.52954E+06, 0.55947E+06, 0.59072E+06, &
         0.62332E+06, 0.65733E+06, 0.69279E+06, 0.72973E+06, 0.76821E+06, &
         0.80827E+06, 0.84996E+06, 0.89332E+06, 0.93840E+06, 0.98526E+06, &
         0.10339E+07, 0.10845E+07, 0.11370E+07, 0.11914E+07, 0.12479E+07, &
         0.13065E+07, 0.13672E+07, 0.14302E+07, 0.14953E+07, 0.15628E+07, &
         0.16327E+07, 0.17050E+07, 0.17798E+07, 0.18571E+07, 0.19371E+07, &
         0.20197E+07, 0.21051E+07, 0.21933E+07, 0.22844E+07, 0.23785E+07, &
         0.24755E+07, 0.25757E+07, 0.26790E+07, 0.27855E+07, 0.28954E+07, &
         0.30086E+07, 0.31253E+07, 0.32455E+07, 0.33693E+07, 0.34967E+07, &
         0.36280E+07, 0.37631E+07, 0.39021E+07, 0.40451E+07, 0.41922E+07, &
         0.43435E+07, 0.44990E+07, 0.46589E+07, 0.48232E+07, 0.49920E+07, &
         0.51654E+07, 0.53436E+07, 0.55265E+07, 0.57143E+07, 0.59071E+07, &
         0.61050E+07/
!...        --       868
      DATA (QOFT(7,J),J=1,119)/ 0.34615E+03, 0.58348E+03, 0.85915E+03, &
         0.11692E+04, 0.15117E+04, 0.18868E+04, 0.22960E+04, 0.27419E+04, &
         0.32278E+04, 0.37579E+04, 0.43366E+04, 0.49686E+04, 0.56591E+04, &
         0.64134E+04, 0.72369E+04, 0.81354E+04, 0.91148E+04, 0.10181E+05, &
         0.11341E+05, 0.12600E+05, 0.13966E+05, 0.15446E+05, 0.17046E+05, &
         0.18775E+05, 0.20640E+05, 0.22649E+05, 0.24810E+05, 0.27132E+05, &
         0.29624E+05, 0.32295E+05, 0.35154E+05, 0.38211E+05, 0.41475E+05, &
         0.44958E+05, 0.48670E+05, 0.52621E+05, 0.56823E+05, 0.61288E+05, &
         0.66026E+05, 0.71052E+05, 0.76376E+05, 0.82011E+05, 0.87972E+05, &
         0.94271E+05, 0.10092E+06, 0.10794E+06, 0.11534E+06, 0.12313E+06, &
         0.13134E+06, 0.13997E+06, 0.14905E+06, 0.15858E+06, 0.16859E+06, &
         0.17909E+06, 0.19010E+06, 0.20164E+06, 0.21373E+06, 0.22638E+06, &
         0.23962E+06, 0.25346E+06, 0.26792E+06, 0.28302E+06, 0.29879E+06, &
         0.31524E+06, 0.33240E+06, 0.35029E+06, 0.36892E+06, 0.38833E+06, &
         0.40853E+06, 0.42956E+06, 0.45142E+06, 0.47416E+06, 0.49778E+06, &
         0.52233E+06, 0.54781E+06, 0.57427E+06, 0.60172E+06, 0.63019E+06, &
         0.65971E+06, 0.69031E+06, 0.72201E+06, 0.75485E+06, 0.78886E+06, &
         0.82405E+06, 0.86048E+06, 0.89815E+06, 0.93711E+06, 0.97739E+06, &
         0.10190E+07, 0.10620E+07, 0.11065E+07, 0.11523E+07, 0.11997E+07, &
         0.12485E+07, 0.12990E+07, 0.13510E+07, 0.14046E+07, 0.14599E+07, &
         0.15169E+07, 0.15756E+07, 0.16361E+07, 0.16984E+07, 0.17626E+07, &
         0.18287E+07, 0.18966E+07, 0.19666E+07, 0.20386E+07, 0.21126E+07, &
         0.21887E+07, 0.22669E+07, 0.23474E+07, 0.24300E+07, 0.25150E+07, &
         0.26022E+07, 0.26919E+07, 0.27839E+07, 0.28784E+07, 0.29753E+07, &
         0.30749E+07/
!...        --       678
      DATA (QOFT(8,J),J=1,119)/ 0.39745E+04, 0.66993E+04, 0.98642E+04, &
         0.13422E+05, 0.17352E+05, 0.21652E+05, 0.26339E+05, 0.31442E+05, &
         0.37000E+05, 0.43058E+05, 0.49669E+05, 0.56885E+05, 0.64766E+05, &
         0.73372E+05, 0.82765E+05, 0.93011E+05, 0.10418E+06, 0.11633E+06, &
         0.12955E+06, 0.14390E+06, 0.15946E+06, 0.17632E+06, 0.19455E+06, &
         0.21424E+06, 0.23547E+06, 0.25835E+06, 0.28296E+06, 0.30939E+06, &
         0.33776E+06, 0.36816E+06, 0.40070E+06, 0.43549E+06, 0.47264E+06, &
         0.51228E+06, 0.55451E+06, 0.59947E+06, 0.64728E+06, 0.69807E+06, &
         0.75198E+06, 0.80915E+06, 0.86971E+06, 0.93381E+06, 0.10016E+07, &
         0.10733E+07, 0.11489E+07, 0.12287E+07, 0.13128E+07, 0.14015E+07, &
         0.14948E+07, 0.15930E+07, 0.16961E+07, 0.18045E+07, 0.19183E+07, &
         0.20378E+07, 0.21629E+07, 0.22942E+07, 0.24316E+07, 0.25754E+07, &
         0.27258E+07, 0.28831E+07, 0.30475E+07, 0.32192E+07, 0.33984E+07, &
         0.35855E+07, 0.37805E+07, 0.39838E+07, 0.41956E+07, 0.44162E+07, &
         0.46458E+07, 0.48847E+07, 0.51332E+07, 0.53916E+07, 0.56601E+07, &
         0.59390E+07, 0.62286E+07, 0.65292E+07, 0.68412E+07, 0.71647E+07, &
         0.75002E+07, 0.78479E+07, 0.82081E+07, 0.85813E+07, 0.89676E+07, &
         0.93676E+07, 0.97814E+07, 0.10209E+08, 0.10652E+08, 0.11110E+08, &
         0.11583E+08, 0.12071E+08, 0.12576E+08, 0.13097E+08, 0.13635E+08, &
         0.14190E+08, 0.14763E+08, 0.15354E+08, 0.15963E+08, 0.16592E+08, &
         0.17239E+08, 0.17906E+08, 0.18593E+08, 0.19301E+08, 0.20030E+08, &
         0.20780E+08, 0.21553E+08, 0.22347E+08, 0.23165E+08, 0.24006E+08, &
         0.24870E+08, 0.25759E+08, 0.26673E+08, 0.27612E+08, 0.28577E+08, &
         0.29568E+08, 0.30585E+08, 0.31631E+08, 0.32704E+08, 0.33805E+08, &
         0.34936E+08/
!...        --       768
      DATA (QOFT(9,J),J=1,119)/ 0.40228E+04, 0.67808E+04, 0.99842E+04, &
         0.13586E+05, 0.17564E+05, 0.21919E+05, 0.26665E+05, 0.31833E+05, &
         0.37461E+05, 0.43596E+05, 0.50286E+05, 0.57589E+05, 0.65562E+05, &
         0.74264E+05, 0.83761E+05, 0.94115E+05, 0.10540E+06, 0.11767E+06, &
         0.13102E+06, 0.14550E+06, 0.16121E+06, 0.17822E+06, 0.19661E+06, &
         0.21646E+06, 0.23788E+06, 0.26094E+06, 0.28574E+06, 0.31239E+06, &
         0.34097E+06, 0.37160E+06, 0.40437E+06, 0.43941E+06, 0.47683E+06, &
         0.51673E+06, 0.55925E+06, 0.60451E+06, 0.65262E+06, 0.70374E+06, &
         0.75799E+06, 0.81550E+06, 0.87643E+06, 0.94092E+06, 0.10091E+07, &
         0.10812E+07, 0.11572E+07, 0.12375E+07, 0.13221E+07, 0.14112E+07, &
         0.15050E+07, 0.16037E+07, 0.17074E+07, 0.18164E+07, 0.19307E+07, &
         0.20507E+07, 0.21765E+07, 0.23084E+07, 0.24464E+07, 0.25909E+07, &
         0.27421E+07, 0.29001E+07, 0.30652E+07, 0.32377E+07, 0.34177E+07, &
         0.36055E+07, 0.38014E+07, 0.40055E+07, 0.42182E+07, 0.44397E+07, &
         0.46703E+07, 0.49102E+07, 0.51597E+07, 0.54191E+07, 0.56886E+07, &
         0.59686E+07, 0.62593E+07, 0.65611E+07, 0.68742E+07, 0.71989E+07, &
         0.75356E+07, 0.78846E+07, 0.82461E+07, 0.86206E+07, 0.90083E+07, &
         0.94097E+07, 0.98249E+07, 0.10254E+08, 0.10699E+08, 0.11158E+08, &
         0.11632E+08, 0.12123E+08, 0.12629E+08, 0.13152E+08, 0.13691E+08, &
         0.14248E+08, 0.14823E+08, 0.15416E+08, 0.16027E+08, 0.16657E+08, &
         0.17307E+08, 0.17976E+08, 0.18665E+08, 0.19375E+08, 0.20106E+08, &
         0.20858E+08, 0.21633E+08, 0.22430E+08, 0.23250E+08, 0.24093E+08, &
         0.24960E+08, 0.25851E+08, 0.26767E+08, 0.27709E+08, 0.28676E+08, &
         0.29670E+08, 0.30691E+08, 0.31739E+08, 0.32815E+08, 0.33919E+08, &
         0.35053E+08/
!...        --       786
      DATA (QOFT(10,J),J=1,119)/ 0.39315E+04, 0.66267E+04, 0.97569E+04, &
         0.13276E+05, 0.17162E+05, 0.21414E+05, 0.26048E+05, 0.31094E+05, &
         0.36590E+05, 0.42581E+05, 0.49120E+05, 0.56260E+05, 0.64061E+05, &
         0.72580E+05, 0.81882E+05, 0.92031E+05, 0.10309E+06, 0.11514E+06, &
         0.12824E+06, 0.14247E+06, 0.15791E+06, 0.17463E+06, 0.19272E+06, &
         0.21226E+06, 0.23333E+06, 0.25604E+06, 0.28047E+06, 0.30673E+06, &
         0.33490E+06, 0.36510E+06, 0.39743E+06, 0.43200E+06, 0.46892E+06, &
         0.50831E+06, 0.55029E+06, 0.59498E+06, 0.64251E+06, 0.69301E+06, &
         0.74662E+06, 0.80347E+06, 0.86370E+06, 0.92747E+06, 0.99491E+06, &
         0.10662E+07, 0.11414E+07, 0.12208E+07, 0.13046E+07, 0.13928E+07, &
         0.14856E+07, 0.15833E+07, 0.16860E+07, 0.17939E+07, 0.19072E+07, &
         0.20261E+07, 0.21508E+07, 0.22814E+07, 0.24182E+07, 0.25614E+07, &
         0.27112E+07, 0.28679E+07, 0.30316E+07, 0.32026E+07, 0.33811E+07, &
         0.35674E+07, 0.37617E+07, 0.39642E+07, 0.41752E+07, 0.43950E+07, &
         0.46237E+07, 0.48618E+07, 0.51094E+07, 0.53668E+07, 0.56343E+07, &
         0.59123E+07, 0.62009E+07, 0.65005E+07, 0.68113E+07, 0.71338E+07, &
         0.74681E+07, 0.78147E+07, 0.81737E+07, 0.85457E+07, 0.89308E+07, &
         0.93295E+07, 0.97420E+07, 0.10169E+08, 0.10610E+08, 0.11066E+08, &
         0.11538E+08, 0.12025E+08, 0.12528E+08, 0.13048E+08, 0.13584E+08, &
         0.14138E+08, 0.14709E+08, 0.15298E+08, 0.15906E+08, 0.16532E+08, &
         0.17178E+08, 0.17843E+08, 0.18528E+08, 0.19234E+08, 0.19961E+08, &
         0.20710E+08, 0.21480E+08, 0.22272E+08, 0.23088E+08, 0.23926E+08, &
         0.24789E+08, 0.25675E+08, 0.26587E+08, 0.27523E+08, 0.28485E+08, &
         0.29474E+08, 0.30489E+08, 0.31532E+08, 0.32603E+08, 0.33701E+08, &
         0.34829E+08/
!...        --       776
      DATA (QOFT(11,J),J=1,119)/ 0.23106E+05, 0.38945E+05, 0.57342E+05, &
         0.78021E+05, 0.10085E+06, 0.12582E+06, 0.15302E+06, 0.18262E+06, &
         0.21482E+06, 0.24989E+06, 0.28812E+06, 0.32983E+06, 0.37535E+06, &
         0.42501E+06, 0.47919E+06, 0.53825E+06, 0.60258E+06, 0.67256E+06, &
         0.74862E+06, 0.83118E+06, 0.92069E+06, 0.10176E+07, 0.11223E+07, &
         0.12354E+07, 0.13574E+07, 0.14887E+07, 0.16299E+07, 0.17816E+07, &
         0.19443E+07, 0.21187E+07, 0.23052E+07, 0.25047E+07, 0.27176E+07, &
         0.29447E+07, 0.31866E+07, 0.34441E+07, 0.37179E+07, 0.40087E+07, &
         0.43173E+07, 0.46444E+07, 0.49910E+07, 0.53578E+07, 0.57456E+07, &
         0.61554E+07, 0.65880E+07, 0.70444E+07, 0.75255E+07, 0.80322E+07, &
         0.85656E+07, 0.91266E+07, 0.97163E+07, 0.10336E+08, 0.10986E+08, &
         0.11668E+08, 0.12383E+08, 0.13133E+08, 0.13918E+08, 0.14739E+08, &
         0.15598E+08, 0.16496E+08, 0.17435E+08, 0.18415E+08, 0.19438E+08, &
         0.20505E+08, 0.21619E+08, 0.22779E+08, 0.23987E+08, 0.25246E+08, &
         0.26556E+08, 0.27920E+08, 0.29337E+08, 0.30811E+08, 0.32343E+08, &
         0.33934E+08, 0.35585E+08, 0.37300E+08, 0.39079E+08, 0.40924E+08, &
         0.42837E+08, 0.44819E+08, 0.46873E+08, 0.49001E+08, 0.51203E+08, &
         0.53483E+08, 0.55842E+08, 0.58282E+08, 0.60805E+08, 0.63414E+08, &
         0.66109E+08, 0.68894E+08, 0.71770E+08, 0.74740E+08, 0.77806E+08, &
         0.80970E+08, 0.84234E+08, 0.87600E+08, 0.91072E+08, 0.94651E+08, &
         0.98339E+08, 0.10214E+09, 0.10605E+09, 0.11009E+09, 0.11424E+09, &
         0.11851E+09, 0.12291E+09, 0.12744E+09, 0.13209E+09, 0.13688E+09, &
         0.14180E+09, 0.14687E+09, 0.15207E+09, 0.15742E+09, 0.16291E+09, &
         0.16855E+09, 0.17435E+09, 0.18030E+09, 0.18641E+09, 0.19268E+09, &
         0.19912E+09/
!...        --       767
      DATA (QOFT(12,J),J=1,119)/ 0.11692E+05, 0.19707E+05, 0.29017E+05, &
         0.39482E+05, 0.51038E+05, 0.63680E+05, 0.77450E+05, 0.92432E+05, &
         0.10873E+06, 0.12649E+06, 0.14584E+06, 0.16694E+06, 0.18996E+06, &
         0.21507E+06, 0.24245E+06, 0.27229E+06, 0.30478E+06, 0.34013E+06, &
         0.37853E+06, 0.42020E+06, 0.46536E+06, 0.51424E+06, 0.56708E+06, &
         0.62411E+06, 0.68559E+06, 0.75178E+06, 0.82296E+06, 0.89939E+06, &
         0.98137E+06, 0.10692E+07, 0.11631E+07, 0.12636E+07, 0.13708E+07, &
         0.14851E+07, 0.16069E+07, 0.17365E+07, 0.18742E+07, 0.20206E+07, &
         0.21758E+07, 0.23404E+07, 0.25148E+07, 0.26992E+07, 0.28943E+07, &
         0.31004E+07, 0.33179E+07, 0.35474E+07, 0.37892E+07, 0.40440E+07, &
         0.43121E+07, 0.45940E+07, 0.48904E+07, 0.52017E+07, 0.55285E+07, &
         0.58713E+07, 0.62306E+07, 0.66071E+07, 0.70014E+07, 0.74140E+07, &
         0.78456E+07, 0.82967E+07, 0.87681E+07, 0.92604E+07, 0.97742E+07, &
         0.10310E+08, 0.10869E+08, 0.11452E+08, 0.12059E+08, 0.12691E+08, &
         0.13348E+08, 0.14033E+08, 0.14745E+08, 0.15484E+08, 0.16253E+08, &
         0.17052E+08, 0.17881E+08, 0.18741E+08, 0.19634E+08, 0.20560E+08, &
         0.21520E+08, 0.22515E+08, 0.23546E+08, 0.24613E+08, 0.25718E+08, &
         0.26862E+08, 0.28046E+08, 0.29270E+08, 0.30536E+08, 0.31845E+08, &
         0.33197E+08, 0.34594E+08, 0.36037E+08, 0.37527E+08, 0.39065E+08, &
         0.40652E+08, 0.42289E+08, 0.43977E+08, 0.45719E+08, 0.47514E+08, &
         0.49363E+08, 0.51270E+08, 0.53233E+08, 0.55255E+08, 0.57337E+08, &
         0.59480E+08, 0.61686E+08, 0.63956E+08, 0.66290E+08, 0.68691E+08, &
         0.71160E+08, 0.73699E+08, 0.76307E+08, 0.78988E+08, 0.81743E+08, &
         0.84572E+08, 0.87478E+08, 0.90462E+08, 0.93525E+08, 0.96669E+08, &
         0.99896E+08/
!...        --       888
      DATA (QOFT(13,J),J=1,119)/ 0.36175E+03, 0.60978E+03, 0.89790E+03, &
         0.12219E+04, 0.15802E+04, 0.19728E+04, 0.24016E+04, 0.28696E+04, &
         0.33807E+04, 0.39394E+04, 0.45506E+04, 0.52196E+04, 0.59521E+04, &
         0.67538E+04, 0.76308E+04, 0.85894E+04, 0.96361E+04, 0.10777E+05, &
         0.12021E+05, 0.13373E+05, 0.14841E+05, 0.16434E+05, 0.18158E+05, &
         0.20023E+05, 0.22037E+05, 0.24208E+05, 0.26547E+05, 0.29061E+05, &
         0.31762E+05, 0.34659E+05, 0.37762E+05, 0.41083E+05, 0.44632E+05, &
         0.48421E+05, 0.52462E+05, 0.56766E+05, 0.61346E+05, 0.66215E+05, &
         0.71386E+05, 0.76873E+05, 0.82688E+05, 0.88848E+05, 0.95365E+05, &
         0.10226E+06, 0.10954E+06, 0.11722E+06, 0.12532E+06, 0.13387E+06, &
         0.14286E+06, 0.15233E+06, 0.16229E+06, 0.17275E+06, 0.18374E+06, &
         0.19528E+06, 0.20737E+06, 0.22006E+06, 0.23335E+06, 0.24726E+06, &
         0.26182E+06, 0.27705E+06, 0.29297E+06, 0.30960E+06, 0.32696E+06, &
         0.34509E+06, 0.36399E+06, 0.38371E+06, 0.40425E+06, 0.42566E+06, &
         0.44794E+06, 0.47114E+06, 0.49527E+06, 0.52036E+06, 0.54644E+06, &
         0.57354E+06, 0.60169E+06, 0.63091E+06, 0.66124E+06, 0.69270E+06, &
         0.72533E+06, 0.75916E+06, 0.79421E+06, 0.83053E+06, 0.86814E+06, &
         0.90708E+06, 0.94737E+06, 0.98907E+06, 0.10322E+07, 0.10768E+07, &
         0.11229E+07, 0.11705E+07, 0.12197E+07, 0.12705E+07, 0.13230E+07, &
         0.13771E+07, 0.14330E+07, 0.14906E+07, 0.15501E+07, 0.16114E+07, &
         0.16745E+07, 0.17397E+07, 0.18067E+07, 0.18759E+07, 0.19470E+07, &
         0.20203E+07, 0.20957E+07, 0.21733E+07, 0.22532E+07, 0.23353E+07, &
         0.24198E+07, 0.25067E+07, 0.25960E+07, 0.26878E+07, 0.27821E+07, &
         0.28790E+07, 0.29785E+07, 0.30807E+07, 0.31857E+07, 0.32934E+07, &
         0.34040E+07/
!...        --       887
      DATA (QOFT(14,J),J=1,119)/ 0.42000E+04, 0.70796E+04, 0.10424E+05, &
         0.14186E+05, 0.18342E+05, 0.22896E+05, 0.27866E+05, 0.33285E+05, &
         0.39199E+05, 0.45659E+05, 0.52720E+05, 0.60444E+05, 0.68895E+05, &
         0.78139E+05, 0.88246E+05, 0.99288E+05, 0.11134E+06, 0.12447E+06, &
         0.13877E+06, 0.15431E+06, 0.17119E+06, 0.18949E+06, 0.20930E+06, &
         0.23071E+06, 0.25383E+06, 0.27875E+06, 0.30558E+06, 0.33442E+06, &
         0.36539E+06, 0.39861E+06, 0.43418E+06, 0.47224E+06, 0.51291E+06, &
         0.55632E+06, 0.60260E+06, 0.65189E+06, 0.70434E+06, 0.76008E+06, &
         0.81927E+06, 0.88206E+06, 0.94862E+06, 0.10191E+07, 0.10937E+07, &
         0.11725E+07, 0.12558E+07, 0.13436E+07, 0.14363E+07, 0.15340E+07, &
         0.16368E+07, 0.17450E+07, 0.18588E+07, 0.19784E+07, 0.21040E+07, &
         0.22358E+07, 0.23741E+07, 0.25190E+07, 0.26708E+07, 0.28297E+07, &
         0.29961E+07, 0.31700E+07, 0.33518E+07, 0.35417E+07, 0.37400E+07, &
         0.39469E+07, 0.41628E+07, 0.43878E+07, 0.46224E+07, 0.48667E+07, &
         0.51210E+07, 0.53858E+07, 0.56611E+07, 0.59475E+07, 0.62451E+07, &
         0.65544E+07, 0.68755E+07, 0.72089E+07, 0.75550E+07, 0.79139E+07, &
         0.82861E+07, 0.86720E+07, 0.90719E+07, 0.94861E+07, 0.99151E+07, &
         0.10359E+08, 0.10819E+08, 0.11294E+08, 0.11786E+08, 0.12294E+08, &
         0.12820E+08, 0.13363E+08, 0.13924E+08, 0.14503E+08, 0.15101E+08, &
         0.15719E+08, 0.16356E+08, 0.17013E+08, 0.17690E+08, 0.18389E+08, &
         0.19109E+08, 0.19851E+08, 0.20616E+08, 0.21404E+08, 0.22215E+08, &
         0.23050E+08, 0.23910E+08, 0.24794E+08, 0.25704E+08, 0.26640E+08, &
         0.27603E+08, 0.28593E+08, 0.29610E+08, 0.30656E+08, 0.31731E+08, &
         0.32835E+08, 0.33969E+08, 0.35133E+08, 0.36329E+08, 0.37556E+08, &
         0.38816E+08/
!...        --       878
      DATA (QOFT(15,J),J=1,119)/ 0.21250E+04, 0.35820E+04, 0.52744E+04, &
         0.71778E+04, 0.92814E+04, 0.11586E+05, 0.14102E+05, 0.16845E+05, &
         0.19839E+05, 0.23108E+05, 0.26680E+05, 0.30588E+05, 0.34861E+05, &
         0.39534E+05, 0.44642E+05, 0.50219E+05, 0.56305E+05, 0.62937E+05, &
         0.70155E+05, 0.78001E+05, 0.86516E+05, 0.95747E+05, 0.10574E+06, &
         0.11653E+06, 0.12819E+06, 0.14075E+06, 0.15427E+06, 0.16881E+06, &
         0.18441E+06, 0.20114E+06, 0.21906E+06, 0.23823E+06, 0.25871E+06, &
         0.28056E+06, 0.30386E+06, 0.32867E+06, 0.35507E+06, 0.38312E+06, &
         0.41291E+06, 0.44450E+06, 0.47799E+06, 0.51344E+06, 0.55095E+06, &
         0.59060E+06, 0.63248E+06, 0.67667E+06, 0.72327E+06, 0.77238E+06, &
         0.82409E+06, 0.87850E+06, 0.93571E+06, 0.99583E+06, 0.10590E+07, &
         0.11252E+07, 0.11947E+07, 0.12675E+07, 0.13438E+07, 0.14237E+07, &
         0.15072E+07, 0.15946E+07, 0.16859E+07, 0.17814E+07, 0.18810E+07, &
         0.19849E+07, 0.20934E+07, 0.22064E+07, 0.23242E+07, 0.24469E+07, &
         0.25747E+07, 0.27076E+07, 0.28459E+07, 0.29897E+07, 0.31391E+07, &
         0.32944E+07, 0.34557E+07, 0.36231E+07, 0.37968E+07, 0.39770E+07, &
         0.41639E+07, 0.43576E+07, 0.45583E+07, 0.47663E+07, 0.49816E+07, &
         0.52045E+07, 0.54352E+07, 0.56739E+07, 0.59207E+07, 0.61759E+07, &
         0.64396E+07, 0.67121E+07, 0.69936E+07, 0.72844E+07, 0.75845E+07, &
         0.78943E+07, 0.82139E+07, 0.85436E+07, 0.88837E+07, 0.92342E+07, &
         0.95956E+07, 0.99680E+07, 0.10352E+08, 0.10747E+08, 0.11154E+08, &
         0.11573E+08, 0.12004E+08, 0.12448E+08, 0.12904E+08, 0.13374E+08, &
         0.13857E+08, 0.14353E+08, 0.14864E+08, 0.15388E+08, 0.15927E+08, &
         0.16481E+08, 0.17050E+08, 0.17634E+08, 0.18234E+08, 0.18849E+08, &
         0.19481E+08/
!...        --       778
      DATA (QOFT(16,J),J=1,119)/ 0.24692E+05, 0.41621E+05, 0.61284E+05, &
         0.83394E+05, 0.10782E+06, 0.13457E+06, 0.16375E+06, 0.19554E+06, &
         0.23020E+06, 0.26801E+06, 0.30930E+06, 0.35443E+06, 0.40375E+06, &
         0.45763E+06, 0.51650E+06, 0.58075E+06, 0.65080E+06, 0.72711E+06, &
         0.81012E+06, 0.90030E+06, 0.99815E+06, 0.11042E+07, 0.12189E+07, &
         0.13428E+07, 0.14765E+07, 0.16206E+07, 0.17757E+07, 0.19423E+07, &
         0.21212E+07, 0.23129E+07, 0.25181E+07, 0.27377E+07, 0.29721E+07, &
         0.32223E+07, 0.34890E+07, 0.37729E+07, 0.40750E+07, 0.43959E+07, &
         0.47365E+07, 0.50978E+07, 0.54807E+07, 0.58860E+07, 0.63147E+07, &
         0.67678E+07, 0.72463E+07, 0.77512E+07, 0.82836E+07, 0.88445E+07, &
         0.94351E+07, 0.10056E+08, 0.10710E+08, 0.11396E+08, 0.12117E+08, &
         0.12873E+08, 0.13666E+08, 0.14497E+08, 0.15367E+08, 0.16279E+08, &
         0.17232E+08, 0.18229E+08, 0.19271E+08, 0.20359E+08, 0.21495E+08, &
         0.22681E+08, 0.23917E+08, 0.25206E+08, 0.26549E+08, 0.27948E+08, &
         0.29404E+08, 0.30920E+08, 0.32496E+08, 0.34135E+08, 0.35838E+08, &
         0.37608E+08, 0.39445E+08, 0.41353E+08, 0.43332E+08, 0.45385E+08, &
         0.47514E+08, 0.49721E+08, 0.52007E+08, 0.54376E+08, 0.56829E+08, &
         0.59367E+08, 0.61995E+08, 0.64712E+08, 0.67523E+08, 0.70429E+08, &
         0.73432E+08, 0.76535E+08, 0.79740E+08, 0.83050E+08, 0.86467E+08, &
         0.89993E+08, 0.93632E+08, 0.97385E+08, 0.10126E+09, 0.10525E+09, &
         0.10936E+09, 0.11360E+09, 0.11796E+09, 0.12246E+09, 0.12709E+09, &
         0.13186E+09, 0.13677E+09, 0.14182E+09, 0.14701E+09, 0.15236E+09, &
         0.15785E+09, 0.16350E+09, 0.16931E+09, 0.17528E+09, 0.18141E+09, &
         0.18771E+09, 0.19418E+09, 0.20082E+09, 0.20764E+09, 0.21465E+09, &
         0.22183E+09/
!...        --       787
      DATA (QOFT(17,J),J=1,119)/ 0.12211E+05, 0.20582E+05, 0.30305E+05, &
         0.41237E+05, 0.53314E+05, 0.66536E+05, 0.80957E+05, 0.96672E+05, &
         0.11380E+06, 0.13250E+06, 0.15292E+06, 0.17524E+06, 0.19965E+06, &
         0.22632E+06, 0.25546E+06, 0.28728E+06, 0.32199E+06, 0.35980E+06, &
         0.40094E+06, 0.44565E+06, 0.49417E+06, 0.54676E+06, 0.60366E+06, &
         0.66516E+06, 0.73152E+06, 0.80305E+06, 0.88002E+06, 0.96276E+06, &
         0.10516E+07, 0.11468E+07, 0.12488E+07, 0.13578E+07, 0.14743E+07, &
         0.15987E+07, 0.17312E+07, 0.18723E+07, 0.20225E+07, 0.21820E+07, &
         0.23514E+07, 0.25310E+07, 0.27214E+07, 0.29230E+07, 0.31362E+07, &
         0.33616E+07, 0.35997E+07, 0.38509E+07, 0.41158E+07, 0.43949E+07, &
         0.46887E+07, 0.49980E+07, 0.53231E+07, 0.56647E+07, 0.60234E+07, &
         0.63998E+07, 0.67946E+07, 0.72084E+07, 0.76418E+07, 0.80955E+07, &
         0.85702E+07, 0.90666E+07, 0.95854E+07, 0.10127E+08, 0.10693E+08, &
         0.11284E+08, 0.11900E+08, 0.12542E+08, 0.13211E+08, 0.13907E+08, &
         0.14633E+08, 0.15388E+08, 0.16173E+08, 0.16990E+08, 0.17838E+08, &
         0.18720E+08, 0.19636E+08, 0.20586E+08, 0.21573E+08, 0.22596E+08, &
         0.23657E+08, 0.24757E+08, 0.25896E+08, 0.27077E+08, 0.28299E+08, &
         0.29565E+08, 0.30874E+08, 0.32229E+08, 0.33630E+08, 0.35079E+08, &
         0.36576E+08, 0.38123E+08, 0.39721E+08, 0.41371E+08, 0.43075E+08, &
         0.44833E+08, 0.46647E+08, 0.48518E+08, 0.50448E+08, 0.52438E+08, &
         0.54489E+08, 0.56603E+08, 0.58780E+08, 0.61023E+08, 0.63332E+08, &
         0.65710E+08, 0.68157E+08, 0.70676E+08, 0.73266E+08, 0.75931E+08, &
         0.78672E+08, 0.81490E+08, 0.84386E+08, 0.87363E+08, 0.90422E+08, &
         0.93564E+08, 0.96791E+08, 0.10011E+09, 0.10351E+09, 0.10700E+09, &
         0.11059E+09/
!...        --       777
      DATA (QOFT(18,J),J=1,119)/ 0.71750E+05, 0.12094E+06, 0.17807E+06, &
         0.24230E+06, 0.31324E+06, 0.39088E+06, 0.47550E+06, 0.56764E+06, &
         0.66800E+06, 0.77740E+06, 0.89677E+06, 0.10271E+07, 0.11694E+07, &
         0.13249E+07, 0.14945E+07, 0.16796E+07, 0.18813E+07, 0.21009E+07, &
         0.23396E+07, 0.25989E+07, 0.28801E+07, 0.31847E+07, 0.35140E+07, &
         0.38698E+07, 0.42535E+07, 0.46669E+07, 0.51115E+07, 0.55893E+07, &
         0.61019E+07, 0.66513E+07, 0.72393E+07, 0.78680E+07, 0.85395E+07, &
         0.92558E+07, 0.10019E+08, 0.10832E+08, 0.11696E+08, 0.12614E+08, &
         0.13588E+08, 0.14621E+08, 0.15716E+08, 0.16875E+08, 0.18100E+08, &
         0.19395E+08, 0.20762E+08, 0.22205E+08, 0.23726E+08, 0.25328E+08, &
         0.27015E+08, 0.28789E+08, 0.30654E+08, 0.32614E+08, 0.34671E+08, &
         0.36830E+08, 0.39093E+08, 0.41465E+08, 0.43949E+08, 0.46549E+08, &
         0.49269E+08, 0.52112E+08, 0.55084E+08, 0.58188E+08, 0.61428E+08, &
         0.64809E+08, 0.68335E+08, 0.72010E+08, 0.75840E+08, 0.79828E+08, &
         0.83979E+08, 0.88299E+08, 0.92792E+08, 0.97463E+08, 0.10232E+09, &
         0.10736E+09, 0.11260E+09, 0.11803E+09, 0.12367E+09, 0.12952E+09, &
         0.13559E+09, 0.14187E+09, 0.14839E+09, 0.15513E+09, 0.16212E+09, &
         0.16935E+09, 0.17683E+09, 0.18457E+09, 0.19257E+09, 0.20085E+09, &
         0.20940E+09, 0.21824E+09, 0.22736E+09, 0.23678E+09, 0.24651E+09, &
         0.25655E+09, 0.26691E+09, 0.27759E+09, 0.28861E+09, 0.29997E+09, &
         0.31167E+09, 0.32374E+09, 0.33616E+09, 0.34896E+09, 0.36214E+09, &
         0.37571E+09, 0.38967E+09, 0.40404E+09, 0.41882E+09, 0.43403E+09, &
         0.44966E+09, 0.46573E+09, 0.48226E+09, 0.49923E+09, 0.51668E+09, &
         0.53460E+09, 0.55301E+09, 0.57191E+09, 0.59131E+09, 0.61123E+09, &
         0.63167E+09/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_O3


!
!     *****************
      SUBROUTINE QT_N2O(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(5) :: XGJ
      REAL(DOUBLE), DIMENSION(5,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 9., 6., 6., 9., 54./
!...      N2O
!...        --       446
      DATA (QOFT(1,J),J=1,119)/ 0.89943E+03, 0.12734E+04, 0.16489E+04, &
         0.20293E+04, 0.24205E+04, 0.28289E+04, 0.32609E+04, 0.37222E+04, &
         0.42180E+04, 0.47529E+04, 0.53312E+04, 0.59572E+04, 0.66348E+04, &
         0.73683E+04, 0.81616E+04, 0.90190E+04, 0.99450E+04, 0.10944E+05, &
         0.12021E+05, 0.13180E+05, 0.14426E+05, 0.15766E+05, 0.17203E+05, &
         0.18745E+05, 0.20396E+05, 0.22162E+05, 0.24051E+05, 0.26069E+05, &
         0.28222E+05, 0.30517E+05, 0.32962E+05, 0.35564E+05, 0.38331E+05, &
         0.41271E+05, 0.44393E+05, 0.47704E+05, 0.51214E+05, 0.54932E+05, &
         0.58868E+05, 0.63030E+05, 0.67429E+05, 0.72075E+05, 0.76979E+05, &
         0.82151E+05, 0.87604E+05, 0.93348E+05, 0.99395E+05, 0.10576E+06, &
         0.11245E+06, 0.11948E+06, 0.12686E+06, 0.13461E+06, 0.14275E+06, &
         0.15128E+06, 0.16021E+06, 0.16958E+06, 0.17938E+06, 0.18964E+06, &
         0.20037E+06, 0.21159E+06, 0.22331E+06, 0.23556E+06, 0.24834E+06, &
         0.26169E+06, 0.27561E+06, 0.29012E+06, 0.30525E+06, 0.32101E+06, &
         0.33743E+06, 0.35452E+06, 0.37230E+06, 0.39080E+06, 0.41004E+06, &
         0.43004E+06, 0.45082E+06, 0.47241E+06, 0.49483E+06, 0.51810E+06, &
         0.54225E+06, 0.56730E+06, 0.59329E+06, 0.62022E+06, 0.64814E+06, &
         0.67707E+06, 0.70703E+06, 0.73806E+06, 0.77018E+06, 0.80342E+06, &
         0.83781E+06, 0.87338E+06, 0.91016E+06, 0.94818E+06, 0.98748E+06, &
         0.10281E+07, 0.10700E+07, 0.11133E+07, 0.11581E+07, 0.12042E+07, &
         0.12519E+07, 0.13010E+07, 0.13517E+07, 0.14040E+07, 0.14579E+07, &
         0.15134E+07, 0.15707E+07, 0.16297E+07, 0.16905E+07, 0.17530E+07, &
         0.18175E+07, 0.18838E+07, 0.19521E+07, 0.20224E+07, 0.20947E+07, &
         0.21690E+07, 0.22455E+07, 0.23242E+07, 0.24050E+07, 0.24881E+07, &
         0.25735E+07/
!...        --       456
      DATA (QOFT(2,J),J=1,119)/ 0.59966E+03, 0.84903E+03, 0.10995E+04, &
         0.13538E+04, 0.16158E+04, 0.18903E+04, 0.21815E+04, 0.24934E+04, &
         0.28295E+04, 0.31927E+04, 0.35862E+04, 0.40128E+04, 0.44752E+04, &
         0.49763E+04, 0.55189E+04, 0.61059E+04, 0.67404E+04, 0.74256E+04, &
         0.81646E+04, 0.89609E+04, 0.98180E+04, 0.10740E+05, 0.11729E+05, &
         0.12791E+05, 0.13930E+05, 0.15149E+05, 0.16453E+05, 0.17847E+05, &
         0.19335E+05, 0.20922E+05, 0.22614E+05, 0.24416E+05, 0.26333E+05, &
         0.28371E+05, 0.30535E+05, 0.32833E+05, 0.35269E+05, 0.37851E+05, &
         0.40585E+05, 0.43478E+05, 0.46537E+05, 0.49769E+05, 0.53182E+05, &
         0.56783E+05, 0.60580E+05, 0.64582E+05, 0.68796E+05, 0.73232E+05, &
         0.77898E+05, 0.82803E+05, 0.87957E+05, 0.93369E+05, 0.99048E+05, &
         0.10501E+06, 0.11125E+06, 0.11780E+06, 0.12465E+06, 0.13182E+06, &
         0.13933E+06, 0.14718E+06, 0.15539E+06, 0.16396E+06, 0.17291E+06, &
         0.18226E+06, 0.19201E+06, 0.20218E+06, 0.21278E+06, 0.22383E+06, &
         0.23534E+06, 0.24733E+06, 0.25980E+06, 0.27278E+06, 0.28628E+06, &
         0.30032E+06, 0.31491E+06, 0.33007E+06, 0.34581E+06, 0.36216E+06, &
         0.37912E+06, 0.39673E+06, 0.41499E+06, 0.43392E+06, 0.45355E+06, &
         0.47389E+06, 0.49496E+06, 0.51678E+06, 0.53937E+06, 0.56276E+06, &
         0.58695E+06, 0.61199E+06, 0.63788E+06, 0.66464E+06, 0.69231E+06, &
         0.72090E+06, 0.75044E+06, 0.78094E+06, 0.81244E+06, 0.84496E+06, &
         0.87853E+06, 0.91316E+06, 0.94889E+06, 0.98573E+06, 0.10237E+07, &
         0.10629E+07, 0.11033E+07, 0.11449E+07, 0.11877E+07, 0.12319E+07, &
         0.12773E+07, 0.13241E+07, 0.13723E+07, 0.14219E+07, 0.14729E+07, &
         0.15254E+07, 0.15793E+07, 0.16349E+07, 0.16919E+07, 0.17506E+07, &
         0.18109E+07/
!...        --       546
      DATA (QOFT(3,J),J=1,119)/ 0.62051E+03, 0.87856E+03, 0.11377E+04, &
         0.14003E+04, 0.16705E+04, 0.19529E+04, 0.22518E+04, 0.25713E+04, &
         0.29149E+04, 0.32859E+04, 0.36873E+04, 0.41220E+04, 0.45929E+04, &
         0.51028E+04, 0.56547E+04, 0.62515E+04, 0.68963E+04, 0.75923E+04, &
         0.83428E+04, 0.91511E+04, 0.10021E+05, 0.10956E+05, 0.11960E+05, &
         0.13036E+05, 0.14190E+05, 0.15425E+05, 0.16746E+05, 0.18158E+05, &
         0.19664E+05, 0.21271E+05, 0.22984E+05, 0.24806E+05, 0.26745E+05, &
         0.28806E+05, 0.30995E+05, 0.33317E+05, 0.35780E+05, 0.38389E+05, &
         0.41151E+05, 0.44073E+05, 0.47162E+05, 0.50425E+05, 0.53871E+05, &
         0.57505E+05, 0.61338E+05, 0.65375E+05, 0.69628E+05, 0.74102E+05, &
         0.78808E+05, 0.83755E+05, 0.88951E+05, 0.94407E+05, 0.10013E+06, &
         0.10614E+06, 0.11243E+06, 0.11902E+06, 0.12593E+06, 0.13316E+06, &
         0.14072E+06, 0.14862E+06, 0.15689E+06, 0.16552E+06, 0.17453E+06, &
         0.18394E+06, 0.19376E+06, 0.20399E+06, 0.21466E+06, 0.22578E+06, &
         0.23737E+06, 0.24942E+06, 0.26198E+06, 0.27503E+06, 0.28861E+06, &
         0.30273E+06, 0.31741E+06, 0.33265E+06, 0.34848E+06, 0.36492E+06, &
         0.38197E+06, 0.39967E+06, 0.41803E+06, 0.43706E+06, 0.45679E+06, &
         0.47723E+06, 0.49840E+06, 0.52033E+06, 0.54303E+06, 0.56653E+06, &
         0.59084E+06, 0.61599E+06, 0.64200E+06, 0.66888E+06, 0.69667E+06, &
         0.72539E+06, 0.75506E+06, 0.78569E+06, 0.81733E+06, 0.84998E+06, &
         0.88369E+06, 0.91846E+06, 0.95433E+06, 0.99132E+06, 0.10295E+07, &
         0.10688E+07, 0.11093E+07, 0.11511E+07, 0.11941E+07, 0.12384E+07, &
         0.12840E+07, 0.13310E+07, 0.13793E+07, 0.14291E+07, 0.14803E+07, &
         0.15329E+07, 0.15871E+07, 0.16428E+07, 0.17000E+07, 0.17589E+07, &
         0.18194E+07/
!...        --       448
      DATA (QOFT(4,J),J=1,119)/ 0.95253E+03, 0.13487E+04, 0.17465E+04, &
         0.21498E+04, 0.25648E+04, 0.29986E+04, 0.34580E+04, 0.39493E+04, &
         0.44779E+04, 0.50488E+04, 0.56669E+04, 0.63366E+04, 0.70625E+04, &
         0.78488E+04, 0.87003E+04, 0.96216E+04, 0.10617E+05, 0.11692E+05, &
         0.12852E+05, 0.14102E+05, 0.15447E+05, 0.16893E+05, 0.18446E+05, &
         0.20112E+05, 0.21898E+05, 0.23811E+05, 0.25856E+05, 0.28042E+05, &
         0.30377E+05, 0.32866E+05, 0.35520E+05, 0.38345E+05, 0.41351E+05, &
         0.44545E+05, 0.47939E+05, 0.51540E+05, 0.55359E+05, 0.59405E+05, &
         0.63689E+05, 0.68222E+05, 0.73015E+05, 0.78078E+05, 0.83424E+05, &
         0.89064E+05, 0.95012E+05, 0.10128E+06, 0.10788E+06, 0.11482E+06, &
         0.12213E+06, 0.12981E+06, 0.13788E+06, 0.14635E+06, 0.15524E+06, &
         0.16456E+06, 0.17433E+06, 0.18457E+06, 0.19530E+06, 0.20652E+06, &
         0.21827E+06, 0.23055E+06, 0.24338E+06, 0.25679E+06, 0.27079E+06, &
         0.28541E+06, 0.30066E+06, 0.31656E+06, 0.33314E+06, 0.35042E+06, &
         0.36841E+06, 0.38715E+06, 0.40666E+06, 0.42695E+06, 0.44805E+06, &
         0.46999E+06, 0.49279E+06, 0.51649E+06, 0.54109E+06, 0.56664E+06, &
         0.59315E+06, 0.62066E+06, 0.64919E+06, 0.67877E+06, 0.70943E+06, &
         0.74121E+06, 0.77413E+06, 0.80822E+06, 0.84351E+06, 0.88004E+06, &
         0.91783E+06, 0.95693E+06, 0.99737E+06, 0.10392E+07, 0.10824E+07, &
         0.11270E+07, 0.11732E+07, 0.12208E+07, 0.12700E+07, 0.13208E+07, &
         0.13732E+07, 0.14272E+07, 0.14830E+07, 0.15405E+07, 0.15999E+07, &
         0.16610E+07, 0.17240E+07, 0.17890E+07, 0.18559E+07, 0.19248E+07, &
         0.19957E+07, 0.20687E+07, 0.21439E+07, 0.22213E+07, 0.23009E+07, &
         0.23828E+07, 0.24671E+07, 0.25537E+07, 0.26428E+07, 0.27343E+07, &
         0.28284E+07/
!...        --       447
      DATA (QOFT(5,J),J=1,119)/ 0.55598E+04, 0.78718E+04, 0.10193E+05, &
         0.12546E+05, 0.14966E+05, 0.17495E+05, 0.20171E+05, 0.23031E+05, &
         0.26106E+05, 0.29426E+05, 0.33018E+05, 0.36908E+05, 0.41121E+05, &
         0.45684E+05, 0.50622E+05, 0.55962E+05, 0.61731E+05, 0.67958E+05, &
         0.74671E+05, 0.81902E+05, 0.89681E+05, 0.98043E+05, 0.10702E+06, &
         0.11665E+06, 0.12697E+06, 0.13801E+06, 0.14983E+06, 0.16244E+06, &
         0.17591E+06, 0.19028E+06, 0.20558E+06, 0.22188E+06, 0.23920E+06, &
         0.25762E+06, 0.27718E+06, 0.29793E+06, 0.31993E+06, 0.34323E+06, &
         0.36791E+06, 0.39401E+06, 0.42160E+06, 0.45074E+06, 0.48151E+06, &
         0.51397E+06, 0.54819E+06, 0.58424E+06, 0.62221E+06, 0.66215E+06, &
         0.70416E+06, 0.74832E+06, 0.79470E+06, 0.84340E+06, 0.89450E+06, &
         0.94808E+06, 0.10042E+07, 0.10631E+07, 0.11247E+07, 0.11892E+07, &
         0.12567E+07, 0.13272E+07, 0.14009E+07, 0.14779E+07, 0.15583E+07, &
         0.16422E+07, 0.17298E+07, 0.18211E+07, 0.19163E+07, 0.20154E+07, &
         0.21187E+07, 0.22263E+07, 0.23382E+07, 0.24546E+07, 0.25757E+07, &
         0.27016E+07, 0.28324E+07, 0.29683E+07, 0.31095E+07, 0.32560E+07, &
         0.34081E+07, 0.35659E+07, 0.37295E+07, 0.38991E+07, 0.40750E+07, &
         0.42572E+07, 0.44459E+07, 0.46414E+07, 0.48437E+07, 0.50531E+07, &
         0.52698E+07, 0.54939E+07, 0.57257E+07, 0.59653E+07, 0.62129E+07, &
         0.64688E+07, 0.67331E+07, 0.70061E+07, 0.72880E+07, 0.75790E+07, &
         0.78792E+07, 0.81891E+07, 0.85086E+07, 0.88382E+07, 0.91780E+07, &
         0.95283E+07, 0.98893E+07, 0.10261E+08, 0.10644E+08, 0.11039E+08, &
         0.11445E+08, 0.11864E+08, 0.12294E+08, 0.12738E+08, 0.13194E+08, &
         0.13663E+08, 0.14145E+08, 0.14641E+08, 0.15151E+08, 0.15675E+08, &
         0.16214E+08/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_N2O


!
!     *****************
      SUBROUTINE QT_CO(T, ISO, GSI, QT)          ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(6) :: XGJ
      REAL(DOUBLE), DIMENSION(6,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1., 2., 1., 6., 2., 12./
!...       CO
!...        --        26
      DATA (QOFT(1,J),J=1,119)/ 0.21942E+02, 0.30949E+02, 0.39960E+02, &
         0.48975E+02, 0.57993E+02, 0.67015E+02, 0.76040E+02, 0.85069E+02, &
         0.94102E+02, 0.10314E+03, 0.11218E+03, 0.12123E+03, 0.13029E+03, &
         0.13936E+03, 0.14845E+03, 0.15756E+03, 0.16669E+03, 0.17585E+03, &
         0.18505E+03, 0.19429E+03, 0.20359E+03, 0.21293E+03, 0.22233E+03, &
         0.23181E+03, 0.24135E+03, 0.25096E+03, 0.26066E+03, 0.27045E+03, &
         0.28032E+03, 0.29030E+03, 0.30037E+03, 0.31053E+03, 0.32081E+03, &
         0.33120E+03, 0.34170E+03, 0.35231E+03, 0.36304E+03, 0.37388E+03, &
         0.38486E+03, 0.39595E+03, 0.40717E+03, 0.41852E+03, 0.42999E+03, &
         0.44160E+03, 0.45334E+03, 0.46522E+03, 0.47723E+03, 0.48937E+03, &
         0.50165E+03, 0.51406E+03, 0.52662E+03, 0.53932E+03, 0.55216E+03, &
         0.56513E+03, 0.57825E+03, 0.59152E+03, 0.60492E+03, 0.61847E+03, &
         0.63217E+03, 0.64601E+03, 0.65999E+03, 0.67412E+03, 0.68840E+03, &
         0.70282E+03, 0.71739E+03, 0.73211E+03, 0.74698E+03, 0.76200E+03, &
         0.77716E+03, 0.79247E+03, 0.80793E+03, 0.82355E+03, 0.83931E+03, &
         0.85522E+03, 0.87128E+03, 0.88749E+03, 0.90386E+03, 0.92037E+03, &
         0.93703E+03, 0.95385E+03, 0.97082E+03, 0.98794E+03, 0.10052E+04, &
         0.10226E+04, 0.10402E+04, 0.10580E+04, 0.10758E+04, 0.10939E+04, &
         0.11121E+04, 0.11304E+04, 0.11489E+04, 0.11675E+04, 0.11864E+04, &
         0.12053E+04, 0.12244E+04, 0.12437E+04, 0.12631E+04, 0.12827E+04, &
         0.13024E+04, 0.13223E+04, 0.13423E+04, 0.13625E+04, 0.13829E+04, &
         0.14034E+04, 0.14240E+04, 0.14448E+04, 0.14658E+04, 0.14870E+04, &
         0.15082E+04, 0.15297E+04, 0.15513E+04, 0.15730E+04, 0.15949E+04, &
         0.16170E+04, 0.16392E+04, 0.16616E+04, 0.16841E+04, 0.17068E+04, &
         0.17296E+04/
!...        --        36
      DATA (QOFT(2,J),J=1,119)/ 0.45875E+02, 0.64721E+02, 0.83574E+02, &
         0.10243E+03, 0.12130E+03, 0.14018E+03, 0.15906E+03, 0.17795E+03, &
         0.19685E+03, 0.21576E+03, 0.23468E+03, 0.25362E+03, 0.27257E+03, &
         0.29156E+03, 0.31059E+03, 0.32966E+03, 0.34879E+03, 0.36799E+03, &
         0.38727E+03, 0.40665E+03, 0.42614E+03, 0.44575E+03, 0.46549E+03, &
         0.48539E+03, 0.50544E+03, 0.52566E+03, 0.54606E+03, 0.56665E+03, &
         0.58744E+03, 0.60843E+03, 0.62965E+03, 0.65108E+03, 0.67275E+03, &
         0.69466E+03, 0.71681E+03, 0.73921E+03, 0.76187E+03, 0.78478E+03, &
         0.80796E+03, 0.83141E+03, 0.85512E+03, 0.87912E+03, 0.90339E+03, &
         0.92795E+03, 0.95279E+03, 0.97792E+03, 0.10033E+04, 0.10291E+04, &
         0.10551E+04, 0.10814E+04, 0.11080E+04, 0.11349E+04, 0.11621E+04, &
         0.11896E+04, 0.12174E+04, 0.12455E+04, 0.12739E+04, 0.13027E+04, &
         0.13317E+04, 0.13611E+04, 0.13908E+04, 0.14208E+04, 0.14510E+04, &
         0.14817E+04, 0.15126E+04, 0.15438E+04, 0.15754E+04, 0.16073E+04, &
         0.16395E+04, 0.16720E+04, 0.17049E+04, 0.17380E+04, 0.17715E+04, &
         0.18053E+04, 0.18394E+04, 0.18739E+04, 0.19086E+04, 0.19437E+04, &
         0.19792E+04, 0.20149E+04, 0.20510E+04, 0.20874E+04, 0.21241E+04, &
         0.21611E+04, 0.21985E+04, 0.22362E+04, 0.22742E+04, 0.23126E+04, &
         0.23513E+04, 0.23903E+04, 0.24296E+04, 0.24693E+04, 0.25093E+04, &
         0.25496E+04, 0.25903E+04, 0.26312E+04, 0.26725E+04, 0.27142E+04, &
         0.27562E+04, 0.27985E+04, 0.28411E+04, 0.28841E+04, 0.29274E+04, &
         0.29710E+04, 0.30150E+04, 0.30593E+04, 0.31039E+04, 0.31489E+04, &
         0.31942E+04, 0.32398E+04, 0.32858E+04, 0.33321E+04, 0.33787E+04, &
         0.34257E+04, 0.34730E+04, 0.35207E+04, 0.35686E+04, 0.36170E+04, &
         0.36656E+04/
!...        --        28
      DATA (QOFT(3,J),J=1,119)/ 0.23024E+02, 0.32483E+02, 0.41946E+02, &
         0.51412E+02, 0.60882E+02, 0.70356E+02, 0.79834E+02, 0.89315E+02, &
         0.98801E+02, 0.10829E+03, 0.11779E+03, 0.12729E+03, 0.13681E+03, &
         0.14634E+03, 0.15589E+03, 0.16546E+03, 0.17506E+03, 0.18470E+03, &
         0.19438E+03, 0.20411E+03, 0.21389E+03, 0.22374E+03, 0.23365E+03, &
         0.24364E+03, 0.25371E+03, 0.26386E+03, 0.27411E+03, 0.28444E+03, &
         0.29489E+03, 0.30543E+03, 0.31608E+03, 0.32685E+03, 0.33773E+03, &
         0.34873E+03, 0.35986E+03, 0.37111E+03, 0.38249E+03, 0.39400E+03, &
         0.40565E+03, 0.41742E+03, 0.42934E+03, 0.44139E+03, 0.45359E+03, &
         0.46592E+03, 0.47841E+03, 0.49103E+03, 0.50380E+03, 0.51672E+03, &
         0.52979E+03, 0.54300E+03, 0.55637E+03, 0.56989E+03, 0.58356E+03, &
         0.59738E+03, 0.61136E+03, 0.62549E+03, 0.63977E+03, 0.65421E+03, &
         0.66881E+03, 0.68357E+03, 0.69847E+03, 0.71354E+03, 0.72877E+03, &
         0.74415E+03, 0.75969E+03, 0.77540E+03, 0.79126E+03, 0.80728E+03, &
         0.82346E+03, 0.83981E+03, 0.85631E+03, 0.87297E+03, 0.88980E+03, &
         0.90679E+03, 0.92394E+03, 0.94125E+03, 0.95873E+03, 0.97636E+03, &
         0.99417E+03, 0.10121E+04, 0.10303E+04, 0.10485E+04, 0.10670E+04, &
         0.10856E+04, 0.11044E+04, 0.11234E+04, 0.11425E+04, 0.11617E+04, &
         0.11812E+04, 0.12008E+04, 0.12206E+04, 0.12405E+04, 0.12606E+04, &
         0.12809E+04, 0.13013E+04, 0.13219E+04, 0.13427E+04, 0.13636E+04, &
         0.13847E+04, 0.14060E+04, 0.14274E+04, 0.14490E+04, 0.14708E+04, &
         0.14927E+04, 0.15148E+04, 0.15371E+04, 0.15595E+04, 0.15821E+04, &
         0.16049E+04, 0.16278E+04, 0.16509E+04, 0.16742E+04, 0.16976E+04, &
         0.17212E+04, 0.17450E+04, 0.17690E+04, 0.17931E+04, 0.18174E+04, &
         0.18418E+04/
!...        --        27
      DATA (QOFT(4,J),J=1,119)/ 0.13502E+03, 0.19046E+03, 0.24593E+03, &
         0.30143E+03, 0.35694E+03, 0.41248E+03, 0.46804E+03, 0.52362E+03, &
         0.57922E+03, 0.63485E+03, 0.69052E+03, 0.74623E+03, 0.80201E+03, &
         0.85786E+03, 0.91382E+03, 0.96991E+03, 0.10262E+04, 0.10826E+04, &
         0.11393E+04, 0.11963E+04, 0.12536E+04, 0.13112E+04, 0.13692E+04, &
         0.14276E+04, 0.14865E+04, 0.15459E+04, 0.16057E+04, 0.16662E+04, &
         0.17272E+04, 0.17888E+04, 0.18510E+04, 0.19139E+04, 0.19774E+04, &
         0.20416E+04, 0.21066E+04, 0.21722E+04, 0.22386E+04, 0.23058E+04, &
         0.23737E+04, 0.24424E+04, 0.25118E+04, 0.25821E+04, 0.26532E+04, &
         0.27251E+04, 0.27978E+04, 0.28714E+04, 0.29458E+04, 0.30211E+04, &
         0.30972E+04, 0.31742E+04, 0.32520E+04, 0.33307E+04, 0.34104E+04, &
         0.34908E+04, 0.35722E+04, 0.36545E+04, 0.37376E+04, 0.38217E+04, &
         0.39066E+04, 0.39925E+04, 0.40793E+04, 0.41670E+04, 0.42556E+04, &
         0.43451E+04, 0.44355E+04, 0.45269E+04, 0.46191E+04, 0.47124E+04, &
         0.48065E+04, 0.49016E+04, 0.49976E+04, 0.50945E+04, 0.51923E+04, &
         0.52912E+04, 0.53909E+04, 0.54916E+04, 0.55932E+04, 0.56957E+04, &
         0.57993E+04, 0.59037E+04, 0.60091E+04, 0.61155E+04, 0.62228E+04, &
         0.63310E+04, 0.64402E+04, 0.65504E+04, 0.66615E+04, 0.67735E+04, &
         0.68866E+04, 0.70005E+04, 0.71154E+04, 0.72313E+04, 0.73481E+04, &
         0.74660E+04, 0.75847E+04, 0.77045E+04, 0.78251E+04, 0.79468E+04, &
         0.80694E+04, 0.81930E+04, 0.83175E+04, 0.84431E+04, 0.85695E+04, &
         0.86970E+04, 0.88254E+04, 0.89548E+04, 0.90852E+04, 0.92164E+04, &
         0.93487E+04, 0.94820E+04, 0.96163E+04, 0.97515E+04, 0.98877E+04, &
         0.10025E+05, 0.10163E+05, 0.10302E+05, 0.10442E+05, 0.10583E+05, &
         0.10725E+05/
!...        --        38
      DATA (QOFT(5,J),J=1,119)/ 0.48251E+02, 0.68086E+02, 0.87930E+02, &
         0.10778E+03, 0.12764E+03, 0.14751E+03, 0.16738E+03, 0.18727E+03, &
         0.20716E+03, 0.22706E+03, 0.24698E+03, 0.26692E+03, 0.28688E+03, &
         0.30687E+03, 0.32691E+03, 0.34701E+03, 0.36717E+03, 0.38742E+03, &
         0.40776E+03, 0.42821E+03, 0.44880E+03, 0.46951E+03, 0.49039E+03, &
         0.51142E+03, 0.53264E+03, 0.55404E+03, 0.57565E+03, 0.59747E+03, &
         0.61951E+03, 0.64179E+03, 0.66430E+03, 0.68707E+03, 0.71008E+03, &
         0.73336E+03, 0.75691E+03, 0.78073E+03, 0.80483E+03, 0.82922E+03, &
         0.85390E+03, 0.87887E+03, 0.90413E+03, 0.92970E+03, 0.95558E+03, &
         0.98176E+03, 0.10082E+04, 0.10351E+04, 0.10622E+04, 0.10896E+04, &
         0.11174E+04, 0.11455E+04, 0.11739E+04, 0.12026E+04, 0.12317E+04, &
         0.12611E+04, 0.12908E+04, 0.13209E+04, 0.13512E+04, 0.13820E+04, &
         0.14130E+04, 0.14444E+04, 0.14762E+04, 0.15082E+04, 0.15407E+04, &
         0.15734E+04, 0.16065E+04, 0.16400E+04, 0.16737E+04, 0.17079E+04, &
         0.17424E+04, 0.17772E+04, 0.18123E+04, 0.18479E+04, 0.18837E+04, &
         0.19199E+04, 0.19565E+04, 0.19934E+04, 0.20306E+04, 0.20683E+04, &
         0.21062E+04, 0.21445E+04, 0.21832E+04, 0.22222E+04, 0.22615E+04, &
         0.23013E+04, 0.23413E+04, 0.23817E+04, 0.24225E+04, 0.24636E+04, &
         0.25051E+04, 0.25470E+04, 0.25892E+04, 0.26317E+04, 0.26746E+04, &
         0.27179E+04, 0.27615E+04, 0.28054E+04, 0.28498E+04, 0.28945E+04, &
         0.29395E+04, 0.29849E+04, 0.30306E+04, 0.30768E+04, 0.31232E+04, &
         0.31701E+04, 0.32173E+04, 0.32648E+04, 0.33127E+04, 0.33610E+04, &
         0.34096E+04, 0.34586E+04, 0.35080E+04, 0.35577E+04, 0.36077E+04, &
         0.36582E+04, 0.37090E+04, 0.37601E+04, 0.38116E+04, 0.38635E+04, &
         0.39158E+04/
!...        --        37
      DATA (QOFT(6,J),J=1,119)/ 0.28263E+03, 0.39877E+03, 0.51497E+03, &
         0.63121E+03, 0.74749E+03, 0.86382E+03, 0.98020E+03, 0.10966E+04, &
         0.12131E+04, 0.13296E+04, 0.14462E+04, 0.15629E+04, 0.16797E+04, &
         0.17966E+04, 0.19138E+04, 0.20313E+04, 0.21490E+04, 0.22672E+04, &
         0.23858E+04, 0.25049E+04, 0.26247E+04, 0.27452E+04, 0.28665E+04, &
         0.29886E+04, 0.31116E+04, 0.32356E+04, 0.33606E+04, 0.34868E+04, &
         0.36141E+04, 0.37427E+04, 0.38725E+04, 0.40036E+04, 0.41361E+04, &
         0.42700E+04, 0.44054E+04, 0.45422E+04, 0.46805E+04, 0.48204E+04, &
         0.49618E+04, 0.51049E+04, 0.52495E+04, 0.53958E+04, 0.55438E+04, &
         0.56934E+04, 0.58448E+04, 0.59979E+04, 0.61527E+04, 0.63092E+04, &
         0.64675E+04, 0.66276E+04, 0.67895E+04, 0.69531E+04, 0.71187E+04, &
         0.72859E+04, 0.74551E+04, 0.76261E+04, 0.77989E+04, 0.79736E+04, &
         0.81501E+04, 0.83285E+04, 0.85088E+04, 0.86910E+04, 0.88751E+04, &
         0.90610E+04, 0.92488E+04, 0.94385E+04, 0.96302E+04, 0.98238E+04, &
         0.10019E+05, 0.10217E+05, 0.10416E+05, 0.10617E+05, 0.10820E+05, &
         0.11025E+05, 0.11233E+05, 0.11441E+05, 0.11652E+05, 0.11865E+05, &
         0.12080E+05, 0.12297E+05, 0.12516E+05, 0.12736E+05, 0.12959E+05, &
         0.13184E+05, 0.13410E+05, 0.13639E+05, 0.13869E+05, 0.14102E+05, &
         0.14336E+05, 0.14573E+05, 0.14811E+05, 0.15051E+05, 0.15294E+05, &
         0.15538E+05, 0.15784E+05, 0.16033E+05, 0.16283E+05, 0.16535E+05, &
         0.16790E+05, 0.17046E+05, 0.17304E+05, 0.17564E+05, 0.17827E+05, &
         0.18091E+05, 0.18357E+05, 0.18626E+05, 0.18896E+05, 0.19168E+05, &
         0.19443E+05, 0.19719E+05, 0.19997E+05, 0.20277E+05, 0.20560E+05, &
         0.20844E+05, 0.21130E+05, 0.21419E+05, 0.21709E+05, 0.22002E+05, &
         0.22296E+05/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_CO


!
!     *****************
      SUBROUTINE QT_CH4(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(3) :: XGJ
      REAL(DOUBLE), DIMENSION(3,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1., 2., 3./
!...      CH4
!...        --       211
!...J-P Chamapion calculations JQSRT
      data (QofT( 1,J),J=1,119)/0.54800E+02,0.91500E+02,0.13410E+03, &
         0.18180E+03,0.23410E+03,0.29070E+03,0.35140E+03,0.41600E+03, &
         0.48450E+03,0.55720E+03,0.63420E+03,0.71600E+03,0.80310E+03, &
         0.89590E+03,0.99520E+03,0.11017E+04,0.12161E+04,0.13393E+04, &
         0.14721E+04,0.16155E+04,0.17706E+04,0.19384E+04,0.21202E+04, &
         0.23172E+04,0.25307E+04,0.27624E+04,0.30137E+04,0.32864E+04, &
         0.35823E+04,0.39034E+04,0.42519E+04,0.46300E+04,0.50402E+04, &
         0.54853E+04,0.59679E+04,0.64913E+04,0.70588E+04,0.76739E+04, &
         0.83404E+04,0.90625E+04,0.98446E+04,0.10691E+05,0.11608E+05, &
         0.12600E+05,0.13674E+05,0.14835E+05,0.16090E+05,0.17447E+05, &
         0.18914E+05,0.20500E+05,0.22212E+05,0.24063E+05,0.26061E+05, &
         0.28218E+05,0.30548E+05,0.33063E+05,0.35778E+05,0.38708E+05, &
         0.41871E+05,0.45284E+05,0.48970E+05,0.52940E+05,0.57230E+05, &
         0.61860E+05,0.66860E+05,0.72250E+05,0.78070E+05,0.84350E+05, &
         0.91130E+05,0.98450E+05,0.10635E+06,0.11488E+06,0.12408E+06, &
         0.13403E+06,0.14480E+06,0.15640E+06,0.16890E+06,0.18240E+06, &
         0.19700E+06,0.21280E+06,0.22980E+06,0.24830E+06,0.26820E+06, &
         0.28970E+06,0.31290E+06,0.33800E+06,0.36520E+06,0.39450E+06, &
         0.42600E+06,0.46000E+06,0.49700E+06,0.53700E+06,0.58100E+06, &
         0.62700E+06,0.67800E+06,0.73300E+06,0.79200E+06,0.85600E+06, &
         0.92500E+06,0.10000E+07,0.10800E+07,0.11670E+07,0.12610E+07, &
         0.13620E+07,0.14720E+07,0.15910E+07,0.17190E+07,0.18600E+07, &
         0.20100E+07,0.21700E+07,0.23400E+07,0.25300E+07,0.27300E+07, &
         0.29500E+07,0.31800E+07,0.34300E+07,0.37000E+07,0.39900E+07, &
         0.42856E+07/
!...        --       311
!...McDowell formula
      DATA (QOFT(2,J),J=1,119)/ 0.10958E+03, 0.18304E+03, 0.26818E+03, &
         0.36356E+03, 0.46820E+03, 0.58141E+03, 0.70270E+03, 0.83186E+03, &
         0.96893E+03, 0.11142E+04, 0.12682E+04, 0.14316E+04, 0.16055E+04, &
         0.17909E+04, 0.19891E+04, 0.22016E+04, 0.24297E+04, 0.26752E+04, &
         0.29399E+04, 0.32255E+04, 0.35342E+04, 0.38680E+04, 0.42294E+04, &
         0.46208E+04, 0.50449E+04, 0.55046E+04, 0.60030E+04, 0.65434E+04, &
         0.71293E+04, 0.77646E+04, 0.84535E+04, 0.92004E+04, 0.10010E+05, &
         0.10888E+05, 0.11838E+05, 0.12869E+05, 0.13984E+05, 0.15193E+05, &
         0.16501E+05, 0.17916E+05, 0.19448E+05, 0.21104E+05, 0.22895E+05, &
         0.24830E+05, 0.26921E+05, 0.29180E+05, 0.31618E+05, 0.34250E+05, &
         0.37090E+05, 0.40152E+05, 0.43454E+05, 0.47012E+05, 0.50845E+05, &
         0.54973E+05, 0.59416E+05, 0.64197E+05, 0.69340E+05, 0.74870E+05, &
         0.80813E+05, 0.87198E+05, 0.94055E+05, 0.10142E+06, 0.10932E+06, &
         0.11779E+06, 0.12688E+06, 0.13662E+06, 0.14706E+06, 0.15824E+06, &
         0.17021E+06, 0.18302E+06, 0.19673E+06, 0.21139E+06, 0.22706E+06, &
         0.24381E+06, 0.26171E+06, 0.28082E+06, 0.30122E+06, 0.32299E+06, &
         0.34621E+06, 0.37097E+06, 0.39737E+06, 0.42551E+06, 0.45548E+06, &
         0.48739E+06, 0.52136E+06, 0.55752E+06, 0.59598E+06, 0.63688E+06, &
         0.68036E+06, 0.72657E+06, 0.77566E+06, 0.82780E+06, 0.88316E+06, &
         0.94191E+06, 0.10043E+07, 0.10704E+07, 0.11405E+07, 0.12148E+07, &
         0.12936E+07, 0.13770E+07, 0.14654E+07, 0.15589E+07, 0.16579E+07, &
         0.17627E+07, 0.18736E+07, 0.19908E+07, 0.21147E+07, 0.22456E+07, &
         0.23840E+07, 0.25301E+07, 0.26844E+07, 0.28474E+07, 0.30193E+07, &
         0.32007E+07, 0.33921E+07, 0.35939E+07, 0.38067E+07, 0.40310E+07, &
         0.42673E+07/
!...        --       212
      DATA (QOFT(3,J),J=1,119)/ 0.44079E+03, 0.73786E+03, 0.10822E+04, &
         0.14679E+04, 0.18912E+04, 0.23493E+04, 0.28400E+04, 0.33627E+04, &
         0.39174E+04, 0.45053E+04, 0.51285E+04, 0.57900E+04, 0.64938E+04, &
         0.72443E+04, 0.80465E+04, 0.89064E+04, 0.98299E+04, 0.10824E+05, &
         0.11895E+05, 0.13051E+05, 0.14300E+05, 0.15652E+05, 0.17115E+05, &
         0.18699E+05, 0.20416E+05, 0.22277E+05, 0.24294E+05, 0.26481E+05, &
         0.28853E+05, 0.31425E+05, 0.34214E+05, 0.37237E+05, 0.40515E+05, &
         0.44067E+05, 0.47916E+05, 0.52087E+05, 0.56604E+05, 0.61495E+05, &
         0.66790E+05, 0.72521E+05, 0.78720E+05, 0.85426E+05, 0.92675E+05, &
         0.10051E+06, 0.10898E+06, 0.11812E+06, 0.12799E+06, 0.13865E+06, &
         0.15014E+06, 0.16254E+06, 0.17591E+06, 0.19031E+06, 0.20583E+06, &
         0.22254E+06, 0.24053E+06, 0.25989E+06, 0.28071E+06, 0.30310E+06, &
         0.32716E+06, 0.35301E+06, 0.38077E+06, 0.41058E+06, 0.44257E+06, &
         0.47688E+06, 0.51367E+06, 0.55311E+06, 0.59537E+06, 0.64064E+06, &
         0.68910E+06, 0.74098E+06, 0.79647E+06, 0.85583E+06, 0.91928E+06, &
         0.98710E+06, 0.10595E+07, 0.11369E+07, 0.12195E+07, 0.13076E+07, &
         0.14017E+07, 0.15019E+07, 0.16088E+07, 0.17227E+07, 0.18441E+07, &
         0.19733E+07, 0.21108E+07, 0.22572E+07, 0.24129E+07, 0.25785E+07, &
         0.27545E+07, 0.29416E+07, 0.31404E+07, 0.33515E+07, 0.35756E+07, &
         0.38135E+07, 0.40659E+07, 0.43336E+07, 0.46174E+07, 0.49183E+07, &
         0.52372E+07, 0.55750E+07, 0.59327E+07, 0.63114E+07, 0.67123E+07, &
         0.71365E+07, 0.75852E+07, 0.80596E+07, 0.85612E+07, 0.90914E+07, &
         0.96515E+07, 0.10243E+08, 0.10868E+08, 0.11527E+08, 0.12224E+08, &
         0.12958E+08, 0.13733E+08, 0.14550E+08, 0.15411E+08, 0.16319E+08, &
         0.17275E+08/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_CH4


!
!     *****************
      SUBROUTINE QT_O2(T, ISO, GSI, QT)          ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(3) :: XGJ
      REAL(DOUBLE), DIMENSION(3,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1., 1., 6./
!...       O2
!...        --        66
      DATA (QOFT(1,J),J=1,119)/ 0.44334E+02, 0.62460E+02, 0.80596E+02, &
         0.98738E+02, 0.11688E+03, 0.13503E+03, 0.15319E+03, 0.17136E+03, &
         0.18954E+03, 0.20775E+03, 0.22600E+03, 0.24431E+03, 0.26270E+03, &
         0.28119E+03, 0.29981E+03, 0.31857E+03, 0.33750E+03, 0.35662E+03, &
         0.37594E+03, 0.39550E+03, 0.41529E+03, 0.43535E+03, 0.45568E+03, &
         0.47630E+03, 0.49722E+03, 0.51844E+03, 0.53998E+03, 0.56185E+03, &
         0.58406E+03, 0.60660E+03, 0.62949E+03, 0.65274E+03, 0.67635E+03, &
         0.70031E+03, 0.72465E+03, 0.74936E+03, 0.77444E+03, 0.79990E+03, &
         0.82574E+03, 0.85197E+03, 0.87858E+03, 0.90558E+03, 0.93297E+03, &
         0.96076E+03, 0.98895E+03, 0.10175E+04, 0.10465E+04, 0.10759E+04, &
         0.11057E+04, 0.11359E+04, 0.11665E+04, 0.11976E+04, 0.12290E+04, &
         0.12609E+04, 0.12931E+04, 0.13258E+04, 0.13590E+04, 0.13925E+04, &
         0.14265E+04, 0.14609E+04, 0.14958E+04, 0.15311E+04, 0.15669E+04, &
         0.16031E+04, 0.16397E+04, 0.16768E+04, 0.17144E+04, 0.17524E+04, &
         0.17909E+04, 0.18298E+04, 0.18692E+04, 0.19091E+04, 0.19495E+04, &
         0.19904E+04, 0.20318E+04, 0.20736E+04, 0.21160E+04, 0.21588E+04, &
         0.22022E+04, 0.22461E+04, 0.22905E+04, 0.23354E+04, 0.23809E+04, &
         0.24268E+04, 0.24734E+04, 0.25204E+04, 0.25680E+04, 0.26162E+04, &
         0.26649E+04, 0.27142E+04, 0.27641E+04, 0.28145E+04, 0.28655E+04, &
         0.29171E+04, 0.29693E+04, 0.30221E+04, 0.30755E+04, 0.31295E+04, &
         0.31841E+04, 0.32393E+04, 0.32951E+04, 0.33516E+04, 0.34087E+04, &
         0.34665E+04, 0.35249E+04, 0.35839E+04, 0.36436E+04, 0.37040E+04, &
         0.37650E+04, 0.38267E+04, 0.38891E+04, 0.39522E+04, 0.40159E+04, &
         0.40804E+04, 0.41455E+04, 0.42114E+04, 0.42780E+04, 0.43452E+04, &
         0.44132E+04/
!...        --        68
      DATA (QOFT(2,J),J=1,119)/ 0.89206E+02, 0.12759E+03, 0.16600E+03, &
         0.20442E+03, 0.24285E+03, 0.28128E+03, 0.31973E+03, 0.35821E+03, &
         0.39672E+03, 0.43530E+03, 0.47398E+03, 0.51281E+03, 0.55183E+03, &
         0.59108E+03, 0.63062E+03, 0.67051E+03, 0.71078E+03, 0.75148E+03, &
         0.79265E+03, 0.83435E+03, 0.87659E+03, 0.91941E+03, 0.96285E+03, &
         0.10069E+04, 0.10517E+04, 0.10971E+04, 0.11432E+04, 0.11901E+04, &
         0.12377E+04, 0.12861E+04, 0.13352E+04, 0.13851E+04, 0.14358E+04, &
         0.14872E+04, 0.15395E+04, 0.15926E+04, 0.16466E+04, 0.17013E+04, &
         0.17569E+04, 0.18134E+04, 0.18706E+04, 0.19288E+04, 0.19877E+04, &
         0.20476E+04, 0.21083E+04, 0.21698E+04, 0.22323E+04, 0.22956E+04, &
         0.23598E+04, 0.24248E+04, 0.24908E+04, 0.25576E+04, 0.26253E+04, &
         0.26940E+04, 0.27635E+04, 0.28339E+04, 0.29052E+04, 0.29775E+04, &
         0.30506E+04, 0.31247E+04, 0.31997E+04, 0.32756E+04, 0.33524E+04, &
         0.34302E+04, 0.35089E+04, 0.35885E+04, 0.36691E+04, 0.37506E+04, &
         0.38331E+04, 0.39166E+04, 0.40010E+04, 0.40864E+04, 0.41727E+04, &
         0.42601E+04, 0.43484E+04, 0.44377E+04, 0.45280E+04, 0.46193E+04, &
         0.47116E+04, 0.48049E+04, 0.48992E+04, 0.49946E+04, 0.50909E+04, &
         0.51883E+04, 0.52868E+04, 0.53863E+04, 0.54868E+04, 0.55884E+04, &
         0.56911E+04, 0.57949E+04, 0.58997E+04, 0.60056E+04, 0.61126E+04, &
         0.62207E+04, 0.63298E+04, 0.64401E+04, 0.65516E+04, 0.66641E+04, &
         0.67778E+04, 0.68926E+04, 0.70085E+04, 0.71256E+04, 0.72439E+04, &
         0.73633E+04, 0.74839E+04, 0.76056E+04, 0.77286E+04, 0.78527E+04, &
         0.79781E+04, 0.81046E+04, 0.82324E+04, 0.83613E+04, 0.84915E+04, &
         0.86229E+04, 0.87556E+04, 0.88895E+04, 0.90247E+04, 0.91611E+04, &
         0.92988E+04/
!...        --        67
      DATA (QOFT(3,J),J=1,119)/ 0.52071E+03, 0.74484E+03, 0.96908E+03, &
         0.11934E+04, 0.14177E+04, 0.16422E+04, 0.18667E+04, 0.20913E+04, &
         0.23161E+04, 0.25413E+04, 0.27671E+04, 0.29936E+04, 0.32212E+04, &
         0.34501E+04, 0.36806E+04, 0.39130E+04, 0.41476E+04, 0.43846E+04, &
         0.46242E+04, 0.48668E+04, 0.51125E+04, 0.53615E+04, 0.56140E+04, &
         0.58701E+04, 0.61300E+04, 0.63938E+04, 0.66617E+04, 0.69337E+04, &
         0.72099E+04, 0.74904E+04, 0.77754E+04, 0.80647E+04, 0.83586E+04, &
         0.86571E+04, 0.89602E+04, 0.92680E+04, 0.95805E+04, 0.98977E+04, &
         0.10220E+05, 0.10547E+05, 0.10878E+05, 0.11215E+05, 0.11556E+05, &
         0.11903E+05, 0.12254E+05, 0.12611E+05, 0.12972E+05, 0.13338E+05, &
         0.13710E+05, 0.14086E+05, 0.14468E+05, 0.14855E+05, 0.15247E+05, &
         0.15644E+05, 0.16046E+05, 0.16453E+05, 0.16866E+05, 0.17283E+05, &
         0.17706E+05, 0.18135E+05, 0.18568E+05, 0.19007E+05, 0.19452E+05, &
         0.19901E+05, 0.20356E+05, 0.20817E+05, 0.21283E+05, 0.21754E+05, &
         0.22231E+05, 0.22713E+05, 0.23201E+05, 0.23695E+05, 0.24194E+05, &
         0.24699E+05, 0.25209E+05, 0.25725E+05, 0.26247E+05, 0.26775E+05, &
         0.27308E+05, 0.27847E+05, 0.28393E+05, 0.28944E+05, 0.29500E+05, &
         0.30063E+05, 0.30632E+05, 0.31207E+05, 0.31788E+05, 0.32375E+05, &
         0.32968E+05, 0.33568E+05, 0.34173E+05, 0.34785E+05, 0.35403E+05, &
         0.36028E+05, 0.36659E+05, 0.37296E+05, 0.37939E+05, 0.38590E+05, &
         0.39246E+05, 0.39909E+05, 0.40579E+05, 0.41256E+05, 0.41939E+05, &
         0.42629E+05, 0.43325E+05, 0.44029E+05, 0.44739E+05, 0.45456E+05, &
         0.46180E+05, 0.46911E+05, 0.47649E+05, 0.48394E+05, 0.49146E+05, &
         0.49905E+05, 0.50671E+05, 0.51445E+05, 0.52226E+05, 0.53014E+05, &
         0.53809E+05/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_O2


!
!     *****************
      SUBROUTINE QT_NO(T, ISO, GSI, QT)          ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(3) :: XGJ
      REAL(DOUBLE), DIMENSION(3,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 3., 2., 3./
!...       NO
!...        --        46
      DATA (QOFT(1,J),J=1,119)/ 0.15840E+03, 0.23971E+03, 0.33080E+03, &
         0.42907E+03, 0.53251E+03, 0.63972E+03, 0.74975E+03, 0.86195E+03, &
         0.97582E+03, 0.10911E+04, 0.12074E+04, 0.13248E+04, 0.14430E+04, &
         0.15621E+04, 0.16820E+04, 0.18027E+04, 0.19243E+04, 0.20468E+04, &
         0.21703E+04, 0.22948E+04, 0.24204E+04, 0.25472E+04, 0.26753E+04, &
         0.28046E+04, 0.29354E+04, 0.30676E+04, 0.32013E+04, 0.33365E+04, &
         0.34734E+04, 0.36120E+04, 0.37522E+04, 0.38942E+04, 0.40379E+04, &
         0.41835E+04, 0.43310E+04, 0.44803E+04, 0.46316E+04, 0.47849E+04, &
         0.49400E+04, 0.50972E+04, 0.52564E+04, 0.54176E+04, 0.55809E+04, &
         0.57462E+04, 0.59137E+04, 0.60832E+04, 0.62548E+04, 0.64286E+04, &
         0.66045E+04, 0.67825E+04, 0.69628E+04, 0.71451E+04, 0.73297E+04, &
         0.75164E+04, 0.77053E+04, 0.78964E+04, 0.80897E+04, 0.82853E+04, &
         0.84830E+04, 0.86830E+04, 0.88852E+04, 0.90896E+04, 0.92963E+04, &
         0.95052E+04, 0.97164E+04, 0.99297E+04, 0.10145E+05, 0.10363E+05, &
         0.10583E+05, 0.10806E+05, 0.11031E+05, 0.11258E+05, 0.11487E+05, &
         0.11718E+05, 0.11952E+05, 0.12188E+05, 0.12426E+05, 0.12667E+05, &
         0.12910E+05, 0.13155E+05, 0.13403E+05, 0.13652E+05, 0.13905E+05, &
         0.14159E+05, 0.14416E+05, 0.14675E+05, 0.14936E+05, 0.15199E+05, &
         0.15465E+05, 0.15733E+05, 0.16004E+05, 0.16277E+05, 0.16552E+05, &
         0.16829E+05, 0.17109E+05, 0.17391E+05, 0.17675E+05, 0.17962E+05, &
         0.18251E+05, 0.18542E+05, 0.18836E+05, 0.19131E+05, 0.19430E+05, &
         0.19730E+05, 0.20033E+05, 0.20338E+05, 0.20646E+05, 0.20955E+05, &
         0.21268E+05, 0.21582E+05, 0.21899E+05, 0.22218E+05, 0.22539E+05, &
         0.22863E+05, 0.23189E+05, 0.23518E+05, 0.23848E+05, 0.24181E+05, &
         0.24517E+05/
!...        --        56
      DATA (QOFT(2,J),J=1,119)/ 0.10942E+03, 0.16560E+03, 0.22856E+03, &
         0.29647E+03, 0.36795E+03, 0.44204E+03, 0.51808E+03, 0.59561E+03, &
         0.67432E+03, 0.75396E+03, 0.83439E+03, 0.91551E+03, 0.99725E+03, &
         0.10796E+04, 0.11625E+04, 0.12460E+04, 0.13302E+04, 0.14150E+04, &
         0.15005E+04, 0.15868E+04, 0.16739E+04, 0.17618E+04, 0.18506E+04, &
         0.19404E+04, 0.20311E+04, 0.21229E+04, 0.22158E+04, 0.23098E+04, &
         0.24050E+04, 0.25013E+04, 0.25989E+04, 0.26976E+04, 0.27977E+04, &
         0.28991E+04, 0.30018E+04, 0.31058E+04, 0.32112E+04, 0.33180E+04, &
         0.34262E+04, 0.35358E+04, 0.36468E+04, 0.37593E+04, 0.38732E+04, &
         0.39885E+04, 0.41054E+04, 0.42237E+04, 0.43436E+04, 0.44649E+04, &
         0.45877E+04, 0.47121E+04, 0.48379E+04, 0.49654E+04, 0.50943E+04, &
         0.52248E+04, 0.53568E+04, 0.54904E+04, 0.56255E+04, 0.57622E+04, &
         0.59004E+04, 0.60403E+04, 0.61816E+04, 0.63246E+04, 0.64692E+04, &
         0.66152E+04, 0.67630E+04, 0.69123E+04, 0.70631E+04, 0.72156E+04, &
         0.73696E+04, 0.75253E+04, 0.76825E+04, 0.78414E+04, 0.80018E+04, &
         0.81638E+04, 0.83275E+04, 0.84927E+04, 0.86596E+04, 0.88280E+04, &
         0.89981E+04, 0.91698E+04, 0.93430E+04, 0.95180E+04, 0.96945E+04, &
         0.98726E+04, 0.10052E+05, 0.10234E+05, 0.10417E+05, 0.10601E+05, &
         0.10788E+05, 0.10975E+05, 0.11165E+05, 0.11356E+05, 0.11549E+05, &
         0.11743E+05, 0.11939E+05, 0.12137E+05, 0.12336E+05, 0.12537E+05, &
         0.12739E+05, 0.12943E+05, 0.13149E+05, 0.13356E+05, 0.13565E+05, &
         0.13776E+05, 0.13988E+05, 0.14202E+05, 0.14418E+05, 0.14635E+05, &
         0.14853E+05, 0.15074E+05, 0.15296E+05, 0.15520E+05, 0.15745E+05, &
         0.15972E+05, 0.16200E+05, 0.16431E+05, 0.16663E+05, 0.16896E+05, &
         0.17131E+05/
!...        --        48
      DATA (QOFT(3,J),J=1,119)/ 0.16695E+03, 0.25269E+03, 0.34876E+03, &
         0.45239E+03, 0.56148E+03, 0.67455E+03, 0.79059E+03, 0.90891E+03, &
         0.10290E+04, 0.11506E+04, 0.12733E+04, 0.13971E+04, 0.15219E+04, &
         0.16476E+04, 0.17742E+04, 0.19017E+04, 0.20302E+04, 0.21598E+04, &
         0.22904E+04, 0.24223E+04, 0.25553E+04, 0.26897E+04, 0.28255E+04, &
         0.29628E+04, 0.31016E+04, 0.32420E+04, 0.33842E+04, 0.35280E+04, &
         0.36736E+04, 0.38211E+04, 0.39704E+04, 0.41217E+04, 0.42750E+04, &
         0.44302E+04, 0.45876E+04, 0.47469E+04, 0.49084E+04, 0.50720E+04, &
         0.52378E+04, 0.54058E+04, 0.55759E+04, 0.57483E+04, 0.59230E+04, &
         0.60999E+04, 0.62791E+04, 0.64605E+04, 0.66443E+04, 0.68304E+04, &
         0.70187E+04, 0.72095E+04, 0.74026E+04, 0.75980E+04, 0.77958E+04, &
         0.79960E+04, 0.81986E+04, 0.84036E+04, 0.86109E+04, 0.88207E+04, &
         0.90328E+04, 0.92474E+04, 0.94644E+04, 0.96839E+04, 0.99057E+04, &
         0.10130E+05, 0.10357E+05, 0.10586E+05, 0.10817E+05, 0.11052E+05, &
         0.11288E+05, 0.11527E+05, 0.11768E+05, 0.12012E+05, 0.12259E+05, &
         0.12507E+05, 0.12759E+05, 0.13012E+05, 0.13269E+05, 0.13527E+05, &
         0.13788E+05, 0.14052E+05, 0.14318E+05, 0.14587E+05, 0.14858E+05, &
         0.15131E+05, 0.15408E+05, 0.15686E+05, 0.15967E+05, 0.16251E+05, &
         0.16537E+05, 0.16825E+05, 0.17116E+05, 0.17410E+05, 0.17706E+05, &
         0.18004E+05, 0.18305E+05, 0.18609E+05, 0.18915E+05, 0.19224E+05, &
         0.19535E+05, 0.19848E+05, 0.20164E+05, 0.20483E+05, 0.20804E+05, &
         0.21127E+05, 0.21453E+05, 0.21782E+05, 0.22113E+05, 0.22447E+05, &
         0.22783E+05, 0.23122E+05, 0.23463E+05, 0.23807E+05, 0.24153E+05, &
         0.24502E+05, 0.24853E+05, 0.25207E+05, 0.25563E+05, 0.25922E+05, &
         0.26283E+05/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_NO


!
!     *****************
      SUBROUTINE QT_SO2(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(2) :: XGJ
      REAL(DOUBLE), DIMENSION(2,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1., 1./
!...      SO2
!...        --       626
      DATA (QOFT(1,J),J=1,119)/ 0.52899E+03, 0.89171E+03, 0.13139E+04, &
         0.17915E+04, 0.23246E+04, 0.29155E+04, 0.35675E+04, 0.42848E+04, &
         0.50723E+04, 0.59352E+04, 0.68794E+04, 0.79109E+04, 0.90366E+04, &
         0.10264E+05, 0.11599E+05, 0.13052E+05, 0.14629E+05, 0.16340E+05, &
         0.18193E+05, 0.20199E+05, 0.22366E+05, 0.24704E+05, 0.27225E+05, &
         0.29938E+05, 0.32855E+05, 0.35987E+05, 0.39346E+05, 0.42944E+05, &
         0.46794E+05, 0.50909E+05, 0.55302E+05, 0.59986E+05, 0.64977E+05, &
         0.70288E+05, 0.75934E+05, 0.81931E+05, 0.88294E+05, 0.95040E+05, &
         0.10219E+06, 0.10975E+06, 0.11774E+06, 0.12619E+06, 0.13511E+06, &
         0.14452E+06, 0.15443E+06, 0.16487E+06, 0.17586E+06, 0.18742E+06, &
         0.19957E+06, 0.21234E+06, 0.22573E+06, 0.23978E+06, 0.25451E+06, &
         0.26995E+06, 0.28611E+06, 0.30302E+06, 0.32071E+06, 0.33920E+06, &
         0.35852E+06, 0.37869E+06, 0.39974E+06, 0.42171E+06, 0.44461E+06, &
         0.46848E+06, 0.49334E+06, 0.51922E+06, 0.54617E+06, 0.57419E+06, &
         0.60334E+06, 0.63363E+06, 0.66511E+06, 0.69780E+06, 0.73174E+06, &
         0.76696E+06, 0.80349E+06, 0.84138E+06, 0.88066E+06, 0.92136E+06, &
         0.96352E+06, 0.10072E+07, 0.10524E+07, 0.10992E+07, 0.11475E+07, &
         0.11976E+07, 0.12493E+07, 0.13028E+07, 0.13580E+07, 0.14151E+07, &
         0.14741E+07, 0.15349E+07, 0.15977E+07, 0.16625E+07, 0.17293E+07, &
         0.17982E+07, 0.18693E+07, 0.19425E+07, 0.20180E+07, 0.20958E+07, &
         0.21758E+07, 0.22583E+07, 0.23432E+07, 0.24305E+07, 0.25204E+07, &
         0.26129E+07, 0.27080E+07, 0.28058E+07, 0.29064E+07, 0.30097E+07, &
         0.31159E+07, 0.32250E+07, 0.33371E+07, 0.34522E+07, 0.35705E+07, &
         0.36918E+07, 0.38164E+07, 0.39442E+07, 0.40754E+07, 0.42099E+07, &
         0.43479E+07/
!...        --       646
      DATA (QOFT(2,J),J=1,119)/ 0.53140E+03, 0.89578E+03, 0.13199E+04, &
         0.17997E+04, 0.23353E+04, 0.29288E+04, 0.35837E+04, 0.43043E+04, &
         0.50953E+04, 0.59621E+04, 0.69104E+04, 0.79465E+04, 0.90772E+04, &
         0.10310E+05, 0.11651E+05, 0.13110E+05, 0.14694E+05, 0.16413E+05, &
         0.18274E+05, 0.20289E+05, 0.22465E+05, 0.24814E+05, 0.27345E+05, &
         0.30070E+05, 0.33000E+05, 0.36145E+05, 0.39519E+05, 0.43133E+05, &
         0.46999E+05, 0.51132E+05, 0.55544E+05, 0.60248E+05, 0.65260E+05, &
         0.70594E+05, 0.76264E+05, 0.82287E+05, 0.88678E+05, 0.95453E+05, &
         0.10263E+06, 0.11022E+06, 0.11825E+06, 0.12674E+06, 0.13569E+06, &
         0.14514E+06, 0.15510E+06, 0.16558E+06, 0.17662E+06, 0.18823E+06, &
         0.20043E+06, 0.21325E+06, 0.22670E+06, 0.24081E+06, 0.25561E+06, &
         0.27111E+06, 0.28733E+06, 0.30432E+06, 0.32208E+06, 0.34065E+06, &
         0.36005E+06, 0.38031E+06, 0.40145E+06, 0.42351E+06, 0.44651E+06, &
         0.47047E+06, 0.49544E+06, 0.52144E+06, 0.54849E+06, 0.57664E+06, &
         0.60591E+06, 0.63633E+06, 0.66794E+06, 0.70077E+06, 0.73485E+06, &
         0.77022E+06, 0.80691E+06, 0.84496E+06, 0.88440E+06, 0.92527E+06, &
         0.96761E+06, 0.10115E+07, 0.10568E+07, 0.11038E+07, 0.11524E+07, &
         0.12027E+07, 0.12546E+07, 0.13083E+07, 0.13638E+07, 0.14211E+07, &
         0.14803E+07, 0.15414E+07, 0.16045E+07, 0.16695E+07, 0.17366E+07, &
         0.18059E+07, 0.18772E+07, 0.19507E+07, 0.20265E+07, 0.21046E+07, &
         0.21850E+07, 0.22678E+07, 0.23531E+07, 0.24408E+07, 0.25310E+07, &
         0.26239E+07, 0.27194E+07, 0.28176E+07, 0.29186E+07, 0.30224E+07, &
         0.31290E+07, 0.32386E+07, 0.33512E+07, 0.34668E+07, 0.35855E+07, &
         0.37074E+07, 0.38324E+07, 0.39608E+07, 0.40925E+07, 0.42276E+07, &
         0.43662E+07/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_SO2


!
!     *****************
      SUBROUTINE QT_NO2(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(1) :: XGJ
      REAL(DOUBLE), DIMENSION(1,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 3./
!...      NO2
!...        --       646
      DATA (QOFT(1,J),J=1,119)/ 0.12046E+04, 0.20297E+04, 0.29875E+04, &
         0.40626E+04, 0.52463E+04, 0.65350E+04, 0.79286E+04, 0.94298E+04, &
         0.11043E+05, 0.12776E+05, 0.14634E+05, 0.16627E+05, 0.18765E+05, &
         0.21056E+05, 0.23511E+05, 0.26143E+05, 0.28961E+05, 0.31979E+05, &
         0.35209E+05, 0.38663E+05, 0.42355E+05, 0.46300E+05, 0.50510E+05, &
         0.55001E+05, 0.59787E+05, 0.64884E+05, 0.70308E+05, 0.76075E+05, &
         0.82201E+05, 0.88704E+05, 0.95602E+05, 0.10291E+06, 0.11065E+06, &
         0.11884E+06, 0.12750E+06, 0.13665E+06, 0.14631E+06, 0.15650E+06, &
         0.16724E+06, 0.17856E+06, 0.19047E+06, 0.20301E+06, 0.21618E+06, &
         0.23002E+06, 0.24456E+06, 0.25981E+06, 0.27580E+06, 0.29256E+06, &
         0.31012E+06, 0.32850E+06, 0.34773E+06, 0.36784E+06, 0.38886E+06, &
         0.41082E+06, 0.43374E+06, 0.45766E+06, 0.48262E+06, 0.50863E+06, &
         0.53574E+06, 0.56398E+06, 0.59339E+06, 0.62398E+06, 0.65581E+06, &
         0.68891E+06, 0.72331E+06, 0.75905E+06, 0.79617E+06, 0.83470E+06, &
         0.87469E+06, 0.91617E+06, 0.95919E+06, 0.10038E+07, 0.10500E+07, &
         0.10979E+07, 0.11474E+07, 0.11988E+07, 0.12519E+07, 0.13068E+07, &
         0.13636E+07, 0.14224E+07, 0.14831E+07, 0.15459E+07, 0.16107E+07, &
         0.16776E+07, 0.17467E+07, 0.18180E+07, 0.18916E+07, 0.19675E+07, &
         0.20458E+07, 0.21265E+07, 0.22097E+07, 0.22954E+07, 0.23837E+07, &
         0.24747E+07, 0.25684E+07, 0.26648E+07, 0.27641E+07, 0.28662E+07, &
         0.29713E+07, 0.30794E+07, 0.31905E+07, 0.33048E+07, 0.34223E+07, &
         0.35430E+07, 0.36670E+07, 0.37944E+07, 0.39253E+07, 0.40597E+07, &
         0.41976E+07, 0.43393E+07, 0.44846E+07, 0.46337E+07, 0.47867E+07, &
         0.49437E+07, 0.51046E+07, 0.52696E+07, 0.54388E+07, 0.56122E+07, &
         0.57900E+07/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_NO2


!
!     *****************
      SUBROUTINE QT_NH3(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(2) :: XGJ
      REAL(DOUBLE), DIMENSION(2,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 3., 2./
!...      NH3
!...        --      4111
      DATA (QOFT(1,J),J=1,119)/ 0.16013E+03, 0.26692E+03, 0.39067E+03, &
         0.52933E+03, 0.68153E+03, 0.84641E+03, 0.10234E+04, 0.12125E+04, &
         0.14136E+04, 0.16272E+04, 0.18537E+04, 0.20937E+04, 0.23481E+04, &
         0.26177E+04, 0.29035E+04, 0.32065E+04, 0.35279E+04, 0.38688E+04, &
         0.42304E+04, 0.46141E+04, 0.50212E+04, 0.54531E+04, 0.59114E+04, &
         0.63976E+04, 0.69133E+04, 0.74602E+04, 0.80401E+04, 0.86549E+04, &
         0.93066E+04, 0.99971E+04, 0.10729E+05, 0.11504E+05, 0.12324E+05, &
         0.13193E+05, 0.14112E+05, 0.15085E+05, 0.16114E+05, 0.17201E+05, &
         0.18352E+05, 0.19567E+05, 0.20851E+05, 0.22208E+05, 0.23640E+05, &
         0.25152E+05, 0.26747E+05, 0.28430E+05, 0.30205E+05, 0.32077E+05, &
         0.34050E+05, 0.36128E+05, 0.38317E+05, 0.40623E+05, 0.43050E+05, &
         0.45605E+05, 0.48292E+05, 0.51119E+05, 0.54091E+05, 0.57215E+05, &
         0.60498E+05, 0.63947E+05, 0.67569E+05, 0.71372E+05, 0.75364E+05, &
         0.79552E+05, 0.83946E+05, 0.88553E+05, 0.93384E+05, 0.98447E+05, &
         0.10375E+06, 0.10931E+06, 0.11513E+06, 0.12122E+06, 0.12760E+06, &
         0.13427E+06, 0.14125E+06, 0.14855E+06, 0.15619E+06, 0.16417E+06, &
         0.17250E+06, 0.18121E+06, 0.19031E+06, 0.19981E+06, 0.20973E+06, &
         0.22008E+06, 0.23088E+06, 0.24215E+06, 0.25390E+06, 0.26615E+06, &
         0.27892E+06, 0.29223E+06, 0.30610E+06, 0.32055E+06, 0.33559E+06, &
         0.35125E+06, 0.36756E+06, 0.38453E+06, 0.40219E+06, 0.42056E+06, &
         0.43967E+06, 0.45953E+06, 0.48019E+06, 0.50165E+06, 0.52396E+06, &
         0.54714E+06, 0.57122E+06, 0.59622E+06, 0.62218E+06, 0.64913E+06, &
         0.67710E+06, 0.70613E+06, 0.73624E+06, 0.76748E+06, 0.79988E+06, &
         0.83347E+06, 0.86829E+06, 0.90439E+06, 0.94180E+06, 0.98056E+06, &
         0.10207E+07/
!...        --      5111
      DATA (QOFT(2,J),J=1,119)/ 0.10697E+03, 0.17832E+03, 0.26100E+03, &
         0.35364E+03, 0.45533E+03, 0.56549E+03, 0.68377E+03, 0.81007E+03, &
         0.94447E+03, 0.10872E+04, 0.12385E+04, 0.13988E+04, 0.15688E+04, &
         0.17490E+04, 0.19399E+04, 0.21424E+04, 0.23571E+04, 0.25848E+04, &
         0.28264E+04, 0.30828E+04, 0.33548E+04, 0.36434E+04, 0.39496E+04, &
         0.42745E+04, 0.46190E+04, 0.49845E+04, 0.53720E+04, 0.57828E+04, &
         0.62182E+04, 0.66796E+04, 0.71684E+04, 0.76862E+04, 0.82344E+04, &
         0.88149E+04, 0.94292E+04, 0.10079E+05, 0.10767E+05, 0.11494E+05, &
         0.12262E+05, 0.13074E+05, 0.13932E+05, 0.14839E+05, 0.15796E+05, &
         0.16806E+05, 0.17872E+05, 0.18997E+05, 0.20183E+05, 0.21434E+05, &
         0.22752E+05, 0.24141E+05, 0.25604E+05, 0.27145E+05, 0.28767E+05, &
         0.30475E+05, 0.32271E+05, 0.34160E+05, 0.36146E+05, 0.38234E+05, &
         0.40428E+05, 0.42733E+05, 0.45154E+05, 0.47696E+05, 0.50364E+05, &
         0.53163E+05, 0.56100E+05, 0.59180E+05, 0.62408E+05, 0.65792E+05, &
         0.69339E+05, 0.73053E+05, 0.76943E+05, 0.81016E+05, 0.85279E+05, &
         0.89740E+05, 0.94406E+05, 0.99287E+05, 0.10439E+06, 0.10972E+06, &
         0.11530E+06, 0.12112E+06, 0.12720E+06, 0.13355E+06, 0.14018E+06, &
         0.14711E+06, 0.15433E+06, 0.16186E+06, 0.16971E+06, 0.17791E+06, &
         0.18645E+06, 0.19534E+06, 0.20462E+06, 0.21428E+06, 0.22434E+06, &
         0.23481E+06, 0.24572E+06, 0.25706E+06, 0.26887E+06, 0.28116E+06, &
         0.29393E+06, 0.30722E+06, 0.32103E+06, 0.33539E+06, 0.35031E+06, &
         0.36581E+06, 0.38191E+06, 0.39864E+06, 0.41600E+06, 0.43403E+06, &
         0.45274E+06, 0.47215E+06, 0.49230E+06, 0.51319E+06, 0.53487E+06, &
         0.55734E+06, 0.58064E+06, 0.60478E+06, 0.62981E+06, 0.65574E+06, &
         0.68260E+06/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_NH3


!
!     *****************
      SUBROUTINE QT_HNO3(T, ISO, GSI, QT)        ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(1) :: XGJ
      REAL(DOUBLE), DIMENSION(1,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 6./
!...     HNO3
!...        --       146
      DATA (QOFT(1,J),J=1,119)/ 0.15010E+05, 0.25316E+05, 0.37374E+05, &
         0.51216E+05, 0.67105E+05, 0.85473E+05, 0.10688E+06, 0.13201E+06, &
         0.16165E+06, 0.19671E+06, 0.23825E+06, 0.28749E+06, 0.34583E+06, &
         0.41490E+06, 0.49657E+06, 0.59302E+06, 0.70673E+06, 0.84054E+06, &
         0.99775E+06, 0.11821E+07, 0.13978E+07, 0.16498E+07, 0.19436E+07, &
         0.22855E+07, 0.26825E+07, 0.31428E+07, 0.36753E+07, 0.42903E+07, &
         0.49993E+07, 0.58151E+07, 0.67523E+07, 0.78269E+07, 0.90572E+07, &
         0.10463E+08, 0.12067E+08, 0.13895E+08, 0.15973E+08, 0.18333E+08, &
         0.21009E+08, 0.24039E+08, 0.27464E+08, 0.31331E+08, 0.35690E+08, &
         0.40597E+08, 0.46115E+08, 0.52310E+08, 0.59257E+08, 0.67037E+08, &
         0.75739E+08, 0.85461E+08, 0.96310E+08, 0.10840E+09, 0.12186E+09, &
         0.13683E+09, 0.15346E+09, 0.17191E+09, 0.19236E+09, 0.21501E+09, &
         0.24006E+09, 0.26774E+09, 0.29830E+09, 0.33200E+09, 0.36914E+09, &
         0.41002E+09, 0.45498E+09, 0.50438E+09, 0.55862E+09, 0.61812E+09, &
         0.68332E+09, 0.75473E+09, 0.83286E+09, 0.91828E+09, 0.10116E+10, &
         0.11134E+10, 0.12245E+10, 0.13456E+10, 0.14775E+10, 0.16210E+10, &
         0.17771E+10, 0.19467E+10, 0.21309E+10, 0.23309E+10, 0.25477E+10, &
         0.27827E+10, 0.30372E+10, 0.33127E+10, 0.36107E+10, 0.39329E+10, &
         0.42809E+10, 0.46567E+10, 0.50623E+10, 0.54997E+10, 0.59711E+10, &
         0.64789E+10, 0.70257E+10, 0.76140E+10, 0.82468E+10, 0.89269E+10, &
         0.96575E+10, 0.10442E+11, 0.11284E+11, 0.12187E+11, 0.13155E+11, &
         0.14193E+11, 0.15304E+11, 0.16494E+11, 0.17767E+11, 0.19129E+11, &
         0.20585E+11, 0.22140E+11, 0.23802E+11, 0.25576E+11, 0.27469E+11, &
         0.29489E+11, 0.31642E+11, 0.33937E+11, 0.36382E+11, 0.38985E+11, &
         0.41757E+11/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_HNO3


!
!     *****************
      SUBROUTINE QT_OH(T, ISO, GSI, QT)          ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(3) :: XGJ
      REAL(DOUBLE), DIMENSION(3,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 2., 2., 3./
!...       OH
!...        --        61
      DATA (QOFT(1,J),J=1,119)/ 0.20066E+02, 0.24774E+02, 0.30309E+02, &
         0.36357E+02, 0.42745E+02, 0.49371E+02, 0.56168E+02, 0.63093E+02, &
         0.70116E+02, 0.77217E+02, 0.84380E+02, 0.91594E+02, 0.98850E+02, &
         0.10614E+03, 0.11346E+03, 0.12081E+03, 0.12818E+03, 0.13557E+03, &
         0.14298E+03, 0.15041E+03, 0.15785E+03, 0.16531E+03, 0.17278E+03, &
         0.18027E+03, 0.18778E+03, 0.19530E+03, 0.20284E+03, 0.21040E+03, &
         0.21797E+03, 0.22556E+03, 0.23318E+03, 0.24082E+03, 0.24848E+03, &
         0.25617E+03, 0.26389E+03, 0.27163E+03, 0.27941E+03, 0.28721E+03, &
         0.29505E+03, 0.30292E+03, 0.31084E+03, 0.31878E+03, 0.32677E+03, &
         0.33480E+03, 0.34287E+03, 0.35099E+03, 0.35915E+03, 0.36736E+03, &
         0.37561E+03, 0.38391E+03, 0.39227E+03, 0.40067E+03, 0.40913E+03, &
         0.41764E+03, 0.42620E+03, 0.43482E+03, 0.44350E+03, 0.45223E+03, &
         0.46102E+03, 0.46987E+03, 0.47878E+03, 0.48775E+03, 0.49679E+03, &
         0.50588E+03, 0.51503E+03, 0.52425E+03, 0.53354E+03, 0.54288E+03, &
         0.55229E+03, 0.56177E+03, 0.57132E+03, 0.58092E+03, 0.59060E+03, &
         0.60035E+03, 0.61016E+03, 0.62004E+03, 0.62999E+03, 0.64001E+03, &
         0.65010E+03, 0.66025E+03, 0.67049E+03, 0.68078E+03, 0.69115E+03, &
         0.70160E+03, 0.71211E+03, 0.72269E+03, 0.73335E+03, 0.74408E+03, &
         0.75488E+03, 0.76576E+03, 0.77671E+03, 0.78773E+03, 0.79883E+03, &
         0.81000E+03, 0.82124E+03, 0.83256E+03, 0.84396E+03, 0.85542E+03, &
         0.86696E+03, 0.87858E+03, 0.89027E+03, 0.90204E+03, 0.91389E+03, &
         0.92580E+03, 0.93781E+03, 0.94988E+03, 0.96203E+03, 0.97425E+03, &
         0.98656E+03, 0.99893E+03, 0.10114E+04, 0.10239E+04, 0.10365E+04, &
         0.10492E+04, 0.10620E+04, 0.10748E+04, 0.10878E+04, 0.11007E+04, &
         0.11138E+04/
!...        --        81
      DATA (QOFT(2,J),J=1,119)/ 0.20124E+02, 0.24876E+02, 0.30457E+02, &
         0.36553E+02, 0.42991E+02, 0.49666E+02, 0.56513E+02, 0.63489E+02, &
         0.70563E+02, 0.77715E+02, 0.84929E+02, 0.92195E+02, 0.99504E+02, &
         0.10685E+03, 0.11423E+03, 0.12164E+03, 0.12907E+03, 0.13654E+03, &
         0.14403E+03, 0.15154E+03, 0.15909E+03, 0.16666E+03, 0.17427E+03, &
         0.18191E+03, 0.18959E+03, 0.19731E+03, 0.20507E+03, 0.21287E+03, &
         0.22073E+03, 0.22863E+03, 0.23658E+03, 0.24459E+03, 0.25266E+03, &
         0.26078E+03, 0.26897E+03, 0.27722E+03, 0.28554E+03, 0.29393E+03, &
         0.30238E+03, 0.31091E+03, 0.31952E+03, 0.32820E+03, 0.33696E+03, &
         0.34579E+03, 0.35471E+03, 0.36371E+03, 0.37279E+03, 0.38196E+03, &
         0.39121E+03, 0.40055E+03, 0.40998E+03, 0.41949E+03, 0.42910E+03, &
         0.43879E+03, 0.44858E+03, 0.45845E+03, 0.46843E+03, 0.47849E+03, &
         0.48865E+03, 0.49890E+03, 0.50924E+03, 0.51969E+03, 0.53022E+03, &
         0.54086E+03, 0.55159E+03, 0.56242E+03, 0.57335E+03, 0.58437E+03, &
         0.59550E+03, 0.60673E+03, 0.61805E+03, 0.62947E+03, 0.64100E+03, &
         0.65263E+03, 0.66435E+03, 0.67618E+03, 0.68811E+03, 0.70014E+03, &
         0.71228E+03, 0.72451E+03, 0.73685E+03, 0.74929E+03, 0.76184E+03, &
         0.77449E+03, 0.78724E+03, 0.80009E+03, 0.81306E+03, 0.82612E+03, &
         0.83929E+03, 0.85256E+03, 0.86594E+03, 0.87942E+03, 0.89301E+03, &
         0.90670E+03, 0.92050E+03, 0.93440E+03, 0.94841E+03, 0.96253E+03, &
         0.97675E+03, 0.99108E+03, 0.10055E+04, 0.10201E+04, 0.10347E+04, &
         0.10495E+04, 0.10643E+04, 0.10793E+04, 0.10944E+04, 0.11096E+04, &
         0.11248E+04, 0.11402E+04, 0.11558E+04, 0.11714E+04, 0.11871E+04, &
         0.12029E+04, 0.12189E+04, 0.12349E+04, 0.12511E+04, 0.12673E+04, &
         0.12837E+04/
!...        --        62
      DATA (QOFT(3,J),J=1,119)/ 0.41032E+02, 0.54704E+02, 0.70201E+02, &
         0.86985E+02, 0.10469E+03, 0.12306E+03, 0.14194E+03, 0.16119E+03, &
         0.18075E+03, 0.20054E+03, 0.22053E+03, 0.24068E+03, 0.26096E+03, &
         0.28135E+03, 0.30183E+03, 0.32241E+03, 0.34305E+03, 0.36376E+03, &
         0.38453E+03, 0.40535E+03, 0.42622E+03, 0.44714E+03, 0.46811E+03, &
         0.48913E+03, 0.51019E+03, 0.53131E+03, 0.55246E+03, 0.57368E+03, &
         0.59495E+03, 0.61627E+03, 0.63766E+03, 0.65912E+03, 0.68064E+03, &
         0.70223E+03, 0.72390E+03, 0.74565E+03, 0.76749E+03, 0.78941E+03, &
         0.81143E+03, 0.83355E+03, 0.85578E+03, 0.87810E+03, 0.90054E+03, &
         0.92310E+03, 0.94577E+03, 0.96857E+03, 0.99149E+03, 0.10145E+04, &
         0.10377E+04, 0.10611E+04, 0.10845E+04, 0.11081E+04, 0.11319E+04, &
         0.11558E+04, 0.11798E+04, 0.12040E+04, 0.12284E+04, 0.12529E+04, &
         0.12776E+04, 0.13025E+04, 0.13275E+04, 0.13527E+04, 0.13781E+04, &
         0.14036E+04, 0.14293E+04, 0.14552E+04, 0.14813E+04, 0.15076E+04, &
         0.15340E+04, 0.15606E+04, 0.15874E+04, 0.16144E+04, 0.16416E+04, &
         0.16690E+04, 0.16965E+04, 0.17243E+04, 0.17522E+04, 0.17804E+04, &
         0.18087E+04, 0.18373E+04, 0.18660E+04, 0.18949E+04, 0.19241E+04, &
         0.19534E+04, 0.19829E+04, 0.20127E+04, 0.20426E+04, 0.20727E+04, &
         0.21031E+04, 0.21336E+04, 0.21644E+04, 0.21954E+04, 0.22266E+04, &
         0.22579E+04, 0.22895E+04, 0.23213E+04, 0.23534E+04, 0.23856E+04, &
         0.24180E+04, 0.24506E+04, 0.24835E+04, 0.25166E+04, 0.25499E+04, &
         0.25834E+04, 0.26171E+04, 0.26510E+04, 0.26852E+04, 0.27195E+04, &
         0.27541E+04, 0.27889E+04, 0.28239E+04, 0.28592E+04, 0.28946E+04, &
         0.29303E+04, 0.29661E+04, 0.30023E+04, 0.30386E+04, 0.30751E+04, &
         0.31119E+04/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_OH


!
!     *****************
      SUBROUTINE QT_HF(T, ISO, GSI, QT)          ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(1) :: XGJ
      REAL(DOUBLE), DIMENSION(1,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 4./
!...       HF
!...        --        19
      DATA (QOFT(1,J),J=1,119)/ 0.95958E+01, 0.12933E+02, 0.16295E+02, &
         0.19666E+02, 0.23043E+02, 0.26425E+02, 0.29809E+02, 0.33195E+02, &
         0.36584E+02, 0.39974E+02, 0.43366E+02, 0.46759E+02, 0.50154E+02, &
         0.53550E+02, 0.56947E+02, 0.60346E+02, 0.63746E+02, 0.67148E+02, &
         0.70550E+02, 0.73955E+02, 0.77361E+02, 0.80769E+02, 0.84179E+02, &
         0.87591E+02, 0.91006E+02, 0.94424E+02, 0.97846E+02, 0.10127E+03, &
         0.10470E+03, 0.10813E+03, 0.11157E+03, 0.11502E+03, 0.11847E+03, &
         0.12193E+03, 0.12540E+03, 0.12888E+03, 0.13236E+03, 0.13586E+03, &
         0.13936E+03, 0.14288E+03, 0.14641E+03, 0.14995E+03, 0.15351E+03, &
         0.15708E+03, 0.16066E+03, 0.16426E+03, 0.16788E+03, 0.17151E+03, &
         0.17516E+03, 0.17882E+03, 0.18251E+03, 0.18621E+03, 0.18994E+03, &
         0.19368E+03, 0.19745E+03, 0.20123E+03, 0.20504E+03, 0.20887E+03, &
         0.21272E+03, 0.21659E+03, 0.22049E+03, 0.22441E+03, 0.22836E+03, &
         0.23233E+03, 0.23632E+03, 0.24034E+03, 0.24439E+03, 0.24846E+03, &
         0.25255E+03, 0.25668E+03, 0.26083E+03, 0.26501E+03, 0.26921E+03, &
         0.27344E+03, 0.27770E+03, 0.28199E+03, 0.28631E+03, 0.29066E+03, &
         0.29503E+03, 0.29944E+03, 0.30387E+03, 0.30833E+03, 0.31282E+03, &
         0.31735E+03, 0.32190E+03, 0.32648E+03, 0.33110E+03, 0.33574E+03, &
         0.34042E+03, 0.34512E+03, 0.34986E+03, 0.35463E+03, 0.35943E+03, &
         0.36426E+03, 0.36913E+03, 0.37402E+03, 0.37895E+03, 0.38391E+03, &
         0.38891E+03, 0.39393E+03, 0.39899E+03, 0.40408E+03, 0.40921E+03, &
         0.41436E+03, 0.41955E+03, 0.42478E+03, 0.43004E+03, 0.43533E+03, &
         0.44065E+03, 0.44601E+03, 0.45140E+03, 0.45683E+03, 0.46229E+03, &
         0.46779E+03, 0.47332E+03, 0.47888E+03, 0.48448E+03, 0.49011E+03, &
         0.49578E+03/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_HF


!
!     *****************
      SUBROUTINE QT_HCL(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(2) :: XGJ
      REAL(DOUBLE), DIMENSION(2,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 8., 8./
!...      HCl
!...        --        15
      DATA (QOFT(1,J),J=1,119)/ 0.34775E+02, 0.48060E+02, 0.61370E+02, &
         0.74692E+02, 0.88024E+02, 0.10136E+03, 0.11471E+03, 0.12806E+03, &
         0.14141E+03, 0.15478E+03, 0.16814E+03, 0.18151E+03, 0.19489E+03, &
         0.20827E+03, 0.22166E+03, 0.23506E+03, 0.24847E+03, 0.26189E+03, &
         0.27533E+03, 0.28878E+03, 0.30225E+03, 0.31575E+03, 0.32928E+03, &
         0.34284E+03, 0.35645E+03, 0.37009E+03, 0.38378E+03, 0.39753E+03, &
         0.41134E+03, 0.42521E+03, 0.43914E+03, 0.45316E+03, 0.46725E+03, &
         0.48142E+03, 0.49568E+03, 0.51003E+03, 0.52448E+03, 0.53902E+03, &
         0.55368E+03, 0.56843E+03, 0.58330E+03, 0.59829E+03, 0.61339E+03, &
         0.62862E+03, 0.64396E+03, 0.65944E+03, 0.67504E+03, 0.69078E+03, &
         0.70665E+03, 0.72265E+03, 0.73880E+03, 0.75508E+03, 0.77151E+03, &
         0.78809E+03, 0.80481E+03, 0.82168E+03, 0.83870E+03, 0.85587E+03, &
         0.87320E+03, 0.89068E+03, 0.90832E+03, 0.92611E+03, 0.94407E+03, &
         0.96218E+03, 0.98046E+03, 0.99889E+03, 0.10175E+04, 0.10363E+04, &
         0.10552E+04, 0.10743E+04, 0.10936E+04, 0.11130E+04, 0.11326E+04, &
         0.11524E+04, 0.11723E+04, 0.11924E+04, 0.12127E+04, 0.12332E+04, &
         0.12538E+04, 0.12746E+04, 0.12956E+04, 0.13168E+04, 0.13381E+04, &
         0.13597E+04, 0.13814E+04, 0.14032E+04, 0.14253E+04, 0.14475E+04, &
         0.14700E+04, 0.14926E+04, 0.15153E+04, 0.15383E+04, 0.15615E+04, &
         0.15848E+04, 0.16083E+04, 0.16320E+04, 0.16559E+04, 0.16800E+04, &
         0.17043E+04, 0.17287E+04, 0.17533E+04, 0.17782E+04, 0.18032E+04, &
         0.18284E+04, 0.18538E+04, 0.18794E+04, 0.19051E+04, 0.19311E+04, &
         0.19573E+04, 0.19836E+04, 0.20102E+04, 0.20369E+04, 0.20638E+04, &
         0.20910E+04, 0.21183E+04, 0.21458E+04, 0.21735E+04, 0.22014E+04, &
         0.22295E+04/
!...        --        17
      DATA (QOFT(2,J),J=1,119)/ 0.34823E+02, 0.48128E+02, 0.61458E+02, &
         0.74801E+02, 0.88152E+02, 0.10151E+03, 0.11488E+03, 0.12825E+03, &
         0.14162E+03, 0.15500E+03, 0.16839E+03, 0.18178E+03, 0.19518E+03, &
         0.20858E+03, 0.22199E+03, 0.23541E+03, 0.24884E+03, 0.26228E+03, &
         0.27574E+03, 0.28921E+03, 0.30270E+03, 0.31622E+03, 0.32977E+03, &
         0.34336E+03, 0.35698E+03, 0.37065E+03, 0.38436E+03, 0.39813E+03, &
         0.41196E+03, 0.42585E+03, 0.43981E+03, 0.45384E+03, 0.46796E+03, &
         0.48215E+03, 0.49644E+03, 0.51081E+03, 0.52528E+03, 0.53986E+03, &
         0.55453E+03, 0.56932E+03, 0.58421E+03, 0.59922E+03, 0.61435E+03, &
         0.62960E+03, 0.64498E+03, 0.66048E+03, 0.67611E+03, 0.69187E+03, &
         0.70777E+03, 0.72381E+03, 0.73998E+03, 0.75630E+03, 0.77276E+03, &
         0.78936E+03, 0.80612E+03, 0.82302E+03, 0.84007E+03, 0.85727E+03, &
         0.87463E+03, 0.89215E+03, 0.90982E+03, 0.92765E+03, 0.94563E+03, &
         0.96378E+03, 0.98209E+03, 0.10006E+04, 0.10192E+04, 0.10380E+04, &
         0.10570E+04, 0.10761E+04, 0.10954E+04, 0.11149E+04, 0.11345E+04, &
         0.11543E+04, 0.11743E+04, 0.11945E+04, 0.12148E+04, 0.12353E+04, &
         0.12560E+04, 0.12768E+04, 0.12979E+04, 0.13191E+04, 0.13405E+04, &
         0.13620E+04, 0.13838E+04, 0.14057E+04, 0.14278E+04, 0.14501E+04, &
         0.14726E+04, 0.14952E+04, 0.15180E+04, 0.15410E+04, 0.15642E+04, &
         0.15876E+04, 0.16112E+04, 0.16349E+04, 0.16589E+04, 0.16830E+04, &
         0.17073E+04, 0.17318E+04, 0.17565E+04, 0.17814E+04, 0.18064E+04, &
         0.18317E+04, 0.18572E+04, 0.18828E+04, 0.19086E+04, 0.19346E+04, &
         0.19609E+04, 0.19873E+04, 0.20139E+04, 0.20406E+04, 0.20676E+04, &
         0.20948E+04, 0.21222E+04, 0.21498E+04, 0.21775E+04, 0.22055E+04, &
         0.22337E+04/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_HCL


!
!     *****************
      SUBROUTINE QT_HBR(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(2) :: XGJ
      REAL(DOUBLE), DIMENSION(2,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 8., 8./
!...      HBr
!...        --        19
      DATA (QOFT(1,J),J=1,119)/ 0.42744E+02, 0.59373E+02, 0.76023E+02, &
         0.92685E+02, 0.10936E+03, 0.12604E+03, 0.14272E+03, 0.15942E+03, &
         0.17612E+03, 0.19282E+03, 0.20954E+03, 0.22626E+03, 0.24299E+03, &
         0.25973E+03, 0.27648E+03, 0.29325E+03, 0.31004E+03, 0.32686E+03, &
         0.34371E+03, 0.36060E+03, 0.37753E+03, 0.39451E+03, 0.41156E+03, &
         0.42868E+03, 0.44587E+03, 0.46314E+03, 0.48051E+03, 0.49798E+03, &
         0.51556E+03, 0.53325E+03, 0.55106E+03, 0.56900E+03, 0.58708E+03, &
         0.60530E+03, 0.62367E+03, 0.64219E+03, 0.66088E+03, 0.67972E+03, &
         0.69874E+03, 0.71793E+03, 0.73730E+03, 0.75685E+03, 0.77659E+03, &
         0.79652E+03, 0.81664E+03, 0.83696E+03, 0.85748E+03, 0.87820E+03, &
         0.89914E+03, 0.92028E+03, 0.94163E+03, 0.96319E+03, 0.98498E+03, &
         0.10070E+04, 0.10292E+04, 0.10516E+04, 0.10743E+04, 0.10972E+04, &
         0.11203E+04, 0.11437E+04, 0.11673E+04, 0.11911E+04, 0.12151E+04, &
         0.12394E+04, 0.12640E+04, 0.12887E+04, 0.13137E+04, 0.13390E+04, &
         0.13645E+04, 0.13902E+04, 0.14162E+04, 0.14424E+04, 0.14689E+04, &
         0.14956E+04, 0.15226E+04, 0.15498E+04, 0.15773E+04, 0.16050E+04, &
         0.16330E+04, 0.16612E+04, 0.16897E+04, 0.17185E+04, 0.17475E+04, &
         0.17767E+04, 0.18062E+04, 0.18360E+04, 0.18660E+04, 0.18963E+04, &
         0.19269E+04, 0.19577E+04, 0.19888E+04, 0.20202E+04, 0.20518E+04, &
         0.20837E+04, 0.21158E+04, 0.21482E+04, 0.21809E+04, 0.22139E+04, &
         0.22471E+04, 0.22806E+04, 0.23143E+04, 0.23484E+04, 0.23827E+04, &
         0.24173E+04, 0.24521E+04, 0.24873E+04, 0.25227E+04, 0.25584E+04, &
         0.25943E+04, 0.26306E+04, 0.26671E+04, 0.27039E+04, 0.27409E+04, &
         0.27783E+04, 0.28159E+04, 0.28538E+04, 0.28920E+04, 0.29305E+04, &
         0.29693E+04/
!...        --        11
      DATA (QOFT(2,J),J=1,119)/ 0.42756E+02, 0.59390E+02, 0.76045E+02, &
         0.92713E+02, 0.10939E+03, 0.12607E+03, 0.14277E+03, 0.15947E+03, &
         0.17617E+03, 0.19288E+03, 0.20960E+03, 0.22633E+03, 0.24306E+03, &
         0.25981E+03, 0.27656E+03, 0.29334E+03, 0.31014E+03, 0.32696E+03, &
         0.34381E+03, 0.36071E+03, 0.37764E+03, 0.39464E+03, 0.41169E+03, &
         0.42881E+03, 0.44601E+03, 0.46329E+03, 0.48066E+03, 0.49813E+03, &
         0.51572E+03, 0.53341E+03, 0.55123E+03, 0.56918E+03, 0.58727E+03, &
         0.60549E+03, 0.62387E+03, 0.64240E+03, 0.66109E+03, 0.67994E+03, &
         0.69896E+03, 0.71816E+03, 0.73754E+03, 0.75710E+03, 0.77684E+03, &
         0.79678E+03, 0.81691E+03, 0.83724E+03, 0.85776E+03, 0.87850E+03, &
         0.89943E+03, 0.92058E+03, 0.94194E+03, 0.96352E+03, 0.98531E+03, &
         0.10073E+04, 0.10295E+04, 0.10520E+04, 0.10747E+04, 0.10976E+04, &
         0.11207E+04, 0.11441E+04, 0.11677E+04, 0.11915E+04, 0.12156E+04, &
         0.12399E+04, 0.12644E+04, 0.12892E+04, 0.13142E+04, 0.13395E+04, &
         0.13650E+04, 0.13907E+04, 0.14167E+04, 0.14429E+04, 0.14694E+04, &
         0.14961E+04, 0.15231E+04, 0.15504E+04, 0.15778E+04, 0.16056E+04, &
         0.16336E+04, 0.16618E+04, 0.16903E+04, 0.17191E+04, 0.17481E+04, &
         0.17773E+04, 0.18069E+04, 0.18367E+04, 0.18667E+04, 0.18970E+04, &
         0.19276E+04, 0.19584E+04, 0.19895E+04, 0.20209E+04, 0.20525E+04, &
         0.20844E+04, 0.21166E+04, 0.21490E+04, 0.21817E+04, 0.22147E+04, &
         0.22479E+04, 0.22814E+04, 0.23152E+04, 0.23492E+04, 0.23835E+04, &
         0.24181E+04, 0.24530E+04, 0.24882E+04, 0.25236E+04, 0.25593E+04, &
         0.25952E+04, 0.26315E+04, 0.26680E+04, 0.27048E+04, 0.27419E+04, &
         0.27793E+04, 0.28169E+04, 0.28549E+04, 0.28931E+04, 0.29316E+04, &
         0.29703E+04/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_HBR


!
!     *****************
      SUBROUTINE QT_HI(T, ISO, GSI, QT)          ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(1) :: XGJ
      REAL(DOUBLE), DIMENSION(1,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 12./
!...       HI
!...        --        17
      DATA (QOFT(1,J),J=1,119)/ 0.82031E+02, 0.11447E+03, 0.14694E+03, &
         0.17943E+03, 0.21194E+03, 0.24445E+03, 0.27699E+03, 0.30953E+03, &
         0.34209E+03, 0.37466E+03, 0.40725E+03, 0.43986E+03, 0.47249E+03, &
         0.50517E+03, 0.53789E+03, 0.57068E+03, 0.60354E+03, 0.63650E+03, &
         0.66957E+03, 0.70278E+03, 0.73614E+03, 0.76967E+03, 0.80340E+03, &
         0.83735E+03, 0.87153E+03, 0.90596E+03, 0.94067E+03, 0.97566E+03, &
         0.10110E+04, 0.10466E+04, 0.10826E+04, 0.11189E+04, 0.11555E+04, &
         0.11926E+04, 0.12300E+04, 0.12679E+04, 0.13061E+04, 0.13448E+04, &
         0.13839E+04, 0.14235E+04, 0.14635E+04, 0.15039E+04, 0.15448E+04, &
         0.15862E+04, 0.16280E+04, 0.16704E+04, 0.17132E+04, 0.17565E+04, &
         0.18003E+04, 0.18446E+04, 0.18894E+04, 0.19347E+04, 0.19806E+04, &
         0.20269E+04, 0.20738E+04, 0.21212E+04, 0.21691E+04, 0.22176E+04, &
         0.22666E+04, 0.23162E+04, 0.23662E+04, 0.24169E+04, 0.24680E+04, &
         0.25198E+04, 0.25720E+04, 0.26249E+04, 0.26783E+04, 0.27322E+04, &
         0.27867E+04, 0.28418E+04, 0.28975E+04, 0.29537E+04, 0.30105E+04, &
         0.30678E+04, 0.31258E+04, 0.31843E+04, 0.32434E+04, 0.33031E+04, &
         0.33633E+04, 0.34242E+04, 0.34856E+04, 0.35477E+04, 0.36103E+04, &
         0.36735E+04, 0.37373E+04, 0.38018E+04, 0.38668E+04, 0.39324E+04, &
         0.39986E+04, 0.40654E+04, 0.41329E+04, 0.42009E+04, 0.42696E+04, &
         0.43388E+04, 0.44087E+04, 0.44792E+04, 0.45503E+04, 0.46221E+04, &
         0.46944E+04, 0.47674E+04, 0.48410E+04, 0.49152E+04, 0.49901E+04, &
         0.50656E+04, 0.51417E+04, 0.52185E+04, 0.52959E+04, 0.53739E+04, &
         0.54526E+04, 0.55319E+04, 0.56118E+04, 0.56924E+04, 0.57736E+04, &
         0.58555E+04, 0.59380E+04, 0.60212E+04, 0.61050E+04, 0.61895E+04, &
         0.62746E+04/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_HI


!
!     *****************
      SUBROUTINE QT_CLO(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(2) :: XGJ
      REAL(DOUBLE), DIMENSION(2,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 4., 4./
!...      ClO
!...        --        56
      DATA (QOFT(1,J),J=1,119)/ 0.53847E+03, 0.76580E+03, 0.10017E+04, &
         0.12511E+04, 0.15168E+04, 0.18001E+04, 0.21014E+04, 0.24206E+04, &
         0.27577E+04, 0.31127E+04, 0.34857E+04, 0.38765E+04, 0.42854E+04, &
         0.47124E+04, 0.51575E+04, 0.56208E+04, 0.61025E+04, 0.66026E+04, &
         0.71211E+04, 0.76582E+04, 0.82138E+04, 0.87882E+04, 0.93813E+04, &
         0.99932E+04, 0.10624E+05, 0.11273E+05, 0.11942E+05, 0.12629E+05, &
         0.13336E+05, 0.14061E+05, 0.14806E+05, 0.15570E+05, 0.16353E+05, &
         0.17155E+05, 0.17976E+05, 0.18816E+05, 0.19676E+05, 0.20555E+05, &
         0.21453E+05, 0.22371E+05, 0.23308E+05, 0.24264E+05, 0.25240E+05, &
         0.26236E+05, 0.27250E+05, 0.28284E+05, 0.29338E+05, 0.30412E+05, &
         0.31505E+05, 0.32617E+05, 0.33749E+05, 0.34901E+05, 0.36072E+05, &
         0.37263E+05, 0.38474E+05, 0.39705E+05, 0.40955E+05, 0.42225E+05, &
         0.43515E+05, 0.44825E+05, 0.46154E+05, 0.47504E+05, 0.48873E+05, &
         0.50262E+05, 0.51672E+05, 0.53101E+05, 0.54549E+05, 0.56019E+05, &
         0.57508E+05, 0.59017E+05, 0.60546E+05, 0.62095E+05, 0.63665E+05, &
         0.65254E+05, 0.66864E+05, 0.68494E+05, 0.70144E+05, 0.71814E+05, &
         0.73504E+05, 0.75215E+05, 0.76946E+05, 0.78698E+05, 0.80470E+05, &
         0.82261E+05, 0.84074E+05, 0.85907E+05, 0.87760E+05, 0.89633E+05, &
         0.91527E+05, 0.93442E+05, 0.95377E+05, 0.97333E+05, 0.99309E+05, &
         0.10131E+06, 0.10332E+06, 0.10536E+06, 0.10742E+06, 0.10950E+06, &
         0.11160E+06, 0.11372E+06, 0.11586E+06, 0.11802E+06, 0.12020E+06, &
         0.12241E+06, 0.12463E+06, 0.12688E+06, 0.12914E+06, 0.13143E+06, &
         0.13374E+06, 0.13607E+06, 0.13842E+06, 0.14079E+06, 0.14318E+06, &
         0.14559E+06, 0.14802E+06, 0.15048E+06, 0.15295E+06, 0.15545E+06, &
         0.15797E+06/
!...        --        76
      DATA (QOFT(2,J),J=1,119)/ 0.54775E+03, 0.77899E+03, 0.10189E+04, &
         0.12726E+04, 0.15430E+04, 0.18313E+04, 0.21378E+04, 0.24627E+04, &
         0.28059E+04, 0.31674E+04, 0.35472E+04, 0.39454E+04, 0.43621E+04, &
         0.47972E+04, 0.52508E+04, 0.57232E+04, 0.62143E+04, 0.67242E+04, &
         0.72531E+04, 0.78010E+04, 0.83678E+04, 0.89537E+04, 0.95589E+04, &
         0.10183E+05, 0.10827E+05, 0.11490E+05, 0.12172E+05, 0.12874E+05, &
         0.13595E+05, 0.14335E+05, 0.15095E+05, 0.15875E+05, 0.16674E+05, &
         0.17493E+05, 0.18332E+05, 0.19190E+05, 0.20068E+05, 0.20965E+05, &
         0.21882E+05, 0.22820E+05, 0.23776E+05, 0.24753E+05, 0.25750E+05, &
         0.26766E+05, 0.27803E+05, 0.28859E+05, 0.29935E+05, 0.31032E+05, &
         0.32148E+05, 0.33284E+05, 0.34441E+05, 0.35617E+05, 0.36814E+05, &
         0.38031E+05, 0.39267E+05, 0.40524E+05, 0.41802E+05, 0.43099E+05, &
         0.44417E+05, 0.45755E+05, 0.47113E+05, 0.48492E+05, 0.49891E+05, &
         0.51310E+05, 0.52750E+05, 0.54210E+05, 0.55690E+05, 0.57191E+05, &
         0.58713E+05, 0.60255E+05, 0.61817E+05, 0.63400E+05, 0.65004E+05, &
         0.66628E+05, 0.68272E+05, 0.69938E+05, 0.71624E+05, 0.73331E+05, &
         0.75058E+05, 0.76806E+05, 0.78575E+05, 0.80364E+05, 0.82175E+05, &
         0.84006E+05, 0.85858E+05, 0.87731E+05, 0.89625E+05, 0.91539E+05, &
         0.93475E+05, 0.95431E+05, 0.97409E+05, 0.99407E+05, 0.10143E+06, &
         0.10347E+06, 0.10553E+06, 0.10761E+06, 0.10972E+06, 0.11184E+06, &
         0.11399E+06, 0.11615E+06, 0.11834E+06, 0.12055E+06, 0.12278E+06, &
         0.12503E+06, 0.12731E+06, 0.12960E+06, 0.13192E+06, 0.13425E+06, &
         0.13661E+06, 0.13899E+06, 0.14139E+06, 0.14382E+06, 0.14626E+06, &
         0.14873E+06, 0.15121E+06, 0.15372E+06, 0.15625E+06, 0.15880E+06, &
         0.16138E+06/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_CLO


!
!     *****************
      SUBROUTINE QT_OCS(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(5) :: XGJ
      REAL(DOUBLE), DIMENSION(5,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1., 1., 2., 4., 1./
!...      OCS
!...        --       622
      DATA (QOFT(1,J),J=1,119)/ 0.20609E+03, 0.29199E+03, 0.37861E+03, &
         0.46737E+03, 0.56024E+03, 0.65929E+03, 0.76649E+03, 0.88361E+03, &
         0.10123E+04, 0.11541E+04, 0.13105E+04, 0.14829E+04, 0.16728E+04, &
         0.18818E+04, 0.21113E+04, 0.23629E+04, 0.26383E+04, 0.29391E+04, &
         0.32672E+04, 0.36245E+04, 0.40128E+04, 0.44343E+04, 0.48911E+04, &
         0.53853E+04, 0.59193E+04, 0.64956E+04, 0.71166E+04, 0.77849E+04, &
         0.85033E+04, 0.92746E+04, 0.10102E+05, 0.10988E+05, 0.11936E+05, &
         0.12949E+05, 0.14032E+05, 0.15186E+05, 0.16416E+05, 0.17726E+05, &
         0.19120E+05, 0.20601E+05, 0.22173E+05, 0.23842E+05, 0.25611E+05, &
         0.27484E+05, 0.29468E+05, 0.31566E+05, 0.33783E+05, 0.36124E+05, &
         0.38595E+05, 0.41202E+05, 0.43949E+05, 0.46842E+05, 0.49888E+05, &
         0.53092E+05, 0.56460E+05, 0.59999E+05, 0.63716E+05, 0.67616E+05, &
         0.71708E+05, 0.75997E+05, 0.80491E+05, 0.85197E+05, 0.90124E+05, &
         0.95278E+05, 0.10067E+06, 0.10630E+06, 0.11219E+06, 0.11833E+06, &
         0.12475E+06, 0.13144E+06, 0.13842E+06, 0.14570E+06, 0.15328E+06, &
         0.16117E+06, 0.16940E+06, 0.17795E+06, 0.18686E+06, 0.19611E+06, &
         0.20574E+06, 0.21574E+06, 0.22613E+06, 0.23692E+06, 0.24813E+06, &
         0.25975E+06, 0.27182E+06, 0.28433E+06, 0.29730E+06, 0.31074E+06, &
         0.32467E+06, 0.33909E+06, 0.35403E+06, 0.36950E+06, 0.38551E+06, &
         0.40207E+06, 0.41920E+06, 0.43691E+06, 0.45522E+06, 0.47415E+06, &
         0.49370E+06, 0.51390E+06, 0.53476E+06, 0.55629E+06, 0.57852E+06, &
         0.60146E+06, 0.62513E+06, 0.64954E+06, 0.67471E+06, 0.70067E+06, &
         0.72742E+06, 0.75499E+06, 0.78339E+06, 0.81265E+06, 0.84279E+06, &
         0.87381E+06, 0.90576E+06, 0.93863E+06, 0.97246E+06, 0.10073E+07, &
         0.10431E+07/
!...        --       624
      DATA (QOFT(2,J),J=1,119)/ 0.21125E+03, 0.29930E+03, 0.38809E+03, &
         0.47911E+03, 0.57437E+03, 0.67603E+03, 0.78610E+03, 0.90643E+03, &
         0.10387E+04, 0.11846E+04, 0.13456E+04, 0.15231E+04, 0.17188E+04, &
         0.19342E+04, 0.21709E+04, 0.24304E+04, 0.27145E+04, 0.30250E+04, &
         0.33638E+04, 0.37328E+04, 0.41339E+04, 0.45694E+04, 0.50415E+04, &
         0.55524E+04, 0.61045E+04, 0.67004E+04, 0.73427E+04, 0.80340E+04, &
         0.87773E+04, 0.95755E+04, 0.10432E+05, 0.11349E+05, 0.12330E+05, &
         0.13380E+05, 0.14500E+05, 0.15696E+05, 0.16970E+05, 0.18327E+05, &
         0.19770E+05, 0.21305E+05, 0.22934E+05, 0.24663E+05, 0.26497E+05, &
         0.28439E+05, 0.30495E+05, 0.32669E+05, 0.34968E+05, 0.37396E+05, &
         0.39958E+05, 0.42661E+05, 0.45510E+05, 0.48511E+05, 0.51669E+05, &
         0.54993E+05, 0.58487E+05, 0.62159E+05, 0.66014E+05, 0.70061E+05, &
         0.74306E+05, 0.78757E+05, 0.83421E+05, 0.88305E+05, 0.93418E+05, &
         0.98767E+05, 0.10436E+06, 0.11021E+06, 0.11632E+06, 0.12270E+06, &
         0.12936E+06, 0.13631E+06, 0.14355E+06, 0.15111E+06, 0.15898E+06, &
         0.16718E+06, 0.17572E+06, 0.18460E+06, 0.19385E+06, 0.20346E+06, &
         0.21346E+06, 0.22385E+06, 0.23464E+06, 0.24585E+06, 0.25748E+06, &
         0.26956E+06, 0.28209E+06, 0.29509E+06, 0.30856E+06, 0.32252E+06, &
         0.33699E+06, 0.35198E+06, 0.36750E+06, 0.38357E+06, 0.40020E+06, &
         0.41741E+06, 0.43521E+06, 0.45362E+06, 0.47264E+06, 0.49231E+06, &
         0.51263E+06, 0.53362E+06, 0.55529E+06, 0.57768E+06, 0.60078E+06, &
         0.62462E+06, 0.64922E+06, 0.67459E+06, 0.70075E+06, 0.72773E+06, &
         0.75554E+06, 0.78419E+06, 0.81372E+06, 0.84413E+06, 0.87546E+06, &
         0.90771E+06, 0.94092E+06, 0.97509E+06, 0.10103E+07, 0.10464E+07, &
         0.10837E+07/
!...        --       632
      DATA (QOFT(3,J),J=1,119)/ 0.41351E+03, 0.58591E+03, 0.76004E+03, &
         0.93907E+03, 0.11273E+04, 0.13289E+04, 0.15481E+04, 0.17884E+04, &
         0.20533E+04, 0.23459E+04, 0.26692E+04, 0.30264E+04, 0.34205E+04, &
         0.38547E+04, 0.43323E+04, 0.48565E+04, 0.54309E+04, 0.60592E+04, &
         0.67451E+04, 0.74928E+04, 0.83064E+04, 0.91903E+04, 0.10149E+05, &
         0.11187E+05, 0.12310E+05, 0.13523E+05, 0.14831E+05, 0.16240E+05, &
         0.17756E+05, 0.19384E+05, 0.21132E+05, 0.23005E+05, 0.25011E+05, &
         0.27157E+05, 0.29449E+05, 0.31896E+05, 0.34506E+05, 0.37286E+05, &
         0.40245E+05, 0.43392E+05, 0.46735E+05, 0.50284E+05, 0.54048E+05, &
         0.58038E+05, 0.62263E+05, 0.66733E+05, 0.71460E+05, 0.76455E+05, &
         0.81728E+05, 0.87292E+05, 0.93159E+05, 0.99341E+05, 0.10585E+06, &
         0.11270E+06, 0.11991E+06, 0.12748E+06, 0.13543E+06, 0.14378E+06, &
         0.15255E+06, 0.16174E+06, 0.17137E+06, 0.18146E+06, 0.19202E+06, &
         0.20308E+06, 0.21465E+06, 0.22674E+06, 0.23937E+06, 0.25257E+06, &
         0.26635E+06, 0.28073E+06, 0.29573E+06, 0.31137E+06, 0.32767E+06, &
         0.34466E+06, 0.36235E+06, 0.38076E+06, 0.39992E+06, 0.41985E+06, &
         0.44057E+06, 0.46211E+06, 0.48450E+06, 0.50775E+06, 0.53189E+06, &
         0.55695E+06, 0.58295E+06, 0.60992E+06, 0.63789E+06, 0.66688E+06, &
         0.69693E+06, 0.72806E+06, 0.76030E+06, 0.79368E+06, 0.82823E+06, &
         0.86399E+06, 0.90097E+06, 0.93923E+06, 0.97878E+06, 0.10197E+07, &
         0.10619E+07, 0.11056E+07, 0.11506E+07, 0.11972E+07, 0.12453E+07, &
         0.12949E+07, 0.13460E+07, 0.13988E+07, 0.14533E+07, 0.15094E+07, &
         0.15673E+07, 0.16270E+07, 0.16884E+07, 0.17518E+07, 0.18170E+07, &
         0.18842E+07, 0.19533E+07, 0.20245E+07, 0.20978E+07, 0.21732E+07, &
         0.22507E+07/
!...        --       623
      DATA (QOFT(4,J),J=1,119)/ 0.83485E+03, 0.11828E+04, 0.15337E+04, &
         0.18934E+04, 0.22697E+04, 0.26712E+04, 0.31059E+04, 0.35809E+04, &
         0.41030E+04, 0.46785E+04, 0.53133E+04, 0.60135E+04, 0.67850E+04, &
         0.76338E+04, 0.85663E+04, 0.95888E+04, 0.10708E+05, 0.11931E+05, &
         0.13265E+05, 0.14718E+05, 0.16298E+05, 0.18012E+05, 0.19870E+05, &
         0.21881E+05, 0.24054E+05, 0.26399E+05, 0.28926E+05, 0.31646E+05, &
         0.34570E+05, 0.37710E+05, 0.41077E+05, 0.44685E+05, 0.48545E+05, &
         0.52672E+05, 0.57078E+05, 0.61780E+05, 0.66790E+05, 0.72125E+05, &
         0.77801E+05, 0.83833E+05, 0.90239E+05, 0.97036E+05, 0.10424E+06, &
         0.11188E+06, 0.11996E+06, 0.12850E+06, 0.13754E+06, 0.14708E+06, &
         0.15715E+06, 0.16777E+06, 0.17896E+06, 0.19076E+06, 0.20317E+06, &
         0.21623E+06, 0.22996E+06, 0.24438E+06, 0.25953E+06, 0.27543E+06, &
         0.29211E+06, 0.30959E+06, 0.32791E+06, 0.34710E+06, 0.36718E+06, &
         0.38820E+06, 0.41017E+06, 0.43314E+06, 0.45713E+06, 0.48219E+06, &
         0.50835E+06, 0.53564E+06, 0.56409E+06, 0.59376E+06, 0.62468E+06, &
         0.65688E+06, 0.69041E+06, 0.72530E+06, 0.76161E+06, 0.79937E+06, &
         0.83862E+06, 0.87941E+06, 0.92179E+06, 0.96581E+06, 0.10115E+07, &
         0.10589E+07, 0.11081E+07, 0.11591E+07, 0.12120E+07, 0.12669E+07, &
         0.13237E+07, 0.13825E+07, 0.14435E+07, 0.15066E+07, 0.15718E+07, &
         0.16394E+07, 0.17093E+07, 0.17815E+07, 0.18562E+07, 0.19334E+07, &
         0.20132E+07, 0.20956E+07, 0.21807E+07, 0.22685E+07, 0.23592E+07, &
         0.24528E+07, 0.25494E+07, 0.26490E+07, 0.27517E+07, 0.28576E+07, &
         0.29667E+07, 0.30792E+07, 0.31951E+07, 0.33145E+07, 0.34374E+07, &
         0.35640E+07, 0.36943E+07, 0.38285E+07, 0.39665E+07, 0.41085E+07, &
         0.42546E+07/
!...        --       822
      DATA (QOFT(5,J),J=1,119)/ 0.21967E+03, 0.31126E+03, 0.40370E+03, &
         0.49862E+03, 0.59823E+03, 0.70481E+03, 0.82050E+03, 0.94724E+03, &
         0.10868E+04, 0.12409E+04, 0.14112E+04, 0.15993E+04, 0.18067E+04, &
         0.20353E+04, 0.22866E+04, 0.25624E+04, 0.28645E+04, 0.31950E+04, &
         0.35558E+04, 0.39490E+04, 0.43767E+04, 0.48413E+04, 0.53452E+04, &
         0.58909E+04, 0.64810E+04, 0.71182E+04, 0.78053E+04, 0.85454E+04, &
         0.93413E+04, 0.10196E+05, 0.11114E+05, 0.12098E+05, 0.13151E+05, &
         0.14277E+05, 0.15480E+05, 0.16764E+05, 0.18133E+05, 0.19592E+05, &
         0.21144E+05, 0.22794E+05, 0.24548E+05, 0.26409E+05, 0.28383E+05, &
         0.30475E+05, 0.32689E+05, 0.35033E+05, 0.37511E+05, 0.40128E+05, &
         0.42892E+05, 0.45808E+05, 0.48882E+05, 0.52121E+05, 0.55532E+05, &
         0.59121E+05, 0.62895E+05, 0.66861E+05, 0.71028E+05, 0.75402E+05, &
         0.79991E+05, 0.84803E+05, 0.89847E+05, 0.95130E+05, 0.10066E+06, &
         0.10645E+06, 0.11251E+06, 0.11883E+06, 0.12545E+06, 0.13236E+06, &
         0.13957E+06, 0.14710E+06, 0.15495E+06, 0.16313E+06, 0.17166E+06, &
         0.18055E+06, 0.18980E+06, 0.19944E+06, 0.20946E+06, 0.21989E+06, &
         0.23073E+06, 0.24200E+06, 0.25371E+06, 0.26587E+06, 0.27850E+06, &
         0.29161E+06, 0.30521E+06, 0.31931E+06, 0.33394E+06, 0.34910E+06, &
         0.36482E+06, 0.38109E+06, 0.39795E+06, 0.41541E+06, 0.43348E+06, &
         0.45217E+06, 0.47151E+06, 0.49151E+06, 0.51219E+06, 0.53356E+06, &
         0.55565E+06, 0.57847E+06, 0.60204E+06, 0.62637E+06, 0.65149E+06, &
         0.67742E+06, 0.70417E+06, 0.73176E+06, 0.76023E+06, 0.78957E+06, &
         0.81982E+06, 0.85100E+06, 0.88313E+06, 0.91622E+06, 0.95031E+06, &
         0.98541E+06, 0.10216E+07, 0.10587E+07, 0.10970E+07, 0.11364E+07, &
         0.11769E+07/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_OCS


!
!     *****************
      SUBROUTINE QT_H2CO(T, ISO, GSI, QT)        ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(3) :: XGJ
      REAL(DOUBLE), DIMENSION(3,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1., 2., 1./
!...     H2CO
!...        --       126
      DATA (QOFT(1,J),J=1,119)/ 0.25934E+03, 0.43623E+03, 0.64143E+03, &
         0.87152E+03, 0.11241E+04, 0.13975E+04, 0.16906E+04, 0.20029E+04, &
         0.23344E+04, 0.26857E+04, 0.30577E+04, 0.34518E+04, 0.38698E+04, &
         0.43138E+04, 0.47860E+04, 0.52890E+04, 0.58256E+04, 0.63985E+04, &
         0.70109E+04, 0.76660E+04, 0.83673E+04, 0.91184E+04, 0.99230E+04, &
         0.10785E+05, 0.11710E+05, 0.12700E+05, 0.13762E+05, 0.14900E+05, &
         0.16119E+05, 0.17425E+05, 0.18823E+05, 0.20320E+05, 0.21923E+05, &
         0.23637E+05, 0.25471E+05, 0.27432E+05, 0.29527E+05, 0.31765E+05, &
         0.34155E+05, 0.36706E+05, 0.39428E+05, 0.42330E+05, 0.45424E+05, &
         0.48720E+05, 0.52231E+05, 0.55968E+05, 0.59945E+05, 0.64175E+05, &
         0.68672E+05, 0.73450E+05, 0.78526E+05, 0.83915E+05, 0.89634E+05, &
         0.95701E+05, 0.10213E+06, 0.10895E+06, 0.11618E+06, 0.12383E+06, &
         0.13193E+06, 0.14049E+06, 0.14956E+06, 0.15914E+06, 0.16927E+06, &
         0.17997E+06, 0.19127E+06, 0.20320E+06, 0.21578E+06, 0.22906E+06, &
         0.24306E+06, 0.25782E+06, 0.27336E+06, 0.28974E+06, 0.30698E+06, &
         0.32513E+06, 0.34422E+06, 0.36430E+06, 0.38542E+06, 0.40761E+06, &
         0.43093E+06, 0.45542E+06, 0.48114E+06, 0.50813E+06, 0.53646E+06, &
         0.56617E+06, 0.59733E+06, 0.63000E+06, 0.66423E+06, 0.70010E+06, &
         0.73767E+06, 0.77701E+06, 0.81818E+06, 0.86127E+06, 0.90635E+06, &
         0.95349E+06, 0.10028E+07, 0.10543E+07, 0.11082E+07, 0.11644E+07, &
         0.12232E+07, 0.12845E+07, 0.13485E+07, 0.14154E+07, 0.14851E+07, &
         0.15578E+07, 0.16337E+07, 0.17127E+07, 0.17952E+07, 0.18810E+07, &
         0.19705E+07, 0.20637E+07, 0.21607E+07, 0.22617E+07, 0.23669E+07, &
         0.24763E+07, 0.25901E+07, 0.27085E+07, 0.28316E+07, 0.29596E+07, &
         0.30926E+07/
!...        --       136
      DATA (QOFT(2,J),J=1,119)/ 0.53173E+03, 0.89447E+03, 0.13153E+04, &
         0.17871E+04, 0.23051E+04, 0.28658E+04, 0.34669E+04, 0.41073E+04, &
         0.47872E+04, 0.55074E+04, 0.62702E+04, 0.70785E+04, 0.79357E+04, &
         0.88462E+04, 0.98147E+04, 0.10846E+05, 0.11946E+05, 0.13121E+05, &
         0.14377E+05, 0.15721E+05, 0.17159E+05, 0.18699E+05, 0.20349E+05, &
         0.22118E+05, 0.24013E+05, 0.26045E+05, 0.28222E+05, 0.30555E+05, &
         0.33055E+05, 0.35733E+05, 0.38601E+05, 0.41671E+05, 0.44958E+05, &
         0.48474E+05, 0.52235E+05, 0.56255E+05, 0.60552E+05, 0.65142E+05, &
         0.70043E+05, 0.75275E+05, 0.80856E+05, 0.86808E+05, 0.93152E+05, &
         0.99913E+05, 0.10711E+06, 0.11478E+06, 0.12293E+06, 0.13161E+06, &
         0.14083E+06, 0.15063E+06, 0.16104E+06, 0.17209E+06, 0.18382E+06, &
         0.19626E+06, 0.20945E+06, 0.22343E+06, 0.23825E+06, 0.25394E+06, &
         0.27054E+06, 0.28812E+06, 0.30671E+06, 0.32636E+06, 0.34713E+06, &
         0.36907E+06, 0.39224E+06, 0.41671E+06, 0.44252E+06, 0.46975E+06, &
         0.49845E+06, 0.52872E+06, 0.56060E+06, 0.59418E+06, 0.62954E+06, &
         0.66676E+06, 0.70591E+06, 0.74710E+06, 0.79040E+06, 0.83591E+06, &
         0.88373E+06, 0.93395E+06, 0.98669E+06, 0.10421E+07, 0.11001E+07, &
         0.11611E+07, 0.12250E+07, 0.12920E+07, 0.13622E+07, 0.14357E+07, &
         0.15128E+07, 0.15934E+07, 0.16779E+07, 0.17662E+07, 0.18587E+07, &
         0.19554E+07, 0.20565E+07, 0.21621E+07, 0.22725E+07, 0.23879E+07, &
         0.25084E+07, 0.26342E+07, 0.27655E+07, 0.29026E+07, 0.30456E+07, &
         0.31947E+07, 0.33502E+07, 0.35124E+07, 0.36814E+07, 0.38575E+07, &
         0.40410E+07, 0.42321E+07, 0.44311E+07, 0.46382E+07, 0.48538E+07, &
         0.50782E+07, 0.53116E+07, 0.55544E+07, 0.58068E+07, 0.60693E+07, &
         0.63421E+07/
!...        --       128
      DATA (QOFT(3,J),J=1,119)/ 0.27198E+03, 0.45755E+03, 0.67282E+03, &
         0.91421E+03, 0.11792E+04, 0.14660E+04, 0.17735E+04, 0.21012E+04, &
         0.24490E+04, 0.28175E+04, 0.32077E+04, 0.36212E+04, 0.40598E+04, &
         0.45256E+04, 0.50211E+04, 0.55488E+04, 0.61116E+04, 0.67127E+04, &
         0.73552E+04, 0.80426E+04, 0.87783E+04, 0.95663E+04, 0.10410E+05, &
         0.11315E+05, 0.12285E+05, 0.13324E+05, 0.14438E+05, 0.15632E+05, &
         0.16911E+05, 0.18281E+05, 0.19748E+05, 0.21319E+05, 0.23000E+05, &
         0.24799E+05, 0.26723E+05, 0.28780E+05, 0.30978E+05, 0.33326E+05, &
         0.35834E+05, 0.38510E+05, 0.41365E+05, 0.44410E+05, 0.47656E+05, &
         0.51115E+05, 0.54798E+05, 0.58719E+05, 0.62891E+05, 0.67329E+05, &
         0.72047E+05, 0.77060E+05, 0.82385E+05, 0.88039E+05, 0.94039E+05, &
         0.10040E+06, 0.10715E+06, 0.11431E+06, 0.12189E+06, 0.12991E+06, &
         0.13841E+06, 0.14740E+06, 0.15691E+06, 0.16696E+06, 0.17759E+06, &
         0.18882E+06, 0.20067E+06, 0.21318E+06, 0.22639E+06, 0.24032E+06, &
         0.25501E+06, 0.27049E+06, 0.28680E+06, 0.30398E+06, 0.32207E+06, &
         0.34111E+06, 0.36114E+06, 0.38221E+06, 0.40436E+06, 0.42765E+06, &
         0.45211E+06, 0.47781E+06, 0.50479E+06, 0.53311E+06, 0.56283E+06, &
         0.59400E+06, 0.62669E+06, 0.66097E+06, 0.69688E+06, 0.73451E+06, &
         0.77393E+06, 0.81520E+06, 0.85840E+06, 0.90360E+06, 0.95090E+06, &
         0.10004E+07, 0.10521E+07, 0.11061E+07, 0.11626E+07, 0.12216E+07, &
         0.12833E+07, 0.13476E+07, 0.14148E+07, 0.14849E+07, 0.15581E+07, &
         0.16344E+07, 0.17140E+07, 0.17969E+07, 0.18834E+07, 0.19735E+07, &
         0.20674E+07, 0.21651E+07, 0.22669E+07, 0.23729E+07, 0.24832E+07, &
         0.25980E+07, 0.27174E+07, 0.28416E+07, 0.29708E+07, 0.31050E+07, &
         0.32446E+07/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_H2CO


!
!     *****************
      SUBROUTINE QT_HOCL(T, ISO, GSI, QT)        ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(2) :: XGJ
      REAL(DOUBLE), DIMENSION(2,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 8., 8./
!...     HOCl
!...        --       165
      DATA (QOFT(1,J),J=1,119)/ 0.17041E+04, 0.28708E+04, 0.42250E+04, &
         0.57456E+04, 0.74211E+04, 0.92470E+04, 0.11225E+05, 0.13359E+05, &
         0.15657E+05, 0.18129E+05, 0.20785E+05, 0.23637E+05, 0.26696E+05, &
         0.29974E+05, 0.33484E+05, 0.37239E+05, 0.41252E+05, 0.45536E+05, &
         0.50105E+05, 0.54973E+05, 0.60152E+05, 0.65659E+05, 0.71507E+05, &
         0.77711E+05, 0.84286E+05, 0.91249E+05, 0.98614E+05, 0.10640E+06, &
         0.11462E+06, 0.12330E+06, 0.13244E+06, 0.14208E+06, 0.15222E+06, &
         0.16289E+06, 0.17411E+06, 0.18589E+06, 0.19825E+06, 0.21123E+06, &
         0.22483E+06, 0.23908E+06, 0.25400E+06, 0.26962E+06, 0.28596E+06, &
         0.30303E+06, 0.32087E+06, 0.33950E+06, 0.35895E+06, 0.37923E+06, &
         0.40038E+06, 0.42243E+06, 0.44539E+06, 0.46930E+06, 0.49419E+06, &
         0.52008E+06, 0.54700E+06, 0.57498E+06, 0.60406E+06, 0.63426E+06, &
         0.66562E+06, 0.69816E+06, 0.73192E+06, 0.76692E+06, 0.80322E+06, &
         0.84083E+06, 0.87979E+06, 0.92014E+06, 0.96192E+06, 0.10052E+07, &
         0.10499E+07, 0.10961E+07, 0.11440E+07, 0.11934E+07, 0.12445E+07, &
         0.12973E+07, 0.13518E+07, 0.14081E+07, 0.14661E+07, 0.15261E+07, &
         0.15879E+07, 0.16516E+07, 0.17174E+07, 0.17851E+07, 0.18550E+07, &
         0.19269E+07, 0.20010E+07, 0.20773E+07, 0.21559E+07, 0.22367E+07, &
         0.23200E+07, 0.24056E+07, 0.24936E+07, 0.25842E+07, 0.26773E+07, &
         0.27730E+07, 0.28714E+07, 0.29724E+07, 0.30763E+07, 0.31829E+07, &
         0.32924E+07, 0.34049E+07, 0.35203E+07, 0.36387E+07, 0.37603E+07, &
         0.38850E+07, 0.40129E+07, 0.41441E+07, 0.42786E+07, 0.44165E+07, &
         0.45579E+07, 0.47028E+07, 0.48512E+07, 0.50033E+07, 0.51592E+07, &
         0.53187E+07, 0.54822E+07, 0.56495E+07, 0.58208E+07, 0.59961E+07, &
         0.61755E+07/
!...        --       167
      DATA (QOFT(2,J),J=1,119)/ 0.17342E+04, 0.29215E+04, 0.42998E+04, &
         0.58473E+04, 0.75524E+04, 0.94107E+04, 0.11423E+05, 0.13595E+05, &
         0.15935E+05, 0.18450E+05, 0.21154E+05, 0.24056E+05, 0.27168E+05, &
         0.30505E+05, 0.34077E+05, 0.37899E+05, 0.41983E+05, 0.46343E+05, &
         0.50993E+05, 0.55947E+05, 0.61218E+05, 0.66822E+05, 0.72774E+05, &
         0.79088E+05, 0.85780E+05, 0.92866E+05, 0.10036E+06, 0.10829E+06, &
         0.11665E+06, 0.12548E+06, 0.13479E+06, 0.14460E+06, 0.15492E+06, &
         0.16578E+06, 0.17719E+06, 0.18918E+06, 0.20177E+06, 0.21497E+06, &
         0.22881E+06, 0.24332E+06, 0.25851E+06, 0.27440E+06, 0.29102E+06, &
         0.30840E+06, 0.32656E+06, 0.34552E+06, 0.36531E+06, 0.38595E+06, &
         0.40748E+06, 0.42991E+06, 0.45328E+06, 0.47762E+06, 0.50295E+06, &
         0.52929E+06, 0.55669E+06, 0.58517E+06, 0.61477E+06, 0.64550E+06, &
         0.67741E+06, 0.71053E+06, 0.74489E+06, 0.78052E+06, 0.81745E+06, &
         0.85573E+06, 0.89539E+06, 0.93645E+06, 0.97897E+06, 0.10230E+07, &
         0.10685E+07, 0.11156E+07, 0.11643E+07, 0.12146E+07, 0.12666E+07, &
         0.13203E+07, 0.13757E+07, 0.14330E+07, 0.14921E+07, 0.15531E+07, &
         0.16160E+07, 0.16809E+07, 0.17478E+07, 0.18168E+07, 0.18878E+07, &
         0.19611E+07, 0.20365E+07, 0.21141E+07, 0.21941E+07, 0.22764E+07, &
         0.23611E+07, 0.24482E+07, 0.25378E+07, 0.26300E+07, 0.27248E+07, &
         0.28222E+07, 0.29223E+07, 0.30251E+07, 0.31308E+07, 0.32393E+07, &
         0.33508E+07, 0.34652E+07, 0.35827E+07, 0.37032E+07, 0.38269E+07, &
         0.39539E+07, 0.40840E+07, 0.42176E+07, 0.43545E+07, 0.44948E+07, &
         0.46387E+07, 0.47861E+07, 0.49372E+07, 0.50920E+07, 0.52506E+07, &
         0.54130E+07, 0.55793E+07, 0.57496E+07, 0.59239E+07, 0.61024E+07, &
         0.62850E+07/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_HOCL


!
!     *****************
      SUBROUTINE QT_N2(T, ISO, GSI, QT)          ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(1) :: XGJ
      REAL(DOUBLE), DIMENSION(1,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1./
!...       N2
!...        --        44
      DATA (QOFT(1,J),J=1,119)/ 0.95487E+02, 0.13466E+03, 0.17386E+03, &
         0.21307E+03, 0.25230E+03, 0.29154E+03, 0.33080E+03, 0.37008E+03, &
         0.40937E+03, 0.44868E+03, 0.48800E+03, 0.52736E+03, 0.56674E+03, &
         0.60616E+03, 0.64562E+03, 0.68515E+03, 0.72475E+03, 0.76445E+03, &
         0.80426E+03, 0.84420E+03, 0.88430E+03, 0.92457E+03, 0.96505E+03, &
         0.10057E+04, 0.10467E+04, 0.10879E+04, 0.11293E+04, 0.11711E+04, &
         0.12132E+04, 0.12556E+04, 0.12984E+04, 0.13416E+04, 0.13851E+04, &
         0.14291E+04, 0.14734E+04, 0.15182E+04, 0.15635E+04, 0.16091E+04, &
         0.16553E+04, 0.17019E+04, 0.17490E+04, 0.17965E+04, 0.18446E+04, &
         0.18932E+04, 0.19422E+04, 0.19918E+04, 0.20419E+04, 0.20926E+04, &
         0.21437E+04, 0.21954E+04, 0.22477E+04, 0.23004E+04, 0.23538E+04, &
         0.24077E+04, 0.24621E+04, 0.25171E+04, 0.25727E+04, 0.26288E+04, &
         0.26856E+04, 0.27428E+04, 0.28007E+04, 0.28591E+04, 0.29181E+04, &
         0.29777E+04, 0.30379E+04, 0.30986E+04, 0.31600E+04, 0.32219E+04, &
         0.32844E+04, 0.33475E+04, 0.34112E+04, 0.34755E+04, 0.35404E+04, &
         0.36059E+04, 0.36720E+04, 0.37387E+04, 0.38060E+04, 0.38739E+04, &
         0.39424E+04, 0.40115E+04, 0.40812E+04, 0.41515E+04, 0.42224E+04, &
         0.42939E+04, 0.43661E+04, 0.44388E+04, 0.45122E+04, 0.45861E+04, &
         0.46607E+04, 0.47359E+04, 0.48117E+04, 0.48882E+04, 0.49652E+04, &
         0.50428E+04, 0.51211E+04, 0.52000E+04, 0.52795E+04, 0.53596E+04, &
         0.54404E+04, 0.55217E+04, 0.56037E+04, 0.56863E+04, 0.57695E+04, &
         0.58533E+04, 0.59378E+04, 0.60229E+04, 0.61086E+04, 0.61950E+04, &
         0.62819E+04, 0.63695E+04, 0.64577E+04, 0.65465E+04, 0.66360E+04, &
         0.67261E+04, 0.68168E+04, 0.69081E+04, 0.70001E+04, 0.70927E+04, &
         0.71859E+04/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_N2


!
!     *****************
      SUBROUTINE QT_HCN(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(3) :: XGJ
      REAL(DOUBLE), DIMENSION(3,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 6., 12., 4./
!...      HCN
!...        --       124
      DATA (QOFT(1,J),J=1,119)/ 0.17143E+03, 0.24209E+03, 0.31285E+03, &
         0.38392E+03, 0.45582E+03, 0.52929E+03, 0.60515E+03, 0.68424E+03, &
         0.76731E+03, 0.85505E+03, 0.94805E+03, 0.10468E+04, 0.11519E+04, &
         0.12637E+04, 0.13826E+04, 0.15090E+04, 0.16435E+04, 0.17863E+04, &
         0.19378E+04, 0.20985E+04, 0.22689E+04, 0.24492E+04, 0.26401E+04, &
         0.28418E+04, 0.30550E+04, 0.32801E+04, 0.35176E+04, 0.37680E+04, &
         0.40318E+04, 0.43097E+04, 0.46021E+04, 0.49097E+04, 0.52330E+04, &
         0.55727E+04, 0.59294E+04, 0.63038E+04, 0.66964E+04, 0.71081E+04, &
         0.75396E+04, 0.79915E+04, 0.84646E+04, 0.89596E+04, 0.94774E+04, &
         0.10019E+05, 0.10585E+05, 0.11176E+05, 0.11793E+05, 0.12437E+05, &
         0.13108E+05, 0.13809E+05, 0.14540E+05, 0.15301E+05, 0.16094E+05, &
         0.16919E+05, 0.17779E+05, 0.18673E+05, 0.19603E+05, 0.20570E+05, &
         0.21575E+05, 0.22619E+05, 0.23704E+05, 0.24831E+05, 0.26000E+05, &
         0.27213E+05, 0.28472E+05, 0.29778E+05, 0.31131E+05, 0.32534E+05, &
         0.33987E+05, 0.35493E+05, 0.37052E+05, 0.38666E+05, 0.40336E+05, &
         0.42064E+05, 0.43852E+05, 0.45701E+05, 0.47612E+05, 0.49587E+05, &
         0.51629E+05, 0.53738E+05, 0.55916E+05, 0.58165E+05, 0.60486E+05, &
         0.62883E+05, 0.65355E+05, 0.67905E+05, 0.70536E+05, 0.73249E+05, &
         0.76045E+05, 0.78927E+05, 0.81897E+05, 0.84957E+05, 0.88108E+05, &
         0.91354E+05, 0.94696E+05, 0.98136E+05, 0.10168E+06, 0.10532E+06, &
         0.10907E+06, 0.11292E+06, 0.11689E+06, 0.12096E+06, 0.12516E+06, &
         0.12946E+06, 0.13389E+06, 0.13844E+06, 0.14311E+06, 0.14791E+06, &
         0.15284E+06, 0.15790E+06, 0.16310E+06, 0.16843E+06, 0.17391E+06, &
         0.17953E+06, 0.18529E+06, 0.19120E+06, 0.19726E+06, 0.20348E+06, &
         0.20986E+06/
!...        --       134
      DATA (QOFT(2,J),J=1,119)/ 0.35186E+03, 0.49693E+03, 0.64221E+03, &
         0.78815E+03, 0.93585E+03, 0.10868E+04, 0.12428E+04, 0.14056E+04, &
         0.15766E+04, 0.17574E+04, 0.19491E+04, 0.21528E+04, 0.23695E+04, &
         0.26002E+04, 0.28457E+04, 0.31068E+04, 0.33845E+04, 0.36795E+04, &
         0.39926E+04, 0.43249E+04, 0.46770E+04, 0.50500E+04, 0.54447E+04, &
         0.58621E+04, 0.63032E+04, 0.67690E+04, 0.72606E+04, 0.77789E+04, &
         0.83252E+04, 0.89005E+04, 0.95062E+04, 0.10143E+05, 0.10813E+05, &
         0.11517E+05, 0.12256E+05, 0.13032E+05, 0.13846E+05, 0.14699E+05, &
         0.15593E+05, 0.16530E+05, 0.17511E+05, 0.18538E+05, 0.19612E+05, &
         0.20734E+05, 0.21908E+05, 0.23134E+05, 0.24414E+05, 0.25750E+05, &
         0.27145E+05, 0.28599E+05, 0.30115E+05, 0.31694E+05, 0.33340E+05, &
         0.35054E+05, 0.36838E+05, 0.38694E+05, 0.40625E+05, 0.42633E+05, &
         0.44720E+05, 0.46889E+05, 0.49142E+05, 0.51481E+05, 0.53910E+05, &
         0.56430E+05, 0.59045E+05, 0.61757E+05, 0.64568E+05, 0.67482E+05, &
         0.70502E+05, 0.73630E+05, 0.76869E+05, 0.80223E+05, 0.83694E+05, &
         0.87285E+05, 0.91000E+05, 0.94843E+05, 0.98815E+05, 0.10292E+06, &
         0.10716E+06, 0.11155E+06, 0.11608E+06, 0.12075E+06, 0.12558E+06, &
         0.13056E+06, 0.13570E+06, 0.14100E+06, 0.14647E+06, 0.15211E+06, &
         0.15793E+06, 0.16392E+06, 0.17009E+06, 0.17646E+06, 0.18301E+06, &
         0.18976E+06, 0.19671E+06, 0.20387E+06, 0.21123E+06, 0.21881E+06, &
         0.22660E+06, 0.23462E+06, 0.24287E+06, 0.25135E+06, 0.26007E+06, &
         0.26903E+06, 0.27824E+06, 0.28771E+06, 0.29743E+06, 0.30742E+06, &
         0.31767E+06, 0.32820E+06, 0.33901E+06, 0.35011E+06, 0.36150E+06, &
         0.37319E+06, 0.38518E+06, 0.39749E+06, 0.41010E+06, 0.42304E+06, &
         0.43631E+06/
!...        --       125
      DATA (QOFT(3,J),J=1,119)/ 0.11863E+03, 0.16755E+03, 0.21653E+03, &
         0.26576E+03, 0.31559E+03, 0.36656E+03, 0.41926E+03, 0.47428E+03, &
         0.53214E+03, 0.59333E+03, 0.65824E+03, 0.72727E+03, 0.80074E+03, &
         0.87898E+03, 0.96227E+03, 0.10509E+04, 0.11452E+04, 0.12454E+04, &
         0.13518E+04, 0.14647E+04, 0.15844E+04, 0.17112E+04, 0.18455E+04, &
         0.19875E+04, 0.21377E+04, 0.22962E+04, 0.24636E+04, 0.26402E+04, &
         0.28263E+04, 0.30224E+04, 0.32289E+04, 0.34461E+04, 0.36745E+04, &
         0.39145E+04, 0.41667E+04, 0.44314E+04, 0.47092E+04, 0.50005E+04, &
         0.53059E+04, 0.56259E+04, 0.59609E+04, 0.63116E+04, 0.66785E+04, &
         0.70622E+04, 0.74633E+04, 0.78823E+04, 0.83200E+04, 0.87769E+04, &
         0.92536E+04, 0.97509E+04, 0.10269E+05, 0.10810E+05, 0.11373E+05, &
         0.11959E+05, 0.12570E+05, 0.13205E+05, 0.13866E+05, 0.14554E+05, &
         0.15268E+05, 0.16011E+05, 0.16782E+05, 0.17583E+05, 0.18415E+05, &
         0.19279E+05, 0.20174E+05, 0.21103E+05, 0.22067E+05, 0.23065E+05, &
         0.24100E+05, 0.25172E+05, 0.26282E+05, 0.27432E+05, 0.28622E+05, &
         0.29853E+05, 0.31127E+05, 0.32445E+05, 0.33807E+05, 0.35215E+05, &
         0.36670E+05, 0.38174E+05, 0.39727E+05, 0.41330E+05, 0.42986E+05, &
         0.44695E+05, 0.46459E+05, 0.48278E+05, 0.50155E+05, 0.52091E+05, &
         0.54086E+05, 0.56143E+05, 0.58263E+05, 0.60447E+05, 0.62696E+05, &
         0.65013E+05, 0.67399E+05, 0.69856E+05, 0.72384E+05, 0.74986E+05, &
         0.77663E+05, 0.80416E+05, 0.83249E+05, 0.86161E+05, 0.89156E+05, &
         0.92233E+05, 0.95397E+05, 0.98648E+05, 0.10199E+06, 0.10542E+06, &
         0.10894E+06, 0.11256E+06, 0.11627E+06, 0.12009E+06, 0.12400E+06, &
         0.12802E+06, 0.13214E+06, 0.13636E+06, 0.14070E+06, 0.14515E+06, &
         0.14971E+06/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_HCN


!
!     *****************
      SUBROUTINE QT_CH3CL(T, ISO, GSI, QT)       ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(2) :: XGJ
      REAL(DOUBLE), DIMENSION(2,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 4., 4./
!...    CH3Cl
!...        --       215
      DATA (QOFT(1,J),J=1,119)/ 0.10106E+05, 0.17025E+05, 0.25055E+05, &
         0.34073E+05, 0.44011E+05, 0.54858E+05, 0.66651E+05, 0.79468E+05, &
         0.93426E+05, 0.10867E+06, 0.12538E+06, 0.14375E+06, 0.16401E+06, &
         0.18641E+06, 0.21121E+06, 0.23871E+06, 0.26925E+06, 0.30317E+06, &
         0.34086E+06, 0.38274E+06, 0.42928E+06, 0.48098E+06, 0.53840E+06, &
         0.60213E+06, 0.67284E+06, 0.75125E+06, 0.83814E+06, 0.93438E+06, &
         0.10409E+07, 0.11587E+07, 0.12890E+07, 0.14328E+07, 0.15916E+07, &
         0.17668E+07, 0.19599E+07, 0.21727E+07, 0.24069E+07, 0.26645E+07, &
         0.29478E+07, 0.32590E+07, 0.36006E+07, 0.39754E+07, 0.43864E+07, &
         0.48367E+07, 0.53297E+07, 0.58692E+07, 0.64591E+07, 0.71038E+07, &
         0.78078E+07, 0.85762E+07, 0.94143E+07, 0.10328E+08, 0.11323E+08, &
         0.12406E+08, 0.13585E+08, 0.14867E+08, 0.16260E+08, 0.17772E+08, &
         0.19414E+08, 0.21195E+08, 0.23126E+08, 0.25218E+08, 0.27484E+08, &
         0.29936E+08, 0.32589E+08, 0.35456E+08, 0.38555E+08, 0.41901E+08, &
         0.45513E+08, 0.49409E+08, 0.53610E+08, 0.58137E+08, 0.63014E+08, &
         0.68264E+08, 0.73913E+08, 0.79989E+08, 0.86521E+08, 0.93539E+08, &
         0.10108E+09, 0.10917E+09, 0.11785E+09, 0.12716E+09, 0.13714E+09, &
         0.14783E+09, 0.15928E+09, 0.17154E+09, 0.18466E+09, 0.19869E+09, &
         0.21369E+09, 0.22972E+09, 0.24685E+09, 0.26513E+09, 0.28465E+09, &
         0.30547E+09, 0.32768E+09, 0.35136E+09, 0.37659E+09, 0.40346E+09, &
         0.43208E+09, 0.46254E+09, 0.49495E+09, 0.52943E+09, 0.56608E+09, &
         0.60504E+09, 0.64643E+09, 0.69039E+09, 0.73707E+09, 0.78660E+09, &
         0.83916E+09, 0.89491E+09, 0.95401E+09, 0.10167E+10, 0.10830E+10, &
         0.11533E+10, 0.12278E+10, 0.13066E+10, 0.13900E+10, 0.14782E+10, &
         0.15715E+10/
!...        --       217
      DATA (QOFT(2,J),J=1,119)/ 0.10265E+05, 0.17294E+05, 0.25452E+05, &
         0.34612E+05, 0.44707E+05, 0.55726E+05, 0.67706E+05, 0.80727E+05, &
         0.94906E+05, 0.11039E+06, 0.12737E+06, 0.14603E+06, 0.16661E+06, &
         0.18936E+06, 0.21456E+06, 0.24250E+06, 0.27352E+06, 0.30798E+06, &
         0.34626E+06, 0.38881E+06, 0.43609E+06, 0.48861E+06, 0.54694E+06, &
         0.61168E+06, 0.68351E+06, 0.76316E+06, 0.85144E+06, 0.94920E+06, &
         0.10574E+07, 0.11771E+07, 0.13094E+07, 0.14556E+07, 0.16169E+07, &
         0.17949E+07, 0.19910E+07, 0.22072E+07, 0.24451E+07, 0.27068E+07, &
         0.29946E+07, 0.33107E+07, 0.36578E+07, 0.40385E+07, 0.44560E+07, &
         0.49134E+07, 0.54143E+07, 0.59624E+07, 0.65617E+07, 0.72166E+07, &
         0.79318E+07, 0.87124E+07, 0.95638E+07, 0.10492E+08, 0.11503E+08, &
         0.12603E+08, 0.13801E+08, 0.15103E+08, 0.16518E+08, 0.18055E+08, &
         0.19723E+08, 0.21532E+08, 0.23494E+08, 0.25619E+08, 0.27921E+08, &
         0.30412E+08, 0.33106E+08, 0.36020E+08, 0.39167E+08, 0.42567E+08, &
         0.46236E+08, 0.50194E+08, 0.54462E+08, 0.59061E+08, 0.64015E+08, &
         0.69349E+08, 0.75088E+08, 0.81261E+08, 0.87896E+08, 0.95026E+08, &
         0.10268E+09, 0.11090E+09, 0.11972E+09, 0.12918E+09, 0.13932E+09, &
         0.15018E+09, 0.16181E+09, 0.17427E+09, 0.18759E+09, 0.20185E+09, &
         0.21709E+09, 0.23337E+09, 0.25077E+09, 0.26935E+09, 0.28918E+09, &
         0.31033E+09, 0.33289E+09, 0.35695E+09, 0.38258E+09, 0.40988E+09, &
         0.43895E+09, 0.46990E+09, 0.50283E+09, 0.53785E+09, 0.57509E+09, &
         0.61467E+09, 0.65672E+09, 0.70138E+09, 0.74880E+09, 0.79912E+09, &
         0.85252E+09, 0.90915E+09, 0.96919E+09, 0.10328E+10, 0.11003E+10, &
         0.11717E+10, 0.12473E+10, 0.13274E+10, 0.14121E+10, 0.15017E+10, &
         0.15965E+10/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_CH3CL


!
!     *****************
      SUBROUTINE QT_H2O2(T, ISO, GSI, QT)        ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(1) :: XGJ
      REAL(DOUBLE), DIMENSION(1,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1./
!...     H2O2
!...        --      1661
      DATA (QOFT(1,J),J=1,119)/ 0.62392E+03, 0.10958E+04, 0.16692E+04, &
         0.23492E+04, 0.31427E+04, 0.40574E+04, 0.51014E+04, 0.62840E+04, &
         0.76157E+04, 0.91085E+04, 0.10776E+05, 0.12633E+05, 0.14696E+05, &
         0.16983E+05, 0.19515E+05, 0.22312E+05, 0.25396E+05, 0.28792E+05, &
         0.32526E+05, 0.36625E+05, 0.41118E+05, 0.46036E+05, 0.51410E+05, &
         0.57275E+05, 0.63667E+05, 0.70623E+05, 0.78185E+05, 0.86394E+05, &
         0.95295E+05, 0.10493E+06, 0.11536E+06, 0.12662E+06, 0.13878E+06, &
         0.15188E+06, 0.16600E+06, 0.18118E+06, 0.19750E+06, 0.21503E+06, &
         0.23383E+06, 0.25398E+06, 0.27556E+06, 0.29864E+06, 0.32333E+06, &
         0.34970E+06, 0.37784E+06, 0.40786E+06, 0.43985E+06, 0.47392E+06, &
         0.51018E+06, 0.54874E+06, 0.58972E+06, 0.63324E+06, 0.67943E+06, &
         0.72843E+06, 0.78037E+06, 0.83540E+06, 0.89366E+06, 0.95530E+06, &
         0.10205E+07, 0.10894E+07, 0.11622E+07, 0.12391E+07, 0.13202E+07, &
         0.14057E+07, 0.14959E+07, 0.15909E+07, 0.16910E+07, 0.17963E+07, &
         0.19072E+07, 0.20237E+07, 0.21463E+07, 0.22750E+07, 0.24102E+07, &
         0.25522E+07, 0.27012E+07, 0.28575E+07, 0.30213E+07, 0.31931E+07, &
         0.33730E+07, 0.35615E+07, 0.37588E+07, 0.39653E+07, 0.41813E+07, &
         0.44072E+07, 0.46433E+07, 0.48901E+07, 0.51479E+07, 0.54171E+07, &
         0.56982E+07, 0.59915E+07, 0.62976E+07, 0.66167E+07, 0.69495E+07, &
         0.72963E+07, 0.76577E+07, 0.80342E+07, 0.84262E+07, 0.88343E+07, &
         0.92591E+07, 0.97011E+07, 0.10161E+08, 0.10639E+08, 0.11136E+08, &
         0.11652E+08, 0.12189E+08, 0.12746E+08, 0.13325E+08, 0.13926E+08, &
         0.14550E+08, 0.15198E+08, 0.15870E+08, 0.16566E+08, 0.17289E+08, &
         0.18038E+08, 0.18814E+08, 0.19619E+08, 0.20452E+08, 0.21315E+08, &
         0.22209E+08/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_H2O2


!
!     *****************
      SUBROUTINE QT_C2H2(T, ISO, GSI, QT)        ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(2) :: XGJ
      REAL(DOUBLE), DIMENSION(2,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1., 8./
!...     C2H2
!...        --      1221
      DATA (QOFT(1,J),J=1,119)/ 0.71617E+02, 0.10121E+03, 0.13092E+03, &
         0.16104E+03, 0.19218E+03, 0.22509E+03, 0.26062E+03, 0.29959E+03, &
         0.34281E+03, 0.39103E+03, 0.44503E+03, 0.50558E+03, 0.57346E+03, &
         0.64950E+03, 0.73457E+03, 0.82960E+03, 0.93557E+03, 0.10535E+04, &
         0.11846E+04, 0.13301E+04, 0.14911E+04, 0.16692E+04, 0.18658E+04, &
         0.20825E+04, 0.23211E+04, 0.25833E+04, 0.28711E+04, 0.31867E+04, &
         0.35323E+04, 0.39102E+04, 0.43230E+04, 0.47735E+04, 0.52645E+04, &
         0.57991E+04, 0.63807E+04, 0.70127E+04, 0.76988E+04, 0.84430E+04, &
         0.92495E+04, 0.10123E+05, 0.11067E+05, 0.12088E+05, 0.13191E+05, &
         0.14381E+05, 0.15664E+05, 0.17047E+05, 0.18536E+05, 0.20137E+05, &
         0.21859E+05, 0.23710E+05, 0.25696E+05, 0.27827E+05, 0.30112E+05, &
         0.32561E+05, 0.35183E+05, 0.37990E+05, 0.40991E+05, 0.44199E+05, &
         0.47626E+05, 0.51285E+05, 0.55189E+05, 0.59353E+05, 0.63791E+05, &
         0.68518E+05, 0.73551E+05, 0.78908E+05, 0.84604E+05, 0.90661E+05, &
         0.97095E+05, 0.10393E+06, 0.11118E+06, 0.11888E+06, 0.12704E+06, &
         0.13569E+06, 0.14486E+06, 0.15457E+06, 0.16485E+06, 0.17572E+06, &
         0.18722E+06, 0.19938E+06, 0.21223E+06, 0.22581E+06, 0.24014E+06, &
         0.25527E+06, 0.27123E+06, 0.28807E+06, 0.30582E+06, 0.32452E+06, &
         0.34423E+06, 0.36498E+06, 0.38683E+06, 0.40982E+06, 0.43401E+06, &
         0.45944E+06, 0.48618E+06, 0.51428E+06, 0.54380E+06, 0.57480E+06, &
         0.60735E+06, 0.64151E+06, 0.67735E+06, 0.71495E+06, 0.75436E+06, &
         0.79568E+06, 0.83898E+06, 0.88434E+06, 0.93184E+06, 0.98158E+06, &
         0.10336E+07, 0.10881E+07, 0.11451E+07, 0.12047E+07, 0.12670E+07, &
         0.13321E+07, 0.14002E+07, 0.14713E+07, 0.15455E+07, 0.16231E+07, &
         0.17040E+07/
!...        --      1231
      DATA (QOFT(2,J),J=1,119)/ 0.28647E+03, 0.40486E+03, 0.52369E+03, &
         0.64419E+03, 0.76874E+03, 0.90040E+03, 0.10425E+04, 0.11984E+04, &
         0.13713E+04, 0.15642E+04, 0.17802E+04, 0.20223E+04, 0.22939E+04, &
         0.25981E+04, 0.29384E+04, 0.33185E+04, 0.37424E+04, 0.42142E+04, &
         0.47386E+04, 0.53203E+04, 0.59646E+04, 0.66769E+04, 0.74634E+04, &
         0.83302E+04, 0.92845E+04, 0.10333E+05, 0.11485E+05, 0.12747E+05, &
         0.14129E+05, 0.15641E+05, 0.17292E+05, 0.19094E+05, 0.21058E+05, &
         0.23197E+05, 0.25523E+05, 0.28051E+05, 0.30796E+05, 0.33773E+05, &
         0.36999E+05, 0.40492E+05, 0.44270E+05, 0.48354E+05, 0.52765E+05, &
         0.57525E+05, 0.62658E+05, 0.68189E+05, 0.74144E+05, 0.80551E+05, &
         0.87439E+05, 0.94840E+05, 0.10279E+06, 0.11131E+06, 0.12045E+06, &
         0.13025E+06, 0.14074E+06, 0.15196E+06, 0.16397E+06, 0.17680E+06, &
         0.19051E+06, 0.20514E+06, 0.22076E+06, 0.23742E+06, 0.25517E+06, &
         0.27408E+06, 0.29421E+06, 0.31564E+06, 0.33842E+06, 0.36265E+06, &
         0.38839E+06, 0.41572E+06, 0.44474E+06, 0.47553E+06, 0.50818E+06, &
         0.54278E+06, 0.57945E+06, 0.61829E+06, 0.65940E+06, 0.70289E+06, &
         0.74890E+06, 0.79754E+06, 0.84894E+06, 0.90324E+06, 0.96057E+06, &
         0.10211E+07, 0.10849E+07, 0.11523E+07, 0.12233E+07, 0.12981E+07, &
         0.13769E+07, 0.14599E+07, 0.15473E+07, 0.16393E+07, 0.17361E+07, &
         0.18378E+07, 0.19447E+07, 0.20571E+07, 0.21752E+07, 0.22992E+07, &
         0.24294E+07, 0.25661E+07, 0.27094E+07, 0.28598E+07, 0.30175E+07, &
         0.31828E+07, 0.33560E+07, 0.35374E+07, 0.37274E+07, 0.39264E+07, &
         0.41346E+07, 0.43525E+07, 0.45805E+07, 0.48188E+07, 0.50681E+07, &
         0.53286E+07, 0.56008E+07, 0.58852E+07, 0.61823E+07, 0.64924E+07, &
         0.68162E+07/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_C2H2


!
!     *****************
      SUBROUTINE QT_C2H6(T, ISO, GSI, QT)        ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

! HIT 08 second isotope
!	add index 2 to xgj
!  duplicate qoft for second isotope

      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(2) :: XGJ
      REAL(DOUBLE), DIMENSION(2,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1., 1./
!...     C2H6
!...        --      1221
      DATA (QOFT(1,J),J=1,119)/ 0.47267E+04, 0.80011E+04, 0.11928E+05, &
         0.16564E+05, 0.21985E+05, 0.28287E+05, 0.35590E+05, 0.44049E+05, &
         0.53862E+05, 0.65277E+05, 0.78597E+05, 0.94191E+05, 0.11250E+06, &
         0.13407E+06, 0.15952E+06, 0.18962E+06, 0.22526E+06, 0.26751E+06, &
         0.31763E+06, 0.37714E+06, 0.44780E+06, 0.53174E+06, 0.63145E+06, &
         0.74989E+06, 0.89056E+06, 0.10576E+07, 0.12559E+07, 0.14912E+07, &
         0.17704E+07, 0.21013E+07, 0.24936E+07, 0.29582E+07, 0.35083E+07, &
         0.41591E+07, 0.49286E+07, 0.58379E+07, 0.69116E+07, 0.81787E+07, &
         0.96728E+07, 0.11433E+08, 0.13506E+08, 0.15945E+08, 0.18812E+08, &
         0.22180E+08, 0.26134E+08, 0.30770E+08, 0.36204E+08, 0.42565E+08, &
         0.50008E+08, 0.58708E+08, 0.68868E+08, 0.80725E+08, 0.94548E+08, &
         0.11065E+09, 0.12940E+09, 0.15119E+09, 0.17652E+09, 0.20593E+09, &
         0.24003E+09, 0.27956E+09, 0.32533E+09, 0.37829E+09, 0.43951E+09, &
         0.51021E+09, 0.59180E+09, 0.68588E+09, 0.79427E+09, 0.91904E+09, &
         0.10625E+10, 0.12275E+10, 0.14168E+10, 0.16341E+10, 0.18831E+10, &
         0.21684E+10, 0.24949E+10, 0.28684E+10, 0.32951E+10, 0.37823E+10, &
         0.43382E+10, 0.49719E+10, 0.56938E+10, 0.65156E+10, 0.74502E+10, &
         0.85125E+10, 0.97190E+10, 0.11088E+11, 0.12641E+11, 0.14401E+11, &
         0.16393E+11, 0.18648E+11, 0.21198E+11, 0.24079E+11, 0.27332E+11, &
         0.31003E+11, 0.35142E+11, 0.39807E+11, 0.45060E+11, 0.50972E+11, &
         0.57620E+11, 0.65091E+11, 0.73483E+11, 0.82902E+11, 0.93467E+11, &
         0.10531E+12, 0.11858E+12, 0.13343E+12, 0.15005E+12, 0.16864E+12, &
         0.18941E+12, 0.21260E+12, 0.23849E+12, 0.26737E+12, 0.29957E+12, &
         0.33545E+12, 0.37541E+12, 0.41987E+12, 0.46934E+12, 0.52432E+12, &
         0.58542E+12/

      DATA (QOFT(2,J),J=1,119)/ 0.47267E+04, 0.80011E+04, 0.11928E+05, &
         0.16564E+05, 0.21985E+05, 0.28287E+05, 0.35590E+05, 0.44049E+05, &
         0.53862E+05, 0.65277E+05, 0.78597E+05, 0.94191E+05, 0.11250E+06, &
         0.13407E+06, 0.15952E+06, 0.18962E+06, 0.22526E+06, 0.26751E+06, &
         0.31763E+06, 0.37714E+06, 0.44780E+06, 0.53174E+06, 0.63145E+06, &
         0.74989E+06, 0.89056E+06, 0.10576E+07, 0.12559E+07, 0.14912E+07, &
         0.17704E+07, 0.21013E+07, 0.24936E+07, 0.29582E+07, 0.35083E+07, &
         0.41591E+07, 0.49286E+07, 0.58379E+07, 0.69116E+07, 0.81787E+07, &
         0.96728E+07, 0.11433E+08, 0.13506E+08, 0.15945E+08, 0.18812E+08, &
         0.22180E+08, 0.26134E+08, 0.30770E+08, 0.36204E+08, 0.42565E+08, &
         0.50008E+08, 0.58708E+08, 0.68868E+08, 0.80725E+08, 0.94548E+08, &
         0.11065E+09, 0.12940E+09, 0.15119E+09, 0.17652E+09, 0.20593E+09, &
         0.24003E+09, 0.27956E+09, 0.32533E+09, 0.37829E+09, 0.43951E+09, &
         0.51021E+09, 0.59180E+09, 0.68588E+09, 0.79427E+09, 0.91904E+09, &
         0.10625E+10, 0.12275E+10, 0.14168E+10, 0.16341E+10, 0.18831E+10, &
         0.21684E+10, 0.24949E+10, 0.28684E+10, 0.32951E+10, 0.37823E+10, &
         0.43382E+10, 0.49719E+10, 0.56938E+10, 0.65156E+10, 0.74502E+10, &
         0.85125E+10, 0.97190E+10, 0.11088E+11, 0.12641E+11, 0.14401E+11, &
         0.16393E+11, 0.18648E+11, 0.21198E+11, 0.24079E+11, 0.27332E+11, &
         0.31003E+11, 0.35142E+11, 0.39807E+11, 0.45060E+11, 0.50972E+11, &
         0.57620E+11, 0.65091E+11, 0.73483E+11, 0.82902E+11, 0.93467E+11, &
         0.10531E+12, 0.11858E+12, 0.13343E+12, 0.15005E+12, 0.16864E+12, &
         0.18941E+12, 0.21260E+12, 0.23849E+12, 0.26737E+12, 0.29957E+12, &
         0.33545E+12, 0.37541E+12, 0.41987E+12, 0.46934E+12, 0.52432E+12, &
         0.58542E+12/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_C2H6


!
!     *****************
      SUBROUTINE QT_PH3(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(1) :: XGJ
      REAL(DOUBLE), DIMENSION(1,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 2./
!...      PH3
!...        --      1111
      DATA (QOFT(1,J),J=1,119)/ 0.29652E+03, 0.49643E+03, 0.72810E+03, &
         0.98777E+03, 0.12729E+04, 0.15820E+04, 0.19145E+04, 0.22708E+04, &
         0.26520E+04, 0.30600E+04, 0.34971E+04, 0.39662E+04, 0.44702E+04, &
         0.50126E+04, 0.55970E+04, 0.62273E+04, 0.69075E+04, 0.76421E+04, &
         0.84357E+04, 0.92933E+04, 0.10220E+05, 0.11222E+05, 0.12304E+05, &
         0.13473E+05, 0.14736E+05, 0.16099E+05, 0.17571E+05, 0.19160E+05, &
         0.20873E+05, 0.22720E+05, 0.24710E+05, 0.26854E+05, 0.29162E+05, &
         0.31646E+05, 0.34317E+05, 0.37188E+05, 0.40273E+05, 0.43585E+05, &
         0.47140E+05, 0.50953E+05, 0.55040E+05, 0.59419E+05, 0.64108E+05, &
         0.69127E+05, 0.74496E+05, 0.80236E+05, 0.86369E+05, 0.92918E+05, &
         0.99909E+05, 0.10737E+06, 0.11532E+06, 0.12380E+06, 0.13282E+06, &
         0.14244E+06, 0.15266E+06, 0.16354E+06, 0.17511E+06, 0.18739E+06, &
         0.20044E+06, 0.21430E+06, 0.22900E+06, 0.24459E+06, 0.26111E+06, &
         0.27862E+06, 0.29716E+06, 0.31680E+06, 0.33757E+06, 0.35954E+06, &
         0.38277E+06, 0.40733E+06, 0.43326E+06, 0.46065E+06, 0.48955E+06, &
         0.52005E+06, 0.55222E+06, 0.58614E+06, 0.62188E+06, 0.65953E+06, &
         0.69917E+06, 0.74091E+06, 0.78483E+06, 0.83103E+06, 0.87960E+06, &
         0.93067E+06, 0.98432E+06, 0.10407E+07, 0.10999E+07, 0.11620E+07, &
         0.12272E+07, 0.12956E+07, 0.13673E+07, 0.14425E+07, 0.15212E+07, &
         0.16038E+07, 0.16902E+07, 0.17808E+07, 0.18755E+07, 0.19746E+07, &
         0.20784E+07, 0.21868E+07, 0.23002E+07, 0.24187E+07, 0.25425E+07, &
         0.26719E+07, 0.28070E+07, 0.29480E+07, 0.30952E+07, 0.32488E+07, &
         0.34091E+07, 0.35762E+07, 0.37504E+07, 0.39320E+07, 0.41213E+07, &
         0.43185E+07, 0.45239E+07, 0.47378E+07, 0.49605E+07, 0.51923E+07, &
         0.54335E+07/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_PH3


!
!     *****************
      SUBROUTINE QT_COF2(T, ISO, GSI, QT)        ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(1) :: XGJ
      REAL(DOUBLE), DIMENSION(1,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1./
!...     COF2
!...        --       269
      DATA (QOFT(1,J),J=1,119)/ 0.54999E+04, 0.92749E+04, 0.13668E+05, &
         0.18643E+05, 0.24224E+05, 0.30487E+05, 0.37547E+05, 0.45543E+05, &
         0.54639E+05, 0.65019E+05, 0.76886E+05, 0.90462E+05, 0.10600E+06, &
         0.12377E+06, 0.14407E+06, 0.16723E+06, 0.19363E+06, 0.22367E+06, &
         0.25780E+06, 0.29650E+06, 0.34031E+06, 0.38982E+06, 0.44568E+06, &
         0.50859E+06, 0.57932E+06, 0.65872E+06, 0.74770E+06, 0.84724E+06, &
         0.95844E+06, 0.10825E+07, 0.12205E+07, 0.13741E+07, 0.15446E+07, &
         0.17336E+07, 0.19428E+07, 0.21742E+07, 0.24296E+07, 0.27113E+07, &
         0.30214E+07, 0.33626E+07, 0.37373E+07, 0.41484E+07, 0.45989E+07, &
         0.50921E+07, 0.56313E+07, 0.62202E+07, 0.68626E+07, 0.75628E+07, &
         0.83251E+07, 0.91542E+07, 0.10055E+08, 0.11033E+08, 0.12093E+08, &
         0.13242E+08, 0.14486E+08, 0.15831E+08, 0.17284E+08, 0.18853E+08, &
         0.20546E+08, 0.22371E+08, 0.24335E+08, 0.26450E+08, 0.28724E+08, &
         0.31167E+08, 0.33790E+08, 0.36605E+08, 0.39623E+08, 0.42856E+08, &
         0.46318E+08, 0.50022E+08, 0.53983E+08, 0.58215E+08, 0.62735E+08, &
         0.67558E+08, 0.72702E+08, 0.78186E+08, 0.84028E+08, 0.90247E+08, &
         0.96865E+08, 0.10390E+09, 0.11138E+09, 0.11933E+09, 0.12777E+09, &
         0.13672E+09, 0.14622E+09, 0.15629E+09, 0.16695E+09, 0.17825E+09, &
         0.19021E+09, 0.20287E+09, 0.21625E+09, 0.23039E+09, 0.24534E+09, &
         0.26113E+09, 0.27779E+09, 0.29538E+09, 0.31392E+09, 0.33348E+09, &
         0.35409E+09, 0.37580E+09, 0.39867E+09, 0.42274E+09, 0.44806E+09, &
         0.47470E+09, 0.50271E+09, 0.53215E+09, 0.56308E+09, 0.59557E+09, &
         0.62968E+09, 0.66548E+09, 0.70304E+09, 0.74243E+09, 0.78374E+09, &
         0.82703E+09, 0.87240E+09, 0.91992E+09, 0.96967E+09, 0.10218E+10, &
         0.10763E+10/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_COF2


!
!     *****************
      SUBROUTINE QT_SF6(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(1) :: XGJ
      REAL(DOUBLE), DIMENSION(1,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1./
!...      SF6
!...        --        29
      DATA (QOFT(1,J),J=1,119)/ 0.46373E+05, 0.78844E+05, 0.11939E+06, &
         0.17183E+06, 0.24247E+06, 0.34059E+06, 0.47963E+06, 0.67906E+06, &
         0.96713E+06, 0.13848E+07, 0.19911E+07, 0.28714E+07, 0.41481E+07, &
         0.59956E+07, 0.86617E+07, 0.12496E+08, 0.17991E+08, 0.25832E+08, &
         0.36971E+08, 0.52724E+08, 0.74895E+08, 0.10595E+09, 0.14923E+09, &
         0.20925E+09, 0.29208E+09, 0.40582E+09, 0.56124E+09, 0.77259E+09, &
         0.10586E+10, 0.14439E+10, 0.19605E+10, 0.26500E+10, 0.35662E+10, &
         0.47781E+10, 0.63747E+10, 0.84689E+10, 0.11205E+11, 0.14765E+11, &
         0.19378E+11, 0.25336E+11, 0.32998E+11, 0.42819E+11, 0.55361E+11, &
         0.71323E+11, 0.91569E+11, 0.11716E+12, 0.14941E+12, 0.18992E+12, &
         0.24065E+12, 0.30398E+12, 0.38283E+12, 0.48069E+12, 0.60182E+12, &
         0.75136E+12, 0.93546E+12, 0.11615E+13, 0.14384E+13, 0.17767E+13, &
         0.21890E+13, 0.26903E+13, 0.32984E+13, 0.40344E+13, 0.49232E+13, &
         0.59942E+13, 0.72819E+13, 0.88272E+13, 0.10678E+14, 0.12889E+14, &
         0.15527E+14, 0.18666E+14, 0.22397E+14, 0.26823E+14, 0.32062E+14, &
         0.38253E+14, 0.45558E+14, 0.54161E+14, 0.64277E+14, 0.76153E+14, &
         0.90072E+14, 0.10636E+15, 0.12539E+15, 0.14759E+15, 0.17345E+15, &
         0.20354E+15, 0.23848E+15, 0.27902E+15, 0.32597E+15, 0.38028E+15, &
         0.44303E+15, 0.51542E+15, 0.59883E+15, 0.69482E+15, 0.80516E+15, &
         0.93182E+15, 0.10770E+16, 0.12434E+16, 0.14336E+16, 0.16511E+16, &
         0.18992E+16, 0.21821E+16, 0.25043E+16, 0.28709E+16, 0.32875E+16, &
         0.37604E+16, 0.42968E+16, 0.49046E+16, 0.55925E+16, 0.63704E+16, &
         0.72492E+16, 0.82411E+16, 0.93596E+16, 0.10620E+17, 0.12038E+17, &
         0.13633E+17, 0.15425E+17, 0.17438E+17, 0.19694E+17, 0.22224E+17, &
         0.25057E+17/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_SF6


!
!     *****************
      SUBROUTINE QT_H2S(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(3) :: XGJ
      REAL(DOUBLE), DIMENSION(3,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1., 1., 4./
!...      H2S
!...        --       121
      DATA (QOFT(1,J),J=1,119)/ 0.47192E+02, 0.78671E+02, 0.11510E+03, &
         0.15589E+03, 0.20061E+03, 0.24896E+03, 0.30070E+03, 0.35571E+03, &
         0.41386E+03, 0.47513E+03, 0.53951E+03, 0.60703E+03, 0.67772E+03, &
         0.75167E+03, 0.82896E+03, 0.90969E+03, 0.99396E+03, 0.10819E+04, &
         0.11736E+04, 0.12692E+04, 0.13689E+04, 0.14727E+04, 0.15809E+04, &
         0.16937E+04, 0.18111E+04, 0.19333E+04, 0.20606E+04, 0.21931E+04, &
         0.23309E+04, 0.24744E+04, 0.26236E+04, 0.27788E+04, 0.29403E+04, &
         0.31081E+04, 0.32825E+04, 0.34638E+04, 0.36522E+04, 0.38478E+04, &
         0.40510E+04, 0.42619E+04, 0.44808E+04, 0.47080E+04, 0.49437E+04, &
         0.51881E+04, 0.54415E+04, 0.57042E+04, 0.59764E+04, 0.62584E+04, &
         0.65505E+04, 0.68529E+04, 0.71660E+04, 0.74899E+04, 0.78251E+04, &
         0.81718E+04, 0.85303E+04, 0.89008E+04, 0.92838E+04, 0.96795E+04, &
         0.10088E+05, 0.10510E+05, 0.10946E+05, 0.11396E+05, 0.11860E+05, &
         0.12339E+05, 0.12833E+05, 0.13342E+05, 0.13867E+05, 0.14408E+05, &
         0.14966E+05, 0.15540E+05, 0.16132E+05, 0.16741E+05, 0.17368E+05, &
         0.18013E+05, 0.18677E+05, 0.19361E+05, 0.20064E+05, 0.20786E+05, &
         0.21529E+05, 0.22293E+05, 0.23078E+05, 0.23885E+05, 0.24714E+05, &
         0.25565E+05, 0.26439E+05, 0.27337E+05, 0.28258E+05, 0.29204E+05, &
         0.30174E+05, 0.31170E+05, 0.32191E+05, 0.33239E+05, 0.34313E+05, &
         0.35414E+05, 0.36543E+05, 0.37700E+05, 0.38886E+05, 0.40101E+05, &
         0.41346E+05, 0.42621E+05, 0.43926E+05, 0.45263E+05, 0.46631E+05, &
         0.48033E+05, 0.49466E+05, 0.50934E+05, 0.52435E+05, 0.53971E+05, &
         0.55542E+05, 0.57149E+05, 0.58792E+05, 0.60472E+05, 0.62190E+05, &
         0.63946E+05, 0.65740E+05, 0.67574E+05, 0.69448E+05, 0.71362E+05, &
         0.73318E+05/
!...        --       141
      DATA (QOFT(2,J),J=1,119)/ 0.47310E+02, 0.78869E+02, 0.11539E+03, &
         0.15628E+03, 0.20112E+03, 0.24959E+03, 0.30147E+03, 0.35661E+03, &
         0.41491E+03, 0.47634E+03, 0.54088E+03, 0.60857E+03, 0.67945E+03, &
         0.75359E+03, 0.83107E+03, 0.91201E+03, 0.99649E+03, 0.10846E+04, &
         0.11766E+04, 0.12724E+04, 0.13724E+04, 0.14765E+04, 0.15850E+04, &
         0.16980E+04, 0.18157E+04, 0.19382E+04, 0.20658E+04, 0.21987E+04, &
         0.23369E+04, 0.24807E+04, 0.26303E+04, 0.27859E+04, 0.29478E+04, &
         0.31160E+04, 0.32909E+04, 0.34727E+04, 0.36615E+04, 0.38576E+04, &
         0.40613E+04, 0.42728E+04, 0.44923E+04, 0.47200E+04, 0.49563E+04, &
         0.52013E+04, 0.54554E+04, 0.57188E+04, 0.59917E+04, 0.62744E+04, &
         0.65672E+04, 0.68704E+04, 0.71843E+04, 0.75090E+04, 0.78451E+04, &
         0.81926E+04, 0.85520E+04, 0.89236E+04, 0.93075E+04, 0.97042E+04, &
         0.10114E+05, 0.10537E+05, 0.10974E+05, 0.11425E+05, 0.11890E+05, &
         0.12370E+05, 0.12866E+05, 0.13376E+05, 0.13903E+05, 0.14445E+05, &
         0.15004E+05, 0.15580E+05, 0.16173E+05, 0.16784E+05, 0.17412E+05, &
         0.18059E+05, 0.18725E+05, 0.19410E+05, 0.20115E+05, 0.20839E+05, &
         0.21584E+05, 0.22350E+05, 0.23137E+05, 0.23946E+05, 0.24777E+05, &
         0.25630E+05, 0.26507E+05, 0.27407E+05, 0.28330E+05, 0.29278E+05, &
         0.30251E+05, 0.31249E+05, 0.32273E+05, 0.33324E+05, 0.34401E+05, &
         0.35505E+05, 0.36637E+05, 0.37797E+05, 0.38985E+05, 0.40204E+05, &
         0.41451E+05, 0.42729E+05, 0.44038E+05, 0.45379E+05, 0.46751E+05, &
         0.48155E+05, 0.49593E+05, 0.51064E+05, 0.52569E+05, 0.54109E+05, &
         0.55684E+05, 0.57295E+05, 0.58943E+05, 0.60627E+05, 0.62349E+05, &
         0.64109E+05, 0.65908E+05, 0.67747E+05, 0.69625E+05, 0.71544E+05, &
         0.73505E+05/
!...        --       131
      DATA (QOFT(3,J),J=1,119)/ 0.18901E+03, 0.31509E+03, 0.46102E+03, &
         0.62437E+03, 0.80349E+03, 0.99713E+03, 0.12044E+04, 0.14247E+04, &
         0.16576E+04, 0.19030E+04, 0.21609E+04, 0.24313E+04, 0.27145E+04, &
         0.30106E+04, 0.33202E+04, 0.36436E+04, 0.39811E+04, 0.43332E+04, &
         0.47005E+04, 0.50835E+04, 0.54827E+04, 0.58987E+04, 0.63321E+04, &
         0.67836E+04, 0.72538E+04, 0.77434E+04, 0.82532E+04, 0.87838E+04, &
         0.93360E+04, 0.99106E+04, 0.10508E+05, 0.11130E+05, 0.11777E+05, &
         0.12449E+05, 0.13147E+05, 0.13874E+05, 0.14628E+05, 0.15412E+05, &
         0.16225E+05, 0.17070E+05, 0.17947E+05, 0.18857E+05, 0.19801E+05, &
         0.20780E+05, 0.21795E+05, 0.22847E+05, 0.23937E+05, 0.25067E+05, &
         0.26236E+05, 0.27448E+05, 0.28702E+05, 0.29999E+05, 0.31342E+05, &
         0.32730E+05, 0.34166E+05, 0.35650E+05, 0.37184E+05, 0.38769E+05, &
         0.40406E+05, 0.42097E+05, 0.43842E+05, 0.45644E+05, 0.47503E+05, &
         0.49421E+05, 0.51399E+05, 0.53439E+05, 0.55542E+05, 0.57709E+05, &
         0.59942E+05, 0.62242E+05, 0.64611E+05, 0.67051E+05, 0.69563E+05, &
         0.72148E+05, 0.74808E+05, 0.77545E+05, 0.80360E+05, 0.83255E+05, &
         0.86232E+05, 0.89291E+05, 0.92435E+05, 0.95667E+05, 0.98986E+05, &
         0.10240E+06, 0.10590E+06, 0.10949E+06, 0.11318E+06, 0.11697E+06, &
         0.12086E+06, 0.12484E+06, 0.12893E+06, 0.13313E+06, 0.13743E+06, &
         0.14184E+06, 0.14637E+06, 0.15100E+06, 0.15575E+06, 0.16062E+06, &
         0.16560E+06, 0.17071E+06, 0.17594E+06, 0.18129E+06, 0.18677E+06, &
         0.19238E+06, 0.19813E+06, 0.20400E+06, 0.21002E+06, 0.21617E+06, &
         0.22246E+06, 0.22890E+06, 0.23548E+06, 0.24221E+06, 0.24909E+06, &
         0.25612E+06, 0.26331E+06, 0.27065E+06, 0.27816E+06, 0.28583E+06, &
         0.29366E+06/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_H2S


!
!     *****************
      SUBROUTINE QT_HCOOH(T, ISO, GSI, QT)       ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(1) :: XGJ
      REAL(DOUBLE), DIMENSION(1,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 4./
!...    HCOOH
!...        --       126
      DATA (QOFT(1,J),J=1,119)/ 0.31899E+04, 0.53773E+04, 0.79205E+04, &
         0.10792E+05, 0.13993E+05, 0.17550E+05, 0.21509E+05, 0.25930E+05, &
         0.30885E+05, 0.36460E+05, 0.42750E+05, 0.49864E+05, 0.57926E+05, &
         0.67071E+05, 0.77453E+05, 0.89243E+05, 0.10263E+06, 0.11783E+06, &
         0.13507E+06, 0.15462E+06, 0.17676E+06, 0.20183E+06, 0.23018E+06, &
         0.26221E+06, 0.29836E+06, 0.33911E+06, 0.38501E+06, 0.43664E+06, &
         0.49467E+06, 0.55981E+06, 0.63286E+06, 0.71470E+06, 0.80628E+06, &
         0.90865E+06, 0.10230E+07, 0.11505E+07, 0.12927E+07, 0.14509E+07, &
         0.16269E+07, 0.18225E+07, 0.20396E+07, 0.22804E+07, 0.25472E+07, &
         0.28425E+07, 0.31692E+07, 0.35301E+07, 0.39285E+07, 0.43681E+07, &
         0.48525E+07, 0.53858E+07, 0.59727E+07, 0.66178E+07, 0.73265E+07, &
         0.81042E+07, 0.89571E+07, 0.98918E+07, 0.10915E+08, 0.12035E+08, &
         0.13259E+08, 0.14597E+08, 0.16057E+08, 0.17650E+08, 0.19387E+08, &
         0.21279E+08, 0.23339E+08, 0.25579E+08, 0.28016E+08, 0.30663E+08, &
         0.33536E+08, 0.36655E+08, 0.40037E+08, 0.43701E+08, 0.47671E+08, &
         0.51967E+08, 0.56614E+08, 0.61639E+08, 0.67068E+08, 0.72930E+08, &
         0.79257E+08, 0.86082E+08, 0.93439E+08, 0.10137E+09, 0.10990E+09, &
         0.11909E+09, 0.12898E+09, 0.13960E+09, 0.15102E+09, 0.16329E+09, &
         0.17646E+09, 0.19059E+09, 0.20575E+09, 0.22200E+09, 0.23941E+09, &
         0.25806E+09, 0.27802E+09, 0.29938E+09, 0.32223E+09, 0.34666E+09, &
         0.37276E+09, 0.40064E+09, 0.43041E+09, 0.46218E+09, 0.49607E+09, &
         0.53221E+09, 0.57074E+09, 0.61179E+09, 0.65551E+09, 0.70206E+09, &
         0.75159E+09, 0.80430E+09, 0.86034E+09, 0.91992E+09, 0.98324E+09, &
         0.10505E+10, 0.11219E+10, 0.11977E+10, 0.12782E+10, 0.13635E+10, &
         0.14540E+10/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_HCOOH


!
!     *****************
      SUBROUTINE QT_HO2(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(1) :: XGJ
      REAL(DOUBLE), DIMENSION(1,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 2./
!...      HO2
!...        --       166
      DATA (QOFT(1,J),J=1,119)/ 0.39277E+03, 0.66062E+03, 0.97123E+03, &
         0.13194E+04, 0.17014E+04, 0.21148E+04, 0.25578E+04, 0.30296E+04, &
         0.35297E+04, 0.40585E+04, 0.46167E+04, 0.52055E+04, 0.58264E+04, &
         0.64809E+04, 0.71707E+04, 0.78978E+04, 0.86641E+04, 0.94715E+04, &
         0.10322E+05, 0.11218E+05, 0.12161E+05, 0.13154E+05, 0.14198E+05, &
         0.15296E+05, 0.16449E+05, 0.17661E+05, 0.18933E+05, 0.20267E+05, &
         0.21666E+05, 0.23133E+05, 0.24669E+05, 0.26277E+05, 0.27960E+05, &
         0.29720E+05, 0.31560E+05, 0.33482E+05, 0.35489E+05, 0.37584E+05, &
         0.39769E+05, 0.42048E+05, 0.44423E+05, 0.46898E+05, 0.49475E+05, &
         0.52157E+05, 0.54948E+05, 0.57850E+05, 0.60868E+05, 0.64003E+05, &
         0.67261E+05, 0.70643E+05, 0.74154E+05, 0.77797E+05, 0.81575E+05, &
         0.85492E+05, 0.89553E+05, 0.93760E+05, 0.98118E+05, 0.10263E+06, &
         0.10730E+06, 0.11213E+06, 0.11713E+06, 0.12230E+06, 0.12765E+06, &
         0.13317E+06, 0.13888E+06, 0.14478E+06, 0.15086E+06, 0.15715E+06, &
         0.16363E+06, 0.17032E+06, 0.17723E+06, 0.18434E+06, 0.19168E+06, &
         0.19924E+06, 0.20704E+06, 0.21506E+06, 0.22333E+06, 0.23185E+06, &
         0.24061E+06, 0.24963E+06, 0.25891E+06, 0.26846E+06, 0.27828E+06, &
         0.28838E+06, 0.29876E+06, 0.30943E+06, 0.32039E+06, 0.33166E+06, &
         0.34323E+06, 0.35512E+06, 0.36732E+06, 0.37985E+06, 0.39271E+06, &
         0.40590E+06, 0.41944E+06, 0.43333E+06, 0.44758E+06, 0.46219E+06, &
         0.47717E+06, 0.49252E+06, 0.50826E+06, 0.52439E+06, 0.54091E+06, &
         0.55784E+06, 0.57518E+06, 0.59293E+06, 0.61112E+06, 0.62973E+06, &
         0.64878E+06, 0.66828E+06, 0.68824E+06, 0.70866E+06, 0.72955E+06, &
         0.75091E+06, 0.77276E+06, 0.79511E+06, 0.81795E+06, 0.84131E+06, &
         0.86518E+06/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_HO2


!
!     *****************
      SUBROUTINE QT_O(T, ISO, GSI, QT)           ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(1) :: XGJ
      REAL(DOUBLE), DIMENSION(1,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 0./
!...        O
!...        --         6
      DATA (QOFT(1,J),J=1,119)/ 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_O


!
!     *****************
      SUBROUTINE QT_CLONO2(T, ISO, GSI, QT)      ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(2) :: XGJ
      REAL(DOUBLE), DIMENSION(2,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 12., 12./
!...   ClONO2
!...        --      5646
      DATA (QOFT(1,J),J=1,119)/ 0.11444E+06, 0.21121E+06, 0.34858E+06, &
         0.53934E+06, 0.80041E+06, 0.11539E+07, 0.16286E+07, 0.22614E+07, &
         0.30992E+07, 0.42015E+07, 0.56426E+07, 0.75152E+07, 0.99344E+07, &
         0.13042E+08, 0.17012E+08, 0.22058E+08, 0.28437E+08, 0.36463E+08, &
         0.46514E+08, 0.59042E+08, 0.74589E+08, 0.93801E+08, 0.11744E+09, &
         0.14643E+09, 0.18181E+09, 0.22486E+09, 0.27705E+09, 0.34009E+09, &
         0.41598E+09, 0.50705E+09, 0.61599E+09, 0.74590E+09, 0.90037E+09, &
         0.10835E+10, 0.13001E+10, 0.15554E+10, 0.18556E+10, 0.22079E+10, &
         0.26200E+10, 0.31012E+10, 0.36615E+10, 0.43126E+10, 0.50675E+10, &
         0.59409E+10, 0.69492E+10, 0.81110E+10, 0.94469E+10, 0.10980E+11, &
         0.12736E+11, 0.14745E+11, 0.17037E+11, 0.19649E+11, 0.22620E+11, &
         0.25994E+11, 0.29819E+11, 0.34150E+11, 0.39044E+11, 0.44568E+11, &
         0.50794E+11, 0.57799E+11, 0.65672E+11, 0.74506E+11, 0.84408E+11, &
         0.95490E+11, 0.10788E+12, 0.12171E+12, 0.13713E+12, 0.15431E+12, &
         0.17342E+12, 0.19465E+12, 0.21822E+12, 0.24435E+12, 0.27329E+12, &
         0.30530E+12, 0.34069E+12, 0.37976E+12, 0.42286E+12, 0.47034E+12, &
         0.52262E+12, 0.58012E+12, 0.64330E+12, 0.71267E+12, 0.78875E+12, &
         0.87214E+12, 0.96344E+12, 0.10633E+13, 0.11725E+13, 0.12918E+13, &
         0.14220E+13, 0.15640E+13, 0.17188E+13, 0.18873E+13, 0.20706E+13, &
         0.22700E+13, 0.24866E+13, 0.27218E+13, 0.29771E+13, 0.32538E+13, &
         0.35537E+13, 0.38784E+13, 0.42299E+13, 0.46100E+13, 0.50208E+13, &
         0.54645E+13, 0.59435E+13, 0.64603E+13, 0.70175E+13, 0.76180E+13, &
         0.82647E+13, 0.89608E+13, 0.97097E+13, 0.10515E+14, 0.11380E+14, &
         0.12310E+14, 0.13307E+14, 0.14378E+14, 0.15526E+14, 0.16756E+14, &
         0.18075E+14/
!...        --      7646
      DATA (QOFT(2,J),J=1,119)/ 0.11735E+06, 0.21659E+06, 0.35745E+06, &
         0.55307E+06, 0.82078E+06, 0.11833E+07, 0.16700E+07, 0.23189E+07, &
         0.31781E+07, 0.43084E+07, 0.57862E+07, 0.77065E+07, 0.10187E+08, &
         0.13374E+08, 0.17445E+08, 0.22619E+08, 0.29161E+08, 0.37391E+08, &
         0.47698E+08, 0.60545E+08, 0.76487E+08, 0.96188E+08, 0.12043E+09, &
         0.15015E+09, 0.18644E+09, 0.23059E+09, 0.28410E+09, 0.34874E+09, &
         0.42657E+09, 0.51995E+09, 0.63167E+09, 0.76489E+09, 0.92329E+09, &
         0.11111E+10, 0.13331E+10, 0.15950E+10, 0.19029E+10, 0.22641E+10, &
         0.26867E+10, 0.31801E+10, 0.37547E+10, 0.44224E+10, 0.51965E+10, &
         0.60921E+10, 0.71261E+10, 0.83174E+10, 0.96873E+10, 0.11260E+11, &
         0.13061E+11, 0.15120E+11, 0.17471E+11, 0.20149E+11, 0.23196E+11, &
         0.26656E+11, 0.30578E+11, 0.35019E+11, 0.40038E+11, 0.45703E+11, &
         0.52087E+11, 0.59270E+11, 0.67343E+11, 0.76403E+11, 0.86556E+11, &
         0.97921E+11, 0.11062E+12, 0.12481E+12, 0.14062E+12, 0.15824E+12, &
         0.17783E+12, 0.19961E+12, 0.22377E+12, 0.25057E+12, 0.28024E+12, &
         0.31308E+12, 0.34936E+12, 0.38943E+12, 0.43362E+12, 0.48232E+12, &
         0.53593E+12, 0.59489E+12, 0.65968E+12, 0.73081E+12, 0.80883E+12, &
         0.89434E+12, 0.98797E+12, 0.10904E+13, 0.12024E+13, 0.13247E+13, &
         0.14582E+13, 0.16038E+13, 0.17625E+13, 0.19353E+13, 0.21233E+13, &
         0.23278E+13, 0.25499E+13, 0.27911E+13, 0.30528E+13, 0.33366E+13, &
         0.36442E+13, 0.39772E+13, 0.43376E+13, 0.47273E+13, 0.51486E+13, &
         0.56036E+13, 0.60948E+13, 0.66248E+13, 0.71962E+13, 0.78119E+13, &
         0.84751E+13, 0.91889E+13, 0.99569E+13, 0.10783E+14, 0.11670E+14, &
         0.12623E+14, 0.13646E+14, 0.14744E+14, 0.15921E+14, 0.17183E+14, &
         0.18535E+14/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_CLONO2


!
!     *****************
      SUBROUTINE QT_NOP(T, ISO, GSI, QT)         ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(1) :: XGJ
      REAL(DOUBLE), DIMENSION(1,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 3./
!...      NO+
!...        --        46
      DATA (QOFT(1,J),J=1,119)/ 0.63956E+02, 0.90185E+02, 0.11642E+03, &
         0.14265E+03, 0.16889E+03, 0.19513E+03, 0.22138E+03, 0.24763E+03, &
         0.27388E+03, 0.30013E+03, 0.32639E+03, 0.35266E+03, 0.37894E+03, &
         0.40523E+03, 0.43155E+03, 0.45790E+03, 0.48429E+03, 0.51074E+03, &
         0.53725E+03, 0.56383E+03, 0.59052E+03, 0.61731E+03, 0.64422E+03, &
         0.67127E+03, 0.69846E+03, 0.72582E+03, 0.75335E+03, 0.78108E+03, &
         0.80901E+03, 0.83715E+03, 0.86552E+03, 0.89413E+03, 0.92298E+03, &
         0.95208E+03, 0.98144E+03, 0.10111E+04, 0.10410E+04, 0.10712E+04, &
         0.11017E+04, 0.11325E+04, 0.11636E+04, 0.11950E+04, 0.12268E+04, &
         0.12588E+04, 0.12912E+04, 0.13239E+04, 0.13570E+04, 0.13903E+04, &
         0.14241E+04, 0.14581E+04, 0.14926E+04, 0.15273E+04, 0.15624E+04, &
         0.15979E+04, 0.16337E+04, 0.16699E+04, 0.17065E+04, 0.17434E+04, &
         0.17806E+04, 0.18183E+04, 0.18563E+04, 0.18947E+04, 0.19334E+04, &
         0.19725E+04, 0.20120E+04, 0.20519E+04, 0.20921E+04, 0.21327E+04, &
         0.21737E+04, 0.22151E+04, 0.22568E+04, 0.22990E+04, 0.23415E+04, &
         0.23844E+04, 0.24276E+04, 0.24713E+04, 0.25153E+04, 0.25598E+04, &
         0.26046E+04, 0.26497E+04, 0.26953E+04, 0.27413E+04, 0.27876E+04, &
         0.28343E+04, 0.28815E+04, 0.29290E+04, 0.29769E+04, 0.30251E+04, &
         0.30738E+04, 0.31229E+04, 0.31723E+04, 0.32222E+04, 0.32724E+04, &
         0.33230E+04, 0.33740E+04, 0.34254E+04, 0.34772E+04, 0.35294E+04, &
         0.35819E+04, 0.36349E+04, 0.36883E+04, 0.37420E+04, 0.37961E+04, &
         0.38507E+04, 0.39056E+04, 0.39609E+04, 0.40166E+04, 0.40727E+04, &
         0.41292E+04, 0.41861E+04, 0.42434E+04, 0.43010E+04, 0.43591E+04, &
         0.44176E+04, 0.44764E+04, 0.45357E+04, 0.45953E+04, 0.46554E+04, &
         0.47158E+04/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_NOP


!
!     *****************
      SUBROUTINE QT_HOBR(T, ISO, GSI, QT)        ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(2) :: XGJ
      REAL(DOUBLE), DIMENSION(2,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 8., 8./
!...     HOBr
!...        --       169
      DATA (QOFT(1,J),J=1,119)/ 0.24445E+04, 0.41206E+04, 0.60683E+04, &
         0.82610E+04, 0.10689E+05, 0.13352E+05, 0.16261E+05, 0.19427E+05, &
         0.22867E+05, 0.26600E+05, 0.30643E+05, 0.35018E+05, 0.39745E+05, &
         0.44844E+05, 0.50338E+05, 0.56249E+05, 0.62599E+05, 0.69410E+05, &
         0.76706E+05, 0.84509E+05, 0.92845E+05, 0.10174E+06, 0.11121E+06, &
         0.12128E+06, 0.13199E+06, 0.14335E+06, 0.15540E+06, 0.16815E+06, &
         0.18165E+06, 0.19591E+06, 0.21096E+06, 0.22684E+06, 0.24358E+06, &
         0.26120E+06, 0.27974E+06, 0.29922E+06, 0.31969E+06, 0.34118E+06, &
         0.36372E+06, 0.38735E+06, 0.41210E+06, 0.43800E+06, 0.46511E+06, &
         0.49345E+06, 0.52307E+06, 0.55400E+06, 0.58628E+06, 0.61997E+06, &
         0.65509E+06, 0.69170E+06, 0.72984E+06, 0.76954E+06, 0.81087E+06, &
         0.85386E+06, 0.89856E+06, 0.94502E+06, 0.99329E+06, 0.10434E+07, &
         0.10955E+07, 0.11495E+07, 0.12055E+07, 0.12636E+07, 0.13238E+07, &
         0.13862E+07, 0.14508E+07, 0.15177E+07, 0.15870E+07, 0.16587E+07, &
         0.17328E+07, 0.18095E+07, 0.18888E+07, 0.19707E+07, 0.20554E+07, &
         0.21428E+07, 0.22331E+07, 0.23263E+07, 0.24225E+07, 0.25217E+07, &
         0.26241E+07, 0.27296E+07, 0.28385E+07, 0.29506E+07, 0.30662E+07, &
         0.31853E+07, 0.33079E+07, 0.34341E+07, 0.35641E+07, 0.36979E+07, &
         0.38355E+07, 0.39771E+07, 0.41228E+07, 0.42725E+07, 0.44265E+07, &
         0.45848E+07, 0.47474E+07, 0.49145E+07, 0.50862E+07, 0.52624E+07, &
         0.54435E+07, 0.56293E+07, 0.58201E+07, 0.60159E+07, 0.62168E+07, &
         0.64229E+07, 0.66343E+07, 0.68511E+07, 0.70734E+07, 0.73013E+07, &
         0.75349E+07, 0.77742E+07, 0.80196E+07, 0.82709E+07, 0.85283E+07, &
         0.87920E+07, 0.90620E+07, 0.93385E+07, 0.96215E+07, 0.99112E+07, &
         0.10208E+08/
!...        --       161
      DATA (QOFT(2,J),J=1,119)/ 0.24350E+04, 0.41047E+04, 0.60448E+04, &
         0.82291E+04, 0.10648E+05, 0.13301E+05, 0.16200E+05, 0.19355E+05, &
         0.22784E+05, 0.26504E+05, 0.30534E+05, 0.34895E+05, 0.39607E+05, &
         0.44691E+05, 0.50169E+05, 0.56063E+05, 0.62394E+05, 0.69186E+05, &
         0.76461E+05, 0.84243E+05, 0.92555E+05, 0.10142E+06, 0.11087E+06, &
         0.12091E+06, 0.13159E+06, 0.14292E+06, 0.15494E+06, 0.16766E+06, &
         0.18112E+06, 0.19534E+06, 0.21036E+06, 0.22620E+06, 0.24289E+06, &
         0.26047E+06, 0.27896E+06, 0.29840E+06, 0.31882E+06, 0.34025E+06, &
         0.36274E+06, 0.38630E+06, 0.41099E+06, 0.43683E+06, 0.46387E+06, &
         0.49215E+06, 0.52169E+06, 0.55255E+06, 0.58475E+06, 0.61836E+06, &
         0.65340E+06, 0.68992E+06, 0.72796E+06, 0.76757E+06, 0.80880E+06, &
         0.85169E+06, 0.89628E+06, 0.94263E+06, 0.99079E+06, 0.10408E+07, &
         0.10927E+07, 0.11466E+07, 0.12025E+07, 0.12605E+07, 0.13205E+07, &
         0.13828E+07, 0.14472E+07, 0.15140E+07, 0.15831E+07, 0.16546E+07, &
         0.17286E+07, 0.18051E+07, 0.18842E+07, 0.19660E+07, 0.20504E+07, &
         0.21377E+07, 0.22277E+07, 0.23207E+07, 0.24167E+07, 0.25157E+07, &
         0.26178E+07, 0.27231E+07, 0.28317E+07, 0.29436E+07, 0.30589E+07, &
         0.31777E+07, 0.33001E+07, 0.34260E+07, 0.35557E+07, 0.36892E+07, &
         0.38265E+07, 0.39678E+07, 0.41131E+07, 0.42626E+07, 0.44162E+07, &
         0.45741E+07, 0.47364E+07, 0.49031E+07, 0.50744E+07, 0.52503E+07, &
         0.54309E+07, 0.56164E+07, 0.58067E+07, 0.60021E+07, 0.62025E+07, &
         0.64081E+07, 0.66191E+07, 0.68354E+07, 0.70572E+07, 0.72846E+07, &
         0.75177E+07, 0.77565E+07, 0.80013E+07, 0.82521E+07, 0.85090E+07, &
         0.87721E+07, 0.90415E+07, 0.93173E+07, 0.95997E+07, 0.98888E+07, &
         0.10185E+08/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_HOBR


!
!     *****************
      SUBROUTINE QT_C2H4(T, ISO, GSI, QT)        ! Total Internal Partition Function
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
!++:  bd-QT
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISO       ! isotope code (HITRAN INDEX)
      REAL(DOUBLE)  :: T                ! temperature in K
      REAL(DOUBLE) , INTENT(OUT) :: GSI ! state independent nuclear
                                        ! degeneracyfactor
      REAL(DOUBLE)  :: QT               ! Total Internal Partition Function
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !, I
      REAL(DOUBLE), DIMENSION(2) :: XGJ
      REAL(DOUBLE), DIMENSION(2,119) :: QOFT
      REAL(DOUBLE), DIMENSION(NT) :: Q
      !REAL(DOUBLE) :: EPS
!-----------------------------------------------
!
      DATA XGJ/ 1., 2./
!...     C2H4
!...        --       221
      DATA (QOFT(1,J),J=1,119)/ 0.95843E+03, 0.16137E+04, 0.23744E+04, &
         0.32285E+04, 0.41694E+04, 0.51963E+04, 0.63143E+04, 0.75337E+04, &
         0.88702E+04, 0.10344E+05, 0.11978E+05, 0.13802E+05, 0.15846E+05, &
         0.18145E+05, 0.20740E+05, 0.23675E+05, 0.27000E+05, 0.30770E+05, &
         0.35048E+05, 0.39905E+05, 0.45420E+05, 0.51680E+05, 0.58786E+05, &
         0.66850E+05, 0.75997E+05, 0.86369E+05, 0.98123E+05, 0.11144E+06, &
         0.12651E+06, 0.14356E+06, 0.16284E+06, 0.18463E+06, 0.20923E+06, &
         0.23699E+06, 0.26831E+06, 0.30360E+06, 0.34334E+06, 0.38808E+06, &
         0.43840E+06, 0.49495E+06, 0.55847E+06, 0.62976E+06, 0.70973E+06, &
         0.79935E+06, 0.89973E+06, 0.10121E+07, 0.11378E+07, 0.12782E+07, &
         0.14351E+07, 0.16102E+07, 0.18055E+07, 0.20231E+07, 0.22656E+07, &
         0.25354E+07, 0.28356E+07, 0.31692E+07, 0.35398E+07, 0.39511E+07, &
         0.44074E+07, 0.49132E+07, 0.54736E+07, 0.60940E+07, 0.67803E+07, &
         0.75392E+07, 0.83776E+07, 0.93035E+07, 0.10325E+08, 0.11452E+08, &
         0.12694E+08, 0.14062E+08, 0.15567E+08, 0.17224E+08, 0.19045E+08, &
         0.21046E+08, 0.23243E+08, 0.25655E+08, 0.28300E+08, 0.31200E+08, &
         0.34377E+08, 0.37856E+08, 0.41662E+08, 0.45826E+08, 0.50378E+08, &
         0.55351E+08, 0.60781E+08, 0.66707E+08, 0.73172E+08, 0.80219E+08, &
         0.87899E+08, 0.96262E+08, 0.10537E+09, 0.11527E+09, 0.12604E+09, &
         0.13775E+09, 0.15047E+09, 0.16428E+09, 0.17927E+09, 0.19553E+09, &
         0.21316E+09, 0.23226E+09, 0.25296E+09, 0.27537E+09, 0.29963E+09, &
         0.32587E+09, 0.35425E+09, 0.38492E+09, 0.41805E+09, 0.45383E+09, &
         0.49246E+09, 0.53413E+09, 0.57908E+09, 0.62754E+09, 0.67977E+09, &
         0.73602E+09, 0.79660E+09, 0.86179E+09, 0.93194E+09, 0.10074E+10, &
         0.10885E+10/
!...        --       231
      DATA (QOFT(2,J),J=1,119)/ 0.39228E+04, 0.66051E+04, 0.97190E+04, &
         0.13215E+05, 0.17066E+05, 0.21270E+05, 0.25846E+05, 0.30838E+05, &
         0.36309E+05, 0.42341E+05, 0.49032E+05, 0.56496E+05, 0.64862E+05, &
         0.74275E+05, 0.84897E+05, 0.96912E+05, 0.11052E+06, 0.12595E+06, &
         0.14347E+06, 0.16335E+06, 0.18592E+06, 0.21155E+06, 0.24064E+06, &
         0.27365E+06, 0.31109E+06, 0.35354E+06, 0.40166E+06, 0.45615E+06, &
         0.51785E+06, 0.58765E+06, 0.66657E+06, 0.75575E+06, 0.85646E+06, &
         0.97011E+06, 0.10983E+07, 0.12428E+07, 0.14055E+07, 0.15886E+07, &
         0.17945E+07, 0.20260E+07, 0.22861E+07, 0.25779E+07, 0.29052E+07, &
         0.32721E+07, 0.36830E+07, 0.41429E+07, 0.46573E+07, 0.52323E+07, &
         0.58744E+07, 0.65912E+07, 0.73906E+07, 0.82816E+07, 0.92740E+07, &
         0.10379E+08, 0.11607E+08, 0.12973E+08, 0.14490E+08, 0.16174E+08, &
         0.18042E+08, 0.20112E+08, 0.22406E+08, 0.24945E+08, 0.27755E+08, &
         0.30861E+08, 0.34293E+08, 0.38083E+08, 0.42266E+08, 0.46878E+08, &
         0.51961E+08, 0.57560E+08, 0.63724E+08, 0.70504E+08, 0.77959E+08, &
         0.86150E+08, 0.95145E+08, 0.10502E+09, 0.11585E+09, 0.12772E+09, &
         0.14072E+09, 0.15496E+09, 0.17054E+09, 0.18759E+09, 0.20622E+09, &
         0.22658E+09, 0.24880E+09, 0.27306E+09, 0.29952E+09, 0.32837E+09, &
         0.35981E+09, 0.39404E+09, 0.43131E+09, 0.47186E+09, 0.51595E+09, &
         0.56387E+09, 0.61594E+09, 0.67247E+09, 0.73382E+09, 0.80038E+09, &
         0.87255E+09, 0.95076E+09, 0.10355E+10, 0.11272E+10, 0.12265E+10, &
         0.13339E+10, 0.14501E+10, 0.15756E+10, 0.17113E+10, 0.18577E+10, &
         0.20159E+10, 0.21865E+10, 0.23705E+10, 0.25688E+10, 0.27826E+10, &
         0.30129E+10, 0.32608E+10, 0.35277E+10, 0.38149E+10, 0.41237E+10, &
         0.44557E+10/

      !EPS = 0.01
!
      GSI = XGJ(ISO)
      Q(:NT) = QOFT(ISO,:NT)

!
!...value depends on temperature range
      IF (T<70. .OR. T>3000.) THEN
         QT = -1.
         WRITE (*, '(A)') '  OUT OF TEMPERATURE RANGE'
         GO TO 99
      ENDIF

      CALL ATOB (T, QT, TDAT, Q, NT)
!
   99 CONTINUE
      RETURN
      END SUBROUTINE QT_C2H4


!
!
!***************************
      SUBROUTINE ATOB(AA, BB, A, B, NPT)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!***************************
!...LaGrange 3- and 4-point interpolation
!...arrays A and B are the npt data points,  given aa, a value of the
!...A variable, the routine will find the corresponding bb value
!
!...input:  aa
!...output: bb
!...Translated by Pacific-Sierra Research 77to90  4.4D  16:34:04   3/24/05
!...Switches: -yf12
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NPT
      REAL(DOUBLE) , INTENT(IN) :: AA
      REAL(DOUBLE) , INTENT(OUT) :: BB
      REAL(DOUBLE) , INTENT(IN) :: A(NT)
      REAL(DOUBLE) , INTENT(IN) :: B(NT)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J
      REAL(DOUBLE) :: A0D1, A0D2, A1D1, A1D2, A2D1, A2D2, A0, A1, A2, A0D3, &
         A1D3, A2D3, A3D1, A3D2, A3D3, A3
!-----------------------------------------------
!
!
!
      DO I = 2, NPT
         IF (A(I) < AA) CYCLE
         IF (I<3 .OR. I==NPT) THEN
!     LaGrange three point interpolation
            J = I
            IF (I < 3) J = 3
            IF (I == NPT) J = NPT
!.....do not devide by zero
            A0D1 = A(J-2) - A(J-1)
            IF (A0D1 == 0.) A0D1 = 0.0001
            A0D2 = A(J-2) - A(J)
            IF (A0D2 == 0.) A0D2 = 0.0001
            A1D1 = A(J-1) - A(J-2)
            IF (A1D1 == 0.) A1D1 = 0.0001
            A1D2 = A(J-1) - A(J)
            IF (A1D2 == 0.) A1D2 = 0.0001
            A2D1 = A(J) - A(J-2)
            IF (A2D1 == 0.) A2D1 = 0.0001
            A2D2 = A(J) - A(J-1)
            IF (A2D2 == 0.) A2D2 = 0.0001
!
            A0 = (AA - A(J-1))*(AA - A(J))/(A0D1*A0D2)
            A1 = (AA - A(J-2))*(AA - A(J))/(A1D1*A1D2)
            A2 = (AA - A(J-2))*(AA - A(J-1))/(A2D1*A2D2)
!
            BB = A0*B(J-2) + A1*B(J-1) + A2*B(J)
!
         ELSE
!     LaGrange four point interpolation
            J = I
!.....do not devide by zero
            A0D1 = A(J-2) - A(J-1)
            IF (A0D1 == 0.) A0D1 = 0.0001
            A0D2 = A(J-2) - A(J)
            IF (A0D2 == 0.) A0D2 = 0.0001
            A0D3 = A(J-2) - A(J+1)
            IF (A0D3 == 0.) A0D3 = 0.0001
!
            A1D1 = A(J-1) - A(J-2)
            IF (A1D1 == 0.) A1D1 = 0.0001
            A1D2 = A(J-1) - A(J)
            IF (A1D2 == 0.) A1D2 = 0.0001
            A1D3 = A(J-1) - A(J+1)
            IF (A1D3 == 0.) A1D3 = 0.0001
!
            A2D1 = A(J) - A(J-2)
            IF (A2D1 == 0.) A2D1 = 0.0001
            A2D2 = A(J) - A(J-1)
            IF (A2D2 == 0.) A2D2 = 0.0001
            A2D3 = A(J) - A(J+1)
            IF (A2D3 == 0.) A2D3 = 0.0001
!
            A3D1 = A(J+1) - A(J-2)
            IF (A3D1 == 0.) A3D1 = 0.0001
            A3D2 = A(J+1) - A(J-1)
            IF (A3D2 == 0.) A3D2 = 0.0001
            A3D3 = A(J+1) - A(J)
            IF (A3D3 == 0.) A3D3 = 0.0001
!
            A0 = (AA - A(J-1))*(AA - A(J))*(AA - A(J+1))
            A0 = A0/(A0D1*A0D2*A0D3)
            A1 = (AA - A(J-2))*(AA - A(J))*(AA - A(J+1))
            A1 = A1/(A1D1*A1D2*A1D3)
            A2 = (AA - A(J-2))*(AA - A(J-1))*(AA - A(J+1))
            A2 = A2/(A2D1*A2D2*A2D3)
            A3 = (AA - A(J-2))*(AA - A(J-1))*(AA - A(J))
            A3 = A3/(A3D1*A3D2*A3D3)
!
            BB = A0*B(J-2) + A1*B(J-1) + A2*B(J) + A3*B(J+1)
         ENDIF
!
         EXIT
      END DO
      RETURN
      END SUBROUTINE ATOB


!

end module tips

