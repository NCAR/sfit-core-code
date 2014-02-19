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

      MODULE VIBFCN

! --- DEFAULT VIBRATION FUNCTION CALCULATION
! --- VERSION 3.90 - TO 99 MOLECULES

      USE PARAMS
      USE RETVPARAM
      USE MOLCPARAM
      !USE WRITEOUT

      IMPLICIT NONE

! --- VIBRATIONAL FUNCTION DATA

! June 2013
! add isotopes to N2O # 6 7 8 / from Toth via ggg linelist

! July 26 2012
! edit for removal of default isotope separation

! 01/09/07
!  Added CH3CHO

! 12/09/05
!  Added PAN

! 10/31/05
!  Expanded max ivib values to 30, and added CH3OH, CH3CN, and C2H6PL

! 7/18/05
!  Changed CFCLO to COCLF
!          CCL2O to COCL2
!          CHFCL2 to CHCL2F

!   Change
! Changes 6/22/04 at DU
!  Filled in missing vibrational values in NMODE and IVIB for:
!                    HONO, CH3F, CH3BR, CH3I, HCOOH, CHFCL2, OCLO
!  Updated vibrational values in IVIB for: C2H6, C2H4
!  Changed COCLF to CFCLO and COCL2 to CCL2O
!

      INTEGER, PARAMETER                        :: MAXVIBVALS = 30
      REAL(DOUBLE), DIMENSION(MAXGAS,LAYMAX+1)  :: QV
      INTEGER, DIMENSION(2,MAXVIBVALS,MOLTOTAL) :: IVIB
      INTEGER, DIMENSION(MOLTOTAL)              :: NMODE
      INTEGER, PRIVATE                          :: J,K,L

      DATA NMODE/ &
        3,  3,  3,  3,  1,  4,  1,  1,  3,  3,                       & !1-10
        4,  9,  1,  2,  1,  1,  1,  1,  3,  6,                       & !11-20
        3,  3,  6,  6, 11,  9,  8,  3,  6,  6,                       & !21-30
        4,  9,  6, 12,  4,  6,  6, 11, 12,  5,                       & !31-40
        1,  9,  6,  6,  6,  9,  3,  9,  1, 15,                       & !41-50
        6,  0,  0,  0,  0,  0,  0,  3,  1,  1,                       & !51-60
       12, 12,  1, 12,  8, 11, 27, 15,  8,  0,                       & !61-70
       18, 33, 27, 27, 21, 30,  0,  0,  0,  0,                       & !71-80
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,                       & !81-90
        0,  0,  0,  0,  0,  0,  0,  0,  0/                             !91-99


      DATA (((IVIB(J,K,L),J=1,2),K=1,MAXVIBVALS),L=1,10)/ &
      3657,1, 1595,1, 3756,1, 54*0,                                  & !H2O
      1333,1,  667,2, 2349,1, 54*0,                                  & !CO2
      1103,1,  701,1, 1042,1, 54*0,                                  & !O3
      1285,1,  589,2, 2224,1, 54*0,                                  & !N2O
      2143,1, 58*0,                                                  & !CO
      2917,1, 1534,2, 3019,3, 1306,3, 52*0,                          & !CH4
      1580,1, 58*0,                                                  & !O2
      1876,1, 58*0,                                                  & !NO
      1152,1,  518,1, 1362,1, 54*0,                                  & !SO2
      1320,1,  750,1, 1617,1, 54*0/                                    !NO2

      DATA (((IVIB(J,K,L),J=1,2),K=1,MAXVIBVALS),L=11,20)/ &
      3337,1,  951,1, 3444,2, 1627,2,  52*0,                         & !NH3
      3550,1, 1710,1, 1303,1, 1325,1,  879,1,  647,1,  579,1,  762,1,& !HNO3
       456,1, 42*0,                                                  &
      3570,1, 58*0,                                                  & !OH
      3962,1, 58*0,                                                  & !HF
      2886,1, 58*0,                                                  & !HCL
      2559,1, 58*0,                                                  & !HBR
      2230,1, 58*0,                                                  & !HI
       844,1, 58*0,                                                  & !CLO
       859,1,  521,2, 2062,1,  54*0,                                 & !OCS
      2782,1, 1746,1, 1500,1, 1167,1, 2843,1, 1249,1,  48*0/           !H2CO

      DATA (((IVIB(J,K,L),J=1,2),K=1,MAXVIBVALS),L=21,30)/ &
      3609,1, 1238,1,  720,1, 54*0,                                  & !HOCL
      3436,1, 1392,1, 1098,1, 54*0,                                  & !HO2
      2869,1, 1435,1, 1408,1,  870,1, 3417,1, 1266,1, 48*0,          & !H2O2
       543,1,  596,1,  790,1, 1263,1, 1700,1, 3591,1, 48*0,          & !HONO
       200,1,  400,1,  500,1,  633,1,  735,1,  803,1,  880,1, 1304,1,& !HO2NO2
      1396,1, 1728,1, 3540,1, 38*0,                                  &
      1728,1,  743,1,  353,1, 1338,1,  614,1,   85,1,  577,1, 1247,1,& !N2O5
       860,1, 42*0,                                                  &
      1735,1, 1292,1,  809,1,  780,1,  560,1,  434,1,  270,1,  711,1,& !CLONO2
      44*0,                                                          &
      2097,1,  712,2, 3312,1, 54*0,                                  & !HCN
      2930,1, 1464,1, 1049,1, 3006,2, 1467,2, 1182,2, 48*0,          & !CH3F
      2966,1, 1355,1,  732,1, 3039,2, 1455,2, 1015,2, 48*0/            !CH3CL

      DATA (((IVIB(J,K,L),J=1,2),K=1,MAXVIBVALS),L=31,40)/ &
      909, 1,  435,2, 1283,3,  632,3, 52*0,                          & !CF4
      1102,1,  667,1,  468,1,  262,1,  321,1,  923,1,  462,1, 1161,1,& !CCL2F2
       437,1, 42*0,                                                  &
      1085,1,  535,1,  350,1,  847,2,  394,2,  241,2, 48*0,          & !CCL3F
       344,1,  526,1, 1010,1, 1386,1, 2951,1,  239,2,  303,2, 724,2, & !CH3CCL3
      1088,2, 1455,2, 3013,2,  239,2, 36*0,                          &
       459,1,  217,2,  799,3,  314,3, 52*0,                          & !CCL4
      1944,1,  963,1,  582,1, 1242,1,  619,1,  774,1, 48*0,          & !COF2
      1868,1, 1095,1,  776,1,  501,1,  415,1,  667,1, 48*0,          & !COCLF
      2985,2, 2969,2, 2954,1, 2896,1, 1472,2, 1468,1, 1388,1, 1379,1,& !C2H6
      1190,2,  995,1,  822,2, 38*0,                                  &
      3026,1, 1623,1, 1342,1, 1023,1, 3103,1, 1236,1,  949,1,        & !C2H4
       943,1, 3106,1,  826,1, 2989,1, 1444,1, 36*0,                  &
      3374,1, 1974,1, 3295,1,  612,2,  729,2, 50*0/                    !C2H2

      DATA (((IVIB(J,K,L),J=1,2),K=1,MAXVIBVALS),L=41,50)/ &
      2359,1, 58*0,                                                  & !N2
      3021,1, 1313,1, 1107,1,  809,1,  596,1,  413,1, 1351,1, 1127,1,& !F22
       366,1, 42*0,                                                  &
      1828,1,  574,1,  303,1,  586,1,  851,1,  445,1, 48*0,          & !COCL2
      2973,1, 1306,1,  611,1, 3056,1, 1443,2,  955,2, 48*0,          & !CH3BR
      2933,1, 1252,1,  533,1, 3060,2, 1436,2,  882,2, 48*0,          & !CH3I
      3570,1, 2943,1, 1770,1, 1387,1, 1229,1, 1105,1,  625,1, 1033,1,& !HCOOH
       638,1, 42*0,                                                  &
      1183,1, 2615,1, 2627,1, 54*0,                                  & !H2S
      3024,1, 1315,1, 1079,1,  744,1,  458,1,  277,1, 1238,1,  807,1,& !CHCL2F
       367,1, 42*0,                                                  &
      1580,1, 58*0,                                                  & !O2CIA
       774,1,  642,2,  948,3,  616,3,  525,3,  347,3,  48*0/           !SF6

      DATA (((IVIB(J,K,L),J=1,2),K=1,MAXVIBVALS),L=51,60)/ &
      1032,1,  647,1,  905,2,  493,2,   52*0,                        & !NF3
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
       945,1,  449,1, 1108,1, 54*0,                                  & !OCLO
      5000,1, 58*0,                                                  & !F134A
      5000,1, 58*0/                                                    !C3H8

      DATA (((IVIB(J,K,L),J=1,2),K=1,MAXVIBVALS),L=61,70)/ &
       344,1,  681,1, 1104,1, 1396,1, 2951,1,  339,2,  303,2,  905,2,& !F142B
      1135,2, 1455,2, 3013,2,  239,2, 36*0,                          &
      1100,1,  807,1,  348,1,  675,1,  372,1, 1250,2,  619,2,  372,2,& !CFC113
       778,2,  271,2,  114,2,   68,2, 36*0,                          &
      5000,1, 58*0,                                                  & !F141B
      3681,1, 3000,1, 2844,1, 1477,1, 1455,1, 1345,1, 1060,1, 1033,1,& !CH3OH
      2960,1, 1477,1, 1165,1,  295,1,   36*0,                        &
      2954,1, 2267,1, 1385,1,  920,1, 3009,2, 1448,2, 1041,2,  362,2,& !CH3CNPL
        44*0,                                                        &
      2985,2, 2969,2, 2954,1, 2896,1, 1472,2, 1468,1, 1388,1, 1379,1,& !C2H6PL
      1190,2,  995,1,  822,2, 38*0,                                  &
      3164,1, 3121,1, 3058,1, 1880,1, 1806,1, 1475,1, 1471,1, 1400,1,& !PAN
      1352,1, 1172,1, 1065,1,  999,1,  984,1,  828,1,  806,1,  736,1,&
       727,1,  616,1,  585,1,  495,1,  373,1,  327,1,  316,1,  100,1,&
        96,1,   82,1,   24,1,  6*0,                                  &
      3005,1, 2917,1, 2822,1, 1743,1, 1441,1, 1400,1, 1352,1, 113,1, & !CH3CHO
       919,1,  509,1, 2967,1, 1420,1,  867,1,  763,1,  150,1,  30*0, &
      2954,1, 2267,1, 1385,1,  920,1, 3009,2, 1448,2, 1041,2,  362,2,& !CH3CN
        44*0,                                                        &
      5000,1, 58*0/                                                    !OTHER

      DATA (((IVIB(J,K,L),J=1,2),K=1,MAXVIBVALS),L=71,80)/ &

      3583,1, 3051,1, 2944,1, 1788,1, 1430,1, 1382,1, 1264,1, 1182,1,& !CH3COOH
       989,1,  847,1,  657,1,  581,1, 2996,1, 1430,1, 1048,1,  642,1,&
       534,1,   93,1,   24*0,                                        &

      3097,1, 3000,1, 2968,1, 1790,1, 1613,1, 1452,1, 1071,1,  991,1,& !C5H8
       906,1,  894,1,   40*0,                                        &

      3392,1, 3104,1, 3020,1, 2972,1, 1714,1, 1623,1, 1400,1, 1247,1,& !MVK
      1180,1,  960,1,   40*0,                                        &

      1700,1, 1000,1,  930,1,  800,1, 52*0,                          & !MACR

      3091,1, 3015,1, 2991,1, 2973,1, 2931,1, 1653,1, 1459,1, 1420,1,& !C3H6
      1378,1, 1298,1, 1174,1,  935,1,  919,1,  428,1, 2953,1, 1442,1,&
      1045,1,  990,1,  912,1,  575,1,  188,1,   18*0,                &

      2975,1,  900,1,  56*0,                                         & !C4H8

      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0/                                                    !OTHER

      DATA (((IVIB(J,K,L),J=1,2),K=1,MAXVIBVALS),L=81,90)/ &
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0/                                                    !OTHER

      DATA (((IVIB(J,K,L),J=1,2),K=1,MAXVIBVALS),L=91,MOLTOTAL)/ &
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0,                                                  & !OTHER
      5000,1, 58*0/                                                    !OTHER

      CONTAINS

!----------------------------------------------------------------------
      SUBROUTINE QVIB( XCS_DETAIL )

!  COMPUTE VIBRATIONAL PARTITION FUNCTION FOR RETRIEVAL AND BACKGROUND
!  GASES FOR EACH LAYER TEMPERATURE AND AT THE REFERENCE TEMPERATURE
!  296 K.  THIS SUBROUTINE IS BASED ON FASCOD1C
!
!     NMODE=NUMBER OF VIBRATIONAL MODES
!           TORSIONAL MODES NOT INCLUDED -- C2H6(#38), CLONO2 (#27)
!     IVIB=ARRAY CONTAINING FUNDAMENTAL FREQUENCIES AND DEGENERACIES

      LOGICAL, INTENT(IN) :: XCS_DETAIL
      INTEGER             :: KLEVEL, I, K, N, IMODE, ISUBS(MAXSPE)
      REAL(DOUBLE)        :: TEMP, FREQ

      KLEVEL = KMAX + 1
      DO I = 1, NGAS
!  --- ICODE=PC GAS CODE (1=H2O,2=CO2,ETC.)
         IF( IFMIX(I) == 0 ) CYCLE
         DO K = 1, KLEVEL
            QV(I,K) = 1.D0
            IF (ICODE(I) - MOLTOTAL > 0) CYCLE
            IF (K <= KMAX) TEMP = T(K)
            IF (K == KLEVEL) TEMP = 296.D0
            N = NMODE(ICODE(I))
            DO IMODE = 1, N
               FREQ = IVIB(1,IMODE,ICODE(I))
               QV(I,K) = QV(I,K)/(1.D0 - EXP((-RCONST2*FREQ/TEMP)))**IVIB(2,IMODE,ICODE(I))
            END DO
         END DO
      END DO

!  --- GET INDEXES OF NGAS THAT WE ARE ACTUALLY USING
      N = 0
      DO I=1, NGAS
         IF( IFMIX(I) == 0 ) CYCLE
          N = N + 1
         ISUBS(N) = I
      END DO

      IF( XCS_DETAIL )THEN
!  --- WRITE OUT VIBRATIONAL PARTITION FUNCTIONS
         WRITE (16, 100)
         WRITE (16, 5, ADVANCE='NO')
         WRITE(16, 6) NAME(ICODE(ISUBS(:N)))
         DO K = 1, KLEVEL
            IF (K <= KMAX) TEMP = T(K)
            IF (K == KLEVEL) TEMP = 296.D0
            WRITE (16, 7, ADVANCE='NO') TEMP
            WRITE(16, 9) QV(ISUBS(:N),K)
         END DO
      ENDIF

      RETURN

    5 FORMAT(' TEMP    ')
    6 FORMAT(7(2X,A7))
    7 FORMAT(F8.2,':')
    9 FORMAT(7(F8.4,1X))
  100 FORMAT(/,' VIBRATIONAL PARTITION FUNCTION FOR GASES'/)

      END SUBROUTINE QVIB

      END MODULE VIBFCN
