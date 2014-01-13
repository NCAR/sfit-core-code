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


      real(double) function lineari (x1, y1, x2, y2, target_x)

! SFIT4 - v003.90 - TO 99 MOLECULES

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------

      USE params, ONLY:  double
!...Translated by Pacific-Sierra Research 77to90  4.4E  13:15:27   1/23/07
!...Switches: -yf12
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double) , intent(in) :: x1
      real(double) , intent(in) :: y1
      real(double) , intent(in) :: x2
      real(double) , intent(in) :: y2
      real(double) , intent(in) :: target_x
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!      write (6, *) 'x1,y1,x2,y2,target_x=', x1, y1, x2, y2, target_x

      lineari = y1 + (y2 - y1)*((target_x - x1)/(x2 - x1))
      return
      end function lineari



      integer function tablelinterp (n, x, y, targetx, value)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  double
!...Translated by Pacific-Sierra Research 77to90  4.4E  13:15:27   1/23/07
!...Switches: -yf12
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      real(double)  :: targetx
      real(double) , intent(out) :: value
      real(double)  :: x(n)
      real(double)  :: y(n)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      real(double) , external :: lineari
!-----------------------------------------------

      if (targetx < x(1)) then
         tablelinterp = -1
         return
      else if (targetx > x(n)) then
         tablelinterp = 1
         return
      else
         i = 2
   10    continue
         if (targetx < x(i)) then
            value = lineari(x(i-1),y(i-1),x(i),y(i),targetx)
!            write (6, *) 'Interpolated value = ', value
            tablelinterp = 0
            return
         endif
         i = i + 1
         if (i <= n) go to 10
         tablelinterp = -2
      endif
      return
      end function tablelinterp



!*****************************************************************************
!
! Purpose:
!
!   Convert beta to a new temperature using Diffusion coefficient rules
!
!   See J.O. Hirschfelder...
!
!
!
! Entry:
!
!   Mol  Molecule ID
!
!   T    temperature to convert to
!
!
!
! (c) GPC
! DRH 20040522.1501
!
!****************************************************************************/
      real(double) function betat (mol, t)

      use datafiles

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------

      USE params
!...Translated by Pacific-Sierra Research 77to90  4.4E  13:15:27   1/23/07
!...Switches: -yf12
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: mol
      real(double) , intent(in) :: t
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: epsilonair = 97.0
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: result, i
      real(double) :: treduced, omegat, omegat0
      real(double), dimension(MOLTOTAL) :: epsilon
      real(double), dimension(82) :: tstar, omegastar
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      integer , external :: tablelinterp
!-----------------------------------------------
!

      data (epsilon(i),i=1,MOLTOTAL)/ &
        230.9, & ! "Hirschfelder, Viscosity (Table 8.6-1)", H2O
        190.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , CO2
         -1.0, & !  none                                  , O3
        220.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , N2O
        110.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , CO
        144.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , CH4
        113.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , O2
        119.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , NO
        252.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , SO2
         -1.0, & !  none                                  , NO2         10
        146.8, & ! "Hirschfelder, Viscosity (Table 8.6-1)", NH3
         -1.0, & ! none                                   , HNO3
         -1.0, & !  none                                  , OH
        330.0, & ! "White,Zerilli,Jones"                  , HF  should be ~1240
        360.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , HCl
         -1.0, & !  none                                  , HBr
        324.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , HI
         -1.0, & !  none                                  , ClO
        335.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , OCS
         -1.0, & !  none                                  , H2CO        20
         -1.0, & !  none                                  , HOCl
         -1.0, & !  none                                  , HO2
         -1.0, & !  none                                  , H2O2
         -1.0, & !  none                                  , HONO
         -1.0, & !  none                                  , HO2N02
         -1.0, & !  none                                  , N205
         -1.0, & !  none                                  , ClONO2
         -1.0, & !  none                                  , HCN
         -1.0, & !  none                                  , CH3F
        855.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , CH3Cl       30
         -1.0, & !  none                                  , CF4
         -1.0, & !  none                                  , CCL2F2
         -1.0, & !  none                                  , CCL3F
         -1.0, & !  none                                  , CH3CCL3
         -1.0, & !  none                                  , CCL4
         -1.0, & !  none                                  , COF2
         -1.0, & !  none                                  , COCLF
        230.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , C2H6
        205.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , C2H4
        185.0, & ! "Hirschfelder, Viscosity (Table I-A)"  , C2H2        40
         91.5, & ! "Hirschfelder, Viscosity (Table I-A)"  , N2
         -1.0, & !  none                                  , CHF2CL
         -1.0, & !  none                                  , COCL2
         -1.0, & !  none                                  , CH3BR
         -1.0, & !  none                                  , CH3I
         -1.0, & !  none                                  , HCOOH
        221.1, & ! "Hirschfelder, Viscosity (Table 8.6-1)", H2S
         -1.0, & !  none                                  , CHCL2F
         -1.0, & !  none                                  , HDO
        200.9, & ! "Hirschfelder, Viscosity (Table I-A)"  , SF6         50
         -1.0, & !  none                                  , NF3
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , CH3D
         -1.0, & !  none                                  , O3668
         -1.0, & !  none                                  , O3686
         -1.0, & !  none                                  , O3667
         -1.0, & !  none                                  , O3676
         -1.0, & !  none                                  , OCLO
         -1.0, & !  none                                  , F134A
         -1.0, & !  none                                  , C3H8        60
         -1.0, & !  none                                  , F142B
         -1.0, & !  none                                  , CFC113
         -1.0, & !  none                                  , F141B
         -1.0, & !  none                                  , CH3OH
         -1.0, & !  none                                  , CH3CNPL
         -1.0, & !  none                                  , C2H6PL
         -1.0, & !  none                                  , PAN
         -1.0, & !  none                                  , CH3CHO
         -1.0, & !  none                                  , CH3CN
         -1.0, & !  none                                  , OTHER       70
         -1.0, & !  none                                  , CH3COOH
         -1.0, & !  none                                  , C5H8
         -1.0, & !  none                                  , MVK
         -1.0, & !  none                                  , MACR
         -1.0, & !  none                                  , C3H6
         -1.0, & !  none                                  , C4H8
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER       80
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER       90
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0, & !  none                                  , OTHER
         -1.0/   !  none                                  , OTHER       99
!
!
! T* values for Diffusion coefficent temperature changes
! See J.O. Hirschfelder...Table I-M
!
      data (tstar(i),i=1,82)/ 0.3000, 0.3500, 0.4000, 0.4500, 0.5000, 0.5500, &
         0.6000, 0.6500, 0.7000, 0.7500, 0.8000, 0.8500, 0.9000, 0.9500, 1.0000&
         , 1.0500, 1.1000, 1.1500, 1.2000, 1.2500, 1.3000, 1.3500, 1.4000, &
         1.4500, 1.5000, 1.5500, 1.6000, 1.6500, 1.7000, 1.7500, 1.8000, 1.8500&
         , 1.9000, 1.9500, 2.0000, 2.1000, 2.2000, 2.3000, 2.4000, 2.5000, &
         2.6000, 2.7000, 2.8000, 2.9000, 3.0000, 3.1000, 3.2000, 3.3000, 3.4000&
         , 3.5000, 3.6000, 3.7000, 3.8000, 3.9000, 4.0000, 4.1000, 4.2000, &
         4.3000, 4.4000, 4.5000, 4.6000, 4.7000, 4.8000, 4.9000, 5.0000, 6.0000&
         , 7.0000, 8.0000, 9.0000, 10.000, 20.000, 30.000, 40.000, 50.000, &
         60.000, 70.000, 80.000, 90.000, 100.00, 200.00, 300.00, 400.00/
!
! Omega* values for Diffusion coefficent temperature changes
! See J.O. Hirschfelder...Table I-M
!
      data (omegastar(i),i=1,82)/ 2.6620, 2.4670, 2.3180, 2.1840, 2.0660, &
         1.9660, 1.8770, 1.7980, 1.7290, 1.6670, 1.6120, 1.5620, 1.5170, 1.4760&
         , 1.4390, 1.4060, 1.3750, 1.3460, 1.3200, 1.2960, 1.2730, 1.2530, &
         1.2330, 1.2150, 1.1980, 1.1820, 1.1670, 1.1530, 1.1400, 1.1280, 1.1160&
         , 1.1050, 1.0940, 1.0840, 1.0750, 1.0570, 1.0410, 1.0260, 1.0120, &
         0.9996, 0.9878, 0.9770, 0.9672, 0.9576, 0.9490, 0.9406, 0.9328, 0.9256&
         , 0.9186, 0.9120, 0.9058, 0.8998, 0.8942, 0.8888, 0.8836, 0.8788, &
         0.8740, 0.8694, 0.8652, 0.8610, 0.8568, 0.8530, 0.8492, 0.8456, 0.8422&
         , 0.8124, 0.7896, 0.7712, 0.7556, 0.7424, 0.6641, 0.6232, 0.5960, &
         0.5756, 0.5596, 0.5464, 0.5352, 0.5256, 0.5130, 0.4644, 0.4360, 0.4170&
         /





!    Get Omega at T0 (ie 296.K)

      if (epsilon(mol) < 0) then
         write (6, *) mol, ' Molecule not supported in GALATRY/BETAT routines'
         CALL SHUTDOWN
         STOP 3
         !stop 'Molecule not supported'
      endif
      treduced = 296./sqrt(epsilon(mol)*epsilonair)
!      write (6, *) 'TReduced = ', treduced

!     Linear interpolate Omega from Tables @ TReduced


      result = tablelinterp(82,tstar,omegastar,treduced,omegat0)
      if (result /= 0) then
         select case (result)
         case (-1)
            write (6, *) 'Temperature ', treduced, ' is lower than table'
            !stop 'Temp too low'
         case (1)
            write (6, *) 'Temperature ', treduced, ' is above table'
            !stop 'Temp too high'
         case (-2)
            write (6, *) 'Didn''t find greater temperature in table'
            !stop 'Temp table error'
         end select
         CALL SHUTDOWN
         STOP 3
      endif



!     Get Omega at T

      treduced = t/sqrt(epsilon(mol)*epsilonair)

!     Linear interpolate Omega from Tables @ TReduced


      result = tablelinterp(82,tstar,omegastar,treduced,omegat)
      if (result /= 0) then
         select case (result)
         case (-1)
            write (6, *) 'Temperature ', treduced, ' is lower than table'
            !stop 'Temp too low'
         case (1)
            write (6, *) 'Temperature ', treduced, ' is above table'
            !stop 'Temp too high'
         case (-2)
            write (6, *) 'Didn''t find greater temperature in table'
            !stop 'Temp table error'
         end select
         CALL SHUTDOWN
         STOP 3
      endif


      betat = omegat/omegat0*sqrt(296./t)

      return
      end function betat

      DOUBLE COMPLEX FUNCTION MYWOFZ (X, Y)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DBLE_COMPLEX, DOUBLE
!...Translated by Pacific-Sierra Research 77to90  4.4E  13:05:37   1/16/07
!...Switches: -yf12
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE)  :: X
      REAL(DOUBLE)  :: Y
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: K, L
!-----------------------------------------------

      CALL HUMLIK (X, Y, K, L)
      MYWOFZ = DCMPLX(K,L)
      RETURN
      END FUNCTION MYWOFZ



      SUBROUTINE HUMLIK(X, Y, K, L)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  DOUBLE
!...Translated by Pacific-Sierra Research 77to90  4.4E  13:05:37   1/16/07
!...Switches: -yf12
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: X  ! Input x
      REAL(DOUBLE) , INTENT(IN) :: Y  ! y value >=0.0
      REAL(DOUBLE) , INTENT(OUT) :: K ! Real (Voigt)
      REAL(DOUBLE) , INTENT(OUT) :: L ! Optional
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
     ! Region boundaries for R=4
      REAL(DOUBLE), PARAMETER :: R0 = 146.7
      REAL(DOUBLE), PARAMETER :: R1 = 14.67
      REAL(DOUBLE), PARAMETER :: RRTPI = 0.56418958 !! 1/SQRT(pi)
     ! for CPF12 algorithm
      REAL(DOUBLE), PARAMETER :: Y0 = 1.5
      REAL(DOUBLE), PARAMETER :: Y0PY0 = Y0 + Y0
      REAL(DOUBLE), PARAMETER :: Y0Q = Y0*Y0
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J !I                 ! Loop variables
      INTEGER :: RG1, RG2, RG3        ! y polynomial flags

      REAL(DOUBLE), DIMENSION(0:5) :: C, S, T

      ! |x|, x^2, y^2, y/SQRT(pi)
      REAL(DOUBLE) :: ABX, XQ, YQ, YRRTPI

      ! |x| on region boundaries
      REAL(DOUBLE) :: XLIM0, XLIM1, XLIM2, XLIM3, XLIM4

      REAL(DOUBLE) :: A0, D0, D2, E0, E2, E4, H0, H2, H4, &
                      H6, P0, P2, P4, P6, P8, Z0, Z2, Z4, &
                    Z6, Z8, B1, F1, F3, F5, Q1, Q3, Q5, Q7

      ! CPF12 temporary valus
      REAL(DOUBLE), DIMENSION(0:5) :: XP, XM, YP, YM, MQ, PQ, MF, PF

      REAL(DOUBLE) :: D, YF, YPY0, YPY0Q, PREVY

!     SAVE preserves values of C, S and T (static) arrays and all variables
!       except I and J between procedure calls
      SAVE C, S, T, RG1, RG2, RG3, ABX, XQ, YQ, YRRTPI, XLIM0, XLIM1, XLIM2, &
         XLIM3, XLIM4, A0, D0, D2, E0, E2, E4, H0, H2, H4, H6, P0, P2, P4, P6, &
         P8, Z0, Z2, Z4, Z6, Z8, B1, F1, F3, F5, Q1, Q3, Q5, Q7, XP, XM, YP, YM&
         , MQ, PQ, MF, PF, D, YF, YPY0, YPY0Q, PREVY
!-----------------------------------------------
!
!     To calculate the Faddeeva function with relative error less than 10^(-R).
!     R0=1.51*EXP(1.144*R) and R1=1.60*EXP(0.554*R) can be set by the the user
!     subject to the constraints 14.88<R0<460.4 and 4.85<R1<25.5
!
!
! Arguments
!
! Constants

      DATA C/ 1.0117281, -0.75197147, 0.012557727, 0.010022008, -0.00024206814&
         , 0.00000050084806/
      DATA S/ 1.393237, 0.23115241, -0.15535147, 0.0062183662, 0.000091908299, &
          - 0.00000062752596/
      DATA T/ 0.31424038, 0.94778839, 1.5976826, 2.2795071, 3.0206370, &
         3.8897249/
!
! Local variables
!
!
!
!
!
!
!
!
!

!**** Start of executable code *****************************************
      IF (ABS(Y - PREVY)<1.E-7 .OR. PREVY==0) THEN
         RG1 = 1                                 ! Set flags
         RG2 = 1
         RG3 = 1
         YQ = Y*Y                                ! y^2
         YRRTPI = Y*RRTPI                        ! y/SQRT(pi)

!      Region boundaries when both K and L are required or when R<>4
         XLIM0 = R0 - Y
         XLIM1 = R1 - Y
         XLIM3 = 3.097*Y - 0.45
!      For speed the following 3 lines should replace the 3 above if R=4 and L is not required   *
!      XLIM0 = 15100.0 + Y*(40.0 + Y*3.6)
!      XLIM1 = 164.0 - Y*(4.3 + Y*1.8)
!      XLIM3 = 5.76*YQ

         XLIM2 = 6.8 - Y
         XLIM4 = 18.1*Y + 1.65
         IF (Y <= 0.000001) THEN                 ! When y<10^-6
            XLIM1 = XLIM0                        ! avoid W4 algorithm
            XLIM2 = XLIM0
         ENDIF
      ENDIF
      PREVY = Y
!.....
      ABX = ABS(X)                               ! |x|
      XQ = ABX*ABX                               ! x^2
      IF (ABX > XLIM0) THEN
         K = YRRTPI/(XQ + YQ)                    ! Region 0 algorithm
         L = K*X/Y

      ELSE IF (ABX > XLIM1) THEN                 ! Humlicek W4 Region 1
         IF (RG1 /= 0) THEN                      ! First point in Region 1
            RG1 = 0
            A0 = YQ + 0.5                        ! Region 1 y-dependents
            D0 = A0*A0
            D2 = YQ + YQ - 1.0
            B1 = YQ - 0.5
         ENDIF
         D = RRTPI/(D0 + XQ*(D2 + XQ))
         K = D*Y*(A0 + XQ)
         L = D*X*(B1 + XQ)

      ELSE IF (ABX > XLIM2) THEN                 ! Humlicek W4 Region 2
         IF (RG2 /= 0) THEN                      ! First point in Region 2
            RG2 = 0
                                                 ! Region 2 y-dependents
            H0 = 0.5625 + YQ*(4.5 + YQ*(10.5 + YQ*(6.0 + YQ)))
            H2 = (-4.5) + YQ*(9.0 + YQ*(6.0 + YQ*4.0))
            H4 = 10.5 - YQ*(6.0 - YQ*6.0)
            H6 = (-6.0) + YQ*4.0
            E0 = 1.875 + YQ*(8.25 + YQ*(5.5 + YQ))
            E2 = 5.25 + YQ*(1.0 + YQ*3.0)
            E4 = 0.75*H6
            F1 = (-1.875) + YQ*(5.25 + YQ*(4.5 + YQ))
            F3 = 8.25 - YQ*(1.0 - YQ*3.0)
            F5 = (-5.5) + YQ*3.0
         ENDIF
         D = RRTPI/(H0 + XQ*(H2 + XQ*(H4 + XQ*(H6 + XQ))))
         K = D*Y*(E0 + XQ*(E2 + XQ*(E4 + XQ)))
         L = D*X*(F1 + XQ*(F3 + XQ*(F5 + XQ)))

      ELSE IF (ABX < XLIM3) THEN                 ! Humlicek W4 Region 3
         IF (RG3 /= 0) THEN                      ! First point in Region 3
            RG3 = 0
                                                 ! Region 3 y-dependents
            Z0 = 272.1014 + Y*(1280.829 + Y*(2802.870 + Y*(3764.966 + Y*(&
               3447.629 + Y*(2256.981 + Y*(1074.409 + Y*(369.1989 + Y*(88.26741&
                + Y*(13.39880 + Y)))))))))
            Z2 = 211.678 + Y*(902.3066 + Y*(1758.336 + Y*(2037.310 + Y*(&
               1549.675 + Y*(793.4273 + Y*(266.2987 + Y*(53.59518 + Y*5.0))))))&
               )
            Z4 = 78.86585 + Y*(308.1852 + Y*(497.3014 + Y*(479.2576 + Y*(&
               269.2916 + Y*(80.39278 + Y*10.0)))))
            Z6 = 22.03523 + Y*(55.02933 + Y*(92.75679 + Y*(53.59518 + Y*10.0)))
            Z8 = 1.496460 + Y*(13.39880 + Y*5.0)
            P0 = 153.5168 + Y*(549.3954 + Y*(919.4955 + Y*(946.8970 + Y*(&
               662.8097 + Y*(328.2151 + Y*(115.3772 + Y*(27.93941 + Y*(4.264678&
                + Y*0.3183291))))))))
            P2 = (-34.16955) + Y*(-1.322256 + Y*(124.5975 + Y*(189.7730 + Y*(&
               139.4665 + Y*(56.81652 + Y*(12.79458 + Y*1.2733163))))))
            P4 = 2.584042 + Y*(10.46332 + Y*(24.01655 + Y*(29.81482 + Y*(&
               12.79568 + Y*1.9099744))))
            P6 = (-0.07272979) + Y*(0.9377051 + Y*(4.266322 + Y*1.273316))
            P8 = 0.0005480304 + Y*0.3183291
            Q1 = 173.2355 + Y*(508.2585 + Y*(685.8378 + Y*(557.5178 + Y*(&
               301.3208 + Y*(111.0528 + Y*(27.62940 + Y*(4.264130 + Y*0.3183291&
               )))))))
            Q3 = 18.97431 + Y*(100.7375 + Y*(160.4013 + Y*(130.8905 + Y*(&
               55.88650 + Y*(12.79239 + Y*1.273316)))))
            Q5 = 7.985877 + Y*(19.83766 + Y*(28.88480 + Y*(12.79239 + Y*&
               1.909974)))
            Q7 = 0.6276985 + Y*(4.264130 + Y*1.273316)
         ENDIF
         D = 1.7724538/(Z0 + XQ*(Z2 + XQ*(Z4 + XQ*(Z6 + XQ*(Z8 + XQ)))))
         K = D*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))))
         L = D*X*(Q1 + XQ*(Q3 + XQ*(Q5 + XQ*(Q7 + XQ*0.3183291))))

      ELSE                                       ! Humlicek CPF12 algorithm
         YPY0 = Y + Y0
         YPY0Q = YPY0*YPY0
         L = 0.0
         DO J = 0, 5
            D = X - T(J)
            MQ(J) = D*D
            MF(J) = 1.0/(MQ(J)+YPY0Q)
            XM(J) = MF(J)*D
            YM(J) = MF(J)*YPY0
            D = X + T(J)
            PQ(J) = D*D
            PF(J) = 1.0/(PQ(J)+YPY0Q)
            XP(J) = PF(J)*D
            YP(J) = PF(J)*YPY0
            L = L + C(J)*(XM(J)+XP(J)) + S(J)*(YM(J)-YP(J))
         END DO

         IF (ABX <= XLIM4) THEN                  ! Humlicek CPF12 Region I
            K = SUM(C(:5)*(YM(:5)+YP(:5))-S(:5)*(XM(:5)-XP(:5)))

         ELSE                                    ! Humlicek CPF12 Region II
            YF = Y + Y0PY0
!            K = K + SUM((C(:5)*(MQ(:5)*MF(:5)-Y0*YM(:5))+S(:5)*YF*XM(:5))/(MQ(:&
            K = SUM((C(:5)*(MQ(:5)*MF(:5)-Y0*YM(:5))+S(:5)*YF*XM(:5))/(MQ(:&
               5)+Y0Q)+(C(:5)*(PQ(:5)*PF(:5)-Y0*YP(:5))-S(:5)*YF*XP(:5))/(PQ(:5&
               )+Y0Q))
            K = Y*K + EXP((-XQ))
         ENDIF
      ENDIF
      RETURN
!.....
      END SUBROUTINE HUMLIK
      real(double) function galatry (x, y, z)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE params, ONLY:  double, dble_complex
!...Translated by Pacific-Sierra Research 77to90  4.4E  13:05:37   1/16/07
!...Switches: -yf12
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double)  :: x
      real(double)  :: y
      real(double) , intent(in) :: z
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: itot, nn, n, i
      real(double) :: rpi, ylim, delta, rn, r1, r2, r3, r4, r5, r6
      complex(dble_complex) :: teta, den, g
      complex(dble_complex), dimension(100) :: b
      complex(dble_complex) :: w, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, &
         t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, &
         t27, t28
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!      real , external :: dreal
      complex(dble_complex) , external :: mywofz
!-----------------------------------------------
!
!
!
!     Returns Galatry profile value at x,y,z where
!
!
!
!     x = sqrt(ln(2.))/gamma_D*(nu-nu0)
!
!     y = sqrt(ln(2.))/gamma_D*gamma_L
!
!     z = sqrt(ln(2.))/gamma_D*beta
!
!
!
!     Must be normalized by multiplying returned value by
!
!     sqrt(ln(2.)/Pi)/gamma_D to have unit area.
!
!     Call for Galatry with z=0 should return Voigt profile.
!
!
!
!     Calls MyWofZ(x,y) which returns the complex error function.
!
!
!
!     Adapted from F. Herbert, JQSRT vol 14,p 943 (1974)
!
!
!
!
!
!
!
!


      rpi = sqrt(3.14159265358979)

      w = dcmplx(y,(-x))

      ylim = 4.*z**0.868

      if (z>0.1D+00 .and. y<0.5D+00 .or. z>0.1D+00 .and. y<ylim .or. z>5.D+00&
          .and. y>12.D+00) then

         go to 2

      else if (z<0.04D+00 .and. y<0.5D+00) then

         w = mywofz(x,y)

         teta = dcmplx(x,y)

         r1 = z**2

         r2 = r1*z

         t3 = teta**2

         t4 = t3**2

         t5 = t4*t3

         r3 = 1/rpi

         r4 = r1**2

         r5 = r4*r1

         r6 = r4*z

         t6 = 2.D0/9.D0*r2*t5*r3 - r5*t4*w/6 + 2.D0/3.D0*z*r3 + w + 83.D0/96.D0&
            *r4*w + 8.D0/105.D0*r6*r3 - 16.D0/15.D0*r2*r3 + r1*w/6 - r5*w/48 + &
            r5*t3*w/6 + 2.D0/45.D0*r5*t5*w
         if (z == 0.0) then
      galatry = DBLE(t6)
                return
         endif

         t7 = t4**2

         t8 = (-29.D0/105.D0*r6*t3*r3) + 8.D0/63.D0*r6*t4*r3 - 4.D0/315.D0*r6*&
            t5*r3 + 13.D0/90.D0*r4*t7*w - 2.D0/9.D0*r1*t5*w - 85.D0/12.D0*r4*t3&
            *w + 29.D0/4.D0*r4*t4*w - 89.D0/45.D0*r4*t5*w + 127.D0/30.D0*r2*t3*&
            r3 - 94.D0/45.D0*r2*t4*r3 - 3.D0/2.D0*r1*t3*w + 4.D0/3.D0*r1*t4*w

         t9 = dcmplx(0.D0,1.D0)

         t10 = t9*r1

         t11 = teta*r3

         t12 = t3*teta

         t13 = t12*r3

         t14 = t9*z

         t15 = t12*w

         t16 = teta*w

         t17 = t9*r4

         t18 = t4*teta

         t19 = t18*r3

         t20 = t9*r2

         t21 = t18*w

         t22 = (-2.D0/3.D0*z*t3*r3) - r5*t7*w/315 + t10*t11 - 11.D0/9.D0*t10*&
            t13 - 2.D0/3.D0*t14*t15 + t14*t16 - 2293.D0/360.D0*t17*t13 + 343.D0&
            /180.D0*t17*t19 - 29.D0/12.D0*t20*t16 + 31.D0/6.D0*t20*t15 - 11.D0/&
            5.D0*t20*t21

         t23 = t4*t12

         t24 = t23*w

         t25 = t9*r6

         t26 = t23*r3

         t27 = t9*r5

         t28 = 2.D0/9.D0*t20*t24 + t25*t16/6 - t25*t15/3 + 2.D0/15.D0*t25*t21&
             - 4.D0/315.D0*t25*t24 + 2.D0/9.D0*t10*t19 + 1121.D0/240.D0*t17*t11&
             - 13.D0/90.D0*t17*t26 - 31.D0/280.D0*t27*t11 + 37.D0/252.D0*t27*&
            t13 - 3.D0/70.D0*t27*t19 + t27*t26/315.D0

!         galatry = dreal(t6 + t8 + t22 + t28)
         galatry = DBLE(t6 + t8 + t22 + t28)

         return

      endif
     ! added nint nov09
      itot = nint( 3. + 37*dexp((-0.6*y)) )

      b(1) = (itot - 1)*z + w

      do i = 2, itot

         b(i) = b(i-1)**(-1)*(itot - i + 1)/2 + (itot - i)*z + w

      end do

      g = 1/(b(itot)*rpi)

!      galatry = dreal(g)
      galatry = DBLE(g)

      return

    2 continue
    ! added nint nov09
      itot = nint( 4. + z**(-1.05)*(1 + 3*exp((-1.1*y))) )

      delta = 1./(2.*z*z)

      teta = dcmplx(delta + y/z,(-x/z))

!      g = (0.E+00,0.E+00)_dble_complex
      g = (0.E+00,0.E+00)

!      den = (1.E+00,0.E+00)_dble_complex
      den = (1.E+00,0.E+00)

      do nn = 1, itot

         n = nn - 1

         rn = dfloat(n)

         den = den*(teta + rn)

         g = g + delta**n/den

      end do

      g = g/(z*rpi)

!      galatry = dreal(g)
      galatry = DBLE(g)

      return

      end function galatry
