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

      REAL(DOUBLE) FUNCTION VOIGT (X, Y)

!      REAL(KIND(0.0D0)) FUNCTION VOIGT (X, Y)

      USE PARAMS

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: X
      REAL(DOUBLE) , INTENT(IN) :: Y
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                     :: I, J, MAXV, MINV, N
      REAL(DOUBLE), DIMENSION(22), SAVE :: B
      REAL(DOUBLE), DIMENSION(15), SAVE :: RI
      REAL(DOUBLE), DIMENSION(15) :: XN, YN
      REAL(DOUBLE), DIMENSION(25), SAVE :: D0, D1, D2, D3, D4, HN
      REAL(DOUBLE), DIMENSION(3)  :: XX, HH
      REAL(DOUBLE), DIMENSION(19) :: NBY2
      REAL(DOUBLE), DIMENSION(21) :: C
      REAL(DOUBLE), SAVE :: CO
      REAL(DOUBLE)       :: H, UU, VV, U, Y2, DX, V
      LOGICAL            :: TRU
!-----------------------------------------------
! ***
! ***   THIS ROUTINE FROM MICHIGAN COMPUTES THE VOIGT FUNCTION ***
! ***      Y/PI*(INTEGRAL (EXP(-T*T)/(Y*Y+(X-T*(X-T))           ***
! ***      IT IS A REVISED VERSION OF DRAYSON'S ORIGINAL      ***
! ***
!      REAL*8 B(22),RI(15),XN(15),YN(15)
!      REAL*8 D0(25),D1(25),D2(25),D3(25),D4(25),HN(25)
!      REAL*8 XX(3),HH(3),NBY2(19),C(21)
      DATA B/ 0.D0, .7093602D-7, 20*0.D+0/
      DATA XN/ 10.D0, 9.D0, 2*8.D0, 7.D0, 6.D0, 5.D0, 4.D0, 7*3.D0/
      DATA YN/ 3*.6D0, .5D0, 2*.4D0, 4*.3D0, 1.D0, .9D0, .8D0, 2*.7D0/
      DATA H/ .201D0/
      DATA XX/ .5246476D0, 1.65068D0, .7071068D0/
      DATA HH/ .2562121D0, .2588268D-1, .2820948D0/
      DATA NBY2/ 9.5D0, 9.D0, 8.5D0, 8.D0, 7.5D0, 7.D0, 6.5D0, 6.D0, 5.5D0, &
         5.D0, 4.5D0, 4.0D0, 3.5D0, 3.0D0, 2.5D0, 2.D0, 1.5D0, 1.D0, .5D0/
      DATA C/ .7093602D-7, -.2518434D-6, .8566874D-6, -.2787638D-5, .866074D-5&
         , -.2565551D-4, .7228775D-4, -.1933631D-3, .4899520D-3, -.1173267D-2, &
         .2648762D-2, -.5623190D-2, .1119601D-1, -.2084976D-1, .3621573D-1, &
         -.5851412D-1, .8770816D-1, -.121664D0, .15584D0, -.184D0, .2D0/
      DATA TRU/ .FALSE./
      IF (.NOT.TRU) THEN
! ***
! ****** REGION I.      COMPUTE DAWSON-S FUNCTION AT MESH POINTS ******
! ***
         TRU = .TRUE.
         DO I = 1, 15
            RI(I) = -I/2.D0
         END DO
         DO I = 1, 25
            HN(I) = H*(I - .5D0)
            CO = 4.D0*HN(I)*HN(I)/25.D0 - 2.D0
            DO J = 2, 21
               B(J+1) = CO*B(J) - B(J-1) + C(J)
            END DO
            D0(I) = HN(I)*(B(22)-B(21))/5.D0
            D1(I) = 1.D0 - 2.D0*HN(I)*D0(I)
            D2(I) = (HN(I)*D1(I)+D0(I))/RI(2)
            D3(I) = (HN(I)*D2(I)+D1(I))/RI(3)
            D4(I) = (HN(I)*D3(I)+D2(I))/RI(4)
         END DO
      ENDIF
      IF (X - 5.D0 < 0.D0) THEN
         IF (Y - 1.D0 > 0.D0) THEN
            IF (X > 1.85D0*(3.6D0 - Y)) GO TO 112
! ***
! ****** REGION II      CONTINUED FRACTION.  COMPUTE # OF TERMS ******
! ***
            IF (Y >= 1.45D0) THEN
               !I = Y + Y
               I = FLOOR(Y + Y)
               !print *, 'a', i, floor(y+y)
            ELSE
               !I = 11.D0*Y
               I = FLOOR(11.D0 * Y)
               !print *, 'b', i, floor(11.d0*y)
            ENDIF
            !J = X + X + 1.85D0
            J = FLOOR(X + X + 1.85D0)
            !print *, 'c', j, floor(X + X + 1.85D0)
            !MAXV = XN(J)*YN(I) + .46D0
            MAXV = FLOOR(XN(J)*YN(I) + .46D0)
            !print *, 'd', maxv, floor(XN(J)*YN(I) + .46D0)
            MINV = MIN0(16,21 - 2*MAXV)
! ***
! ******   EVALUATE CONTINED FRACTION
! ***
            UU = Y
            VV = X
            DO J = MINV, 19
               U = NBY2(J)/(UU*UU + VV*VV)
               UU = Y + U*UU
               VV = X - U*VV
            END DO
            VOIGT = UU/(UU*UU + VV*VV)/1.772454D0
            RETURN
         ENDIF
         Y2 = Y*Y
         IF (X + Y > 5.D0) GO TO 113
! ***
! ****** REGION I.   COMPUTE DAWSON-S FUNCTION AT X FRON TAYLOR SERIES.
! ***
         !N = X/H
         N = FLOOR(X/H)
         !print *, 'e', n, floor(X/H)
         DX = X - HN(N+1)
         U = (((D4(N+1)*DX+D3(N+1))*DX+D2(N+1))*DX+D1(N+1))*DX + D0(N+1)
         V = 1.D0 - 2.D0*X*U
! ******
! ****** TAYLOR SERIES EXPANSION ABOUT Y=0.0   ******
! ******
!      VV=EXP(Y2-X*X)*COS(2.*X*Y)/1.128379D0-Y*V
         VV = EXP(Y2 - X*X)*COS(2.D0*X*Y)/1.128379D0 - Y*V
         UU = -Y
         !MAXV = 5.D0 + (12.5D0 - X)*.8D0*Y
         MAXV = FLOOR(5.D0 + (12.5D0 - X)*.8D0*Y)
         !print *, 'f', maxv, floor(5.D0 + (12.5D0 - X)*.8D0*Y)
         DO I = 2, MAXV, 2
            U = (X*V + U)/RI(I)
            V = (X*U + V)/RI(I+1)
            UU = -UU*Y2
            VV = VV + V*UU
         END DO
         VOIGT = 1.128379D0*VV
         RETURN
      ENDIF
  112 CONTINUE
      Y2 = Y*Y
      IF (Y >= 11.D0 - .6875D0*X) THEN
! ***
! ****** REGION IIIB.  2-POINT GAUSS-HERMITE QUADRATURE.  ******
! ***
         U = X - XX(3)
         V = X + XX(3)
         VOIGT = Y*(HH(3)/(Y2+U*U)+HH(3)/(Y2+V*V))
         RETURN
! ***
! ****** REGION IIIA.  4-POINT GAUS-HERMITE QUADRATURE   ******
! ***
      ENDIF
  113 CONTINUE
      U = X - XX(1)
      V = X + XX(1)
      UU = X - XX(2)
      VV = X + XX(2)
      VOIGT = Y*(HH(1)/(Y2+U*U)+HH(1)/(Y2+V*V)+HH(2)/(Y2+UU*UU)+HH(2)/(Y2+VV*VV&
         ))
      RETURN
      END FUNCTION VOIGT
