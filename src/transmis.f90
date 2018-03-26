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

      MODULE TRANSMIS

      USE params
      USE bandparam
      USE retvparam
      USE vibfcn
      USE xsections
      USE molcparam
      USE lineparam

      IMPLICIT NONE

! --- TCONV and TCALC now allocated in setup
      COMPLEX(DBLE_COMPLEX), DIMENSION(:), ALLOCATABLE :: TCONV
      REAL(DOUBLE), DIMENSION(:,:),   ALLOCATABLE :: TCALC
      REAL(DOUBLE), DIMENSION(:,:),   ALLOCATABLE :: TCALC_I !mp
      REAL(DOUBLE), DIMENSION(:,:,:), ALLOCATABLE :: TCALC_E !mp
      REAL(DOUBLE), DIMENSION(:,:,:), ALLOCATABLE :: TCALC_S !mp
      REAL(DOUBLE), DIMENSION(:,:,:), ALLOCATABLE :: CROSS_FACMAS !mp

      CONTAINS

!----------------------------------------------------------------------
      SUBROUTINE MASSPATH( K )

      INTEGER, INTENT(IN) :: K
      INTEGER             :: KK

!dt = T(k) scale factor
! T is already perturbed from torg
      if( k .ne. -1 )then
         !T(K) = TORG(K) * DT
         !do i=1, nspec

!print*, k, nspec, ccc(:nspec,k), T(K), TORG(k), CORG(:nspec,k)

            !CCC(i,K) = CCC(i,K) * (TORG(K) / DT * TORG(K) )

            ! reset  ccc
!            CCC(:nspec,:Kmax) = CCC(:nspec,:Kmax)
            ! update only layer k
!            CCC(:nspec,K) = CORG(:nspec,K) * (TORG(K) / T(K) )
            ! reset to last iteration and update k

            do kk = 1, kmax
               CCC(:nspec,kk) = CORG(:nspec,kk) * (TORG(kk) / T(kk) )
            enddo

!print*, k, nspec, ccc(:nspec,k), T(K), TORG(k), CORG(:nspec,k)

!print*, k, ccc(:nspec,k)
         !enddo
      endif

      DO KK = 1, KMAX
         PMASMX(KK) = 0.D0
         PMASMX(KK) = DMAX1(MAXVAL(CCC(:NSPEC,KK)),PMASMX(KK))
         !print *,  PMASMX(KK)
      END DO

      END SUBROUTINE MASSPATH


!----------------------------------------------------------------------
      SUBROUTINE TALL

!  --- MAKE APPROPRIATE CALLS TO TRANS SUBROUTINE TO COMPUTE ALL
!  --- MONCHROMATIC TRANSMITTANCES (IPARM=1 CALL)

      INTEGER :: MONONE, MXONE, IBAND, N

      MONONE = 1
      MXONE  = 1
!  --- COMPUTE MONOCHROMATIC TRANSMITTANCES FOR ALL SCANS
      DO IBAND = 1, NBAND

         N = NSCAN(IBAND)

         IF (N == 0) CYCLE
         CALL NTRAN (IBAND, 1, 1, MONONE, MXONE)

         MONONE = MONONE + NM(IBAND)*NSCAN(IBAND)
         MXONE  = MXONE  + NM(IBAND)

      END DO

!  --- COPY TRANSMISSION ARRAY

      TCALC  (2,:NMONSM)         = TCALC  (1,:NMONSM)
      TCALC_I(2,:NMONSM)         = TCALC_I(1,:NMONSM)
      TCALC_S(2,:NMONSM,:KMAX)   = TCALC_S(1,:NMONSM,:KMAX)
      TCALC_E(2,:NMONSM,:KMAX+1) = TCALC_E(1,:NMONSM,:KMAX+1)

      RETURN

      END SUBROUTINE TALL


!----------------------------------------------------------------------
      SUBROUTINE NTRAN(IBAND, JMIN, IPOINT, MONONE, MXONE)

!     COMPUTE MONOCHROMATIC TRANSMITTANCES AT EACH WAVELENGTH
!     FOR BANDPASS IBAND.  TRANMSITTANCES ARE COMPUTED FOR SPECTRA
!     BETWEEN SPECTRUM ISCAN(IBAND,JMIN) AND ISCAN(IBAND,NSCAN(IBAND)).
!
!          IBAND=BAND PASS NUMBER
!          JSCAN=ISCAN(IBAND,J)=SPECTRUM NUMBER FOR CALCULATIONS
!          IPOINT,MONONE-INDICES IN TCALC ARRAY FOR CALCULATED TRANSMITTANCES

      INTEGER, INTENT(IN) :: IBAND, JMIN, IPOINT, MONONE, MXONE

      INTEGER :: NMON, NSCANS, INDXX, KSMAX2, K, JSCAN, IR, ICINDX, ICINDX2
      INTEGER :: MSTOR, MXMAX, J, MADD, I, ALT
      REAL(DOUBLE) :: FACMAS, WAVE_NR

!  --- NMON=NUMBER OF MONOCHROMATIC POINTS FOR THE BANDPASS CALCULATION
      NMON   = NM(IBAND)
      NSCANS = NSCAN(IBAND)

!  --- ZERO APPROPRIATE TRANSMISSION ARRAY ELEMENTS FOR CROSS SECTION
!  ---  CALCULATIONS
      MADD = MONONE
      DO INDXX = 1, NSCANS
         IF (INDXX >= JMIN) THEN
            TCALC  (IPOINT,MADD:NMON-1+MADD)         = 0.D0
            TCALC_I(IPOINT,MADD:NMON-1+MADD)         = 0.D0
            TCALC_E(IPOINT,MADD:NMON-1+MADD,:KMAX+1) = 0.D0
            TCALC_S(IPOINT,MADD:NMON-1+MADD,:KMAX)   = 0.D0
         ENDIF
         MADD = MADD + NM(IBAND)
      END DO

!  --- MAXIMUM LAYER FOR SUMMING CROSS SECTIONS CALCULATIONS

      KSMAX2 = KZTAN(ISCAN(IBAND,NSCANS))

!                   ------------LOOP OVER LAYERS
      DO K = 1, KSMAX2
         MADD = MONONE

         ! ------------LOOP OVER SPECTRA
         DO INDXX = 1, NSCANS

            IF (INDXX >= JMIN) THEN
               JSCAN = ISCAN(IBAND,INDXX) ! jscan picks out mass paths in eg ccc not kscan

               IF (K <= KZTAN(JSCAN)) THEN

                  FACMAS = CCC(JSCAN,K)/PMASMX(K)

                  ! ------------LOOP OVER FREQUENCIES
                  MXMAX = MXONE + NMON - 1
                  DO J = 1, NMON

                     ICINDX = MXONE + J - 1
                     MSTOR  = MADD  + J - 1

                     ! --- DON'T APPLY SHIFT TO FIRST POINT
                     CROSS_FACMAS(1,K,MSTOR) = CROSS(1,K,ICINDX) * FACMAS

                     TCALC(IPOINT,MSTOR) = TCALC(IPOINT,MSTOR) + (X(1,K)/XORG(1,K)) * CROSS_FACMAS(1,K,MSTOR)

                     IF (IEMISSION/=0) THEN
                        ! Transmission calculated below the layer ALT, needed
                        ! for calculation of contribution to emission from
                        ! layer ALT to the ground
                        DO ALT=1,KSMAX2
                           IF (ZBAR(ALT) > ZBAR(K)) THEN
                              TCALC_E(IPOINT,MSTOR,ALT) = &
                              TCALC_E(IPOINT,MSTOR,ALT) + (X(1,K)/XORG(1,K))*CROSS_FACMAS(1,K,MSTOR)
                           ENDIF
                        ENDDO
                     ENDIF

                     IF (IFDIFF) THEN
                        ! ------------LOOP OVER RETRIEVAL GASES
                        DO IR = 2, NRET
                           ! --- APPLY DIFFERENTIAL WAVENUMBER SHIFT
                           ICINDX2 = ICINDX + ISHIFT(IR-1)
                           ! --- FIXUP AT ENDPOINTS
                           ICINDX2 = MAX0(MXONE,ICINDX2)
                           ICINDX2 = MIN0(MXMAX,ICINDX2)
                           CROSS_FACMAS(IR,K,MSTOR) = CROSS(IR,K,ICINDX2)*FACMAS

                           TCALC(IPOINT,MSTOR) = &
                                TCALC(IPOINT,MSTOR) + (X(IR,K)/XORG(IR,K))*CROSS_FACMAS(IR,K,MSTOR)
                           IF (IEMISSION/=0) THEN
                              TCALC_E(IPOINT,MSTOR,KSMAX2) = 1.0D0
                              DO ALT=1,KSMAX2
                 ! Transmission calculated below the layer ALT, needed
                 ! for calculation of contribution to emission from
                 ! layer ALT to the ground
                                 IF (ZBAR(ALT) > ZBAR(K)) THEN
                                    TCALC_E(IPOINT,MSTOR,ALT) = &
                                    TCALC_E(IPOINT,MSTOR,ALT) + (X(IR,K)/XORG(IR,K)) * CROSS_FACMAS(IR,K,MSTOR)
                                 ENDIF
                              ENDDO
                           ENDIF
                        ENDDO

                        ! ------------BACKGROUND GASES
                        CROSS_FACMAS(NRET+1,K,MSTOR) = CROSS(NRET+1,K,ICINDX2)*FACMAS
                        TCALC(IPOINT,MSTOR) = TCALC(IPOINT,MSTOR) + CROSS_FACMAS(NRET+1,K,MSTOR)

                        IF (IEMISSION/=0) THEN
                           DO ALT=1,KSMAX2
                              IF (ZBAR(ALT) > ZBAR(K)) THEN
                                 TCALC_E(IPOINT,MSTOR,ALT) = &
                                 TCALC_E(IPOINT,MSTOR,ALT) + CROSS_FACMAS(NRET+1,K,MSTOR)
                              ENDIF
                           ENDDO
                        ENDIF
                     ELSE
                        ! ------------LOOP OVER RETRIEVAL GASES
                        DO IR = 2, NRET
                           CROSS_FACMAS(IR,K,MSTOR) = CROSS(IR,K,ICINDX)*FACMAS
                           TCALC(IPOINT,MSTOR) = &
                           TCALC(IPOINT,MSTOR) + (X(IR,K)/XORG(IR,K)) * CROSS_FACMAS(IR,K,MSTOR)
                           IF (IEMISSION/=0) THEN
                 ! Transmission calculated below the layer ALT, needed
                 ! for calculation of contribution to emission from
                 ! layer ALT to the ground
                              DO ALT=1,KSMAX2
                                 IF (ZBAR(ALT) > ZBAR(K)) THEN
                                    TCALC_E(IPOINT,MSTOR,ALT) = &
                                    TCALC_E(IPOINT,MSTOR,ALT) + (X(IR,K)/XORG(IR,K)) * CROSS_FACMAS(IR,K,MSTOR)
                                 END IF
                              END DO
                           END IF
                        END DO

                        ! ------------BACKGROUND GASES
                        CROSS_FACMAS(NRET+1,K,MSTOR) = CROSS(NRET+1,K,ICINDX)*FACMAS
                        TCALC(IPOINT,MSTOR) = TCALC(IPOINT,MSTOR) + CROSS_FACMAS(NRET+1,K,MSTOR)
                        IF (IEMISSION/=0) THEN
                          DO ALT=1,KSMAX2
                 ! Transmission calculated below the layer ALT, needed
                 ! for calculation of contribution to emission from
                 ! layer ALT to the ground
                             IF (ZBAR(ALT) > ZBAR(K)) THEN
                                TCALC_E(IPOINT,MSTOR,ALT) = &
                                TCALC_E(IPOINT,MSTOR,ALT) + CROSS_FACMAS(NRET+1,K,MSTOR)
                             END IF
                          END DO
                       END IF
                    ENDIF
                 END DO
              ENDIF
           ENDIF
           MADD = MADD + NM(IBAND)
        END DO
     END DO
     !  --- COMPUTE MONOCHROMATIC TRANSMITTANCES FROM CROSS SECTION SUMS
     MADD = MONONE
     IF (IEMISSION/=0) THEN
        !--- COMPUTE MONOCHROMATIC RADIATION CROSS SECTIONS FOR EMISSION
        MADD = MONONE
        DO INDXX = 1, NSCANS
           IF (INDXX >= JMIN) THEN
              DO I = 1, NM(IBAND)
                 WAVE_NR = WSTART(IBAND) + (I-1)*DN(IBAND)

                 IF (EMISSION_OBJECT.EQ.'M') THEN
                    TCALC(IPOINT, MADD+I-1) &
                         = (PLANCK(WAVE_NR,EMISSION_T_BACK) + 1.D-6 *PLANCK(WAVE_NR ,6000.0D0)) &
                         * EXP((-TCALC(IPOINT,MADD+I-1)))
                 ELSE
                    TCALC(IPOINT, MADD+I-1) &
                         = PLANCK(WAVE_NR,EMISSION_T_BACK) &
                         * EXP((-TCALC(IPOINT,MADD+I-1)))
                 END IF
                 DO K=1, KSMAX2
                    IF( ABS( TCALC_E(IPOINT,MADD+I-1,K)) .GT. 664.0 ) THEN
                       ! LIMIT SO ONLY GET EXPONENT < 300
                       TCALC_E(IPOINT,MADD+I-1,K) = 664.0D0
                    ENDIF
                    ! Transmission from altitude K to the ground
                    TCALC_E(IPOINT,MADD+I-1,K) = EXP(-TCALC_E(IPOINT,MADD+I-1,K))
                 END DO
                 ! KSMAX2+1 MEANS INFINITY NOT UNDERGROUND !!!!
                 TCALC_E(IPOINT,MADD+I-1,KMAX+1) = TCALC(IPOINT, MADD+I-1)
                 TCALC_S(IPOINT, MADD+I-1, 1)=0.D0
                 DO K=2,KSMAX2
                    ! Calculates the spectrum, Planck is the emission of a black body at the temperature
                    ! of T(K) in layer K,
                    ! TCALC_E(IPOINT,MADD+I-1,K) - TCALC_E(IPOINT,MADD+I-1,K-1) is the absorption
                    ! of the layer K alone
                    TCALC_S(IPOINT, MADD+I-1, K) = PLANCK(WAVE_NR,T(K))&
                         *(TCALC_E(IPOINT,MADD+I-1,K) - TCALC_E(IPOINT,MADD+I-1,K-1))
                    ! TCALC contains the spectrum
                    TCALC(IPOINT, MADD+I-1) = TCALC(IPOINT, MADD+I-1) &
                         + TCALC_S(IPOINT, MADD+I-1, K)
                    ! TCALC_S WILL BE USED FOR COMPUTING THE DERIVATIVES
                    TCALC_S(IPOINT, MADD+I-1, K) = TCALC_S(IPOINT, MADD+I-1, K)&
                         + TCALC_S(IPOINT, MADD+I-1, K-1)
                 END DO
                 ! Contribution to spectrum from layer K is kept for calculating the derivatives
                 DO K=2,KSMAX2
                    TCALC_E(IPOINT,MADD+I-1,K) = PLANCK(WAVE_NR,T(K))*(TCALC_E(IPOINT,MADD+I-1,K) - TCALC_E(IPOINT,MADD+I-1,K-1))
                 END DO
              END DO
           ENDIF
           MADD = MADD + NM(IBAND)
        END DO
     ELSE
        DO INDXX = 1, NSCANS
           IF (INDXX >= JMIN) THEN
              DO I = 1, NMON
                 IF( ABS( TCALC(IPOINT,MADD+I-1)) .GT. 664.0 ) THEN
                    ! LIMIT SO ONLY GET EXPONENT < 300
                    TCALC(IPOINT,MADD+I-1) = 664.0D0
                 ENDIF
                 TCALC(IPOINT,MADD+I-1) = EXP((-TCALC(IPOINT,MADD+I-1)))
              END DO
           ENDIF
           MADD = MADD + NM(IBAND)
        END DO
     END IF
     RETURN
   END SUBROUTINE NTRAN


!--------------------------------------------------------------------------------
   SUBROUTINE GASNTRAN( IR, IBAND, JMIN, IPOINT, MONONE, MXONE)

!     COMPUTE MONOCHROMATIC TRANSMITTANCES AT EACH WAVELENGTH
!     FOR BANDPASS IBAND FOR ONE GAS.  TRANMSITTANCES ARE COMPUTED FOR SPECTRA
!     BETWEEN SPECTRUM ISCAN(IBAND,JMIN) AND ISCAN(IBAND,NSCAN(IBAND)).
!
!          IR= GAS NUMBER
!          IBAND=BAND PASS NUMBER
!          JSCAN=ISCAN(IBAND,J)=SPECTRUM NUMBER FOR CALCULATIONS
!          IPOINT,MONONE-INDICES IN TCALC ARRAY FOR CALCULATED TRANSMITTANCES

      INTEGER :: IR
      INTEGER :: IBAND
      INTEGER :: JMIN
      INTEGER :: IPOINT
      INTEGER :: MONONE
      INTEGER :: MXONE

      INTEGER :: NMON, NSCANS, KSMAX2, K, JSCAN,ICINDX, MSTOR, MXMAX, J, MADD, I, ALT
      REAL(DOUBLE) :: FACMAS, XFAC, WAVE_NR

! MP update gas spectrum for emission eg wave_nr...in progress


!  --- NMON=NUMBER OF MONOCHROMATIC POINTS FOR THE BANDPASS CALCULATION
      NMON   = NM(IBAND)
      NSCANS = NSCAN(IBAND)

!  --- ZERO APPROPRIATE TRANSMISSION ARRAY ELEMENTS FOR CROSS SECTION
!  ---  CALCULATIONS
      MADD = MONONE
      TCALC(IPOINT,MADD:NMON-1+MADD) = 0.D0
!  --- MAXIMUM LAYER FOR SUMMING CROSS SECTIONS CALCULATIONS
      KSMAX2 = KZTAN(ISCAN(IBAND,NSCANS))
!                   ------------LOOP OVER LAYERS
      DO K = 1, KSMAX2
         MADD = MONONE
               JSCAN = ISCAN(IBAND,JMIN)
               IF (K <= KZTAN(JSCAN)) THEN
                  FACMAS = CCC(JSCAN,K)/PMASMX(K)
!                   ------------LOOP OVER FREQUENCIES
                  DO J = 1, NMON
                        XFAC = X(IR,K)/XORG(IR,K)
                        ICINDX = MXONE + J - 1
                        MSTOR = MADD + J - 1
                        WAVE_NR = WSTART(IBAND) + (J-1)*DN(IBAND)
                        IF (IR/=1 .AND. IFDIFF) THEN
!  --- APPLY DIFFERENTIAL WAVENUMBER SHIFT
                           ICINDX = ICINDX + ISHIFT(IR-1)
                           MXMAX = MXONE + NMON - 1
!  --- FIXUP AT ENDPOINTS
                           ICINDX = MAX0(MXONE,ICINDX)
                           ICINDX = MIN0(MXMAX,ICINDX)
                        ENDIF
                        TCALC(IPOINT,MSTOR) = TCALC(IPOINT,MSTOR) + XFAC*CROSS(&
                           IR,K,ICINDX)*FACMAS
                        IF (IEMISSION/=0) THEN
                           ! Transmission calculated below the layer ALT, needed
                           ! for calculation of contribution to emission from
                           ! layer ALT to the ground
                           DO ALT=1,KSMAX2
                              IF (ZBAR(ALT) > ZBAR(K)) THEN
                                 TCALC_E(IPOINT,MSTOR,ALT) = &
                                      TCALC_E(IPOINT,MSTOR,ALT) + (X(IR,K)/XORG(IR,K)) * CROSS_FACMAS(IR,K,MSTOR)
                              END IF
                           END DO
                        END IF
                     END DO
                  ENDIF
               END DO
!  --- COMPUTE MONOCHROMATIC TRANSMITTANCES FROM CROSS SECTION SUMS
      MADD = MONONE
            DO I = 1, NMON
               if( ABS( TCALC(IPOINT,MADD+I-1)) .GT. 664.0 ) THEN
                  ! LIMIT SO ONLY GET EXPONENT < 300
                  TCALC(IPOINT,MADD+I-1) = 664.0d0
               ENDIF
               if (IEMISSION.EQ.1) then
                  ! Background
                  TCALC(IPOINT, MADD+I-1) &
                       = PLANCK(WAVE_NR,EMISSION_T_BACK) &
                       * EXP((-TCALC(IPOINT,MADD+I-1)))
                  TCALC_E(IPOINT,MADD+I-1,KMAX+1) = TCALC(IPOINT, MADD+I-1)
                  TCALC_S(IPOINT, MADD+I-1, 1)=0.D0
                  DO K=2,KSMAX2
                     ! Calculates the spectrum, Planck is the emission of a black body at the temperature
                     ! of T(K) in layer K,
                     ! TCALC_E(IPOINT,MADD+I-1,K) - TCALC_E(IPOINT,MADD+I-1,K-1) is the absorption
                     ! of the layer K alone
                     TCALC_S(IPOINT, MADD+I-1, K) = PLANCK(WAVE_NR,T(K))&
                          *(TCALC_E(IPOINT,MADD+I-1,K) - TCALC_E(IPOINT,MADD+I-1,K-1))
                     ! TCALC contains the spectrum
                     TCALC(IPOINT, MADD+I-1) = TCALC(IPOINT, MADD+I-1) &
                          + TCALC_S(IPOINT, MADD+I-1, K)
                  end DO
               end if
               TCALC(IPOINT,MADD+I-1) = EXP((-TCALC(IPOINT,MADD+I-1)))
            endDO
            RETURN

          END SUBROUTINE GASNTRAN

!---------------------------------------------------------------------------
      SUBROUTINE ZERONTRAN(IBAND, IPOINT, MONONE)

!     COMPUTE MONOCHROMATIC TRANSMITTANCES AT EACH WAVELENGTH
!     FOR BANDPASS IBAND SET TO 1.0.  TRANMSITTANCES ARE COMPUTED FOR SPECTRA
!     BETWEEN SPECTRUM ISCAN(IBAND,JMIN) AND ISCAN(IBAND,NSCAN(IBAND)).
!
!          IBAND=BAND PASS NUMBER
!          JSCAN=ISCAN(IBAND,J)=SPECTRUM NUMBER FOR CALCULATIONS
!          IPOINT,MONONE-INDICES IN TCALC ARRAY FOR CALCULATED TRANSMITTANCES

      INTEGER, INTENT(IN) :: IBAND, IPOINT, MONONE

      INTEGER :: NMON, MADD, I
      !REAL(DOUBLE) :: FACMAS

!  --- NMON=NUMBER OF MONOCHROMATIC POINTS FOR THE BANDPASS CALCULATION
      NMON = NM(IBAND)
      !NSCANS = NSCAN(IBAND)
!  --- ZERO APPROPRIATE TRANSMISSION ARRAY ELEMENTS FOR CROSS SECTION
!  ---  CALCULATIONS
      MADD = MONONE
      TCALC(IPOINT,MADD:NMON-1+MADD) = 0.D0
!  --- COMPUTE MONOCHROMATIC TRANSMITTANCES FROM CROSS SECTION SUMS
      MADD = MONONE
            DO I = 1, NMON
               if( ABS( TCALC(IPOINT,MADD+I-1)) .GT. 664.0 ) THEN
                  ! LIMIT SO ONLY GET EXPONENT < 300
                  TCALC(IPOINT,MADD+I-1) = 664.0d0
               ENDIF
               TCALC(IPOINT,MADD+I-1) = EXP((-TCALC(IPOINT,MADD+I-1)))
            END DO

      RETURN

      END SUBROUTINE ZERONTRAN


!----------------------------------------------------------------------
      REAL (DOUBLE) FUNCTION PLANCK (F, T)
! --- CALCULATES THE PLANCK FUNCTION [W(M**2 SR CM**-1)] FOR &
! --- A GIVEN FREQUENCY F [CM^-1] AND TEMPERATURE T [K]
! --- MATHIAS PALM 2007

        REAL(DOUBLE) :: F, T
! --- CALCULATE CONSTANTS FOR PLANCK FUNCTION TO SPEED UP
        PLANCK = PLANCK_C1 * F**3 / (EXP(PLANCK_C2 * F / T ) - 1.0D0)

      END FUNCTION PLANCK

      END MODULE TRANSMIS
