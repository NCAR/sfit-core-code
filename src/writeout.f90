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

      MODULE WRITEOUT

      USE PARAMS
      USE DATAFILES
      USE RETVPARAM
      USE MOLCPARAM
      USE ISOTOPE
      USE BANDPARAM

      IMPLICIT NONE

      LOGICAL, PARAMETER :: F_WRTDETAIL = .TRUE.

      LOGICAL :: F_WRTBIN            = .FALSE.
      LOGICAL :: F_WRTSMTH           = .FALSE.
      LOGICAL :: F_WRTK              = .FALSE.
      LOGICAL :: F_WRTKOUT           = .FALSE.
      LOGICAL :: F_WRTAK             = .FALSE.
      LOGICAL :: F_WRTKB             = .FALSE.
      LOGICAL :: F_WRTAB             = .FALSE.
      LOGICAL :: F_WRTG              = .FALSE.
      LOGICAL :: F_WRTSTV            = .FALSE.
      LOGICAL :: F_WRTPARM           = .FALSE.
      LOGICAL :: F_WRTRPRF           = .FALSE.
      LOGICAL :: F_WRTAPRF           = .FALSE.
      LOGICAL :: F_WRTGASSPC         = .FALSE.
      LOGICAL :: F_WRTSA             = .FALSE.
      LOGICAL :: F_WRTSHAT           = .FALSE.
      LOGICAL :: F_WRTSAINV          = .FALSE.
      LOGICAL :: F_WRTSEINV          = .FALSE.
      LOGICAL :: F_WRTSMEAS          = .FALSE.
      LOGICAL :: F_WRTSUMRY          = .FALSE.
      LOGICAL :: F_WRTPBP            = .FALSE.
!      LOGICAL :: F_WRTPBP_KB         = .FALSE.
      LOGICAL :: F_WRTCHANNEL        = .FALSE.
      LOGICAL :: F_WRTRAYTC          = .FALSE.
      LOGICAL :: F_WRTSOLSPEC        = .FALSE.
      LOGICAL :: F_WRTLM             = .FALSE.
      LOGICAL :: XSC_DETAIL          = .FALSE.

      INTEGER :: OUTPUTLEVL    = 0
      INTEGER :: GASOUTTYPE    = 1     ! 1: FINAL ITERATION ONLY
                                       ! 2: PLUS EACH ITERATION
                                       ! 3: PLUS EACH LAYER
      INTEGER :: RAYOUTTYPE    = 1     ! 1: SELECTION OF SA'S ONLY
                                       ! 2: PLUS DETAILED RAYTRACING
                                       ! 3: PLUS OLD MIX, MS, PT FILE

      CONTAINS

     SUBROUTINE INIT_WRITEOUT()

     ! DEFINE THE VERBOSITY OF THE OUTPUT
     ! 0 - ONLY OPERATIONAL (DEFAULT)
     ! 1 - SOME MORE OUTPUT FOR EXPERIMENTS
     ! 2 - MORE DIAGNOSTICS
     ! 3 - EVERYTHING ELSE IS PUT OUT, MAINLY FOR DEBUGGING

     IMPLICIT NONE

     IF (OUTPUTLEVL.GE.0) THEN
        F_WRTBIN            = .TRUE.
     END IF
     IF (OUTPUTLEVL.GE.1) THEN
        F_WRTSTV            = .TRUE.
        F_WRTRPRF           = .TRUE.
        F_WRTAPRF           = .TRUE.
        F_WRTPBP            = .TRUE.
        F_WRTSUMRY          = .TRUE.
        F_WRTK              = .TRUE.
        F_WRTAK             = .TRUE.
        F_WRTSA             = .TRUE.
     END IF
     IF (OUTPUTLEVL.GE.2) THEN
        F_WRTAB             = .TRUE.
        F_WRTKB             = .TRUE.
        F_WRTSMEAS          = .TRUE.
        F_WRTSAINV          = .TRUE.
        F_WRTSEINV          = .TRUE.
        F_WRTSHAT           = .TRUE.
        F_WRTSMTH           = .TRUE.
        F_WRTPARM           = .TRUE.
        F_WRTGASSPC         = .TRUE.
     END IF
     IF (OUTPUTLEVL.EQ.3) THEN
        F_WRTCHANNEL        = .TRUE.
        F_WRTRAYTC          = .TRUE.
        F_WRTSOLSPEC        = .TRUE.
        F_WRTLM             = .TRUE.
!        F_WRTPBP_KB         = .TRUE.
        XSC_DETAIL          = .TRUE.
     END IF
     IF (OUTPUTLEVL.GT.3) THEN
        WRITE(16,*) 'WRITEOUT:INIT_WRITEOUT: OUTPUT LEVEL CAN ONLY BE 0, 1, 2 OR 3 : ', OUTPUTLEVL
        WRITE( 0,*) 'WRITEOUT:INIT_WRITEOUT: OUTPUT LEVEL CAN ONLY BE 0, 1, 2 OR 3 : ', OUTPUTLEVL
        CALL SHUTDOWN
        STOP '2'
     END IF

     END SUBROUTINE INIT_WRITEOUT


!----------------------------------------------------------------------
      SUBROUTINE MIXOUT(KZMAXLAY, KVERT)

!  --- CREATES FILE OF REVISED MIXING RATIOS - REVISED FEB 19, 1992

      IMPLICIT NONE

      INTEGER, INTENT(IN)             :: KZMAXLAY, KVERT
      INTEGER                         :: I, K, J
      REAL(DOUBLE), DIMENSION(LAYMAX) :: GASMIX
      CHARACTER (LEN=7)               :: GNAME

      CALL FILEOPEN( 17, 1 )
      CALL FILEOPEN( 22, 1 )

      DO I = 1, MOLTOTAL

         GNAME = NAME(I)

!  --- SEE IF THIS GAS ID HAS BEEN USED FOR A ISOTOPE - IF SO CHANGE THE NAME
         DO K = 1, NISOSEP
            IF( I .EQ. NEWID(K) ) GNAME = NEWNAME(K)
         ENDDO

!  --- SEE IF THIS GAS IS A RETRIEVAL GAS
         DO J = 1, NRET
            IF (GNAME == NAME(IGAS(J))) GOTO 11
         END DO

!  --- NOT A RETRIEVAL GAS
         GASMIX(:KZMAXLAY) = XGAS(J,:KZMAXLAY)
         GO TO 13

!  --- RETRIEVAL GAS, UPDATE MIXING RATIOS
   11    CONTINUE
         GASMIX(:KZMAXLAY) = X(J,:KZMAXLAY)

!  --- OUTPUT OF MIXING RATIOS & PARTIAL COLUMNS
   13    CONTINUE

         WRITE (17, 501) GNAME
         WRITE (17, 883) (GASMIX(K),K=1,KMAX)
         WRITE (22, 501) GNAME
         WRITE (22, 883) (GASMIX(K)*CCC(KVERT,K),K=1,KMAX)

      END DO

      CALL FILECLOSE( 17, 1 )
      CALL FILECLOSE( 22, 1 )

      RETURN

!  111 FORMAT(/,' OUTPUT FILE OPEN ERROR IN MIXOUT-UNIT 17'/,' Filename: "',A,'"')
!  222 FORMAT(/,' INPUT FILE OPEN ERROR IN MIXOUT-UNIT 12'/,' Filename: "',A,'"')
  501 FORMAT(A7)
!  503 FORMAT(8E10.4)
!     --- Use next format statement to write fixed format like
!         from fastcode
! 883 FORMAT(8(E10.4))
! 883 FORMAT(8(E10.4,1x))
  883 FORMAT(6(ES12.4))

      RETURN

      END SUBROUTINE MIXOUT


!  --- WRITE OUT A SUMMARY OF RETRIEVAL PARAMETERS
      SUBROUTINE WRTSMRY( DOF, ITER, CHI_2_Y, FOVDIA, RMS, NLEV, VOSUM, VERSUM )

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)                           :: DOF(3), CHI_2_Y, RMS
      REAL(DOUBLE), DIMENSION(MAXBND), INTENT(IN)        :: FOVDIA
      REAL(DOUBLE), DIMENSION(MOLMAX,LAYMAX), INTENT(IN) :: VERSUM, VOSUM
      INTEGER, INTENT(IN)                                :: ITER, NLEV
      INTEGER I, J

      CALL FILEOPEN( 20, 1 )

      WRITE(20,'(A,A)') TRIM(TAG), ' RETRIEVAL SUMMARY '

      WRITE(20,101) NFITS
      DO I=1, NFITS
         WRITE(20,*) STITLE(I)
      ENDDO

      WRITE(20,101) NRET
      WRITE(20,'(A)') ' IRET   GAS_NAME  IFPRF       APR_COLUMN    RET_COLUMN'
      DO I=1,NRET
         WRITE(20,102) I, ADJUSTR(NAME(IGAS(I))), IFPRF(I), VOSUM(I,NLEV), VERSUM(I,NLEV)
      END DO

      WRITE(20,101) NBAND
      WRITE(20,'(A,A)') 'IBAND       NUSTART        NUSTOP         SPACE     NPTSB     PMAX    FOVDIA     ', &
                       'MEAN_SNR  NSCAN  JSCAN     INIT_SNR     CALC_SNR'
      DO I=1,NBAND
         WRITE(20,103) I, WAVE3(I), WAVE4(I), SPAC(I), NPRIM(I), PMAX(I), FOVDIA(I), &
                           SUM(SNR_CLC(I,1:NSCAN(I)))/DBLE(NSCAN(I)), NSCAN(I)
         DO J=1,NSCAN(I)
            WRITE(20,104) J, SCNSNR(I,J), SNR_CLC(I,J)
         ENDDO
      ENDDO

      WRITE(20,100) '       FITRMS       CHI_2_Y      DOFS_ALL      DOFS_TRG      DOFS_TPR      ITER  MAX_ITER CONVERGED   DIVWARN'
      WRITE(20,105) RMS, CHI_2_Y, DOF(1), DOF(2), DOF(3), ITER, ITRMAX, CONVERGE, DIVWARN

      CALL FILECLOSE( 20, 1 )

      RETURN

 100 FORMAT( /, A )
 101 FORMAT( /, I10 )
 102 FORMAT( I5, 4X, A7, L7, 3X, 2ES14.5 )
 103 FORMAT( I5, 2F14.5, 2X,F12.9, I10, F9.2, F10.6, F13.6, I7, 3F13.6 )
 104 FORMAT( I103, 3F13.6 )
 105 FORMAT( F13.6, 1x, F13.6, 3F14.3, 2I10, 2L10 )

       END SUBROUTINE WRTSMRY


!  --- WRITE OUT TABLE OF PROFILES APRIORI ATMOSPHERE & VMR
      SUBROUTINE WRTAPRF( NRET, NLEV, KVERT )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NRET, NLEV, KVERT
      INTEGER             :: K, KK

      CALL FILEOPEN( 87, 2 )

      WRITE(87,*) TRIM(TAG), ' APRIORI Z, P, T, AIRMASS & PROFILES'
      WRITE(87,*) NMOL, NLEV, NRET, GAS(:NRET)
      WRITE(87,408)0,0,0,0,0,(K,K=1,NMOL)
      WRITE(87,409) (ADJUSTR(TRIM(TRIM(NAME(K)))), K=1,NMOL)
      DO KK = 1, NLEV
         WRITE(87,407) Z(KK), ZBAR(KK), TORG(KK), PMBORG(KK), CORG(KVERT, KK), (FXORG(K,KK), K=1,NMOL)
      END DO

      CALL FILECLOSE( 88, 1 )

      RETURN

  409 FORMAT('      Z   ZBAR  TEMPERATURE       PRESSURE        AIRMASS', 99(A15))
  408 FORMAT(2(I8),I13,110(I15))
  407 FORMAT(2(F8.3),F13.3,255(ES15.4))

      END SUBROUTINE WRTAPRF


!  --- WRITE OUT TABLE OF PROFILES RETRIEVED ATMOSPHERE & VMR
      SUBROUTINE WRTRPRF( NRET, NLEV, KVERT )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NRET, NLEV, KVERT
      INTEGER             :: K, KK

      DO KK = 1, NRET
         DO K = 1, NLEV
            IF( X(KK,K) .LT. 0.0D0 )THEN
               PRINT*, 'NEGATIVE MIXING RATIO VALUES FOUND FOR : ', GAS(KK)
               WRITE(16,*) 'NEGATIVE MIXING RATIO VALUES FOUND FOR : ', GAS(KK)
               EXIT
            ENDIF
         ENDDO
      ENDDO

      DO K=1, NMOL
         DO KK = 1, NRET
            IF( TRIM(NAME(K)) .EQ. TRIM(GAS(KK)) )FXORG(K,:NLEV) = X(KK,:NLEV)
         ENDDO
      ENDDO
      CALL FILEOPEN( 88, 2 )

      WRITE(88,*) TRIM(TAG), ' RETRIEVED Z, P, T, AIRMASS & PROFILES'
      WRITE(88,*) NMOL, NLEV, NRET, GAS(:NRET)
      WRITE(88,408)0,0,0,0,0,(K,K=1,NMOL)
      WRITE(88,409) (ADJUSTR(TRIM(TRIM(NAME(K)))), K=1,NMOL)
      DO KK = 1, NLEV
         WRITE(88,407) Z(KK), ZBAR(KK), T(KK), PMBORG(KK), CCC(KVERT, KK), (FXORG(K,KK), K=1,NMOL)
      END DO

      CALL FILECLOSE( 88, 1 )

      DO K=1, NMOL
         DO KK = 1, NRET
            IF( TRIM(NAME(K)) .EQ. TRIM(NAME(IGAS(KK))) )FXORG(K,:NLEV) = XORG(KK,:NLEV)
         ENDDO
      ENDDO

      RETURN

  409 FORMAT('      Z   ZBAR  TEMPERATURE       PRESSURE        AIRMASS', 99(A15))
  408 FORMAT(2(I8),I13,110(I15))
  407 FORMAT(2(F8.3),F13.3,255(ES15.4))

      END SUBROUTINE WRTRPRF


!  --- WRITE OUT OBSERVED, CALCULATED, AND DIFFERENCES - pbpfile
      SUBROUTINE WRTPBP( TOBS, YHAT )

      IMPLICIT NONE

      REAL(DOUBLE), DIMENSION(MMAX), INTENT(IN) :: TOBS, YHAT
      INTEGER                      :: IA, IB, N, J, K, IBAND, JSCAN, INDXX, IOUT
      REAL(DOUBLE)                 :: WAVE = 0.0D0
      REAL(DOUBLE), DIMENSION(12)  :: FX = 0.0D0

      CALL FILEOPEN( 8, 2 )
      WRITE(8,*) TRIM(TAG), ' OBSERVED, FITTED AND DIFFERENCE SPECTRA'
      WRITE(8,*) NFITS, NBAND
      K = 0
      DO IBAND = 1, NBAND
         N = NSCAN(IBAND)
         IF (N == 0) CYCLE
         DO J = 1, N
            K = K + 1
            JSCAN = ISCAN(IBAND,J)
            WRITE(8,'(A80)') STITLE(K)
            WRITE(8, 37905) ISPEC(JSCAN), SPAC(IBAND), NPRIM(IBAND), WSTART(IBAND), WSTOP(IBAND), &
                            Z(KZTAN(JSCAN)), IBAND, J, NRETB(IBAND)
            IB = 0
            INDXX = ISCNDX(1,IBAND,J) -1
 3790       CONTINUE
            IA = IB + 1
            IB = IA + 11
            IB = MIN0(NPRIM(IBAND),IB)
            FX(:IB-IA+1) = TOBS(IA+INDXX:IB+INDXX) - YHAT(IA+INDXX:IB+INDXX)
            WAVE = WSTART(IBAND) + (IA - 1)*SPAC(IBAND)
            WRITE (8, 3794) WAVE, (TOBS(IOUT+INDXX),IOUT=IA,IB)
            WRITE (8, 3795)       (YHAT(IOUT+INDXX),IOUT=IA,IB)
            WRITE (8, 3795)       (FX  (IOUT-IA+1), IOUT=IA,IB)
            IF (IB < NPRIM(IBAND)) GO TO 3790
         END DO
      END DO
      CALL FILECLOSE( 8, 1 )

       RETURN

37905 FORMAT(I12, 1P, E25.15, I12, 0P, F25.15, 2F25.15, 3i5)
 3794 FORMAT(F14.6,12ES15.5E4)
 3795 FORMAT(14X,12ES15.5E4)

      END SUBROUTINE WRTPBP

!  --- WRITE OUT STATE VECTOR
      SUBROUTINE WRTSTV( NLEV, ITER, ISMIX, VERSUM, VOSUM, PNAME, XHAT, XAPR )

      IMPLICIT NONE

      CHARACTER (LEN=14), DIMENSION(NMAX), INTENT(IN)    :: PNAME
      REAL(DOUBLE),DIMENSION(NMAX), INTENT(IN)           :: XHAT, XAPR
      REAL(DOUBLE), DIMENSION(MOLMAX,LAYMAX), INTENT(IN) :: VERSUM, VOSUM
      INTEGER, INTENT(IN)      :: NLEV, ITER, ISMIX
      INTEGER                  :: I, J

!  --- FILE OF FINAL MIXING RATIOS AND COLUMNS TO TAPE18 - statevec
      CALL FILEOPEN( 18, 2 )
      WRITE (18, *) TRIM(TAG), ' STATE VECTOR'
      WRITE (18, 504) NLEV, ITER, ITRMAX, IFTEMP, CONVERGE, DIVWARN
      WRITE (18, *) 'z (Mid-points, Find layer boundaries in station.layers)'
      !      WRITE(18,*)'z'
      WRITE (18, 508) (ZBAR(I),I=1,NLEV)
      WRITE (18, *) 'p'
      WRITE (18, 506) (PMB(I),I=1,NLEV)
      WRITE (18, *) 't'
      WRITE (18, 506) (TORG(I),I=1,NLEV)
      IF( IFTEMP )THEN
         WRITE (18, *) 'tr'
         WRITE (18, 506) (T(I),I=1,NLEV)
      ENDIF
      WRITE (18, *)
      WRITE (18, *) NRET
      DO I = 1, NRET
         WRITE (18, 507) 'A Priori', NAME(IGAS(I))
         WRITE (18, 506) VOSUM(I,NLEV)
         WRITE (18, 506) (XORG(I,J),J=1,NLEV)
         WRITE (18, 507) 'Retrieved', NAME(IGAS(I))
         WRITE (18, 506) VERSUM(I,NLEV)
         WRITE (18, 506) (X(I,J),J=1,NLEV)
      END DO
      WRITE (18, *)
      WRITE (18, *) ISMIX
      WRITE (18, 502) (PNAME(I),I=1,ISMIX)
      WRITE (18, 506) (XAPR(I),I=1,ISMIX)
      WRITE (18, 506) (XHAT(I),I=1,ISMIX)
      CALL FILECLOSE( 18, 1 )

      RETURN

  502 FORMAT(5(1X,A14))
  504 FORMAT(3I5,3L5)
  506 FORMAT(5(1X,1P,E14.4))
  507 FORMAT(A20)
  508 FORMAT(5(1X,F14.4))

      END SUBROUTINE WRTSTV

      END MODULE WRITEOUT
