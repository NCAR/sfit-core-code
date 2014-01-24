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

    MODULE SFIT4

      USE PARAMS
      USE VIBFCN
      USE RETVPARAM
      USE FRWDMDL
      USE BANDPARAM
      USE MOLCPARAM
      USE DATAFILES
      USE SOLAR
      USE OPT
      USE DIAGNOSTIC
      USE SYNSPEC
      USE LINEPARAM
      USE INITIALIZE
      USE READIN
      USE WRITEOUT
      USE BINPUT_4_0
      USE RAYTRACE

      IMPLICIT NONE

      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: SED, XAPR, XHAT, YHAT, Y

      CONTAINS

      SUBROUTINE SFIT( )

! SFIT4 VERSION 003.82 2011
! SEE FILE COMMENTS.SAVED FOR OLDER COMMENTS AND HISTORY


      IMPLICIT NONE

      INTEGER :: NLEV = 0, NAERR = 0 !, SVAR=1
      INTEGER :: KZMAXLAY, INDXX, KVERT, ITER=0, NEGFLAG
      INTEGER :: I, IBAND, IND, J, K, KK, N
      INTEGER, DIMENSION(MAXSPE) :: NLAY

      REAL(DOUBLE) :: YMAX = 0.0D0
      REAL(DOUBLE) :: SERR = 0.0D0
      REAL(DOUBLE) :: SIGB = 0.0D0
      REAL(DOUBLE) :: SIGC = 0.0D0
!      REAL(DOUBLE) :: SIGZ = 0.0D0
      REAL(DOUBLE) :: AIRCOL = 0.0D0
      REAL(DOUBLE) :: SIGMA


!      REAL(DOUBLE), DIMENSION (MOLMAX) :: DELTA_XM

      REAL(DOUBLE), DIMENSION(12)            :: FX = 0.0D0
      REAL(DOUBLE), DIMENSION(MOLMAX,LAYMAX) :: VERSUM = 0.0D0, VOSUM = 0.0D0


! ------------------------------------------------------------------------------

! --- SET FLAGS FOR OUTPUT
      CALL INIT_WRITEOUT()

! --- READ IN RETRIEVAL LAYERING FROM STATION.LAYERS FILE
      CALL READLAYRS( NLEV )
      KMAX = NLEV

! --- PRINT INVERSION PARAMETERS NOW THAT WE HAVE NLEV
!     ISOTOPE
      CALL READCK1( NLEV, NEGFLAG )

! --- PRINT OUT FM, RT PARAMETERS
      CALL READCK2( CPNAM )

! --- CHECK WE ONLY FIT PHASE AND MODULATION FUNCTIONS TYPE = 2
      IF( IEAP .NE.2 .AND. F_RTAPOD )  GOTO 667
      IF( IEPHS.NE.2 .AND. F_RTPHASE ) GOTO 668

! --- CHECK THAT BACKGROUND AND SHIFT PARMS ARE IN RANGE
      IF( NBACK .LT.1 .OR. NBACK.GT.3 )  GOTO 664
      IF( ISPARM.LT.0 .OR. ISPARM.GT.3 ) GOTO 665

! --- PRINT OUT BAND PARAMETERS
      CALL READCK3( )

! --- READ IN SPECTRA
      CALL SETUP1( )

! --- ALLOCATE MASS PATH VECTORS BEFORE RAYTRACE CALL
      ALLOCATE( CCC(NSPEC+1,NLEV), CORG(NSPEC+1,NLEV) )

! --- PERFORM RAYTRACE FOR ALL SZA'S (NSPEC)
      IF( RAYOUTTYPE .GE. 2 )THEN
         CALL FILEOPEN( 73, 2 )
         WRITE(73,*) TRIM(TAG), ' RAYTRACE DETAIL FILE'
      ENDIF
      CALL LBLATM( 0, NLEV )

      IF( RAYTONLY )THEN
         IF( RAYOUTTYPE .GE. 2 )CALL FILECLOSE( 73, 1 )
         WRITE(16,*) ': COMPUTING RAYTRACE ONLY.'
         WRITE(00,*) ': COMPUTING RAYTRACE ONLY.'
         CALL SHUTDOWN
         STOP '3'
      ENDIF
      IF( NLEV .LT. KMAX ) KMAX = NLEV

! --- ZENITH AIRMASS VECTOR IN CCC
      KVERT = NSPEC +1

!  --- FIND MAXIMUM MASS PATH FOR EACH LAYER
      CALL MASSPATH( -1 )

! --- MORE SETUP GET HITRAN, CHECK MODULATION EFF & PHASE, SOLAR SPECTRUM...
      CALL SETUP2

! --- FINISH SETUP CALCULATE CROSS SECTIONS, -1 means, crosssections for all levels.
      CALL SETUP3( XSC_DETAIL, -1 )

! --- IF EMISSION, SNR IS THE NOISE ON THE MEASURMENT GIVEN IN (MW/(CM^2*SR*CM-1)) MP
! MUSTFIX
! COMMENTED OUT JWH JUN 2012
! check initialize:filse
!      DO I=1, NBAND
!         IF (IEMISSION/=1.AND.(.NOT.IENORM(i)/=1)) SNR = 1/SNR
!      ENDDO

!      DO K = 1, NSPEC +1
!         print*, astang(k), appang(k), ispec(k), rearth(k), reflat(k), xvb(k)
!      enddo

!      do j=1, nband
!      do k=1, NSCAN(j)
!         print*, j, k, ISCAN(j,k), appang(ISCAN(j,k))
!      enddo
 !     enddo

! --- THIS LOOP DOES NOT DO MUCH BUT REMINDS OF WHERE TO SETUP VARIABLE LAYER SCHEMES
! --- KZTAN STILL USED
      DO K = 1, NSPEC
         NLAY(K)  = KMAX
         KZTAN(K) = NLAY(K)   ! # LAYERS IN THIS MASSPATH
         KZMAXLAY = KZTAN(K)
      END DO

! ---- INITIALIZE STATE VECTOR
      CALL INIT_PARM()

! --- SETUP ARRAYS FOR OPTIMAL ESTIMATION SUBROUTINE
! --- FILL MATRIX SED COVARIANCE OF MEASURED SPECTRUM - OFF DIAGONAL ELEMENTS ASSUMED TO BE ZERO

! --- ALLOCATE SE
      ALLOCATE( SE(NFIT), SED(NFIT), STAT=NAERR )
      IF (NAERR /= 0) THEN
         WRITE(16, *) 'COULD NOT ALLOCATE SE ARRAY'
         WRITE(16, *) 'ERROR NUMBER = ', NAERR
         WRITE(00, *) 'COULD NOT ALLOCATE SE ARRAY'
         WRITE(00, *) 'ERROR NUMBER = ', NAERR
         CALL SHUTDOWN
         STOP '4'
      ENDIF
      SE(:)    = 0.0D0
      SED(:)   = 0.0D0

      WRITE (*, *) ' INITIALIZING VARIANCE VECTOR...'
      CALL FILSE (SED, NFIT)
      SE = SED(:NFIT)

! --- ALLOCATE COVARIANCE ARRAYS
      ALLOCATE( SA(NVAR,NVAR), SAINV(NVAR,NVAR), SHAT(NVAR,NVAR), STAT=NAERR )
      IF (NAERR /= 0) THEN
         WRITE(16, *) 'COULD NOT ALLOCATE SA ARRAY'
         WRITE(16, *) 'ERROR NUMBER = ', NAERR
         WRITE(00, *) 'COULD NOT ALLOCATE SA ARRAY'
         WRITE(00, *) 'ERROR NUMBER = ', NAERR
         CALL SHUTDOWN
         STOP '4'
      ENDIF
      SA(:,:)    = 0.0D0
      SHAT(:,:)  = 0.0D0
      SAINV(:,:) = 0.0D0

      CALL FILSA( SA )

! --- ALLOCATE WORKING AND RESULT ARRAYS
      ALLOCATE( XHAT(NVAR), YHAT(NFIT), XAPR(NVAR), Y(NFIT), STAT=NAERR )
      IF (NAERR /= 0) THEN
         WRITE(16, *) 'COULD NOT ALLOCATE XHAT ARRAY'
         WRITE(16, *) 'ERROR NUMBER = ', NAERR
         WRITE(00, *) 'COULD NOT ALLOCATE XHAT ARRAY'
         WRITE(00, *) 'ERROR NUMBER = ', NAERR
         CALL SHUTDOWN
         STOP '4'
      ENDIF
      XHAT(:) = 0.0D0
      YHAT(:) = 0.0D0
      XAPR(:) = 0.0D0
      Y(:)    = 0.0D0

! --- STORE OBSERVED SPECTRA AS Y FOR OE
      Y(:NFIT)    = TOBS(:NFIT)
      XAPR(:NVAR) = PARM(:NVAR)

      IF( F_WRTPARM )THEN
         CALL FILEOPEN( 89, 1 )
         WRITE(89,*) TRIM(TAG), ' STATE VECTOR FACTORS BY ITERATION N VECTOR'
         WRITE(89,*) NVAR
         WRITE(89,263) (I,I=1,NVAR)
         WRITE(89,262) ADJUSTR(PNAME(:NVAR))
      ENDIF

!  --- CALL OPTIMAL ESTIMATION SUBROUTINE
      WRITE (16, 3695) NFIT, NVAR
      WRITE (16, 3696)

!  --- SAVE NVAR KEY TO DETAIL
      WRITE(16,502) (PNAME(I),I=1,NVAR)
      WRITE(16,*)''

      WRITE (*, *) ' BEGIN ITERATIVE RETRIEVAL LOOP...'
      CALL OPT_3(Y, PARM, XHAT, YHAT, NFIT, NVAR, CONVERGE, ITRMAX, TOL, RETFLG, DIVWARN, ITER, ISMIX, NLEV )

      CALL FILECLOSE( 91, 1 )

!  --- CLOSE DETAILED STATE VECTOR
      IF( F_WRTPARM )CALL FILECLOSE( 89, 1 )

!  --- IF FM HAD ABNORMAL TERMINATION, THEN JUST LEAVE
      IF( Tflg )GOTO 700

      IF (.NOT.RETFLG) THEN
         CALL FILEOPEN( 19, 2 )
         WRITE(19,*) TRIM(TAG), ' SYNTHETIC SPECTRUM FROM FM - NO FIT'
         IND = 0
         DO I = 1, NBAND
            DO K = 1, NSPEC
               WRITE (19, *) 'Simulation from opt_fm'
               WRITE (19, *) WSTART(I), WSTOP(I), SPAC(I)
               WRITE (19, *) NPRIM(I)
               DO J = 1, NPRIM(I)
                  WRITE (19, *) YHAT(IND+J)
               END DO
               IND = IND + NPRIM(I)
            END DO
         END DO
         CALL FILECLOSE( 19, 1 )
      ENDIF

! --- RENORMALIZATION OF FINAL SPECTRA
      TOBS_ORIG(:NFIT) = TOBS(:NFIT)
      IND = 0
      DO I = 1, NBAND
         ! Normalized if emission is switched off or the emission spectra are normalized
         IF (IEMISSION.EQ.0.OR.IENORM(I).EQ.1) THEN
            N = NSCAN(I)
            DO K = 1, N
               YMAX = 0.D0
               DO J = 1, NPRIM(I)
                  YMAX = MAX(YMAX,YHAT(IND+J))
               END DO
               YHAT(IND+1:NPRIM(I)+IND) = YHAT(IND+1:NPRIM(I)+IND)/YMAX
               TOBS(IND+1:NPRIM(I)+IND) = TOBS(IND+1:NPRIM(I)+IND)/YMAX
               IND = IND + NPRIM(I)
            END DO
         ENDIF
      END DO

!  --- FINAL UPDATE OF UNCERTAINTIES OF MIXING RATIOS WHEN APPROPRIATE
      INDXX = ISMIX
      DO KK = 1, NRET
         IF( IFPRF(KK) )THEN
            N = NLAYERS
            DO K = 1, N
               INDXX = INDXX + 1
               IF (ILOGRETRIEVAL(KK)/=0) THEN !MP
                  X(KK,K) = EXP(XHAT(INDXX))
               ELSE
                  X(KK,K) = XHAT(INDXX)*XORG(KK,K)
               END IF
               !PRINT *, 'VMR',INDXX,SHAT(INDXX,INDXX)
               !SIG(K,KK) = SQRT(ABS(SM(INDXX,INDXX)))*X(KK,K)
            END DO
         ELSE
            INDXX = INDXX + 1
            X(KK,:KMAX) = XHAT(INDXX)*XORG(KK,:KMAX)
              !print *, 'vmr',indxx,shat(indxx,indxx)
            !SERR = SQRT(ABS(SHAT(INDXX,INDXX)))
            WRITE (16, 350) NAME(IGAS(KK)), XHAT(INDXX) !, SERR
         ENDIF
      END DO

! --- PRINT OUT FINAL BACKGROUND PARAMETERS
      IF (NBACK > 1) THEN
         WRITE (16, 3006)
         DO I = 1, NFITS
            SELECT CASE (NBACK)
            CASE (2)
               J = I*(NBACK - 1)
              !print *, 'back',j,shat(j,j)
               SIGB = SQRT(ABS(SHAT(J,J)))
               WRITE (16, 3007) I, XHAT(J), -999. !SIGB
            CASE (3)
               J = I*(NBACK - 2)
              !print *, 'back',j,shat(j,j)
               SIGB = SQRT(ABS(SHAT(J,J)))
               SIGC = SQRT(ABS(SHAT(J+1,J+1)))
               WRITE (16, 3007) I, XHAT(J), -999., XHAT(J+1), -999. !SIGC
            END SELECT
         END DO
      ENDIF

!  --- PRINT OUT FINAL WAVENUMBER SHIFTS
!  --- CONVERT PARAMETER TO WAVENUMBERS
      !DO IBAND=1, NBAND
      IF( ISPARM > 0 ) THEN
         WRITE( 16, 3008)
         SELECT CASE (ISPARM)
         CASE (1)       ! SINGLE SHIFT FOR ALL BANDS (FITS = BANDS * SPECS)
            WRITE( 16, 3011)
            IBAND = 1
            J = NBKFIT+1
            WSHFT = 0.5D0*(WAVE3(1)+WAVE4(1))*((WAVFAC(IBAND) + XHAT(J)) - 1.D0)
            SWSHFT = SQRT(ABS(SHAT(J,J)))
            SWSHFT = 0.5D0*(WAVE3(1)+WAVE4(1))*((WAVFAC(IBAND) + SWSHFT) - 1.D0)
            WRITE (16, 3009) WSHFT, SWSHFT
         CASE (2)       ! INDEPENDENT SHIFT FOR EACH BAND
            WRITE( 16, 3012)
            DO I=1, NBAND
               J = NBKFIT+I
               WSHFT = 0.5D0*(WAVE3(I)+WAVE4(I))*((WAVFAC(I) + XHAT(J)) - 1.D0)
               SWSHFT = SQRT(ABS(SHAT(J,J)))
               SWSHFT = 0.5D0*(WAVE3(I)+WAVE4(I))*((WAVFAC(I) + SWSHFT) - 1.D0)
               WRITE (16, 3014) I, WSHFT, SWSHFT
            ENDDO
         CASE (3)        ! INDEPENDENT SHIFT FOR EACH FIT (BAND * SPEC)
            WRITE( 16, 3013)
            DO IBAND=1, NBAND
               N = NSCAN(IBAND)
               DO I=1, N
                  J = NBKFIT + IBAND + I - 1
                  WSHFT = 0.5D0*(WAVE3(I)+WAVE4(I))*((WAVFAC(IBAND) + XHAT(J)) - 1.D0)
                  SWSHFT = SQRT(ABS(SHAT(J,J)))
                  SWSHFT = 0.5D0*(WAVE3(I)+WAVE4(I))*((WAVFAC(IBAND) + SWSHFT) - 1.D0)
                  WRITE (16, 3015) IBAND, I, WSHFT, SWSHFT
               ENDDO
            ENDDO
         END SELECT
      ENDIF
      !ENDDO

!  --- PRINT OUT ZERO LEVEL OFFSETS - RETRIEVED OR NOT
      IF (NZERO /= 0) THEN
         K = NBKFIT + NSHIFT
         WRITE (16, 3001)
         DO I = 1, NBAND
            !IF (IZERO(I) == 0) CYCLE
            N = NSCAN(I)
            DO J = 1, N
               IF( IZERO(I) == 1 ) THEN
                  !print *, 'zero',j,shat(j,j)
                  K = K + 1
                  !SIGZ = SQRT(ABS(SHAT(K,K)))
                  WRITE (16, 3002) I, J, XHAT(K), -999. !SIGZ
               ELSE
                  WRITE (16, 3010) I, J, ZSHIFT(I,J), "      NOT FIT"
               ENDIF
            END DO
         END DO
      ENDIF

!  --- PRINT OUT SOLAR CO PARAMETERS
!  --- PRINTS IF THEY ARE USED REGARDLESS OF WHETHER THEY ARE FIT
      IF( IFCO )THEN
!         CPARM(4) = CPARM(4) - 1.D0
!         CPARM(5) = CPARM(5) - 1.D0
         WRITE (16, 3003)
         WRITE (16, 3004) (CPNAM(I),CPARM(I),I=1,5)
      ENDIF

!  --- PRINT OUT FINAL CHANNEL PARMS
      CALL PRINT_CHANNEL_PARMS( 16 )

!  --- SUM UP COLUMNS
      DO I = 1, NRET
         VERSUM(I,1) = X(I,1)*CCC(KVERT,1)
         VOSUM(I,1)  = XORG(I,1)*CCC(KVERT,1)
         IF( IFPRF(I) ) FX(I) = (SIG(1,I)*CCC(KVERT,1))**2
         DO K = 2, KMAX
            VERSUM(I,K) = VERSUM(I,K-1) + X(I,K)*CCC(KVERT,K)
            VOSUM(I,K)  = VOSUM(I,K-1)  + XORG(I,K)*CCC(KVERT,K)
            IF( IFPRF(I) ) FX(I) = FX(I) + (SIG(K,I)*CCC(KVERT,K))**2
         END DO
         IF( IFPRF(I) ) FX(I) = SQRT(ABS(FX(I)))
         IF( I == 1 ) SERR = 100.0D0*FX(I)/VERSUM(I,KMAX)
      END DO
      AIRCOL = 0.0D0
      DO K = 1, KMAX
         AIRCOL = AIRCOL + CCC(KVERT,K)
      ENDDO

!  --- CALCULATE DEGREES OF FREEDOM FOR SIGNAL USING APOSTERIORI SOLUTION
!  --- ONLY IF REALLY RETRIEVED, SOME MATRICES ARE NOT CALCULATED
      IF ( RETFLG ) CALL DOFS(NFIT,NVAR,ISMIX,NLEV)

      INDXX = ISMIX
      DO KK = 1, NRET
         IF( IFPRF(KK) )THEN
            N = NLAYERS
            DO K = 1, N
               INDXX = INDXX + 1
               IF (ILOGRETRIEVAL(KK)/=0) THEN !MP
                  X(KK,K) = EXP(XHAT(INDXX))
               ELSE
                  X(KK,K) = XHAT(INDXX)*XORG(KK,K)
               END IF
               !PRINT *, 'VMR',INDXX,SHAT(INDXX,INDXX)
               !SIG(KK,K) = SQRT(ABS(SM(INDXX,INDXX)))*XORG(KK,K)
            END DO
         end IF
      end DO

!  --- PRINT OUT RETRIEVED  PROFILES
      PRINT *,''
      IF ( RETFLG ) THEN
         INDXX = ISMIX
         DO K = NRET, 1, -1
            IF( .NOT. IFPRF(K) )CYCLE
            N = NLAYERS
            !PRINT 408, NAME(IGAS(K)), VERSUM(K,N), FX(K), 100.0D0*FX(K)/VERSUM(K,N)
            !WRITE(16,408) NAME(IGAS(K)), VERSUM(K,N), FX(K), 100.0D0*FX(K)/VERSUM(K,N)
            WRITE(16,406) !NAME(IGAS(K))
            DO KK = 1, N
               INDXX = INDXX + 1
               SIGMA = SQRT(ABS(SM(INDXX,INDXX)))*XORG(K,KK)
               !            print 0,SQRT(ABS(SM(INDXX,INDXX)))
               WRITE (16, 407) Z(KK), ZBAR(KK), XORG(K,KK), X(K,KK), SIGMA, VOSUM(K,KK), VERSUM(K,KK)
            END DO
         END DO

         WRITE(16,253)
      END IF

!  --- WRITE OUT TABLE OF PROFILES APRIORI ATMOSPHERE & VMR
      IF( F_WRTAPRF )CALL WRTAPRF( NRET, NLEV, KVERT )


!  --- WRITE OUT TABLE OF PROFILES RETRIEVED ATMOSPHERE & VMR
      IF( F_WRTRPRF )CALL WRTRPRF( NRET, NLEV, KVERT )


!  --- WRITE OUT OBSERVED, CALCULATED, AND DIFFERENCES - pbpfile
      IF( F_WRTPBP )CALL WRTPBP( TOBS, YHAT )


!  --- WRITE OUT STATE VECTOR
      IF( F_WRTSTV )CALL WRTSTV( NLEV, ITER, ISMIX, VERSUM, VOSUM, PNAME, XHAT, XAPR )


!  --- WRITE OUT A SUMMARY OF RETRIEVAL PARAMETERS
      print *, RETFLG
      IF( RETFLG .AND. F_WRTSUMRY ) CALL WRTSMRY( DOF, ITER, CHI_2_Y, FOVDIA, RMS, NLEV, VOSUM, VERSUM )


!  --- PRINT A SHORT SUMMARY TO THE CONSOLE
      IF( NRET .NE. 0 )THEN
         PRINT *,''
         PRINT 460, ( NAME(IGAS(I)), IFPRF(I), I=1, NRET )
         PRINT 461, (  VOSUM(I,NLEV), I=1, NRET )
         PRINT 461, ( VERSUM(I,NLEV), I=1, NRET )
         !PRINT 462, ITER, ITRMAX, RMS, NVAR, CONVERGE, DIVWARN, DOF, SERR, SNR, AIRCOL
         PRINT 463, ITER, ITRMAX, RMS, NVAR, CONVERGE, DIVWARN, DOF(2), SNR, CHI_2_Y!, AIRCOL
         PRINT *,''
      ENDIF

! --- UNCOMMENT NEXT LINE TO ACTIVATE OUTPUT OF RETRIEVED MIX FILE
!      CALL MIXOUT (KZMAXLAY, KVERT)

  700 CONTINUE


! --- KB MATRIX CALCULATED?

      IF( F_KB )CALL KBCALCULATE( NLEV )


      IF( F_WRTRAYTC )CALL FILECLOSE( 73, 1 )

! --- DEALLOCATE ARRAYS
      CALL RELEASE_MEM_INT
      CALL RELEASE_MEM_DIA
      CALL RELEASE_MEM_OPT
      CALL RELEASE_MEM_LP
      CALL RELEASE_MEM_RTP
      CALL RELEASE_MEM_SFT

      IF( IFCO )CALL SOLARFH ( 2 )
      RETURN

  664 CONTINUE
      WRITE (16, 591)
      RETURN
  665 CONTINUE
      WRITE (16, 592)
      RETURN
  667 CONTINUE
      WRITE (16, 250)
      RETURN
  668 CONTINUE
      WRITE (16, 251)
      RETURN

  250 FORMAT(/,' CAN ONLY RETRIEVE APODIZATION PARAMETERS FOR POLYNOMIAL FW.APOD_FCN.TYPE=2')
  251 FORMAT(/,' CAN ONLY RETRIEVE PHASE PARAMETERS FOR POLYNOMIAL FW.PHASE_FCN.TYPE=2')
  253 FORMAT(/,  '.END OF RETRIEVAL.')
  262 FORMAT( 5X,2000( 12X, A14 ))
  263 FORMAT( 2000I26 )

  350 FORMAT(/,' MOLECULE = ',A7,' PROFILE SCALE FACTOR =',F7.4,' +/-',F7.4)
  406 FORMAT(/,'  RETRIEVED VERTICAL PROFILE:',/,&
      '  Z[km] ZBAR[km] APRIORI_VMR RETRIEVED_VMR  SIGMA_VMR  APRIORI_COL  RETR&
      &IEVE_COL')
      !408 FORMAT(/' GAS: ', A7, ' COLUMN: ', ES12.4, ' +/- ',ES12.4,1X,F6.3,'%')
  407 FORMAT(2(F7.2),255(ES13.4))
  460 FORMAT(10(2X,A7,':',L2))
  461 FORMAT(10(1P,E12.4))
      ! divergance warning to abbr with Dvrg seems misleading
!  462 FORMAT(/"Itr/Mx:",I2.2,"/",I2.2," %RMS=",F5.3," FitPrm=",I3,' CVRG:',L1, &
!         ' DIVW:',L1," DOFS=",F5.3, " SNR=",F6.0, " CHI_2_Y=",ES12.4, " AIRCOL=",ES12.4)
  463 FORMAT(/"Itr/Mx:",I2.2,"/",I2.2," %RMS=",F5.3," FitPrm=",I3,' CVRG:',L1, &
         ' DIVW:',L1," DOFS=",F5.3, " SNR=",F6.0, " CHI_2_Y=",F9.4 )!, " AIRCOL=",ES12.4)

         !' Dvrg:',L1," DOFS=",F5.3, " %CERR=",F4.1, " SNR=",F5.0, " AIRCOL=",ES12.4)
!  464 FORMAT('pcol ', A7, 2(F6.3,1X), 2(F6.1,1X), 7ES12.4)
!  464 FORMAT('pcol ', A7, 2(F6.3,1X),  (F6.1,1X), 7ES12.4)

  502 FORMAT(5(1X,A14))
!  509 FORMAT(F9.3, F9.3, I5, I5)
  591 FORMAT(/,' ABORT NBACK OUT OF RANGE (1-3 VALID)')
  592 FORMAT('/ ABORT ISPARM OUT OF RANGE (0-3 VALID)')

 3001 FORMAT(/,' RETRIEVED ZERO LEVEL OFFSETS',/,&
               'BAND    ISCAN     ZERO     SIGMA')
 3002 FORMAT(I3,I9,2ES13.5)
 3003 FORMAT(/,' RETRIEVED SOLAR SIMULATION PARAMETERS:')
 3004 FORMAT(2(1X,A14,'=',1P,E13.4,1X))
! 3005 FORMAT(/,' RETRIEVED WAVENUMBER SCALE FACTOR =',F12.9)
 3006 FORMAT(/' RETRIEVED BACKGROUND PARAMETERS:'/,&
         ' BAND   SLOPE       SIGMA     CURVATURE     SIGMA')
 3007 FORMAT(I3,4(1P,E12.4))
 3008 FORMAT(/,' RETRIEVED WAVENUMBER SHIFTS:')
 3009 FORMAT(2ES13.5)
 3010 FORMAT(I3,I9,ES13.5,A)
 3011 FORMAT('  SHIFT          SIGMA')
 3012 FORMAT(' BAND     SHIFT          SIGMA')
 3013 FORMAT(' BAND   SCAN      SHIFT          SIGMA')
 3014 FORMAT(I3, 3X, 2ES13.5)
 3015 FORMAT(I5, I5, 3X, 2ES13.5)

 3695 FORMAT(/,' NFIT =',I5,' NVAR =',I3)
 3696 FORMAT(/,' PRINT OUT OF PARAMETERS FOR EACH ITERATION:',/)

      END SUBROUTINE SFIT


      SUBROUTINE KBCALCULATE( NLEV )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NLEV
      CHARACTER (LEN=255) :: VAL
      LOGICAL             :: HFLG, IFPRF_1_ORIG
      INTEGER             :: I, J, K, L1, L2, L3, ORIG_NVAR, POS, NL = 1

      WRITE(16,254)

! ---  SET DETAILED OUTPUT FILES TO FALSE
      HFLG         = .FALSE.
      F_WRTSTV     = .FALSE.
      F_WRTK       = .FALSE.
      F_WRTGASSPC  = .FALSE.
      F_WRTCHANNEL = .FALSE.
      F_WRTPARM    = .FALSE.
      F_WRTRAYTC   = .FALSE.
      XSC_DETAIL   = .FALSE.

      ifprf_1_orig = ifprf(1)

      ! Define new statevector for calculating KB-matrix
      if (f_kb_profile) then
         ! is the first retrieval gas already retrieved by column?
         val = adjustl(trim(s_kb_profile_gases))
         nrlgas = 0
         pos = index(adjustl(val),' ')
         do
            !                 print *,val
            if (len_trim(val).eq.0) exit
            nrprfgas = nrprfgas + 1
            if (pos.gt.0) then
               read(val(1:pos),*) s_kb_prf_gas(nrprfgas)
            else
               read(val(1:len_trim(adjustl(val))),*) s_kb_prf_gas(nrprfgas)
               exit
            end if
            val = adjustl(val(pos+1:len(val)))
            pos = index(trim(adjustl(val)),' ')
         end do
         do k = 1,nrprfgas
            do j = 1,nret
               if (trim(adjustl(s_kb_prf_gas(k))).eq.trim(adjustl(gas(j)))) then
                  ! only calculated if originally it was not a profile
                  if(.not.ifprf(j)) ifprf_kb(j) = .true.
                  ! but now it needs to be set to profile in order to setup correctly
                  ifprf(j) = .true.
               end if
            end do
         end do
      end if

! ---  DEFINE NEW STATEVECTOR FOR CALCULATING KB-MATRIX
      IF( F_KB_SLOPE .AND. NBACK.LT.2 ) then
         F_BACKG = .TRUE.
         NBACK = 2
      END IF
      IF( F_KB_CURVATURE .AND. NBACK.LT.3 )  then
         NBACK = 3
         F_BACKG = .TRUE.
      END IF
      IF( F_KB_SOLSHFT ) THEN
         F_RTSOL(4) = .TRUE.
         IFCO = .TRUE.
      END IF
      IF( F_KB_SOLSTRNTH ) THEN
         F_RTSOL(5) = .TRUE.
         IFCO = .TRUE.
      END IF
      ! ERROR FOR SIMPLE PHASE ONLY WHEN NO ERROR FOR EMPIRICAL PHASE IS NOT CALCULATED
      IF( F_KB_PHASE .AND..NOT. F_KB_EPHS)   IFPHASE = .TRUE.
      IF( F_KB_TEMP )                        IFTEMP = .TRUE.
      IF( F_KB_IFDIFF )                      IFDIFF = .TRUE.
      IF( F_KB_EAP.AND..NOT.F_EAPOD ) then
         F_RTAPOD = .TRUE.
         F_EAPOD  = .TRUE.
         IEAP = 2
         NEAP = 3
         EAPF(:NEAP) = 1.0D0
         EAPPAR = 1.0D0
      end IF
      IF( F_KB_EPHS.AND..NOT.F_EPHASE ) then
         F_RTPHASE = .TRUE.
         F_EPHASE = .TRUE.
         IEPHS = 2
         NEPHS = 3
         EPHSF(:NEPHS) = 1.0D0
         EPHSPAR = 1.0D0
      end IF
      IF( F_KB_ZSHIFT )  THEN
         IZERO(:NBAND) = 1
         F_ZSHIFT(:NBAND) = .true.
      END IF
      IF( F_KB_WSHIFT )THEN
         ISPARM = 3
         F_WSHIFT = .TRUE.
      ENDIF
      IF( F_KB_SZA )THEN
         IFSZA = 1
         DO I = 1, NSPEC
            ASTANG0(I) = ASTANG(I)
         ENDDO
      ENDIF
      DO I = 1, NBAND
         IF( F_KB_FOV )THEN
            IFFOV = 1
            OMEGA0(I) = OMEGA(I)
         ENDIF
         IF( F_KB_OPD )THEN
            IFOPD = 1
            PMAX0(I) = PMAX(I)
         ENDIF
      ENDDO

      if( f_kb_line )then
         ifline = 1
         ! --- find for which gases are Kb for line parameters are calculated
         select case (trim(adjustl(s_kb_line_gases)))
         case ('target')
            s_kb_line_gas(1) = trim(adjustl(gas(1)))
            niline = 1
            npline = 1
            ntline = 1
            nrlgas = 1
         case ('retrieval')
            do k = 1,nret
               s_kb_line_gas(k) = trim(adjustl(gas(k)))
            end do
            niline = nret
            npline = nret
            ntline = nret
            nrlgas = nret
         case default
            val = adjustl(trim(s_kb_line_gases))
            nrlgas = 0
            pos = index(adjustl(val),' ')
            do
!                 print *,val
               if (len_trim(val).eq.0) exit
               nrlgas = nrlgas + 1
               if (pos.gt.0) then
                  read(val(1:pos),*) s_kb_line_gas(nrlgas)
               else
                  read(val(1:len_trim(adjustl(val))),*) s_kb_line_gas(nrlgas)
                  exit
               end if
               val = adjustl(val(pos+1:len(val)))
               pos = index(trim(adjustl(val)),' ')
            end do
            niline = nrlgas
            npline = nrlgas
            ntline = nrlgas
         end select
      end if

! --- SET THE ENTRIES OF THE STATEVECTOR AS A PRIORI IN THE NEW PARM VECTOR TO MAKE
!     SURE THE KB MATRICES ARE CALCULATED AS DEVIATIONS FROM THE RETRIEVED STATE

! --- SETUP NEW PARM ARRAY
      ORIG_PNAME(:NVAR) = PNAME(:NVAR)
      ORIG_NVAR = NVAR
      RETFLG = .FALSE.
      CALL INIT_PARM()

      IS_IN_KB(:NVAR) = .TRUE.
      do k=1, NVAR
         do i=1, ORIG_NVAR
            ! --- DWNUMSHIFT IS SET TO RETRIEVED PARAMETER OF IWNUMSHIFT
            if( ORIG_PNAME(i)(:9) .eq. 'IWNumShft' .and. PNAME(k)(:9) .eq. 'DWNumShft' )then
               print *, PNAME(k)
               parm(k) = xhat(i)
               IS_IN_KB(k) = .false.
               !               exit
            end if
            ! --- IWNUMSHIFT GETS SET TO RETRIEVED VALUE OF SWNUMSHIFT
            if( ORIG_PNAME(i)(:9) .eq. 'SWNumShft' .and. PNAME(k)(:9) .eq. 'IWNumShft' )then
               ORIG_PNAME(i) = ''
               parm(k) = xhat(i)
               IS_IN_KB(k) = .false.
               exit
            end if
! --- DON'T COMPUTE K FOR RETRIEVED B
            if (ORIG_PNAME(i).eq.PNAME(k)) then
               ORIG_PNAME(i) = ''
               parm(k) = xhat(i)
               IS_IN_KB(k) = .false.
               ! CHECK IF THE ORIGINALLY RETRIEVED GAS IS A COLUMN
               DO J = 1,NRET
                  ! IF SO, CALCULATE A KB ENTRY FOR THIS GAS AS A PROFILE
                  IF ((PNAME(K).EQ.GAS(J)).AND.(IFPRF_KB(J).AND.IFPRF(J))) THEN
                     IS_IN_KB(K) = .TRUE.
                     PARM(K:K+NLEV) = XHAT(I)
                  END IF
               END DO
               exit
            end if
         end do
      end do


!         PRINT *, PNAME
!         print *, XHAT(:ORIG_NVAR)
!         print *, PARM(:NVAR)
!         print *, IS_IN_KMATRIX(:NVAR)


      WRITE(16,250)
      DO I=1, NVAR
         WRITE(16,251) I, PNAME(I), PARM(I), SPARM(I), IS_IN_KB(I)
      ENDDO

      IF( ALLOCATED(KHAT) )DEALLOCATE( KHAT )
      ALLOCATE(KHAT(NFIT,NVAR))
      TOBS(:NFIT) = TOBS_ORIG(:NFIT)

      !write(100,*) YHAT(:NFIT)
      CALL FM(PARM, YHAT, KHAT, NFIT, NVAR, .TRUE., -1, HFLG )
      !write(100,*) YHAT(:NFIT)
      !close(100)

      ! APPEND NAMES ILINE ENTRIES IN KB MATRIX WITH GASNAMES
      ! MIGTH BE BETTER IN THE SUBFUNCTION INIT_PARM, BUT ALL IMPLICATIONS HAVE TO BE CHECKED!!!

      L1 = 1
      L2 = 1
      L3 = 1
      DO I = 1,NVAR
         SELECT CASE (PNAME(i))
         CASE ('LineInt')
            PNAME(I) = 'LineInt'//'_'//trim(s_kb_line_gas(L1))
            L1 = L1 + 1
         CASE ('LinePAir')
            PNAME(I) = 'LinePAir'//'_'//trim(s_kb_line_gas(L2))
            L2 = L2 + 1
         CASE ('LineTAir')
            PNAME(I) = 'LineTAir'//'_'//trim(s_kb_line_gas(L3))
            L3 = L3 + 1
         END SELECT
      END DO

         

      CALL FILEOPEN(90,1)
      WRITE(90,*) TRIM(TAG), ' KB VECTORS FOR MODEL PARAMETERS BI'
      WRITE(90,*) NFIT, COUNT(IS_IN_KB(:NVAR),1), -1, -1
      WRITE(90,260) ADJUSTR( PACK( PNAME(:NVAR), IS_IN_KB(:NVAR) ))
      DO J=1, NFIT
         WRITE(90,261) PACK(KHAT(J,:), IS_IN_KB(:NVAR))
      END DO
      CALL FILECLOSE( 90, 1 )

! --- WRITE OUT AB (G*KB) IN FRACTIONS OF A PRIORI, THIS CORRESPONDS TO FORMULA
!     3.16 ON PAGE 48 IN RODGERS BOOK AND CAN DIRECTLY BE USED FOR THE ERROR CALCULATION
      IF( F_WRTAB )THEN
         IF( ALLOCATED(A) )DEALLOCATE(A)
         ALLOCATE(A(ORIG_NVAR, NVAR))
         CALL MULT( G, KHAT, A, ORIG_NVAR, NFIT, NVAR )

         CALL FILEOPEN(92,1)
         WRITE(92,*) TRIM(TAG), ' DY#KB=AB MATRIX FOR MODEL PARAMETERS BI'
         IF( IFPRF_1_ORIG )NL = NLEV
         WRITE(92,*) NLEV, COUNT(IS_IN_KB(:NVAR),1), -1, -1
         WRITE(92,260) ADJUSTR( PACK( PNAME(:NVAR), IS_IN_KB(:NVAR) ))
         DO J=1, NLEV
            WRITE(92,261) PACK(A(J+ISMIX, :), IS_IN_KB(:NVAR))
         ENDDO
         CALL FILECLOSE( 92, 1 )
      ENDIF

      RETURN

  251 FORMAT(I5, 3X, A10, 2F15.7, 4X, L1 )
  250 FORMAT(/, '    I   PARAMETER   POSTERORI_VALUE     SIGMA     COMPUTED_IN_KB')
!  252 FORMAT(' COMPUTING KB FOR PARAMTETERS :',255(/,3X,A14))
  254 FORMAT(/, 'BEGIN KB CALCULATIONS:',/)
  260 FORMAT( 2000( 12X, A17 ))
  261 FORMAT( 2000ES29.18E3 )


      END SUBROUTINE KBCALCULATE

      SUBROUTINE RELEASE_MEM_SFT

      IF( ALLOCATED( XHAT ))DEALLOCATE( XHAT )
      IF( ALLOCATED( YHAT ))DEALLOCATE( YHAT )
      IF( ALLOCATED( XAPR ))DEALLOCATE( XAPR )
      IF( ALLOCATED( Y )   )DEALLOCATE( Y )
      IF( ALLOCATED( SED ) )DEALLOCATE( SED )

      END SUBROUTINE RELEASE_MEM_SFT

      END MODULE SFIT4
