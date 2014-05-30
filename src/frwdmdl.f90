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

      MODULE FRWDMDL

      USE PARAMS
      USE TRANSMIS
      USE DATAFILES
      USE MOLCPARAM
      USE SOLAR
      USE SYNSPEC
      USE CHANNEL
      USE BANDPARAM
      USE RETVPARAM
      USE WRITEOUT
      USE INITIALIZE

      IMPLICIT NONE

      !CHARACTER, DIMENSION(NMAX)           :: PNAME*14
      CHARACTER (LEN=7), DIMENSION(MOLMAX) :: SDV_GAS
      CHARACTER (LEN=7), DIMENSION(MOLMAX) :: LM_GAS

      CONTAINS

!------------------------------------------------------------------------------
      SUBROUTINE FM(XN, YN, KN, NFIT, NVAR, KFLG, ITER, TFLG )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NFIT !, M - COLUMNS - 1 SPECTRUM
      INTEGER, INTENT(IN) :: NVAR !, N - ROWS    - FOR EACH FIT PARAMETER
      INTEGER, INTENT(IN) :: ITER
      LOGICAL, INTENT(IN) :: KFLG
      LOGICAL, INTENT(INOUT) :: TFLG
      REAL(DOUBLE), INTENT(IN) :: XN(NVAR)
      REAL(DOUBLE), INTENT(OUT) :: YN(NFIT)
      REAL(DOUBLE), INTENT(OUT) :: KN(NFIT,NVAR)

      LOGICAL :: BUG1 = .FALSE., IFCOSAVE=.FALSE.

      CHARACTER :: GASFNAME*(IFLNMSZ)
      CHARACTER :: TITLE*(80)
      INTEGER :: III, NGB
      INTEGER :: NVAR1, KFIT, KFIT2, KZERO, KPHASE, JATMOS, IPARM, I, NCOUNT, &
         KK, K, MXONE, IBAND, N, JSCAN, MONONE, N1, N2, N3, J, MSHIFT, NS, NR, &
         NS1, NS2
      LOGICAL :: XRET, TRET, FLINE, FSZA

      REAL(DOUBLE), DIMENSION(3)      :: B
      REAL(DOUBLE), DIMENSION(NMAX)   :: PARM
      REAL(DOUBLE), DIMENSION(MMAX)   :: YC
      REAL(DOUBLE), DIMENSION(NFIT)   :: AMP_Y
      REAL(DOUBLE), DIMENSION(NMONSM) :: Y_INFTY, DELTA_Y
      REAL(DOUBLE), DIMENSION(MAXSPE) :: ZSHIFTSAV
      REAL(DOUBLE), DIMENSION(NFIT)   :: WAVE_X
      REAL(DOUBLE) :: DEL, SUMSQ, WSCALE, DWAVE, DSHIFT, FRACS, PHI, SMM, YS, &
         BKGND, YCAVE, FX, TEMPP, YCMAX !, STDEV
      REAL(DOUBLE) , DIMENSION(:,:), allocatable    :: store_line

      COMPLEX(DBLE_COMPLEX) :: TCALL, TCALH, TCALI

!  --- PARAMETER INCREMENT FOR PARTIALS

      COMPLEX(DBLE_COMPLEX), DIMENSION(:), ALLOCATABLE :: TCONVSAV
      INTEGER :: NAERR

! -----------------------------------------------------------
!     COMPUTES MONOCHROMATIC SPECTRUM, CONVOLVES WITH INSTRUMENTAL
!     PROFILE, AND FINDS OBSERVED MINUS CALCULATED AMPLITUDES
!     AT ALL POINTS -
!     COMPUTE NEW SPECTRUM (IPARM=1) WITH UPDATED XN
!     THEN COMPUTE KN (NVAR X NFIT) FOR NEXT ITERATION
!     OUTER LOOP - NVAR  INNER LOOP NFIT
! -----------------------------------------------------------

!nbkfit nshift nzero nsolar neaprt nephsrt ndiff nphase  nret*(kmax||1) ntemp channels
!                    do fft                              ismix

      TFLG = .FALSE.
      BUG1 = .FALSE. !.TRUE.
      ! if line parameters are disturbed, get some space to store original ones
      if (nrlgas /= 0 .and. .not. allocated(store_line)) then
         allocate(store_line(4,LNMAX)) 
         store_line(:,:) = 0.0d0
      end if
      IF (KFLG) THEN
         NVAR1 = NVAR + 1
      ELSE
         NVAR1 = 1
      ENDIF

!      IF( F_WRTCHANNEL )THEN            !FOR_CHANNEL_VALUE_FILES
!         IF (ITER == 1) CALL FILEOPEN(30, 1 )
!         WRITE(30,'(A,A,A)' )TRIM(TAG), ' CHANNEL SPECTRUM ITER 1'
!         IF (ITER == 2) CALL FILEOPEN(40, 1 )
!         WRITE(40,'(A,A,A)' )TRIM(TAG), ' CHANNEL SPECTRUM ITER 2'
!      ENDIF


      PARAM: DO ICOUNT = 1, NVAR1
         KFIT   = 0
         KFIT2  = 0
         KZERO  = 0
         KPHASE = 0
         JATMOS = 0
         SUMSQ  = 0.D0
         IPARM  = 0
         DEL = 0.1D-05

         IF( ICOUNT .GT. 1 )IPARM = ICOUNT -1

         IF( BUG1 )PRINT *,''
         IF( BUG1 )WRITE(0,202) 'TOP :  ', ITER, ICOUNT, IPARM, NTEMP1, NVAR1

! --- RESET PARM TO PERTURB EACH INDIVIDUALLY
         PARM(:NVAR) = XN
         IF (ICOUNT .GT. 1) THEN
            PARM(IPARM) = PARM(IPARM) + DEL
            IF( BUG1 )WRITE(0,203) '     PARM: ', IPARM, PNAME(IPARM), PARM(IPARM)+del, PARM(IPARM)
         ENDIF

!  --- ADJUST THESE PARAMETERS IN BAND/SPEC LOOPS BELOW
         NCOUNT = NBKFIT + NSHIFT + NZERO

!  --- SOLAR SPECTRUM PARAMETERS
         IF( IFCO )THEN
            DO I = 1, 5
               IF( F_RTSOL(I) )THEN
                  NCOUNT = NCOUNT + 1
                  CPARM(I) = PARM(NCOUNT)!*CIPARM(I)
                  !print *, ncount, nsolar1, parm(ncount)
               ENDIF
            END DO
         ENDIF

!  --- EMPIRICAL APODIZATION
         IF( F_RTAPOD )THEN
            IF (NEAPRT > 0) THEN
               EAPF(:NEAPRT) = EAPF0(:NEAPRT)*PARM(NCOUNT+1:NEAPRT+NCOUNT)
               NCOUNT = NEAPRT + NCOUNT
            ENDIF
         ENDIF

!  --- EMPIRICAL PHASE FUNCTION
         IF( F_RTPHASE )THEN
            IF (NEPHSRT > 0) THEN
               EPHSF(:NEPHSRT) = EPHSF0(:NEPHSRT)*PARM(NCOUNT+1:NEPHSRT+NCOUNT)
               NCOUNT = NEPHSRT + NCOUNT
            ENDIF
         ENDIF

!  --- UPDATE DIFFERENTIAL SHIFT
         IF( IFDIFF )THEN
            IF (NDIFF > 0) THEN
               ISHIFT(:NDIFF) = NINT(PARM(NCOUNT+1:NDIFF+NCOUNT)/.6D-06)
               NCOUNT = NDIFF + NCOUNT
            ENDIF
         ENDIF

!  ---  UPDATE NCOUNT INDEX TO INCLUDE PHASE ERROR PARAMETERS
         NCOUNT = NCOUNT + NPHASE

!  --- UPDATE THE LINE PARAMETERS ACCORDING TO THE SETUP
         ! --- Intensities
         if ( nrlgas /= 0 ) then
            K = ICOUNT - NCOUNT -1
            KK = NILINE + NPLINE + NTLINE + 2
            FLINE = .FALSE.
            ! setup3 one time more than perturbation, the last time with the original line parameters again.
            if ((K.gt.0).and.(k.lt.kk)) then
                  do k = 1, nrlgas
                     DO I=LINE1(1), LINE2(NBAND)
                        IF( TRIM(s_kb_line_gas(k)) .EQ.  TRIM(NAME(ICODE(LGAS(I)))))THEN
                           !                        print *, i, k, AZERO(i), ST296(i), ICODE(LGAS(I)), trim(NAME(ICODE(LGAS(I)))), ' ', trim(s_kb_line_gas(k))
                           store_line(2,i) = ST296(i)
                           store_line(3,i) = AAA(i)
                           store_line(4,i) = TDLIN(i)
                           if (niline /= 0) ST296(i) = ST296(i)* (1.0d0 + parm(ncount+k))
                           if (npline /= 0) AAA(i) = AAA(i)* (1.0d0 + parm(ncount+niline+k))
                           if (ntline /= 0) TDLIN(i) = TDLIN(i)* (1.0d0 + parm(ncount+niline+npline+k))
                        end IF
                     end do
                  end DO
                  CALL SETUP3( XSC_DETAIL, -1 )
                  ! set back line parameters
                  do k = 1, nrlgas
                     DO I=LINE1(1), LINE2(NBAND)
                        IF( TRIM(s_kb_line_gas(k)) .EQ.  TRIM(NAME(ICODE(LGAS(I)))))THEN
                           !                        print *, i, k, AZERO(i), ST296(i), ICODE(LGAS(I)), trim(NAME(ICODE(LGAS(I)))), ' ', trim(s_kb_line_gas(k))
                           ST296(i) = store_line(2,i)
                           AAA(i)   = store_line(3,i)
                           TDLIN(i) = store_line(4,i)
                        end IF
                     end do
                  end DO
                  FLINE = .true.
               end if
               NCOUNT = NCOUNT + NILINE + NPLINE + NTLINE
            end if

         FSZA = .false.
         if (ifsza /= 0) then
            ! setup2 and setup3 must run one time more than perturbation sza in order to get the old state again
            k = ICOUNT - NCOUNT -1
            do kk = 1,nspec
               astang(kk) = astang0(kk)*(1.0d0+parm(ncount+kk))
            end do
            if (k.gt. 0 .and. k.lt.nspec+2) THEN
               CALL LBLATM( 0, KMAX )
               CALL SETUP3( XSC_DETAIL, -1 )
               FSZA = .true.
            end if
            NCOUNT = NCOUNT + NSPEC
         end if

         do k = 1,nband
            ! Error in Field of View
            if (iffov /= 0) then
               ncount = ncount + 1
               OMEGA(k) = OMEGA0(k)*(1.0d0 + parm(ncount))
            end if
            ! Error in Field of MaxOPD
            if (ifopd /= 0) then
               ncount = ncount + 1
               if (ICOUNT == NCOUNT+1) then
                  ! usual DEL = 0.1D-5 is to small to cause any KB.
                  parm(ncount) = parm(ncount) - DEL + 0.1
                  DEL = 0.1
               end if
               PMAX(k) = PMAX0(k) * (1.0d0 + parm(ncount))
            end if
         end do


!  ---  UPDATE VMRS OF RETRIEVAL GASES
         DELTA_Y(:NFIT) = 0.0D0
         DELTA_Y(:NMONSM) = 0.0D0
         XRET = .FALSE.
         DO KK = 1, NRET
            IF( IFPRF(KK) )THEN
               K = ICOUNT-NCOUNT-1
               IF( ANALYTIC_K .AND. K .GT. 0 .AND. K .LT. KMAX+1 )THEN
                  ! K-MATRICES ARE CALCULATED SEMI-ANALYTICALLY.
                  ! THIS IS WAY FASTER.
                  XRET=.TRUE.
                  IF (IEMISSION.EQ.1) THEN
                     IF (K.GT.1) THEN
                        DELTA_Y(:NMONSM) = &
                          - CROSS_FACMAS(KK,K,:NMONSM) * (TCALC_E(2,:NMONSM, KMAX+1) &
                          + TCALC_S(2,:NMONSM,K-1) - TCALC_E(2,:NMONSM,K))
                     ELSE
                       DELTA_Y(:NMONSM) = -CROSS_FACMAS(KK,K,:NMONSM)* TCALC_E(2,:NMONSM, KMAX+1)
                     ENDIF
                  ELSE
                     DELTA_Y(:NMONSM) = -CROSS_FACMAS(KK,K,:NMONSM)*Y_INFTY(:NMONSM)
                  ENDIF
                  IF( ILOGRETRIEVAL(KK)/=0 )THEN
                     ! THE KMATRIX IS ANALYTICALLY CALCULATED, NEED THE ORIGINAL STATE VECTOR
                     PARM(IPARM) = XN(IPARM)
                     X(KK,:KMAX) = EXP(PARM(NCOUNT+1:KMAX+NCOUNT))
                     DELTA_Y(:NMONSM) = DELTA_Y(:NMONSM) / XORG(KK,K) * X(KK,K)
                  ENDIF
               ELSEIF( ILOGRETRIEVAL(KK) /= 0 )THEN
                  X(KK,:KMAX) = EXP(PARM(NCOUNT+1:KMAX+NCOUNT))
               ELSE
                  X(KK,:KMAX) = PARM(NCOUNT+1:KMAX+NCOUNT)*XORG(KK,:KMAX)
               ENDIF
               NCOUNT = NCOUNT + KMAX
            ELSE

! --- SCALING VERTICAL DISTRIBUTION
               NCOUNT = NCOUNT + 1
               X(KK,:KMAX) = PARM(NCOUNT)*XORG(KK,:KMAX)
            ENDIF
         END DO

! --- TEMPERATURE RETRIEVAL
         IF( IFTEMP ) THEN
            !IF( BUG1 )PRINT *, IFTEMP, IPARM, NCOUNT, NTEMP1, NTEMP, PARM(NCOUNT+1:NCOUNT+1)
            TRET = .FALSE.
            ! --- ONLY CONSIDERING PROFILE FIT
            K = IPARM - NCOUNT 
            ! KMAX + 1 passes to un-perturb final temperature
            IF( K .GE. 1 .AND. K .LE. KMAX+1)THEN
               !IF( NCOUNT+1 .GE. NTEMP1 .AND. NCOUNT+1 .LT. NTEMP1 + NTEMP )THEN
               TRET = .TRUE.
               !if(ntemp1 .eq. ncount+1) print*, k, t(k), torg(k)
               !print*, ncount+1, ntemp1, k
               T(:KMAX) = PARM(NCOUNT+1:NCOUNT+KMAX) * TORG(:KMAX)
               !print*,PARM(NCOUNT+1:NCOUNT+KMAX)
               !if(ntemp1 .eq. ncount+1) print*, k, t(k), torg(k), ITER, KMAX
               NCOUNT = NCOUNT + KMAX
                  !CALL LBLATM( ITER, KMAX )
                  IF (K .GT. KMAX) K = KMAX
                  CALL MASSPATH( K )
                  CALL SETUP3( XSC_DETAIL, k )
            ENDIF ! K
            !write(0,'(2f14.5)') (t(kk),torg(kk), kk=1,kmax)
         ENDIF ! IFTEMP


! --- UPDATE TO SOLAR SPECTRAL CALCULATIONS - ALL BANDS AT ONCE
         IF (NSOLAR /= 0) THEN
            IF (IPARM == 0) GO TO 6
            IF ((IPARM .LT. NSOLAR1) .OR. (IPARM .GT. NSOLAR1+NSOLAR)) GOTO 7
    6       CONTINUE
            IF( BUG1 )PRINT *,'    SOLSARFH', IPARM, NSOLAR1, NSOLAR
            CALL SOLARFH( 1 )
         ENDIF

!  --- SKIP OVER MONOCHROMATIC CALCULATIONS IF NOT REQUIRED
    7    CONTINUE

         IF (IPARM == 0) GO TO 8
!  --- DO NOT CALCULATE SPECTRUM IF NOT NECESSARY.
         NCOUNT = NBKFIT + NSHIFT + NZERO + NSOLAR + NEAPRT  + NEPHSRT +  NPHASE! + NDIFF
         IF( BUG1 )PRINT *, '    IPARM, NCOUNT : ', IPARM, NCOUNT
         IF (IPARM .LT. NCOUNT) GO TO 9

!  --- COMPUTE MONOCHROMATIC TRANSMITTANCES
!  --- ANAYLITC K-MATICES MAY BE CHOSEN IN PARAM_M.F90 MP
    8    CONTINUE

         IF ((.NOT.ANALYTIC_K).OR.(.NOT.XRET).OR.(TRET).OR.(ICOUNT.EQ.1).or.FLINE.or.FSZA) THEN
            CALL TALL
            IF( BUG1 )PRINT*, '    TALL', IPARM
            !print*, nmonsm, TCALC(1,:100)
            !stop
         ELSE
            IF( BUG1 )PRINT*, '    TALL/DIFF', IPARM
            ! THERE COME SOME MORE OPERATIONS ON THE NEW SPECTRUM.
            ! NOT YET INVESTIGATED WHICH ARE LINEAR.
            ! DECREASE NEW SPECTRUM BY DEL, WORKS MORE EXACT WITH FURTHER
            ! OPERATIONS. MAY BE REMOVED LATER ON. MP
            TCALC(1,:NMONSM) = Y_INFTY(:NMONSM) + DELTA_Y(:NMONSM)*DEL
            TCALC(2,:NMONSM) = TCALC(1,:NMONSM)
         END IF

         IF (ICOUNT.EQ.1) Y_INFTY(:NMONSM) = TCALC(2,:NMONSM)

    9    CONTINUE
         MONONE = 1
         MXONE = 1

!  --- RETRIEVE CHANNEL PARMS FROM STATE VECTOR-----------------------------!PWJ

         CALL RETRIEVE_CHANNEL_PARMS (PARM)

!  --- LOOP OVER BANDPASSES ----------------------------------------------------
         BAND: DO IBAND = 1, NBAND
            N = NSCAN(IBAND)
            IF (N == 0) CYCLE

!  --- LOOP OVER SPECTRA -------------------------------------------------------
            SPEC: DO JSCAN = 1, N

!  --- DETERMINE CURRENT BACKGROUND PARAMETERS CORRESPONDING
!  --- TO SPECTRUM AND BANDPASS
               B(1) = 1.0D0
               IF (NBACK <= 1) THEN
                  B(2) = 0.D0
                  B(3) = 0.D0
               ELSE
                  KFIT = KFIT + 1
                  B(2) = PARM(KFIT)
                  IF (NBACK <= 2) THEN
                     B(3) = 0.D0
                  ELSE
                     KFIT = KFIT + 1
                     B(3) = PARM(KFIT)
                  ENDIF
               ENDIF

!  --- DETERMINE CURRENT WAVENUMBER SCALE MULTIPLIER
               IF (ISPARM /= 2) THEN
                  IF (ISPARM == 0) GO TO 3
                  IF (ISPARM == 3) GO TO 131

!  --- SINGLE PARAMETER FOR ALL BANDPASSES
                  WSCALE = PARM(NBKFIT+1)
                  GO TO 14
               ENDIF

!  --- INDEPENDENT PARAMETER FOR EACH BANDPASS
               WSCALE = PARM(NBKFIT+IBAND)
               GO TO 14

!  --- INDEPENDENT PARAMETER FOR EACH FIT
  131          CONTINUE
               KFIT2 = KFIT2 + 1
               WSCALE = PARM(NBKFIT+KFIT2)
               GO TO 14

!  --- NO WAVENUMBER SHIFT
    3          CONTINUE
               WSCALE = 0.0D0

!  --- CALCULATE SHIFT IN WAVENUMBERS
   14          CONTINUE

               DWAVE = 0.5D0*(WAVE3(IBAND)+WAVE4(IBAND))*((WAVFAC(IBAND) + WSCALE) - 1.D0)

!  --- CALCULATE NUMBER OF MONOCHROMATIC POINTS TO SHIFT
               DSHIFT = DWAVE/DN(IBAND)
               MSHIFT = NINT(DSHIFT)
               FRACS = DSHIFT - MSHIFT

!  --- DETERMINE ZERO LEVEL OFFSET TO APPLY
               IF (IZERO(IBAND) == 1) THEN
                  ! kzero increments for each band
                  KZERO = KZERO + 1
                  ZSHIFT(IBAND,JSCAN) = PARM(NBKFIT+NSHIFT+KZERO)
                  ZSHIFTSAV(JSCAN) = ZSHIFT(IBAND,JSCAN)
               ELSE IF (IZERO(IBAND) == 2 ) THEN
                  ! if we're not calculating it then use shift from band from this spec that we are fitting
                  ZSHIFT(IBAND,JSCAN) = ZSHIFTSAV(JSCAN)
               ENDIF

!  --- DETERMINE PHASE ERROR TO APPLY
               PHI = 0.D0
               IF( IFPHASE )THEN
                  KPHASE = KPHASE + 1
                  PHI = PARM(NBKFIT+NSHIFT+NZERO+NEAPRT+NEPHSRT+NSOLAR+NDIFF+KPHASE)
!                  NCOUNT = NBKFIT+NSHIFT+NZERO+NEAPRT+NEPHSRT+NSOLAR+NDIFF+KPHASE
               ENDIF

!  --- COMPUTE FFTS IF REQUIRED
               IF (IPARM == 0) GO TO 15
               NCOUNT = NBKFIT + NSHIFT + NZERO
               IF (IPARM <= NCOUNT) GO TO 16
   15          CONTINUE
               IF( BUG1 )PRINT *, '    FSPEC...', IPARM
               CALL FSPEC1 (IBAND, MONONE, MXONE)
               CALL FSPEC2 (IBAND, MONONE, PHI)

!  --- COMPUTE RESIDUALS
   16          CONTINUE

               N1 = NSTART(IBAND) + MSHIFT + MONONE - 1
               N2 = N1 + (NPRIM(IBAND)-1)*NSPAC(IBAND)

!  --- CHECK FOR TCALC OVERFLOW OF BAND SPACE
               IF (N1<MONONE .OR. N2>=MONONE+NM(IBAND)) GO TO 17
               N3 = NPRIM(IBAND)

               SMM = 0.D0
               DO J = 1, N3

                  I = N1 + (J - 1)*NSPAC(IBAND)
                  JATMOS = JATMOS + 1
                  TCALL = TCONV(I)
                  TCALH = TCONV(I+1)
                  TCALI = TCALL + FRACS*(TCALH - TCALL)
                  YS = (J - 1)*SPAC(IBAND)
                  WAVE_X(JATMOS) = YS + WSTART(IBAND)
                  BKGND = B(1)*(1.0D0 + B(2)*YS+B(3)*YS*YS)
                  BKGND = BKGND*(1.0D0/(1.0D0 + ZSHIFT(IBAND,JSCAN)))

!-- FIT CHANNEL PARMS IF NEEDED ----------------------------------------!PWJ

                  IF (NBEAM_OF_BAND(IBAND) /= 0) THEN
                     YC(JATMOS) = FIT_CHANNEL_PARMS(IBAND,JSCAN,J,DBLE(BKGND),TCALI,TEMPP)
                     IF( F_WRTCHANNEL )THEN
                        AMP_Y(JATMOS) = TEMPP    !for_channel_value_files
                     ENDIF
                  ELSE
                     YC(JATMOS) = BKGND*(DBLE(TCALI) + ZSHIFT(IBAND,JSCAN))
                  ENDIF
!print *,jatmos, yc(jatmos)
                  SMM = SMM + YC(JATMOS)
               END DO

               ! -- normalization of spectra only when absorption
               ! spectra only or normalization is explicitely
               ! required for emission spectra. mp
               IF (IEMISSION.EQ.0 .OR. IENORM(IBAND).NE.0) THEN
                  YCAVE = SMM/N3
                  YC(JATMOS-N3+1:JATMOS) = YC(JATMOS-N3+1:JATMOS)/YCAVE
               ELSE
                  YCAVE = 1.0D0
               END IF

!  --- WRITE SPECTRA BY GAS, BAND, SCAN & ITERATION
               IF( F_WRTGASSPC .AND. (ICOUNT .EQ. NVAR1) &
                 .AND.( (GASOUTTYPE .EQ. 1) .AND. (ITER .EQ. -1)  &
                 .OR.   (GASOUTTYPE .EQ. 2) ))THEN
!  --- SAVE TCALC
                  ALLOCATE (TCONVSAV(NMONSM), STAT=NAERR)
                  IF (NAERR /= 0) THEN
                       WRITE (16, *) 'FWRDMDL: COULD NOT ALLOCATE TCONVSAV ARRAY, ERROR NUMBER = ', NAERR
                       WRITE ( 0, *) 'FWRDMDL: COULD NOT ALLOCATE TCONVSAV ARRAY, ERROR NUMBER = ', NAERR
                       CALL SHUTDOWN
                       STOP 2
                  ENDIF
                  TCONVSAV(:NMONSM) = TCONV(:NMONSM)
                  TCALC(1,:NMONSM) = TCALC(2,:NMONSM)

!  --- TEMPORARILY TURN OF FLAG SO FSPEC1 WILL NOT APPLY TCO
                  IFCOSAVE = IFCO
                  IFCO = .FALSE.

!  --- ALL SPECTRA
                  IF( GASOUTTYPE .EQ. 1 .AND. ITER .EQ. -1 )THEN
                     WRITE(GASFNAME,610)IBAND,JSCAN
                  ELSEIF( GASOUTTYPE .EQ. 2 )THEN
                     IF (ITER == -1) THEN
                        WRITE(GASFNAME,610)IBAND,JSCAN
                     ELSE
                        WRITE(GASFNAME,620)IBAND,JSCAN,ITER
                     ENDIF
                  ENDIF

                  WRITE(TITLE,710) 'ALL', IBAND, JSCAN, ITER
                  OPEN(UNIT=80, FILE=GASFNAME, STATUS='REPLACE', ERR=555)
                  WRITE (80, 640) TITLE
                  WRITE (80, *) WSTART(IBAND), WSTOP(IBAND), SPAC(IBAND), N3
                  YCMAX = maxval(YC(JATMOS-N3+1:JATMOS))
                  DO III=JATMOS-N3+1,JATMOS
                     WRITE (80, *) YC(III)/YCMAX
                  ENDDO
                  CLOSE (80)

!  --- LOOP OVER GASES IN BAND
                  DO NR = 1, NRETB(IBAND)
                     !print*, iband, nr, nretb(iband), igasb(iband,nr),  icount, NGASB(iband,nr), GASB(IBAND,NR)
                     NGB = NGASB(IBAND,NR)
                     CALL GASNTRAN(NGB,IBAND,JSCAN,2,MONONE,MXONE)
!  --- COMPUTE FFTS
                     CALL FSPEC1 (IBAND, MONONE, MXONE)
                     CALL FSPEC2 (IBAND, MONONE, PHI)

                     IF( GASOUTTYPE .EQ. 1 .AND. ITER .EQ. -1 )THEN
                        WRITE(GASFNAME,690) TRIM(GASB(IBAND,NR)),IBAND,JSCAN
                     ELSEIF( GASOUTTYPE .EQ. 2 )THEN
                        IF (ITER == -1) THEN
                           WRITE(GASFNAME,690) TRIM(GASB(IBAND,NR)),IBAND,JSCAN
                        ELSE
                           WRITE(GASFNAME,700) TRIM(GASB(IBAND,NR)),IBAND,JSCAN,ITER
                        ENDIF
                     ENDIF
                     WRITE(TITLE,710) TRIM(GASB(IBAND,NR)), IBAND, JSCAN, ITER

                     OPEN(UNIT=80, FILE=GASFNAME, STATUS='REPLACE', ERR=555)
                     WRITE (80, 640) TITLE
                     WRITE (80, *) WSTART(IBAND), WSTOP(IBAND), SPAC(IBAND), N3
                     DO J = 1, N3
                        I = N1 + (J - 1)*NSPAC(IBAND)
                        WRITE (80, *) DBLE(TCONV(I))
                     ENDDO
                     CLOSE (80)
                  ENDDO

!  --- FINALLY SOLAR SPECTRUM
                  IFCO = IFCOSAVE
                  IF( IFCO )THEN
                     CALL ZERONTRAN( IBAND, 2, MONONE )
                     !  --- COMPUTE FFTS
                     CALL FSPEC1 (IBAND, MONONE, MXONE)
                     CALL FSPEC2 (IBAND, MONONE, PHI)

                     IF( GASOUTTYPE .EQ. 1 .AND. ITER .EQ. -1 )THEN
                        WRITE(GASFNAME,730)IBAND,JSCAN
                     ELSEIF( GASOUTTYPE .EQ. 2 )THEN
                        IF (ITER == -1 ) THEN
                           WRITE(GASFNAME,730)IBAND,JSCAN
                        ELSE
                           WRITE(GASFNAME,740)IBAND,JSCAN,ITER
                        ENDIF
                     ENDIF
                     WRITE(TITLE,710) 'SOLAR', IBAND, JSCAN, ITER

                     OPEN(UNIT=80, FILE=GASFNAME, STATUS='REPLACE', ERR=555)
                     WRITE (80, 640) TITLE
                     WRITE (80, *) WSTART(IBAND), WSTOP(IBAND), SPAC(IBAND), N3
                     DO J = 1, N3
                        I = N1 + (J - 1)*NSPAC(IBAND)
                        WRITE (80, *) DBLE(TCONV(I))
                     ENDDO
                     CLOSE (80)
                  ENDIF

                  GOTO 557

                  ! File Open Error handler
555               CONTINUE
                  WRITE (6, 556)

557               CONTINUE
                  IFCO=IFCOSAVE
                  ! RESTORE TCALC
                  TCALC(2,:NMONSM) = TCALC(1,:NMONSM)
                  TCONV(:NMONSM) = TCONVSAV(:NMONSM)
                  DEALLOCATE (TCONVSAV)

               ENDIF ! WRITE GASOUT FILES

               MONONE = MONONE + NM(IBAND) ! UP BAND WIDTH FOR EACH SCAN
            END DO SPEC

            MXONE = MXONE + NM(IBAND)   ! INDEX IN TCO AND CROSS ARRAYS AS START OF CURRENT BAND
         END DO BAND

         DO I = 1, NFIT
            FX = TOBS(I) - YC(I)
            SUMSQ = SUMSQ + FX*FX
         END DO

         IF (ICOUNT .EQ. 1) THEN
            RMS = 100.D0*SQRT(SUMSQ/NFIT)
            IF( ITER .LT. 0 ) THEN
               WRITE (16, 27) SNR, RMS, NVAR, NFIT
               WRITE (*, 27)  SNR, RMS, NVAR, NFIT
            ELSE
               WRITE (16, 26) ITER, SNR, RMS, NVAR, NFIT
!               WRITE (*, 28) ITER, RMS
            ENDIF

            WRITE (16, 888) (PARM(I),I=1,NVAR)

            IF( F_WRTCHANNEL )THEN
               IF( ITER == 1 )THEN       !for_channel_value_files
                  CALL FILEOPEN(30, 1 )
                  WRITE(30,'(A,A,A)' )TRIM(TAG), ' CHANNEL SPECTRUM ITER 1'
                  WRITE(30, 3000) (WAVE_X(I),AMP_Y(I),I=1,NFIT)
                  CALL FILECLOSE( 30, 1)
               ENDIF
               IF (ITER == 2) THEN       !for_channel_value_files
                  CALL FILEOPEN(40, 1 )
                  WRITE(40,'(A,A,A)' )TRIM(TAG), ' CHANNEL SPECTRUM ITER 2'
                  WRITE(40, 3000) (WAVE_X(I),AMP_Y(I),I=1,NFIT)
                  CALL FILECLOSE( 40, 1 )
               ENDIF
            ENDIF

!  --- UPDATE CALCULATED SPECTRUM ARRAY ON FIRST TRIP THROUGH UPDATE SPECTRA
            YN = YC(:NFIT)

         ELSE

!print*, yc(:5)

!  --- UPDATE ARRAY OF PARTIAL DERIVATIVES
            !IF (ANALYTIC_K.AND.XRET) THEN
            !   KN(:NFIT,IPARM) = (YC(:NFIT)-YN)/DEL
            !ELSE
               KN(:NFIT,IPARM) = (YC(:NFIT)-YN)/DEL
            !END IF
            IF( BUG1 ) &
            WRITE(0,204) '   KN: ', tret, ICOUNT, IPARM, PARM(IPARM), SUM(KN(:NFIT,IPARM))/REAL(NFIT,8), &
                         SQRT(SUM(KN(:NFIT,IPARM)**2))
            !write(0,'(10(e11.4,1x))'), kn(:nfit,iparm)
            !write(0,'(10(e11.4,1x))'),
            !if(iparm .eq. 4)write(0,'(4d22.14)') (yc(kk), yn(kk), yc(kk)-yn(kk), (yc(kk)-yn(kk))/del, kk=1,nfit)
         ENDIF

      ENDDO PARAM

! --- ZERO K MATRIX FOR MOLECULES NOT INCLUDED IN FIT OF A BAND
      BAND1: DO IBAND = 1, NBAND
         NS = NSCAN(IBAND)
         IF (NS == 0) CYCLE
         IF( IBAND > 1 ) THEN
            NS1 = SUM(NPRIM(1:IBAND-1))+1
            NS2 = SUM(NPRIM(1:IBAND))
         ELSE
            NS1 = 1
            NS2 = SUM(NPRIM(1:IBAND))
         ENDIF

         SPEC1: DO JSCAN = 1, NS
            NR = NRETB(IBAND)

            RET1: DO KK = 1, NRET
               IF( NGIDX(KK,0,IBAND) == 0 ) THEN
                 KN( NS1:NS2 , NGIDX(KK,1,0): NGIDX(KK,2,0) ) = 0.0D0
               ELSE
               ENDIF

            END DO RET1
         END DO SPEC1
      END DO BAND1

 !  --- PRINT OUT PARM ARRAY BY ITERATION
      IF( F_WRTPARM )WRITE(89,261) ITER, PARM(:NVAR)
      IF (ALLOCATED(STORE_LINE)) DEALLOCATE(STORE_LINE)


      RETURN

   17 CONTINUE
      WRITE (16, 18) N1, N2, IBAND, NSTART(IBAND), MSHIFT, MONONE, NPRIM(IBAND), NSPAC(IBAND)
      WRITE (16,*) "WAVENUMBER SHIFT OUT OF SPECTRAL RANGE."
      WRITE ( 0,*) "WAVENUMBER SHIFT OUT OF SPECTRAL RANGE."
      TFLG=.TRUE.
      RETURN

 18   FORMAT(/,' !!! ABORT !!! TCALC ARRAY OVERFLOW : ',/,' N1    =',I10, &
         ' N2    =',I10,' IBAND =',I6,/,' NSTART=',I6,' MSHIFT=',I10, &
         ' MONONE=',I6,/,' NPRIM =',I6,' NSPAC =',I6)
 26   FORMAT(/,' ITER=',I2,' AVGSNR=',F12.4,' RMS(%)=',F10.7,' NVAR=',I3,' NFIT=',I6)
 27   FORMAT(/,' FINAL:   AVGSNR=',F12.4,' RMS(%)=',F10.7,' NVAR=',I3,' NFIT=',I6)
! 28   FORMAT(/,/,' ITER=',I2,' RMS(%)=',F10.7)

 !162  FORMAT(/,' EFFECTIVE APODIZATION PARAMETER =',F8.3)
 !163  FORMAT(F8.3)
  261 FORMAT( I5, 2000ES26.18 )
 !201  FORMAT(A, 2I5, 2X, E14.6, 3F15.7)
 202  FORMAT(A, 5I5, 2X, A10, E14.6, 3F14.8)
 203  FORMAT(A, I5, 2X, A, 3F14.6)
 204  FORMAT(A, l1, 2I5, 2X, E14.6, 3D16.8)
 556  FORMAT(/,' COULD NOT CREATE INDIVIDUAL GAS FILE')
 610  FORMAT('spc.all.',I2.2,'.',I2.2,'.final')
 620  FORMAT('spc.all.',I2.2,'.',I2.2,'.',I2.2)
 !630  FORMAT('SFIT4 ALLGASES file')
 640  FORMAT(A80)
 690  FORMAT('spc.',a,'.',I2.2,'.',I2.2,'.final')
 700  FORMAT('spc.',a,'.',I2.2,'.',I2.2,'.',I2.2)
 710  FORMAT('GAS ',a7,' BAND ', I2, ' SCAN ', I2, ' ITER ', I3)
 730  FORMAT('spc.sol.',I2.2,'.',I2.2,'.final')
 740  FORMAT('spc.sol.',I2.2,'.',I2.2,'.',I2.2)
! 750  FORMAT('GAS SOLAR',' BAND ', I2, ' SCAN ', I2, ' ITER ', I3)

  888 FORMAT(5(1P,E14.7,1X))
 3000 FORMAT(2(E14.6))

      END SUBROUTINE FM

      END MODULE FRWDMDL
