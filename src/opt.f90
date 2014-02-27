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

      MODULE OPT

      USE PARAMS
      USE DATAFILES
      USE FRWDMDL
      USE SYNSPEC
      USE MATRIX
      USE WRITEOUT
      USE BANDPARAM

      IMPLICIT NONE

      LOGICAL  :: TFLG
      LOGICAL  :: ALL_SPEC_OUT

      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: SE
      REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: SA, SAINV, SHAT
      REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: KS, KSK
      REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: KHAT
      REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: G_LM

      LOGICAL        :: F_LM = .FALSE.          ! TRUE FOR LEVENBERG MARQUARDT
      REAL (DOUBLE)  :: GAMMA_START, GAMMA_INC, GAMMA_DEC, STOP_CRITERION, CHI_2_Y
      REAL (DOUBLE)  :: CONVERGENCE

      CONTAINS

      SUBROUTINE OPT_3(Y, XA, XHAT, YHAT, M, N, CONVERGE, MAXITER, TOL, RETFLG, DIVWARN, ITER, ISMIX, NLEV  )

! 17SEP02
!    - OUTPUT FORMAT OF KFILE 2000E16.8 - COMPATIBLE W/ IDL
!      32676 CHARS PER LINE

!
!
! VERSION 2; 1/17/91; BJC
! VERSION X; 6/2/95; BJC
! VERSION 3; 10/12/95; BJC
! RETRIEVAL/FORWARD MODEL SWITCH (RETFLG) ADDED
! VERSION 4; 8/13/98; PSM
! DIVWARN AND ITER PASSED OUT, & WRTK OPTION ADDED
!
! KFLG TELLS FM WHETHER TO CALCULATE WEIGHTING FUNCTIONS - BJC 12/96
!
! CHANGE FROM VERSION 1: SE IS NOW A VECTOR; A NEW MATRIX ARITHMETIC
! SUBROUTINE (MULTDIAG) HAS BEEN ADDED

      LOGICAL, INTENT(INOUT)      :: CONVERGE, RETFLG, DIVWARN
      INTEGER, INTENT(IN)         :: N, M, ISMIX, MAXITER, NLEV
      INTEGER, INTENT(INOUT)      :: ITER
      REAL(DOUBLE), INTENT(IN)    :: TOL
      REAL(DOUBLE), INTENT(INOUT) :: XA(N) ! PARM ON INPUT
      REAL(DOUBLE), INTENT(INOUT) :: XHAT(N)
      REAL(DOUBLE), INTENT(INOUT) :: Y(M)
      REAL(DOUBLE), INTENT(INOUT) :: YHAT(M)

      LOGICAL        :: KFLG, FILOPEN, PRTFLG, OPT_DEBUG, DEBUG
      REAL(DOUBLE)   :: YN_OLD(M), SEINVDY(M), SEINVDY_OLD_SE(M)
      REAL(DOUBLE)   :: XN_OLD(N), SAINVDX(N)
      REAL(DOUBLE)   :: GAMMA, RED_GAMMA, INC_GAMMA
      REAL(DOUBLE)   :: CHI_2, CHI_2_X, CHI_2_Y_OLD_SE
      REAL(DOUBLE)   :: CHI_2_OLD, D_CHI_2, CHI_2_OLD_SE, D_CHI_2_OLD_SE

      INTEGER      :: I, J, ONE, SGN, KK, NS, IYDX1, IYDX2, JSCAN, IBAND, NAERR
      REAL(DOUBLE) :: RMSDELY, RMSDY, RMSIM1, CHGY, UNCY, SQRMS, VQRMS, RMSKDX

      REAL(DOUBLE), DIMENSION(:), ALLOCATABLE   :: SEINV, SEINV_OLD, DY, DELY, DYMKDX
      REAL(DOUBLE), DIMENSION(:), ALLOCATABLE   :: XNP1, XN, DX, DX_OLD, KSDYMKDX, DELX

      ! NEED TWO EXTRA VECTORS OF SIZE N FOR LEVENBERG MARQUARDT -- MP
      REAL(DOUBLE), DIMENSION(:), ALLOCATABLE   :: KSDYMKDX_LM, GSAINVDX, G
      REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: GSAINV, SPKSKINV

      ! FOR CALCULATING AVK WHEN LM
      REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: IDNN, T2, T3, T4
      REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: T1, T5

      REAL(DOUBLE), DIMENSION(:),  ALLOCATABLE :: KN_OLD, KN

      ! EQUIVALANCED
      REAL(DOUBLE), DIMENSION(MMAX)            :: YN, KDX
      REAL(DOUBLE), DIMENSION(NMAX*NMAX)       :: SPKSK, SHATINV
      REAL(DOUBLE), DIMENSION(NMAX*MMAX)       :: KNT, KHATT

      EQUIVALENCE (KNT, KHATT), (SPKSK, SHATINV), (YN, KDX)
      ! REMOVED (DX, KSDYMKDX), BECAUSE NEED BOTH AT THE SAME TIME FOR
      ! LEVENBERG MARQUARDT, (DX, KSDYMKDX)

      DATA ONE/ 1 /

      ALLOCATE( KSDYMKDX_LM(N), GSAINVDX(N), G(N) )
      ALLOCATE( GSAINV(N,N), IDNN(N,N), T2(N,N), T3(N,N), T4(N,N) )
      ALLOCATE( T1(N,M), T5(N,M), KN_OLD(M*N), KN(M*N) )
      ALLOCATE( SEINV(M), SEINV_OLD(M), DY(M), DELY(M), DYMKDX(M) )
      ALLOCATE( XNP1(N), XN(N), DX(N), DX_OLD(N), KSDYMKDX(N), DELX(N) )
      ALLOCATE( KS(N,M), KSK(N,N), SPKSKINV(N,N), STAT=NAERR )
      ALLOCATE( KHAT(M,N), G_LM(N,M), STAT=NAERR )
      ALLOCATE( SNR_CLC(NBAND, MAXVAL(NSCAN(:NBAND))) )

      YN = 0.0D0
      DY = 0.0D0
      XN = 0.0D0
      DX = 0.0D0
      KN = 0.0D0
      DELX = 0.0D0
      XNP1 = 0.0D0
      XN_OLD = 0.0D0
      FILOPEN   = .FALSE.
      DEBUG     = .FALSE.
      PRTFLG    = .FALSE.
      DIVWARN   = .FALSE.
      CONVERGE  = .FALSE.
      OPT_DEBUG = .FALSE.

      G_LM = 0.0
      IDNN = 0.0
      DO I=1,N
         IDNN(I,I) = 1.0
      ENDDO

! --- LEVENBERG MARQUARDT METHOD MP
      IF( F_LM )THEN
         GAMMA = GAMMA_START
         RED_GAMMA = GAMMA_DEC
         INC_GAMMA = GAMMA_INC
      ELSE
         GAMMA = 0.0
         RED_GAMMA = 0.0
         INC_GAMMA = 0.0
      ENDIF

      IF( DEBUG )WRITE (*, *) ' M,N:', M, N
      KFLG = .TRUE.

      RMSDELY = 0.D0
      DO I = 1, M
         DELY(I) = SQRT(SE(I))
         RMSDELY = RMSDELY + SE(I)
      END DO

      RMSDELY = SQRT(RMSDELY/M)
      IF (DEBUG) WRITE (0, *) 'RMSDELY= ', RMSDELY

! --- XN WORKING VERSION OF XA (PARM)
      XN(:N) = XA(:N)
      CALL INVRT( SA, SAINV, N )

! --- SUBSTITUTE SAINV VALUES FOR ANY PROFILE RETRIEVE GAS FROM FILE "SA.INPUT"
! --- LOOP OVER RETRIEVAL GASES BUT CALL ONCE IF NEEDED
      DO KK = 1, NRET
! --- PICK OUT GASES WITH SET FLAG (IFOFF=5)
         IF ( CORRELATE(KK) .AND. IFPRF(KK) .AND. ( IFOFF(KK) .EQ. 5 )) THEN
            CALL GETSAINV( ISMIX )
            EXIT
         ENDIF
      ENDDO

      SEINV(:M) = 1.D0/SE(:M)
      ITER = 0
      RMSDY = 0.D0

! --- FOR FORWARD MODEL CALC. ONLY:
      IF(( .NOT. RETFLG ) .OR. MAXITER .EQ. 0 )THEN

         IF( .NOT. F_WRTK )THEN
            KFLG = .FALSE.
         ENDIF
         ITER = -1
         XHAT(:N) = XA(:N)
         SHAT(:N, :N) = SA(:N, :N)
         CALL FM (XHAT, YHAT, KHAT, M, N, KFLG, ITER, TFLG )
         WRITE (0, *) 'OPT RETURNING EARLY...'

         IF ( ALL_SPEC_OUT ) THEN
            WRITE(91,*) Y(:M)
            WRITE(91,*) YHAT(:M)
            WRITE(91,*)
         END IF

         IF( F_WRTK )THEN
            CALL FILEOPEN( 66, 2 )
            WRITE(66,*) TRIM(TAG), ' K MATRIX M SPECTRA ROWS X N PARAM COLUMNS'
            WRITE(66,*) M, N, ISMIX, NLEV
            WRITE(66,260) ADJUSTR(PNAME(:N))
            DO I = 1, M
               WRITE(66,261) (KHAT(I,J),J=1,N)
            END DO
            CALL FILECLOSE( 66, 1 )
         ENDIF

         IF( F_WRTSEINV )THEN
            CALL FILEOPEN( 67, 2 )
            WRITE(67,*) TRIM(TAG), ' SEINV (DIAGONAL) M X 1'
            WRITE(67,*) M, 1
            WRITE(67,261) (SEINV(I),I=1,M)
            CALL FILECLOSE(67,1)
         ENDIF

         IF( F_WRTSAINV )THEN
            CALL FILEOPEN( 69, 2 )
            WRITE(69, *) TRIM(TAG), ' SAINV N X N (BLOCK DIAGONAL) '
            WRITE(69,*) N, N
            WRITE(69,260) ADJUSTR(PNAME(:N))
            DO I = 1, N
               WRITE (69,261) (SAINV(I,J), J=1, N)
            END DO
            CALL FILECLOSE( 69, 1 )
         ENDIF

         RETURN

      ENDIF

      IF( F_WRTLM .AND. F_LM )THEN
         CALL FILEOPEN( 70, 2 )
         WRITE( 70,* ) TRIM(TAG),' LM OUTPUT DETAIL'
         WRITE( 70,* ) "GAMMA                 = ", GAMMA
         WRITE( 70,* ) "RED GAMMA             = ", RED_GAMMA
         WRITE( 70,* ) "INC GAMMA             = ", INC_GAMMA
         WRITE(70, '(7(A10,1X))') 'CHI_2_X', 'CHI_2_Y', 'CHI_2', 'D_CHI_2', &
                    'CHI_2_LIN', 'DIFF_LIN_EXACT', 'RATIO_D_CHI'
      END IF

      CHI_2_OLD = HUGE(CHI_2_OLD)
      SEINV_OLD = SEINV
!
! ITERATION LOOP
!
      WRITE(16,26) SNR, N, M
      WRITE(6,26) SNR, N, M

      WRITE(6,313) 'ITER', 'RMS', 'GAMMA','CHI_2_X','CHI_2_Y','CHI_2','CHI_2_OLD','D_CHI_2'

   10 CONTINUE
      ITER = ITER + 1
      IF (DEBUG) WRITE (*, *) ' ITERATION:', ITER

! --- CALC YN, KN
      CALL FM (XN, YN, KN, M, N, KFLG, ITER, TFLG )
      !PRINT *,'RTN FM'
      IF ( ALL_SPEC_OUT ) THEN
         WRITE(91,*) Y(:M)
         WRITE(91,*) YN(:M)
         WRITE(91,*)
      END IF

      IF (TFLG) RETURN

      IF (DEBUG) THEN
         WRITE (0, *) 'B- RET FROM FM'
         WRITE (0, '(6(ES12.4))') (YN(I),I=1,M)
      ENDIF

! --- CALC X(N+1)
      SGN = -1
! --- NEEDED SOMEWHERE ELSE, CALCULATED HERE TO SAVE SOME TIME
      CALL ADDSUB( Y, YN, DY, M, ONE, SGN )
      CALL TRNSPS( KN, KNT, N, M )
      CALL ADDSUB( XN, XN_OLD, DX_OLD, N, ONE, SGN )
      CALL ADDSUB( XA, XN, DX, N, ONE, SGN )

      CALL MULT( SAINV(:N,:N), DX, SAINVDX, N, N, ONE )
      !CALL MULT( DX, SAINVDX, CHI_2_X, ONE, N, ONE )
      CHI_2_X = DOT_PRODUCT( DX(:N), SAINVDX(:N) )
      ! THE EXPECTED VALUE OF THE COST FUNCTION IS M AT THE MINIMUM
      ! (RODGERS, 2000, P89)
      CHI_2_X = CHI_2_X / M

      SEINVDY(:M) = SEINV(:M)*DY(:M)
      !CALL MULT( DY, SEINVDY, CHI_2_Y, ONE, M, ONE)
      CHI_2_Y = DOT_PRODUCT( DY(:M), SEINVDY(:M) )
      CHI_2_Y = CHI_2_Y / M

      SEINVDY_OLD_SE(:M) = SEINV_OLD(:M)*DY(:M)
      !CALL MULT( DY, SEINVDY, CHI_2_Y, ONE, M, ONE)
      CHI_2_Y_OLD_SE = DOT_PRODUCT( DY(:M), SEINVDY_OLD_SE(:M) )
      CHI_2_Y_OLD_SE = CHI_2_Y_OLD_SE / M

      CHI_2 = CHI_2_X + CHI_2_Y
      CHI_2_OLD_SE = CHI_2_X + CHI_2_Y_OLD_SE

      IF (ITER.GT.1)THEN
         D_CHI_2 = CHI_2_OLD - CHI_2
         D_CHI_2_OLD_SE = CHI_2_OLD - CHI_2_OLD_SE
         DO IBAND=1,NBAND
            IF( IFCALCSE ) THEN
               WRITE(*,314) ITER, RMS, GAMMA, CHI_2_X, CHI_2_Y,        CHI_2,        CHI_2_OLD, D_CHI_2
               WRITE(*,315)                 CHI_2_Y_OLD_SE, CHI_2_OLD_SE,            D_CHI_2_OLD_SE
               PRTFLG = .TRUE.
               EXIT
            ENDIF
         ENDDO
         IF (.NOT.PRTFLG) THEN
            WRITE(*,314) ITER, RMS, GAMMA, CHI_2_X, CHI_2_Y, CHI_2, CHI_2_OLD, D_CHI_2
         ENDIF
      ELSE
         WRITE(*,314) ITER, RMS, GAMMA, CHI_2_X, CHI_2_Y
      ENDIF

      WRITE(16, '(A)') 'COST FUNCTION'
      WRITE(16, '(4(A,1X))') 'CHI_2_X', 'CHI_2_Y', 'CHI_2', 'CHI_2_OLD-CHI_2'
      WRITE(16, '(4ES11.3)') CHI_2_X, CHI_2_Y, CHI_2, D_CHI_2

      IF( ITER.GT.1 &
         .AND. (CONVERGENCE .GT. 0.0) &
         .AND. (D_CHI_2_OLD_SE .LT. CONVERGENCE) &
         .AND. ((CHI_2_OLD.GT.CHI_2_OLD_SE) .OR. (ABS(D_CHI_2_OLD_SE).LT.1.0E-5) .OR. (ABS(D_CHI_2).LT. 1.0E-5))) THEN
         ! CONVERGED
         !PRINT*, "CONVERGE = .TRUE."
         CONVERGE = .TRUE.
         GOTO 20
      END IF

      IF( F_LM .AND. ( ITER .GT. 1 ))THEN

         IF( F_WRTLM )WRITE(70, '(4ES11.3)') CHI_2_X, CHI_2_Y, CHI_2, D_CHI_2

         IF( CHI_2_OLD .GT. CHI_2_OLD_SE )THEN
            GAMMA = GAMMA / RED_GAMMA
            SEINV_OLD(:M) = SEINV(:M)
         ELSE
            GAMMA     = GAMMA * INC_GAMMA
            XN(:N)    = XN_OLD(:N)
            YN(:M)    = YN_OLD(:M)
            KN(:M*N)  = KN_OLD(:M*N)
            CHI_2     = CHI_2_OLD
            SEINV(:M) = SEINV_OLD(:M)
            ! CHI_2_LIN = CHI_2_LIN_OLD;
! --- THOSE HAVE TO BE RECALULATED, CAN SURELY BE SHORTCUT A LITTLE BIT.
            SGN = -1
            CALL ADDSUB (Y, YN, DY, M, ONE, SGN)
! --- HAS BEEN SOMEWHERE ELSE, JUST SHIFTED, NEED KNT LATER
            CALL TRNSPS (KN, KNT, N, M)
            CALL ADDSUB (XA, XN, DX, N, ONE, SGN)
         END IF
      ELSE
         SGN = -1
         SEINV_OLD(:M) = SEINV(:M)
      END IF

! --- END CALCULATION COST FUNCTION
! --- KEEP OLD STATE INFORMATION -- MP
      XN_OLD(:N)   = XN(:N)
      YN_OLD(:M)   = YN(:M)
      KN_OLD(:M*N) = KN(:M*N)
      CHI_2_OLD    = CHI_2;
     ! CHI_2_LIN_OLD = CHI_2_LIN;

      RMSIM1 = RMSDY
      RMSDY = DOT_PRODUCT(DY(:M),DY(:M))
      RMSDY = SQRT(RMSDY/M)
      IF(( ITER .GT. 1 ) .AND.(( RMSDY-RMSIM1 ) > 2*RMSDELY ))THEN
         WRITE(6, *) RMSDY, RMSIM1, RMSDELY
         WRITE(6, *) 'THREAT OF DIVERGENCE AFTER', ITER, ' ITERATIONS'
         DIVWARN = .TRUE.
      ENDIF

      IF( DEBUG )THEN
         WRITE (0, *) 'XA:'
         WRITE (0, '(6(ES12.4))') XA
         WRITE (0, *) 'XN:'
         WRITE (0, '(6(ES12.4))') XN
         WRITE (0, *) 'DX:'
         WRITE (0, '(6(ES12.4))') DX
         WRITE (0, *) 'YN:'
         WRITE (0, '(6(ES12.4))') YN
         WRITE (0, *) 'DY:'
         WRITE (0, '(6(ES12.4))') DY
      ENDIF

      CALL MULT (KN, DX, KDX, M, N, ONE)

      CALL ADDSUB (DY, KDX, DYMKDX, M, ONE, SGN)

      IF( IFCALCSE )THEN
         BND: DO IBAND = 1, NBAND
            NS = NSCAN(IBAND)
            IF (NS == 0) CYCLE
            SPC: DO JSCAN = 1, NS
               IYDX1 = ISCNDX(1,IBAND,JSCAN)
               IYDX2 = ISCNDX(2,IBAND,JSCAN)
               SQRMS = (IYDX2-IYDX1)/DOT_PRODUCT(DY(IYDX1:IYDX2),DY(IYDX1:IYDX2))
               VQRMS = SQRT(1.D0/SQRMS)
               SEINV(IYDX1:IYDX2) = SQRMS
               DELY(IYDX1:IYDX2)  = VQRMS
               WRITE(*,305) IBAND, JSCAN, SQRT(SQRMS)
            ENDDO SPC
         ENDDO BND
      ENDIF
      !RMSDELY = 1.D0/SUM(SEINV(:M))
      !RMSDELY = SQRT(RMSDELY/M)
      !WRITE(*,300)' RMS:DELY = ',RMSDELY, 'AVGSNR = ',SNR,' TOL = ',TOL

!203   CONTINUE

      CALL MULTDIAG (KNT, SEINV, KS(:N,:M), N, M)
      CALL MULT (KS(:N,:M), DYMKDX, KSDYMKDX, N, M, ONE)
      CALL MULT (KS(:N,:M), KN, KSK, N, M, N)

      IF( DEBUG )THEN
!          WRITE(*,'(6(ES12.4))')'KS',KS
!          WRITE(*,'(6(ES12.4))')'KSDYMKDX:',KSDYMKDX
!          WRITE(*,*)'KSK:',KSK
      ENDIF

      SGN = 1
! --- USE LEVENBERG MARQUARDT METHOD -- MP
      IF( F_LM )THEN
         G(:N) = GAMMA
         CALL MULTDIAG ( SAINV(:N,:N), G, GSAINV(:N,:N), N, N )
         CALL MULT ( GSAINV(:N,:N), DX, GSAINVDX, N, N, ONE )
         CALL ADDSUB ( KSDYMKDX, GSAINVDX, KSDYMKDX_LM, N, ONE, -1 )
         KSDYMKDX(:N) = KSDYMKDX_LM(:N)

         G(:N) = 1.0D0 + GAMMA
         CALL MULTDIAG ( SAINV(:N,:N), G, GSAINV(:N,:N), N, N )
         CALL ADDSUB ( GSAINV(:N,:N), KSK, SPKSK, N, N, SGN )
         CALL INVRT (SPKSK, SPKSKINV, N)

         ! CALCULATE T FOR AVK
         ! (CECCHERINI & RIDOLFI, ACP, 10, 3131-3139, 2009)
         CALL MULT(SPKSKINV, KS, T1,N,N,M)
         CALL MULT(SPKSKINV, KSK, T2,N,N,N)
         CALL MULT(SPKSKINV, SAINV, T3,N,N,N)
         CALL ADDSUB(IDNN,T2,T4,N,N,-1)
         CALL ADDSUB(T4,T3,T2,N,N,-1)         ! T2 IS RECYCLED
         CALL MULT(T2,G_LM,T5,N,N,M)
         CALL ADDSUB(T1,T5,G_LM,N,M,1)           ! T IS UPDATED
      ELSE
         CALL ADDSUB( SAINV(:N,:N), KSK, SPKSK, N, N, SGN )
         CALL INVRT ( SPKSK, SPKSKINV, N )
      END IF

      CALL MULT (SPKSKINV, KSDYMKDX, DELX, N, N, ONE)
      CALL ADDSUB (XA, DELX, XNP1, N, ONE, SGN)

      IF( DEBUG )THEN
!          WRITE(*,*)'SPKSK:',SPKSK
!          WRITE(*,*)'SPKSKINV',SPKSKINV
         WRITE (0, *) 'DELX:'
         WRITE (0, '(6(ES12.4))') DELX
         WRITE (0, *) ' XNP1: '
         WRITE (0, '(6(ES12.4))') XNP1(:N)
      ENDIF

! --- CHECK CONVERGENCE
      SGN = -1
      IF( CONVERGENCE .GT. 0.0 )THEN
         CONVERGE = .FALSE.
      ELSE
         CALL ADDSUB (XNP1, XN, DX, N, ONE, SGN)
         IF( DEBUG )THEN
            WRITE (0, *) ' DX: '
            WRITE (0, '(6(ES12.4))') DX(:N)
         ENDIF
         CALL MULT (KN, DX, KDX, M, N, ONE)
         IF( DEBUG )THEN
            WRITE (0, *) 'TOL:'
            WRITE (0, '(6(ES12.4))') TOL
         ENDIF

         RMSKDX = 0.0D0
         DO I=1,M
            RMSKDX = RMSKDX + KDX(I)*KDX(I)
         ENDDO
         RMSKDX = SQRT(RMSKDX/M)
         WRITE(*,300)' RMS:KDX = ', RMSKDX

         IF( DEBUG )THEN
            DO I = 1, M
               WRITE (0, *) 'KDX(I),DELY(I):', KDX(I), DELY(I)
               CHGY = ABS(KDX(I))
               UNCY = ABS(TOL*DELY(I))
               IF( CHGY .GT. UNCY )THEN
                  CONVERGE = .FALSE.
                  EXIT
               ENDIF
               CONVERGE = .TRUE.
               GOTO 20 
            END DO
         ELSE
            DO I = 1, M
               CHGY = ABS(KDX(I))
               UNCY = ABS(TOL*DELY(I))
               IF( CHGY .GT. UNCY )THEN
                  CONVERGE = .FALSE.
                  EXIT
               ENDIF
               CONVERGE = .TRUE.
               GOTO 20
            END DO
         ENDIF
      END IF

! --- CHECK IF ITERATIONS EXCEEDED WITH NO CONVERGERNCE
      IF( ITER .GE. MAXITER )THEN
         IF( F_LM )THEN
            XNP1(:N) = XN_OLD(:N)
         END IF
         WRITE(6,307) ITER
         SNR = SQRT(SUM(SEINV(:M))/M)
         CONVERGE = .FALSE.
         GO TO 20
      ELSE
         XN(:N) = XNP1(:N)
         GOTO 10
      ENDIF

!  --- CONVERGED

   20 CONTINUE

      !PRINT*, D_CHI_2_OLD_SE, CHI_2_OLD, CHI_2_OLD_SE, ABS(D_CHI_2_OLD_SE), ABS(D_CHI_2)


      XHAT(:N) = XNP1(:N)
      CALL FM (XHAT, YHAT, KHAT, M, N, KFLG, -1, TFLG )

      IF( F_WRTK )THEN
         CALL FILEOPEN( 66, 2 )
         WRITE(66,*) TRIM(TAG), ' K MATRIX M SPECTRA ROWS X N PARAM COLUMNS'
         WRITE(66,*) M, N, ISMIX, NLEV
         WRITE(66,260) ADJUSTR(PNAME(:N))
         DO I = 1, M
            WRITE(66,261) (KHAT(I,J),J=1,N)
         END DO
         CALL FILECLOSE( 66, 1 )
      ENDIF

      IF( F_WRTSEINV )THEN
         CALL FILEOPEN( 67, 2 )
         WRITE(67,*)  TRIM(TAG), ' SEINV (DIAGONAL) M X 1'
         WRITE(67,*) M, 1
         WRITE(67,261) (SEINV(I),I=1,M)
         CALL FILECLOSE(67,1)
      ENDIF

      IF( F_WRTSAINV )THEN
         CALL FILEOPEN( 69, 2 )
         WRITE(69,*)  TRIM(TAG), ' SAINV N X N (BLOCK DIAGONAL)'
         WRITE(69,*) N, N
         WRITE(69,260) ADJUSTR(PNAME(:N))
         DO I = 1, N
            WRITE(69,261) (SAINV(I,J), J=1, N)
         END DO
         CALL FILECLOSE( 69, 1 )
      ENDIF

      IF( DEBUG )THEN
         WRITE (0, *) 'YHAT:'
         WRITE (0, '(6(ES12.4))') YHAT(1:M)
         WRITE (0, *) 'XHAT:'
         WRITE (0, '(6(ES12.4))') XHAT(1:N)
      ENDIF

      CALL TRNSPS (KHAT, KHATT, N, M)
      CALL MULTDIAG (KHATT, SEINV, KS(:N,:M), N, M)
      CALL MULT (KS(:N,:M), KHAT, KSK, N, M, N)
      CALL ADDSUB( Y, YN, DY, M, ONE, -1 )
      CHI_2_Y = DOT_PRODUCT( DY(:M), SEINVDY(:M) )
      CHI_2_Y = CHI_2_Y / M

!  --- CALCULATE THE SNR FOR THE RESULT PER BAND AND SCAN
      WRITE(6,303)
      DO IBAND = 1, NBAND
         NS = NSCAN(IBAND)
         IF (NS == 0) CYCLE
         DO JSCAN = 1, NS
            IYDX1 = ISCNDX(1,IBAND,JSCAN)
            IYDX2 = ISCNDX(2,IBAND,JSCAN)
            SQRMS = (IYDX2-IYDX1)/DOT_PRODUCT(DY(IYDX1:IYDX2),DY(IYDX1:IYDX2))
            SNR_CLC(IBAND,JSCAN) = SQRT(SQRMS)
            WRITE(*,305) IBAND, JSCAN, SNR_CLC(IBAND,JSCAN)
         ENDDO
      ENDDO

      SGN = 1
      IF( F_LM )THEN
         ! CALCULATE T (Gain matrix G) FOR AVK ITERATIVELY
         ! (CECCHERINI & RIDOLFI, ACP, 10, 3131-3139, 2009)
         CALL MULTDIAG ( SAINV(:N,:N), G, GSAINV(:N,:N), N, N )
         CALL ADDSUB ( GSAINV(:N,:N), KSK, SPKSK, N, N, SGN )
         CALL INVRT (SPKSK, SPKSKINV, N)
         CALL MULT ( SPKSKINV, KS, T1,N,N,M)
         CALL MULT ( SPKSKINV, KSK, T2,N,N,N)
         CALL MULT ( SPKSKINV, SAINV, T3,N,N,N)
         CALL ADDSUB ( IDNN,T2,T4,N,N,-1)
         CALL ADDSUB ( T4,T3,T2,N,N,-1)         ! T2 IS RECYCLED
         CALL MULT ( T2,G_LM,T5,N,N,M)
         CALL ADDSUB ( T1,T5,G_LM,N,M,1)           ! T IS UPDATED
      END IF
      CALL ADDSUB( SAINV(:N,:N), KSK, SHATINV, N, N, SGN)
      CALL INVRT( SHATINV, SHAT(:N,:N), N)

      IF( F_WRTSHAT )THEN
         CALL FILEOPEN( 64, 1 )
         WRITE(64,*)  TRIM(TAG), ' SHAT N X N (BLOCK DIAGONAL)'
         WRITE(64,*) N, N
         WRITE(64,260) ADJUSTR(PNAME(:N))
         DO I=1,N
            WRITE(64,261) (SHAT(I,J), J=1, N)
         END DO
         CALL FILECLOSE( 64, 1 )
      ENDIF

      IF( F_WRTLM .AND. F_LM )CALL FILECLOSE( 70, 1 )

!  --- DEALLOCATE LOCAL ARRAYS
      DEALLOCATE( KSDYMKDX_LM, GSAINVDX, G )
      DEALLOCATE( GSAINV, IDNN, T2, T3, T4 )
      DEALLOCATE( T1, T5, KN_OLD, KN )
      DEALLOCATE( SEINV, SEINV_OLD, DY, DELY, DYMKDX )
      DEALLOCATE( XNP1, XN, DX, DX_OLD, KSDYMKDX, DELX )
      DEALLOCATE( SPKSKINV )

      RETURN

 26   FORMAT(/,' AVGSNR=',F12.4,' NVAR=',I3,' NFIT=',I6,/)

 260 FORMAT( 2000( 12X, A14 ))
 261 FORMAT( 2000ES26.18 )
 300  FORMAT( 3(A16, ES11.4 ))
 303 FORMAT(/,'   BAND   SCAN   CALCULATED RMSSNR')
! 304  FORMAT( "    BAND    SCAN      SEINV         DELY         SNR" )
 305  FORMAT( 2I7,F20.2 )
! 306  FORMAT( A20, 2ES12.4 )
 307  FORMAT(/, ' NO CONVERGENCE AFTER', I4, ' ITERATIONS')
! 308  FORMAT( A20, ES12.4 )
 313  FORMAT(A4,1X,14(A9,1X))
 314  FORMAT(I4,1X,F9.4,1X,ES9.2,1X,12(F9.3,1X))
 315  FORMAT(2(F12.6,1X)13x,9F12.6)

      END SUBROUTINE OPT_3


SUBROUTINE GETSAINV( ISMIX )

      REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: SAINP
      REAL(DOUBLE), DIMENSION(:),   ALLOCATABLE :: SAROOT
      LOGICAL ::FILOPEN
      INTEGER :: ISMIX, KK, I, J, INDXX, NAERR

      FILOPEN   = .FALSE.

      INDXX = ISMIX
      DO KK = 1, NRET
! --- PICK OUT GASES WITH SET FLAG (IFOFF=5)
         IF (( IFPRF(KK) ) .AND. ( IFOFF(KK) == 5 )) THEN
            ALLOCATE( SAINP(LAYMAX,LAYMAX), SAROOT(LAYMAX), STAT=NAERR )
            IF (NAERR /= 0) THEN
               WRITE (16, *) 'OPT : COULD NOT ALLOCATE SAINP ARRAY, ERROR NUMBER = ', NAERR
               WRITE ( 0, *) 'OPT : COULD NOT ALLOCATE SAINP ARRAY, ERROR NUMBER = ', NAERR
               CALL SHUTDOWN
               STOP 4
            ENDIF
            IF ( .NOT. FILOPEN ) THEN
               CALL FILEOPEN( 62, 3 )
               WRITE(16,301) TRIM(TFILE(62))
               FILOPEN = .TRUE.
            ENDIF
            WRITE(16,302) KK, IFOFF(KK)
! --- GET DIAGONAL FRACTIONAL VARIANCES FROM BINPUT && ALREADY SQUARED
            DO I = 1, NLAYERS
               SAROOT(I) = SQRT( SA(I+INDXX,I+INDXX) )
            END DO
! --- READ IN PUT NEW SAINV VALUES
            READ(62,*) SAINP(:NLAYERS, :NLAYERS )
! --- SCALE AND STORE THEM IN SAINV
            DO I = 1, NLAYERS
               DO J = 1, NLAYERS
                  SAINV(I+INDXX,J+INDXX) = SAINP(I,J)*( 1.0D0/(SAROOT(I)*SAROOT(J)))
               END DO
            END DO
            DEALLOCATE( SAINP, SAROOT )
         ELSE
            WRITE(16,303) KK, IFOFF(KK)
         ENDIF
! --- BUMP UP INDEX OF MIXING RATIO BLOCK IN SA(INV) MATRIX
         INDXX = INDXX + NLAYERS
      ENDDO
      IF( FILOPEN ) CALL FILECLOSE( 62, 2 )

      RETURN

 301  FORMAT(/," FILE: ", A, " OPENED IN OPT" )
 302  FORMAT("  RETRIEVAL GAS # :",I3," HAS IFOFF FLAG:", I3, &
             " ...READING IN NEW SAINV VALUES." )
 303  FORMAT("  RETRIEVAL GAS # :",I3," HAS IFOFF FLAG:", I3, " ...SKIPPING." )

      END SUBROUTINE GETSAINV

      SUBROUTINE RELEASE_MEM_OPT

! --- DEALLOCATE PUBLIC ARRAYS

      IF( ALLOCATED( SA )) DEALLOCATE( SA )
      IF( ALLOCATED( SHAT )) DEALLOCATE( SHAT )
      IF( ALLOCATED( KS )) DEALLOCATE( KS )
      IF( ALLOCATED( SE )) DEALLOCATE( SE )
      IF( ALLOCATED( G_LM )) DEALLOCATE( G_LM )
      IF( ALLOCATED( KHAT )) DEALLOCATE( KHAT )
      IF( ALLOCATED( KSK )) DEALLOCATE( KSK )
      IF( ALLOCATED( SAINV )) DEALLOCATE( SAINV )
      IF( ALLOCATED( SNR_CLC )) DEALLOCATE( SNR_CLC )
      IF( ALLOCATED( KS )) DEALLOCATE( KS )


      END SUBROUTINE RELEASE_MEM_OPT

      END MODULE OPT
